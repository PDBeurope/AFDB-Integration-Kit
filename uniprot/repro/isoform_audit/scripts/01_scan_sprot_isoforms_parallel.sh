#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPRO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${REPRO_DIR}/../../.." && pwd)"

SHARDS_ROOT="/mnt/disks/toolkit-data/uniprot_shards/2025_04/sprot"
JOBS=40
RELEASE="2025_04"
OUT_DIR="${REPRO_DIR}/outputs"
LOG_DIR="${REPRO_DIR}/logs"
FOCUS_ACCESSIONS="P42167 O43687 Q61029"

if [[ -x "${REPO_ROOT}/.venv/bin/python" ]]; then
  PYTHON_CMD="${REPO_ROOT}/.venv/bin/python"
else
  PYTHON_CMD="python3"
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Parallel isoform audit for Swiss-Prot shards.

Options:
  --shards-root PATH     Swiss-Prot shards dir (default: ${SHARDS_ROOT})
  -j, --jobs N           Max parallel shard jobs (default: ${JOBS})
  --release TAG          Release label passed to parser (default: ${RELEASE})
  --out-dir PATH         Output dir (default: ${OUT_DIR})
  --log-dir PATH         Log dir (default: ${LOG_DIR})
  --focus "A B C"        Space-separated primary accessions to report in detail
                         (default: "${FOCUS_ACCESSIONS}")
  --python CMD           Python executable (default: ${PYTHON_CMD})
  -h, --help             Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --shards-root) SHARDS_ROOT="$2"; shift 2 ;;
    -j|--jobs) JOBS="$2"; shift 2 ;;
    --release) RELEASE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --log-dir) LOG_DIR="$2"; shift 2 ;;
    --focus) FOCUS_ACCESSIONS="$2"; shift 2 ;;
    --python) PYTHON_CMD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${OUT_DIR}/per_shard" "${LOG_DIR}"

TASKS_FILE="${LOG_DIR}/sprot_isoform_tasks.txt"
find "${SHARDS_ROOT}" -maxdepth 1 -type f -name 'sprot-shard-*.dat.gz' | sort > "${TASKS_FILE}"

if [[ ! -s "${TASKS_FILE}" ]]; then
  echo "[error] No shard files found in ${SHARDS_ROOT}" >&2
  exit 1
fi

echo "[info] Running isoform audit across $(wc -l < "${TASKS_FILE}") Swiss-Prot shards (max parallel: ${JOBS})"
cat "${TASKS_FILE}" | xargs -r -P "${JOBS}" -I {} bash -c '
  set -euo pipefail
  shard_file="{}"
  shard_name="$(basename "${shard_file}" .dat.gz)"
  out_json="'"${OUT_DIR}"'/per_shard/${shard_name}.json"
  log_file="'"${LOG_DIR}"'/${shard_name}.log"
  echo "[info] ${shard_name}" > "${log_file}"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    '"${PYTHON_CMD}"' "'"${SCRIPT_DIR}"'/scan_isoform_shard.py" \
    --input "${shard_file}" \
    --output "${out_json}" \
    --release "'"${RELEASE}"'" \
    --focus '"${FOCUS_ACCESSIONS}"' >> "${log_file}" 2>&1
'

SUMMARY_JSON="${OUT_DIR}/summary.json"
SUMMARY_TSV="${OUT_DIR}/summary.tsv"

"${PYTHON_CMD}" - "${OUT_DIR}/per_shard" "${SUMMARY_JSON}" "${SUMMARY_TSV}" <<'PY'
import glob
import json
import os
import sys

per_shard_dir, out_json, out_tsv = sys.argv[1:4]
files = sorted(glob.glob(os.path.join(per_shard_dir, "*.json")))
if not files:
    raise SystemExit("No per-shard JSON outputs found.")

totals = {
    "records_scanned": 0,
    "records_with_alt_products": 0,
    "declared_isoforms": 0,
    "emittable_isoforms": 0,
    "emitted_isoforms": 0,
    "missing_emittable_isoforms": 0,
    "non_local_isoforms": 0,
    "displayed_not_dash1_isoforms": 0,
}
focus_examples = []
missing_examples = []
rows = []

for path in files:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    shard_name = os.path.basename(path).replace(".json", "")
    rows.append(
        {
            "shard": shard_name,
            "records_with_alt_products": data["records_with_alt_products"],
            "declared_isoforms": data["declared_isoforms"],
            "emittable_isoforms": data["emittable_isoforms"],
            "emitted_isoforms": data["emitted_isoforms"],
            "missing_emittable_isoforms": data["missing_emittable_isoforms"],
        }
    )
    for key in totals:
        totals[key] += data.get(key, 0)
    focus_examples.extend(data.get("focus_examples", []))
    missing_examples.extend(data.get("missing_examples", []))

summary = {
    "totals": totals,
    "focus_examples": focus_examples,
    "missing_examples_sample": missing_examples[:100],
}
with open(out_json, "w", encoding="utf-8") as handle:
    json.dump(summary, handle, indent=2)
    handle.write("\n")

with open(out_tsv, "w", encoding="utf-8") as handle:
    handle.write(
        "shard\trecords_with_alt_products\tdeclared_isoforms\temittable_isoforms\t"
        "emitted_isoforms\tmissing_emittable_isoforms\n"
    )
    for row in rows:
        handle.write(
            f"{row['shard']}\t{row['records_with_alt_products']}\t{row['declared_isoforms']}\t"
            f"{row['emittable_isoforms']}\t{row['emitted_isoforms']}\t{row['missing_emittable_isoforms']}\n"
        )
PY

echo "[info] Done."
echo "[info] Summary: ${SUMMARY_JSON}"
echo "[info] Per-shard table: ${SUMMARY_TSV}"
