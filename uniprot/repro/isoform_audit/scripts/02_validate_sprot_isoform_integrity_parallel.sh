#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPRO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${REPRO_DIR}/../../.." && pwd)"

SHARDS_ROOT="/mnt/disks/toolkit-data/uniprot_shards/2025_04/sprot"
JOBS=40
RELEASE="2025_04"
OUT_DIR="${REPRO_DIR}/outputs_integrity"
LOG_DIR="${REPRO_DIR}/logs_integrity"

if [[ -x "${REPO_ROOT}/.venv/bin/python" ]]; then
  PYTHON_CMD="${REPO_ROOT}/.venv/bin/python"
else
  PYTHON_CMD="python3"
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Validate isoform sequence/length integrity across Swiss-Prot shards.

Options:
  --shards-root PATH     Swiss-Prot shards dir (default: ${SHARDS_ROOT})
  -j, --jobs N           Max parallel shard jobs (default: ${JOBS})
  --release TAG          Release label passed to parser (default: ${RELEASE})
  --out-dir PATH         Output dir (default: ${OUT_DIR})
  --log-dir PATH         Log dir (default: ${LOG_DIR})
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
    --python) PYTHON_CMD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${OUT_DIR}/per_shard" "${LOG_DIR}"
TASK_FILE="${LOG_DIR}/sprot_integrity_tasks.txt"
find "${SHARDS_ROOT}" -maxdepth 1 -type f -name 'sprot-shard-*.dat.gz' | sort > "${TASK_FILE}"

if [[ ! -s "${TASK_FILE}" ]]; then
  echo "[error] No shard files found in ${SHARDS_ROOT}" >&2
  exit 1
fi

echo "[info] Validating isoform integrity across $(wc -l < "${TASK_FILE}") Swiss-Prot shards (max parallel: ${JOBS})"
cat "${TASK_FILE}" | xargs -r -P "${JOBS}" -I {} bash -c '
  set -euo pipefail
  shard_file="{}"
  shard_name="$(basename "${shard_file}" .dat.gz)"
  out_json="'"${OUT_DIR}"'/per_shard/${shard_name}.json"
  log_file="'"${LOG_DIR}"'/${shard_name}.log"
  echo "[info] ${shard_name}" > "${log_file}"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    '"${PYTHON_CMD}"' "'"${SCRIPT_DIR}"'/validate_isoform_integrity_shard.py" \
      --input "${shard_file}" \
      --output "${out_json}" \
      --release "'"${RELEASE}"'" >> "${log_file}" 2>&1
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
    raise SystemExit("No per-shard outputs found.")

totals = {}
issue_sample = []
rows = []

for path in files:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    shard = os.path.basename(path).replace(".json", "")
    stats = data.get("stats", {})
    rows.append(
        {
            "shard": shard,
            "records_with_alt_products": stats.get("records_with_alt_products", 0),
            "declared_isoforms": stats.get("declared_isoforms", 0),
            "local_candidate_isoforms": stats.get("local_candidate_isoforms", 0),
            "local_emitted": stats.get("local_emitted", 0),
            "local_length_mismatch": stats.get("local_length_mismatch", 0),
            "local_seq_mismatch_reconstructed": stats.get("local_seq_mismatch_reconstructed", 0),
        }
    )
    for key, value in stats.items():
        totals[key] = totals.get(key, 0) + value
    if len(issue_sample) < 200:
        issue_sample.extend(data.get("issues_sample", []))
        issue_sample = issue_sample[:200]

with open(out_json, "w", encoding="utf-8") as handle:
    json.dump({"totals": totals, "issues_sample": issue_sample}, handle, indent=2)
    handle.write("\n")

with open(out_tsv, "w", encoding="utf-8") as handle:
    handle.write(
        "shard\trecords_with_alt_products\tdeclared_isoforms\tlocal_candidate_isoforms\t"
        "local_emitted\tlocal_length_mismatch\tlocal_seq_mismatch_reconstructed\n"
    )
    for row in rows:
        handle.write(
            f"{row['shard']}\t{row['records_with_alt_products']}\t{row['declared_isoforms']}\t"
            f"{row['local_candidate_isoforms']}\t{row['local_emitted']}\t"
            f"{row['local_length_mismatch']}\t{row['local_seq_mismatch_reconstructed']}\n"
        )
PY

echo "[info] Validation complete."
echo "[info] Summary: ${SUMMARY_JSON}"
echo "[info] Per-shard table: ${SUMMARY_TSV}"
