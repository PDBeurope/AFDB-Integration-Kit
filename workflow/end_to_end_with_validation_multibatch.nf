import java.nio.file.Paths

nextflow.enable.dsl=2

/*
 * End-to-end orchestration for ColabFold outputs -> AFDB ingest artefacts.
 * Steps:
 *  - Convert ColabFold score/PDB to AFDB pLDDT/PAE JSON + per-model manifests
 *  - Merge manifests for downstream UniProt tooling
 *  - Export per-model/chain metadata JSONs and batch them
 *  - Emit ModelCIF generator input JSONs from the merged manifests/DB/template
 *  - Generate ModelCIF, annotate with DSSP, add headers back into PDB, and emit BCIF
 */

// Helper to supply defaults without triggering undefined-param warnings
def paramOrDefault(name, defaultValue) {
    def v = params.containsKey(name) ? params.get(name) : null
    if( v == null )
        return defaultValue
    def text = v.toString()
    return text.trim() ? v : defaultValue
}

// Base parameter defaults (override with --param)
POSIX_ROOT          = paramOrDefault('posix_root', '/mnt/disks/data/first_set/data').toString()
MAPPING_FILE        = paramOrDefault('mapping_file', '/mnt/disks/data/first_set/config/sample_af_mapping.tsv').toString()
SAMPLE_SIZE         = paramOrDefault('sample_size', 50)
COLABFOLD_MANIFEST  = paramOrDefault('colabfold_manifest', "${projectDir}/examples/colabfold/config/colabfold_manifest.csv").toString()
DATASET_DIR         = paramOrDefault('dataset_dir', '/mnt/disks/data/first_set').toString()
DATASET_CONFIG      = paramOrDefault('dataset_config', "${DATASET_DIR}/config/dataset_config.json").toString()
ARCHIVE_ROOT        = paramOrDefault('archive_root', '').toString()
CONVERTER_DUCKDB    = paramOrDefault('converter_duckdb', '').toString()
UNIPROT_RELEASE     = paramOrDefault('uniprot_release', '2025_04').toString()
UNIPROT_DB          = paramOrDefault('uniprot_db', file(DATASET_DIR).resolve("db/uniprot_${UNIPROT_RELEASE}.duckdb").toString()).toString()
MODEL_BATCH_PREFIX  = paramOrDefault('model_batch_prefix', 'AF-metadata').toString()
CHAIN_BATCH_PREFIX  = paramOrDefault('chain_batch_prefix', 'AF-chain-metadata').toString()
BATCH_SIZE          = paramOrDefault('batch_size', 10_000)
WORKERS             = paramOrDefault('workers', 10)  // Parallel workers for batch processing
PYTHON_CMD          = paramOrDefault('python_cmd', 'python3').toString()
MODEL_VERSION       = paramOrDefault('model_version', 'v1').toString()
BATCH_LABEL         = paramOrDefault('batch_label', '').toString().trim()
VALIDATION_OUTPUT_DIR = paramOrDefault('validation_output_dir', file(DATASET_DIR).resolve(BATCH_LABEL ? "validation/${BATCH_LABEL}" : 'validation').toString()).toString()

datasetDir = file(DATASET_DIR)
def repoDir = projectDir.parent ?: projectDir

// Allow grouping batch-specific outputs under a unique label to avoid collisions across concurrent runs
def batchSubdirPath(basePath) {
    def base = file(basePath)
    BATCH_LABEL ? base.resolve(BATCH_LABEL).toString() : base.toString()
}

// Helpers for POSIX-sharded per-accession layout
def shardChunks(String afId) {
    def m = (afId =~ /AF-(\d{16})/)
    if (!m) {
        throw new IllegalArgumentException("AF ID ${afId} does not contain a 16-digit numeric segment.")
    }
    def digits = m[0][1]
    (0..3).collect { digits.substring(it * 4, (it + 1) * 4) }
}

def shardPathForAf(String afId)     { file(Paths.get(POSIX_ROOT, *shardChunks(afId)).toString()) }

def derivedBaseDir(String afId)      { shardPathForAf(afId).resolve('derived') }
def derivedScoresDir(String afId)    { derivedBaseDir(afId).resolve('scores') }
def derivedManifestsDir(String afId) { derivedBaseDir(afId).resolve('manifests') }
def derivedMetadataDir(String afId)  { derivedBaseDir(afId).resolve('metadata') }
def derivedStructuresDir(String afId){ derivedBaseDir(afId).resolve('structures') }

// Derived paths
PER_MODEL_MANIFEST_DIR = paramOrDefault('per_model_manifest_dir', batchSubdirPath(datasetDir.resolve('per_accession/manifests/models'))).toString()
PER_CHAIN_MANIFEST_DIR = paramOrDefault('per_chain_manifest_dir', batchSubdirPath(datasetDir.resolve('per_accession/manifests/chains'))).toString()
MERGED_MANIFEST_DIR    = paramOrDefault('merged_manifest_dir', batchSubdirPath(datasetDir.resolve('config'))).toString()
PER_MODEL_JSON_DIR     = paramOrDefault('per_model_json_dir', batchSubdirPath(datasetDir.resolve('per_accession/models'))).toString()
PER_CHAIN_JSON_DIR     = paramOrDefault('per_chain_json_dir', batchSubdirPath(datasetDir.resolve('per_accession/chains'))).toString()
BATCH_MODEL_OUTDIR     = paramOrDefault('batch_model_outdir', batchSubdirPath(datasetDir.resolve('batches/models'))).toString()
BATCH_CHAIN_OUTDIR     = paramOrDefault('batch_chain_outdir', batchSubdirPath(datasetDir.resolve('batches/chains'))).toString()
MODELCIF_INPUT_DIR     = paramOrDefault('modelcif_outdir', batchSubdirPath(datasetDir.resolve('modelcif_inputs'))).toString()
MMICIF_OUTDIR          = paramOrDefault('mmcif_outdir', batchSubdirPath(datasetDir.resolve('mmcif'))).toString()
BCIF_OUTDIR            = paramOrDefault('bcif_outdir', batchSubdirPath(datasetDir.resolve('bcif'))).toString()
FINAL_PDB_OUTDIR       = paramOrDefault('pdb_outdir', batchSubdirPath(datasetDir.resolve('pdb'))).toString()
PROVIDER_JSON          = paramOrDefault('provider_json', datasetDir.resolve('config/provider.json').toString()).toString()

// projectDir points to the workflow/ dir; templates live in repo root by default.
MODELCIF_TEMPLATE   = paramOrDefault('modelcif_template', "${repoDir}/uniprot/templates/colabfold_example_modelcif_metadata.json").toString()
TOOLKIT_CMD         = paramOrDefault('toolkit_cmd', "${PYTHON_CMD} ${repoDir}/main.py").toString()

/*
 * Read up to N rows from the mapping file (tab-delimited). If N <= 0, read all.
 */
def loadMappingRowsLimited(String path, int limit) {
    def rows = []
    def f = file(path)
    if( !f.exists() )
        throw new IllegalArgumentException("Mapping file not found: ${path}")
    def reader = f.newReader("UTF-8")
    try {
        String line
        int count = 0
        while( (line = reader.readLine()) != null ) {
            if( limit > 0 && count >= limit )
                break
            if( !line.trim() )
                continue
            rows << line.split("\t", -1)
            count++
        }
    } finally {
        reader?.close()
    }
    Channel.from(rows)
}

/*
 * Build channel of AF IDs with meta_v1 JSON + model_v1 PDB from the POSIX-sharded layout.
 * Limits to SAMPLE_SIZE rows from the mapping file (default 50 for test runs).
 */
def buildAfInputChannel() {
    def sampleLimit = (SAMPLE_SIZE as int)
    loadMappingRowsLimited(MAPPING_FILE, sampleLimit)
        .map { row ->
            def afDisplay = row[0]?.toString()?.trim()
            def afNumeric = (row.size() > 1 ? row[1] : afDisplay)?.toString()?.trim()
            def afId = afNumeric ?: afDisplay
            def shardDir = shardPathForAf(afId)
            def flatMeta    = ARCHIVE_ROOT ? file(Paths.get(ARCHIVE_ROOT, "${afId}-meta_v1.json").toString()) : ''
            def flatPdb     = ARCHIVE_ROOT ? file(Paths.get(ARCHIVE_ROOT, "${afId}-model_v1.pdb").toString()) : ''
            tuple(afId, shardDir, flatMeta, flatPdb)
        }
}

process PREPARE_AF_ASSETS {
    tag { model_id }

    input:
    tuple val(model_id), val(shard_dir), val(flat_meta_src), val(flat_pdb_src)

    output:
    tuple val(model_id), val(shard_dir), path("${model_id}-meta_v1.json"), path("${model_id}-model_v1.pdb"), emit: prepared_assets

    script:
    def shardPath = file(shard_dir)
    def metaOnDisk = shardPath.resolve("${model_id}-meta_v1.json")
    def pdbOnDisk  = shardPath.resolve("${model_id}-model_v1.pdb")
    def flatMeta = flat_meta_src ? file(flat_meta_src) : ''
    def flatPdb  = flat_pdb_src ? file(flat_pdb_src) : ''
    """
    set -euo pipefail
    shard_dir=${shardPath}
    meta_path=${metaOnDisk}
    pdb_path=${pdbOnDisk}
    flat_meta=${flatMeta}
    flat_pdb=${flatPdb}
    mkdir -p "${shard_dir}"

    if [[ ! -f "\${meta_path}" || ! -f "\${pdb_path}" ]]; then
        if [[ ! -f "\${flat_meta}" || ! -f "\${flat_pdb}" ]]; then
            echo "Missing ${model_id} inputs in archive_root (${ARCHIVE_ROOT})." >&2
            exit 1
        fi
        if gzip -t "\${flat_meta}" >/dev/null 2>&1; then
            gzip -dc "\${flat_meta}" > "\${meta_path}"
        else
            cp "\${flat_meta}" "\${meta_path}"
        fi
        if gzip -t "\${flat_pdb}" >/dev/null 2>&1; then
            gzip -dc "\${flat_pdb}" > "\${pdb_path}"
        else
            cp "\${flat_pdb}" "\${pdb_path}"
        fi
    fi

    cp "\${meta_path}" "${model_id}-meta_v1.json"
    cp "\${pdb_path}" "${model_id}-model_v1.pdb"
    """
}

// Batch version: prepare all AF assets in a single process using parallel I/O
process BATCH_PREPARE_AF_ASSETS {
    tag "batch-prepare"

    input:
    path mapping_file
    val posix_root
    val archive_root

    output:
    path "prepared/*-meta_v1.json", emit: meta_files
    path "prepared/*-model_v1.pdb", emit: pdb_files
    path "prepared_manifest.tsv", emit: manifest

    script:
    def archiveRootStr = archive_root ?: 'NO_ARCHIVE'
    """
    set -euo pipefail
    mkdir -p prepared

    python3 - <<'PYTHON_SCRIPT'
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

POSIX_ROOT = Path("${posix_root}")
ARCHIVE_ROOT = Path("${archiveRootStr}") if "${archiveRootStr}" != "NO_ARCHIVE" else None
MAPPING_FILE = Path("${mapping_file}")
OUTPUT_DIR = Path("prepared")

def shard_path(af_id: str) -> Path:
    digits = af_id.replace("AF-", "")[:16]
    chunks = [digits[i:i+4] for i in range(0, 16, 4)]
    return POSIX_ROOT.joinpath(*chunks)

def prepare_model(af_id: str) -> tuple:
    shard_dir = shard_path(af_id)
    meta_on_disk = shard_dir / f"{af_id}-meta_v1.json"
    pdb_on_disk = shard_dir / f"{af_id}-model_v1.pdb"

    # Try shard location first
    if meta_on_disk.exists() and pdb_on_disk.exists():
        meta_src, pdb_src = meta_on_disk, pdb_on_disk
    elif ARCHIVE_ROOT:
        # Fall back to archive root
        flat_meta = ARCHIVE_ROOT / f"{af_id}-meta_v1.json"
        flat_pdb = ARCHIVE_ROOT / f"{af_id}-model_v1.pdb"
        flat_meta_gz = ARCHIVE_ROOT / f"{af_id}-meta_v1.json.gz"
        flat_pdb_gz = ARCHIVE_ROOT / f"{af_id}-model_v1.pdb.gz"

        # Ensure shard dir exists
        shard_dir.mkdir(parents=True, exist_ok=True)

        # Handle meta file (possibly gzipped)
        if flat_meta_gz.exists():
            with gzip.open(flat_meta_gz, 'rb') as f_in:
                meta_on_disk.write_bytes(f_in.read())
        elif flat_meta.exists():
            shutil.copy(flat_meta, meta_on_disk)
        else:
            return af_id, shard_dir, "error", f"Missing meta: {flat_meta}"

        # Handle PDB file (possibly gzipped)
        if flat_pdb_gz.exists():
            with gzip.open(flat_pdb_gz, 'rb') as f_in:
                pdb_on_disk.write_bytes(f_in.read())
        elif flat_pdb.exists():
            shutil.copy(flat_pdb, pdb_on_disk)
        else:
            return af_id, shard_dir, "error", f"Missing pdb: {flat_pdb}"

        meta_src, pdb_src = meta_on_disk, pdb_on_disk
    else:
        return af_id, shard_dir, "error", "Files not found and no archive_root"

    # Copy to output dir
    out_meta = OUTPUT_DIR / f"{af_id}-meta_v1.json"
    out_pdb = OUTPUT_DIR / f"{af_id}-model_v1.pdb"
    shutil.copy(meta_src, out_meta)
    shutil.copy(pdb_src, out_pdb)

    return af_id, shard_dir, "ok", ""

# Read model IDs from mapping file
model_ids = []
with open(MAPPING_FILE) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.strip().split("\\t")
        af_id = parts[1] if len(parts) > 1 and parts[1] else parts[0]
        model_ids.append(af_id.strip())

# Process in parallel
results = []
with ThreadPoolExecutor(max_workers=${WORKERS}) as executor:
    results = list(executor.map(prepare_model, model_ids))

# Write manifest
with open("prepared_manifest.tsv", "w") as f:
    f.write("model_id\\tshard_dir\\tstatus\\tmessage\\n")
    for af_id, shard_dir, status, msg in results:
        f.write(f"{af_id}\\t{shard_dir}\\t{status}\\t{msg}\\n")

ok_count = sum(1 for r in results if r[2] == "ok")
err_count = sum(1 for r in results if r[2] == "error")
print(f"Prepared {ok_count} models, {err_count} errors")
PYTHON_SCRIPT
    """
}

process VALIDATE_ASSETS {
    tag { model_id }

    input:
    tuple val(model_id), path(meta_json), path(pdb_file)

    output:
    path "${model_id}.tsv", emit: validation_tsv

    script:
    """
    set -euo pipefail
    python - <<'PY'
import math
from pathlib import Path

import gemmi
import orjson

meta_path = Path("${meta_json}")
pdb_path = Path("${pdb_file}")

def load_json(path: Path):
    return orjson.loads(path.read_bytes())

def plddt_len(meta: dict):
    for key in ("plddt", "confidence", "plddt_scores"):
        if key in meta and isinstance(meta[key], list):
            return len(meta[key])
    return 0

def pae_dim(meta: dict):
    for key in ("predicted_aligned_error", "pae"):
        val = meta.get(key)
        if isinstance(val, list) and val:
            rows = len(val)
            cols = len(val[0]) if isinstance(val[0], list) else 0
            return rows if rows == cols else f"{rows}x{cols}"
    return "0"

def pdb_residue_count_gemmi(path: Path):
    # Fast PDB residue counting using gemmi
    structure = gemmi.read_structure(str(path))
    if not structure:
        return 0
    model = structure[0]
    residues = set()
    for chain in model:
        for residue in chain:
            if residue.name == "HOH":
                continue
            # Check for NaN coordinates in first atom
            for atom in residue:
                if math.isnan(atom.pos.x) or math.isnan(atom.pos.y) or math.isnan(atom.pos.z):
                    raise ValueError(f"NaN coordinates in {path}")
                break  # Only check first atom per residue
            residues.add((chain.name, residue.seqid.num))
    return len(residues)

def pdb_residue_count_fallback(path: Path):
    # Fallback line-by-line PDB parsing
    residues = set()
    with path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                if math.isnan(x) or math.isnan(y) or math.isnan(z):
                    raise ValueError(f"NaN coordinates in {path}")
                chain_id = line[21].strip()
                try:
                    res_seq = int(line[22:26])
                except ValueError:
                    continue
                residues.add((chain_id, res_seq))
    return len(residues)

def pdb_residue_count(path: Path):
    # Try gemmi first, fall back to line parsing
    try:
        return pdb_residue_count_gemmi(path)
    except Exception:
        return pdb_residue_count_fallback(path)

status = "ok"
reason = ""
plddt_n = 0
pae_d = "0"
pdb_n = 0

try:
    meta = load_json(meta_path)
    plddt_n = plddt_len(meta)
    pae_d = pae_dim(meta)
    pdb_n = pdb_residue_count(pdb_path)

    if isinstance(pae_d, str) and "x" in pae_d:
        status = "mismatch"
        reason = f"PAE not square ({pae_d})"
    elif pae_d == "0":
        status = "mismatch"
        reason = "PAE missing"
    elif plddt_n != pdb_n or (isinstance(pae_d, int) and pae_d != plddt_n):
        status = "mismatch"
        reason = f"pLDDT={plddt_n}, PDB residues={pdb_n}, PAE={pae_d}"
except Exception as exc:
    status = "mismatch"
    reason = f"json_or_parse_error: {exc}"

out = Path("${model_id}.tsv")
out.write_text("\\t".join(map(str, ["${model_id}", status, reason, plddt_n, pdb_n, pae_d])) + "\\n")
PY
    """
}

// Batch validation - processes all models in a single process (100 processes -> 1)
process BATCH_VALIDATE_ASSETS {
    tag "batch-validate"

    input:
    path model_ids_file
    path input_dir

    output:
    path "validation_results.tsv", emit: validation_tsv

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/batch_validate_assets.py \\
      --ids ${model_ids_file} \\
      --input-dir ${input_dir} \\
      --output validation_results.tsv \\
      --workers ${WORKERS}
    """
}

process CONVERT_COLABFOLD {
    tag { model_id }
    errorStrategy 'ignore'
    maxRetries 0

    // Publish per-model outputs into the derived POSIX layout and aggregate manifest dirs
    publishDir { derivedScoresDir(model_id) }, mode: 'copy', overwrite: true, pattern: '*_v1.json'
    publishDir { derivedManifestsDir(model_id) }, mode: 'copy', overwrite: true, pattern: 'chains/*_afid_mapping.csv'
    publishDir { derivedManifestsDir(model_id) }, mode: 'copy', overwrite: true, pattern: 'models/*_model_metadata.csv'
    publishDir PER_CHAIN_MANIFEST_DIR, mode: 'copy', overwrite: true, pattern: 'chains/*_afid_mapping.csv'
    publishDir PER_MODEL_MANIFEST_DIR, mode: 'copy', overwrite: true, pattern: 'models/*_model_metadata.csv'

    input:
    tuple val(model_id), path(meta_json), path(pdb_file)

    output:
    path "chains/${model_id}_afid_mapping.csv", emit: chain_manifest
    path "models/${model_id}_model_metadata.csv", emit: model_manifest
    path "${model_id}-confidence_v1.json", emit: plddt_json
    path "${model_id}-predicted_aligned_error_v1.json", emit: pae_json

    script:
    def duckdbArg = CONVERTER_DUCKDB ? "--duckdb ${file(CONVERTER_DUCKDB)}" : ""
    """
    set -euo pipefail
    mkdir -p chains models
    afdb-colabfold-convert \\
      ${meta_json} \\
      ${pdb_file} \\
      --manifest ${file(COLABFOLD_MANIFEST)} \\
      --model-entity-id ${model_id} \\
      ${duckdbArg} \\
      --chain-manifest-dir chains \\
      --model-manifest-dir models \\
      --outdir .
    """
}

// Write prepared model IDs to file for batch validation
process WRITE_PREPARED_IDS {
    tag "write-prepared-ids"

    input:
    val model_ids

    output:
    path "prepared_ids.txt", emit: ids_file

    script:
    def idsText = model_ids.join('\n')
    """
    cat > prepared_ids.txt <<'EOF'
${idsText}
EOF
    """
}

// Write validated model IDs to file for batch convert
process WRITE_VALIDATED_IDS {
    tag "write-validated-ids"

    input:
    val model_ids

    output:
    path "validated_ids.txt", emit: ids_file

    script:
    def idsText = model_ids.join('\n')
    """
    cat > validated_ids.txt <<'EOF'
${idsText}
EOF
    """
}

// Batch version of CONVERT_COLABFOLD - processes all models in a single process
process BATCH_CONVERT_COLABFOLD {
    tag "batch-convert-colabfold"
    errorStrategy 'ignore'

    publishDir PER_CHAIN_MANIFEST_DIR, mode: 'copy', overwrite: true, pattern: 'chains/*.csv'
    publishDir PER_MODEL_MANIFEST_DIR, mode: 'copy', overwrite: true, pattern: 'models/*.csv'

    input:
    tuple path(model_ids_file), path(input_dir), path(colabfold_manifest), path(db_file)

    output:
    path "chains/*_afid_mapping.csv", emit: chain_manifests
    path "models/*_model_metadata.csv", emit: model_manifests
    path "scores/*-confidence_v1.json", emit: plddt_jsons
    path "scores/*-predicted_aligned_error_v1.json", emit: pae_jsons

    script:
    def duckdbArg = db_file.name != 'NO_DB' ? "--duckdb ${db_file}" : ""
    """
    set -euo pipefail
    mkdir -p chains models scores
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/batch_convert_colabfold.py \\
      --model-ids-file ${model_ids_file} \\
      --input-dir ${input_dir} \\
      --manifest ${colabfold_manifest} \\
      ${duckdbArg} \\
      --output-dir scores \\
      --chain-manifest-dir chains \\
      --model-manifest-dir models \\
      --workers ${WORKERS}
    """
}

process MERGE_MANIFESTS {
    tag "merge-manifests"

    publishDir MERGED_MANIFEST_DIR, mode: 'copy', overwrite: true

    input:
    val chain_files
    val model_files

    output:
    tuple path("uniprot_afid_mapping.csv"), path("uniprot_model_metadata.csv"), emit: merged_manifests

    script:
    """
    set -euo pipefail
    # Use awk for efficient CSV merging (skip header except first file)
    awk 'FNR==1 && NR!=1 {next} {print}' ${chain_files.join(' ')} > uniprot_afid_mapping.csv
    awk 'FNR==1 && NR!=1 {next} {print}' ${model_files.join(' ')} > uniprot_model_metadata.csv
    """
}

// Legacy per-model export processes (kept for reference, replaced by batch versions)
// process EXPORT_MODEL_METADATA { ... }
// process EXPORT_CHAIN_METADATA { ... }

process WRITE_MODEL_IDS {
    tag "write-model-ids"

    input:
    val model_ids_list

    output:
    path "model_ids.txt", emit: model_ids_file

    script:
    def idsLines = model_ids_list.collect { it.toString().trim() }.findAll { it }.join('\n')
    """
    cat > model_ids.txt <<'EOF'
${idsLines}
EOF
    """
}

process BATCH_EXPORT_MODEL_METADATA {
    tag "batch-model-export"

    publishDir PER_MODEL_JSON_DIR, mode: 'copy', overwrite: true

    input:
    tuple path(model_ids_file), path(db_file), path(chain_manifest), path(model_manifest)

    output:
    path "models/*.json", emit: model_jsons

    script:
    """
    set -euo pipefail
    mkdir -p models
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/batch_export_metadata.py \\
      --model-ids ${model_ids_file} \\
      --db ${db_file} \\
      --config ${file(DATASET_CONFIG)} \\
      --mapping ${chain_manifest} \\
      --model-manifest ${model_manifest} \\
      --output-dir models \\
      --export-type model \\
      --workers ${WORKERS}
    """
}

process BATCH_EXPORT_CHAIN_METADATA {
    tag "batch-chain-export"

    publishDir PER_CHAIN_JSON_DIR, mode: 'copy', overwrite: true

    input:
    tuple path(model_ids_file), path(db_file), path(chain_manifest), path(model_manifest)

    output:
    path "chains/*.json", emit: chain_jsons

    script:
    """
    set -euo pipefail
    mkdir -p chains
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/batch_export_metadata.py \\
      --model-ids ${model_ids_file} \\
      --db ${db_file} \\
      --config ${file(DATASET_CONFIG)} \\
      --mapping ${chain_manifest} \\
      --model-manifest ${model_manifest} \\
      --output-dir chains \\
      --export-type chain \\
      --workers ${WORKERS}
    """
}

process COMBINE_MODEL_METADATA {
    tag "combine-model-json"

    publishDir BATCH_MODEL_OUTDIR, mode: 'copy', overwrite: true

    input:
    path json_files

    output:
    path "${MODEL_BATCH_PREFIX}-*.json"

    script:
    """
    set -euo pipefail
    mkdir -p input_jsons
    cp ${json_files} input_jsons/
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/combine_metadata.py \\
      --input-dir input_jsons \\
      --output-dir ./ \\
      --output-prefix ${MODEL_BATCH_PREFIX} \\
      --chunk-size ${BATCH_SIZE}
    """
}

process COMBINE_CHAIN_METADATA {
    tag "combine-chain-json"

    publishDir BATCH_CHAIN_OUTDIR, mode: 'copy', overwrite: true

    input:
    path json_files

    output:
    path "${CHAIN_BATCH_PREFIX}-*.json"

    script:
    """
    set -euo pipefail
    mkdir -p input_jsons
    cp ${json_files} input_jsons/
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/combine_metadata.py \\
      --input-dir input_jsons \\
      --output-dir ./ \\
      --output-prefix ${CHAIN_BATCH_PREFIX} \\
      --chunk-size ${BATCH_SIZE}
    """
}

// Legacy per-model export process (kept for reference, replaced by batch version)
// process EXPORT_MODELCIF_INPUT { ... }

process BATCH_EXPORT_MODELCIF_INPUT {
    tag "batch-modelcif-export"

    publishDir MODELCIF_INPUT_DIR, mode: 'copy', overwrite: true

    input:
    tuple path(model_ids_file), path(db_file), path(merged_chain_manifest)

    output:
    path "modelcif/*.json", emit: modelcif_jsons

    script:
    """
    set -euo pipefail
    mkdir -p modelcif
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/batch_export_modelcif_input.py \\
      --model-ids ${model_ids_file} \\
      --manifest ${merged_chain_manifest} \\
      --db ${db_file} \\
      --template ${file(MODELCIF_TEMPLATE)} \\
      --output-dir modelcif \\
      --workers ${WORKERS}
    """
}

// Legacy per-model processes (kept for reference)
process RUN_MODELCIF_GEN {
    tag { model_id }

    publishDir { derivedStructuresDir(model_id) }, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(pdb_file), path(modelcif_json)

    output:
    tuple val(model_id), path("${model_id}-model-${MODEL_VERSION}.cif"), emit: mmcif_file

    script:
    """
    set -euo pipefail
    ${TOOLKIT_CMD} run-modelcif-gen \\
      -p ${pdb_file} \\
      -m ${modelcif_json} \\
      -o ${model_id}-model-${MODEL_VERSION}.cif
    """
}

process RUN_DSSP {
    tag { model_id }

    publishDir { derivedStructuresDir(model_id) }, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(mmcif_file)

    output:
    tuple val(model_id), path("${model_id}-model-${MODEL_VERSION}.cif"), emit: dssp_mmcif

    script:
    """
    set -euo pipefail
    ${TOOLKIT_CMD} run-dssp \\
      -i ${mmcif_file} \\
      -o ${model_id}-model-${MODEL_VERSION}.cif
    """
}

process RUN_MODELPDB {
    tag { model_id }

    publishDir { derivedStructuresDir(model_id) }, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(mmcif_file), path(pdb_file)

    output:
    path "${model_id}-model-${MODEL_VERSION}.pdb", emit: enriched_pdb

    script:
    """
    set -euo pipefail
    ${TOOLKIT_CMD} run-modelpdb-gen \\
      -c ${mmcif_file} \\
      -p ${pdb_file} \\
      -r ${file(PROVIDER_JSON)} \\
      -o ${model_id}-model-${MODEL_VERSION}.pdb
    """
}

process RUN_CIF2BCIF {
    tag { model_id }

    publishDir { derivedStructuresDir(model_id) }, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(mmcif_file)

    output:
    path "${model_id}-model-${MODEL_VERSION}.bcif", emit: bcif_file

    script:
    """
    set -euo pipefail
    export PATH="/home/mitsenkov/.npm-global/bin:${PATH}"
    ${TOOLKIT_CMD} run-cif2bcif \\
      -i ${mmcif_file} \\
      -o ${model_id}-model-${MODEL_VERSION}.bcif
    """
}

process BATCH_RUN_MODELCIF_GEN {
    tag "batch-modelcif-gen"

    input:
    path pdb_files
    path metadata_files

    output:
    path "output/*.cif", emit: mmcif_files

    script:
    """
    set -euo pipefail
    mkdir -p output pdb_dir metadata_dir

    # Use symlinks instead of copies (faster)
    for f in *.pdb; do
        [ -f "\$f" ] && ln -s "\$(pwd)/\$f" pdb_dir/
    done
    for f in *.json; do
        [ -f "\$f" ] && ln -s "\$(pwd)/\$f" metadata_dir/
    done

    ${TOOLKIT_CMD} run-batch-modelcif-gen \\
      --input-dir pdb_dir \\
      --metadata-dir metadata_dir \\
      --output-dir output \\
      --model-version ${MODEL_VERSION} \\
      --skip-validation \\
      --skip-alignment
    """
}

process BATCH_RUN_DSSP {
    tag "batch-dssp"

    input:
    path cif_files

    output:
    path "output/*.cif", emit: dssp_files

    script:
    """
    set -euo pipefail
    mkdir -p output

    # Files are already staged by Nextflow - process directly
    ${TOOLKIT_CMD} batch-dssp --input-dir . --output-dir output --workers ${WORKERS}
    """
}

process BATCH_RUN_MODELPDB {
    tag "batch-modelpdb"

    input:
    path cif_files
    path pdb_files
    path provider_json

    output:
    path "output/*.pdb", emit: enriched_pdbs

    script:
    """
    set -euo pipefail
    mkdir -p output cif_dir pdb_dir

    # Stage CIF files
    for f in *-model-*.cif; do
        [ -f "\$f" ] && cp "\$f" cif_dir/
    done

    # Stage PDB files
    for f in *-model_v1.pdb; do
        [ -f "\$f" ] && cp "\$f" pdb_dir/
    done

    ${TOOLKIT_CMD} run-batch-modelpdb-gen \\
      --cif-dir cif_dir \\
      --pdb-dir pdb_dir \\
      --provider ${provider_json} \\
      --output-dir output \\
      --model-version ${MODEL_VERSION}
    """
}

process BATCH_RUN_CIF2BCIF {
    tag "batch-cif2bcif"

    input:
    path cif_files

    output:
    path "output/*.bcif", emit: bcif_files

    script:
    """
    set -euo pipefail
    mkdir -p input_cifs output

    # Stage CIF files to input directory
    for f in *.cif; do
        [ -f "\$f" ] && cp "\$f" input_cifs/
    done

    # Use gemmi-based batch CIF2BCIF (faster than molstar subprocess)
    ${TOOLKIT_CMD} batch-cif2bcif --input-dir input_cifs --output-dir output --workers ${WORKERS}
    """
}

process CLEANUP_META_JSON {
    tag { model_id }

    input:
    tuple val(model_id), val(meta_path), val(ready)

    output:
    val model_id

    script:
    """
    set -euo pipefail
    if [[ -f "${meta_path}" ]]; then
        rm -f "${meta_path}"
    fi
    """
}

process CLEANUP_PDB_FILE {
    tag { model_id }

    input:
    tuple val(model_id), val(pdb_path), val(ready)

    output:
    val model_id

    script:
    """
    set -euo pipefail
    if [[ -f "${pdb_path}" ]]; then
        rm -f "${pdb_path}"
    fi
    """
}

// Batch cleanup - replaces CLEANUP_META_JSON and CLEANUP_PDB_FILE (200 processes -> 1)
process BATCH_CLEANUP_FILES {
    tag "batch-cleanup"

    input:
    val meta_paths
    val pdb_paths
    val ready

    output:
    val true, emit: done

    script:
    def metaList = meta_paths.join(' ')
    def pdbList = pdb_paths.join(' ')
    """
    set -euo pipefail
    for f in ${metaList}; do
        [ -f "\$f" ] && rm -f "\$f" || true
    done
    for f in ${pdbList}; do
        [ -f "\$f" ] && rm -f "\$f" || true
    done
    """
}

process CLEANUP_PER_AF_MANIFESTS {
    tag "cleanup-per-af-manifests"

    input:
    tuple val(model_ids), val(ready)

    output:
    path "cleanup_manifests.done"

    script:
    def idsLines = model_ids.collect { it.toString().trim() }.findAll { it }.join('\n')
    """
    set -euo pipefail
    export POSIX_ROOT_ENV="${POSIX_ROOT}"
    export PER_MODEL_MANIFEST_DIR_ENV="${PER_MODEL_MANIFEST_DIR}"
    export PER_CHAIN_MANIFEST_DIR_ENV="${PER_CHAIN_MANIFEST_DIR}"
    cat > ids.txt <<'EOF'
${idsLines}
EOF
    python - <<'PY'
import os
import pathlib

ids = (pathlib.Path("ids.txt").read_text().splitlines() if pathlib.Path("ids.txt").exists() else [])
posix_root = pathlib.Path(os.environ["POSIX_ROOT_ENV"])
per_model_dir = pathlib.Path(os.environ["PER_MODEL_MANIFEST_DIR_ENV"])
per_chain_dir = pathlib.Path(os.environ["PER_CHAIN_MANIFEST_DIR_ENV"])

for acc in ids:
    digits = acc.replace("AF-","")[:16]
    shard = posix_root.joinpath(*[digits[i:i+4] for i in range(0,16,4)])
    derived_manifest_root = shard / "derived" / "manifests"
    targets = [
        derived_manifest_root / "models" / f"{acc}_model_metadata.csv",
        derived_manifest_root / "chains" / f"{acc}_afid_mapping.csv",
        per_model_dir / f"{acc}_model_metadata.csv",
        per_chain_dir / f"{acc}_afid_mapping.csv",
    ]
    for target in targets:
        try:
            target.unlink()
        except FileNotFoundError:
            pass
PY
    echo "ok" > cleanup_manifests.done
    """
}

process CLEANUP_METADATA_JSONS {
    tag "cleanup-metadata-jsons"

    input:
    tuple val(model_ids), val(ready)

    output:
    path "cleanup_metadata.done"

    script:
    def idsLines = model_ids.collect { it.toString().trim() }.findAll { it }.join('\n')
    """
    set -euo pipefail
    export POSIX_ROOT_ENV="${POSIX_ROOT}"
    cat > ids.txt <<'EOF'
${idsLines}
EOF
    python - <<'PY'
import pathlib, sys
import os
from pathlib import Path
posix_root = pathlib.Path(os.environ["POSIX_ROOT_ENV"])
ids = (Path("ids.txt").read_text().splitlines() if Path("ids.txt").exists() else [])
for acc in ids:
    digits = acc.replace("AF-","")[:16]
    shard = posix_root.joinpath(*[digits[i:i+4] for i in range(0,16,4)])
    target = shard / "derived" / "metadata" / f"{acc}.json"
    try:
        target.unlink()
    except FileNotFoundError:
        pass
PY
    echo "ok" > cleanup_metadata.done
    """
}

workflow {
    main:
        if( !UNIPROT_DB )
            throw new IllegalArgumentException("Parameter 'uniprot_db' is required and must point to a prebuilt DuckDB file.")

        // Batch prepare: single process prepares all assets
        def mapping_file = file(MAPPING_FILE)
        def posix_root_str = POSIX_ROOT
        def posix_root_path = file(POSIX_ROOT)  // For processes that need path input
        def archive_root_str = ARCHIVE_ROOT ?: ''

        def batch_prepared = BATCH_PREPARE_AF_ASSETS(mapping_file, posix_root_str, archive_root_str)

        // Build prepared_assets channel by joining manifest with output files
        def meta_files_by_id = batch_prepared.meta_files.flatten().map { f ->
            def name = f.name
            def model_id = name.replace('-meta_v1.json', '')
            tuple(model_id, f)
        }
        def pdb_files_by_id = batch_prepared.pdb_files.flatten().map { f ->
            def name = f.name
            def model_id = name.replace('-model_v1.pdb', '')
            tuple(model_id, f)
        }

        // Join manifest info (for shard_dir) with actual output files
        def manifest_rows = batch_prepared.manifest
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.status == 'ok' }
            .map { row -> tuple(row.model_id, file(row.shard_dir)) }

        def prepared_assets = manifest_rows
            .join(meta_files_by_id)
            .join(pdb_files_by_id)
            .map { model_id, shard_dir, meta_json, pdb_file ->
                tuple(model_id, shard_dir, meta_json, pdb_file)
            }

        // Batch validation: collect all prepared IDs, validate in one process
        def prepared_ids_list = prepared_assets.map { t -> t[0] }.collect()
        def prepared_ids_file = WRITE_PREPARED_IDS(prepared_ids_list).ids_file
        def batch_validation = BATCH_VALIDATE_ASSETS(prepared_ids_file, posix_root_path)

    def validation_rows = batch_validation.validation_tsv.splitCsv(header: false, sep: '\t')

    def mismatches = validation_rows
        .filter { row -> row[1] != 'ok' }
        .map { row -> "${row[0]}\t${row[2]}\t${row[3]}\t${row[4]}\t${row[5]}" }
    def mismatches_with_header = Channel.of("af_id\treason\tplddt_len\tpdb_residues\tpae_dim").concat(mismatches)
    mismatches_with_header.collectFile(
        name: 'mismatches.tsv',
        storeDir: VALIDATION_OUTPUT_DIR,
        cache: false,
        newLine: true,
        sort: false
    )

    def ok_ids = validation_rows
        .filter { row -> row[1] == 'ok' }
        .map { row -> row[0] }
    def ok_with_header = Channel.of("af_id").concat(ok_ids)
    ok_with_header.collectFile(
        name: 'ok_ids.tsv',
        storeDir: VALIDATION_OUTPUT_DIR,
        cache: false,
        newLine: true,
        sort: false
    )

    def ok_flags = ok_ids.map { id -> tuple(id, true) }
    def prepared_by_id = prepared_assets.map { t -> tuple(t[0], t) }
    def ok_prepared = prepared_by_id.join(ok_flags).map { id, payload, flag -> payload }

    def assets_for_structures = ok_prepared.map { t -> tuple(t[0], t[3]) }
    def assets_for_meta_cleanup = ok_prepared
    def assets_for_pdb_cleanup = ok_prepared

    // Collect validated IDs for batch processing
    def validated_ids_list = ok_prepared.map { t -> t[0] }.collect()
    def validated_ids_file = WRITE_VALIDATED_IDS(validated_ids_list).ids_file

    // Batch convert ColabFold outputs (1 process instead of N)
    def posix_root_file = file(POSIX_ROOT)
    def colabfold_manifest_file = file(COLABFOLD_MANIFEST)
    def converter_db_file = CONVERTER_DUCKDB ? file(CONVERTER_DUCKDB) : file("NO_DB")
    def batch_convert_input = validated_ids_file.map { ids_file ->
        tuple(ids_file, posix_root_file, colabfold_manifest_file, converter_db_file)
    }
    def converted = BATCH_CONVERT_COLABFOLD(batch_convert_input)

    def chain_files = converted.chain_manifests.flatten().collect()
    def model_files = converted.model_manifests.flatten().collect()
    def merged = MERGE_MANIFESTS(chain_files, model_files).merged_manifests
    def merge_done = merged.map { true }.last()

    def model_ids = merged
        .map { chain_manifest, model_manifest -> chain_manifest }
        .splitCsv(header: true)
        .map { row -> (row.model_entity_id ?: row.model_id)?.toString()?.trim() }
        .filter { it }
        .unique()
    def model_ids_list = model_ids.collect()

    def chain_manifest_path = file("${MERGED_MANIFEST_DIR}/uniprot_afid_mapping.csv")
    def model_manifest_path = file("${MERGED_MANIFEST_DIR}/uniprot_model_metadata.csv")
    def db_file = file(UNIPROT_DB)

    // Write model IDs to file for batch processing
    def model_ids_file = WRITE_MODEL_IDS(model_ids_list).model_ids_file

    // Batch export metadata (processes all models in single invocations)
    def batch_model_input = model_ids_file.map { ids_file -> tuple(ids_file, db_file, chain_manifest_path, model_manifest_path) }
    def batch_chain_input = model_ids_file.map { ids_file -> tuple(ids_file, db_file, chain_manifest_path, model_manifest_path) }

    def model_jsons = BATCH_EXPORT_MODEL_METADATA(batch_model_input).model_jsons.flatten()
    def chain_jsons = BATCH_EXPORT_CHAIN_METADATA(batch_chain_input).chain_jsons.flatten()

    // Collect all JSON files and pass to combine processes
    def model_jsons_collected = model_jsons.collect()
    def chain_jsons_collected = chain_jsons.collect()

    def model_batches = COMBINE_MODEL_METADATA(model_jsons_collected)
    def chain_batches = COMBINE_CHAIN_METADATA(chain_jsons_collected)

    // Batch export ModelCIF inputs
    def modelcif_batch_input = model_ids_file.map { ids_file -> tuple(ids_file, db_file, chain_manifest_path) }
    def modelcif_jsons_flat = BATCH_EXPORT_MODELCIF_INPUT(modelcif_batch_input).modelcif_jsons.flatten()

    // Collect PDB files and metadata JSONs for batch processing
    def pdb_channel = assets_for_structures.map { t -> tuple(t[0], t[1]) }
    def pdb_files_collected = pdb_channel.map { t -> t[1] }.collect()
    def metadata_files_collected = modelcif_jsons_flat.collect()

    // Batch ModelCIF generation (1 process instead of N)
    def mmcif_files_flat = BATCH_RUN_MODELCIF_GEN(pdb_files_collected, metadata_files_collected).mmcif_files.flatten()

    // Batch DSSP (1 process instead of N)
    def mmcif_files_collected = mmcif_files_flat.collect()
    def dssp_files_flat = BATCH_RUN_DSSP(mmcif_files_collected).dssp_files.flatten()

    // Batch ModelPDB enrichment (1 process instead of N)
    def dssp_files_collected = dssp_files_flat.collect()
    def provider_json_file = file(PROVIDER_JSON)
    def modelpdb_files_flat = BATCH_RUN_MODELPDB(dssp_files_collected, pdb_files_collected, provider_json_file).enriched_pdbs.flatten()

    // Batch CIF2BCIF conversion (1 process instead of N)
    def bcif_files = BATCH_RUN_CIF2BCIF(dssp_files_collected).bcif_files.flatten()

        // Collect file paths for batch cleanup
        def meta_paths = assets_for_meta_cleanup
            .map { t ->
                def id = t[0]
                shardPathForAf(id).resolve("${id}-meta_v1.json").toString()
            }
            .collect()

        def pdb_paths = assets_for_pdb_cleanup
            .map { t ->
                def id = t[0]
                shardPathForAf(id).resolve("${id}-model_v1.pdb").toString()
            }
            .collect()

        // Wait for modelpdb to finish before cleanup
        def modelpdb_done = modelpdb_files_flat.map { true }.last()

        // Batch cleanup (1 process instead of 200)
        BATCH_CLEANUP_FILES(meta_paths, pdb_paths, modelpdb_done)

        def manifest_cleanup_ready = merge_done.combine(model_ids_list).map { vals ->
            def ready = vals[0]
            def ids = vals.size() > 1 ? vals[1..-1] : []
            tuple(ids, ready)
        }
        CLEANUP_PER_AF_MANIFESTS(manifest_cleanup_ready)

        def model_batches_done = model_batches.map { true }.last()
        def chain_batches_done = chain_batches.map { true }.last()
        def cleanup_ready = model_batches_done.combine(chain_batches_done).map { a,b -> true }
        // combine() flattens value channels, so pull the ready flag from index 0 and the IDs from the tail
        def cleanup_inputs = cleanup_ready.combine(model_ids_list).map { vals ->
            def ready = vals[0]
            def ids = vals.size() > 1 ? vals[1..-1] : []
            tuple(ids, ready)
        }
        CLEANUP_METADATA_JSONS(cleanup_inputs)
}
