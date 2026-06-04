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

process VALIDATE_ASSETS {
    tag { model_id }

    input:
    tuple val(model_id), path(meta_json), path(pdb_file)

    output:
    path "${model_id}.tsv", emit: validation_tsv

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} - <<'PY'
import json
import math
from pathlib import Path

meta_path = Path("${meta_json}")
pdb_path = Path("${pdb_file}")

def load_json(path: Path):
    with path.open() as fh:
        return json.load(fh)

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

def pdb_residue_count(path: Path):
    residues = set()
    with path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                # Fail fast if any coordinate is NaN
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
    cd ${repoDir} && ${PYTHON_CMD} -m afdb_integration_kit.colabfold.converter \\
      \${OLDPWD}/${meta_json} \\
      \${OLDPWD}/${pdb_file} \\
      --manifest ${file(COLABFOLD_MANIFEST)} \\
      --model-entity-id ${model_id} \\
      ${duckdbArg} \\
      --chain-manifest-dir \${OLDPWD}/chains \\
      --model-manifest-dir \${OLDPWD}/models \\
      --outdir \${OLDPWD}
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
    # Avoid "argument list too long" by streaming files one-by-one instead of passing all to tail at once
    {
      first=1
      for f in ${chain_files.join(' ')}; do
        if [[ \$first -eq 1 ]]; then
          head -n1 "\$f"
          first=0
        fi
        tail -n +2 "\$f"
      done
    } > uniprot_afid_mapping.csv

    {
      first=1
      for f in ${model_files.join(' ')}; do
        if [[ \$first -eq 1 ]]; then
          head -n1 "\$f"
          first=0
        fi
        tail -n +2 "\$f"
      done
    } > uniprot_model_metadata.csv
    """
}

process EXPORT_MODEL_METADATA {
    tag { model_id }

    publishDir { derivedMetadataDir(model_id) }, mode: 'copy', overwrite: true
    publishDir PER_MODEL_JSON_DIR, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(db_file), path(merged_chain_manifest), path(merged_model_manifest)

    output:
    path "${model_id}.json"

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/export_model_metadata.py \\
      --model-entity-id ${model_id} \\
      --db ${db_file} \\
      --config ${file(DATASET_CONFIG)} \\
      --mapping ${merged_chain_manifest} \\
      --model-manifest ${merged_model_manifest} \\
      --out ${model_id}.json
    """
}

process EXPORT_CHAIN_METADATA {
    tag { model_id }

    publishDir { derivedMetadataDir(model_id) }, mode: 'copy', overwrite: true
    publishDir PER_CHAIN_JSON_DIR, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(db_file), path(merged_chain_manifest), path(merged_model_manifest)

    output:
    path "${model_id}.json"

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/export_chain_metadata.py \\
      --model-entity-id ${model_id} \\
      --db ${db_file} \\
      --config ${file(DATASET_CONFIG)} \\
      --mapping ${merged_chain_manifest} \\
      --model-manifest ${merged_model_manifest} \\
      --out ${model_id}.json
    """
}

process COMBINE_MODEL_METADATA {
    tag "combine-model-json"

    publishDir BATCH_MODEL_OUTDIR, mode: 'copy', overwrite: true

    input:
    val ready

    output:
    path "${MODEL_BATCH_PREFIX}-*.json"

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/combine_metadata.py \\
      --input-dir ${file(PER_MODEL_JSON_DIR)} \\
      --output-dir ./ \\
      --output-prefix ${MODEL_BATCH_PREFIX} \\
      --chunk-size ${BATCH_SIZE}
    """
}

process COMBINE_CHAIN_METADATA {
    tag "combine-chain-json"

    publishDir BATCH_CHAIN_OUTDIR, mode: 'copy', overwrite: true

    input:
    val ready

    output:
    path "${CHAIN_BATCH_PREFIX}-*.json"

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/combine_metadata.py \\
      --input-dir ${file(PER_CHAIN_JSON_DIR)} \\
      --output-dir ./ \\
      --output-prefix ${CHAIN_BATCH_PREFIX} \\
      --chunk-size ${BATCH_SIZE}
    """
}

process EXPORT_MODELCIF_INPUT {
    tag { model_id }

    publishDir MODELCIF_INPUT_DIR, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), path(db_file), path(merged_chain_manifest)

    output:
    tuple val(model_id), path("${model_id}.json"), emit: modelcif_input

    script:
    """
    set -euo pipefail
    ${PYTHON_CMD} ${repoDir}/uniprot/scripts/export_modelcif_input.py \\
      --model-id ${model_id} \\
      --manifest ${merged_chain_manifest} \\
      --db ${db_file} \\
      --template ${file(MODELCIF_TEMPLATE)} \\
      --out ${model_id}.json
    """
}

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
    ${PYTHON_CMD} - <<'PY'
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
    ${PYTHON_CMD} - <<'PY'
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

        def af_manifest = buildAfInputChannel()
        def prepared_assets = PREPARE_AF_ASSETS(af_manifest).prepared_assets

        def validation_inputs = prepared_assets.map { t -> tuple(t[0], t[2], t[3]) }
        def validation = VALIDATE_ASSETS(validation_inputs)

    def validation_rows = validation.validation_tsv.splitCsv(header: false, sep: '\t')

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

    def assets_for_convert = ok_prepared.map { t -> tuple(t[0], t[2], t[3]) }
    def assets_for_structures = ok_prepared.map { t -> tuple(t[0], t[3]) }
    def assets_for_meta_cleanup = ok_prepared
    def assets_for_pdb_cleanup = ok_prepared

    def converted = CONVERT_COLABFOLD(assets_for_convert)

    def chain_files = converted.chain_manifest.collect()
    def model_files = converted.model_manifest.collect()
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

    def metadata_inputs = model_ids.map { model_id -> tuple(model_id, db_file, chain_manifest_path, model_manifest_path) }

    def ids_for_model_meta = metadata_inputs
    def ids_for_chain_meta = metadata_inputs

    def model_jsons = EXPORT_MODEL_METADATA(ids_for_model_meta)
    def chain_jsons = EXPORT_CHAIN_METADATA(ids_for_chain_meta)

    def model_jsons_done = model_jsons.map { true }.last()
    def chain_jsons_done = chain_jsons.map { true }.last()

    def model_batches = COMBINE_MODEL_METADATA(model_jsons_done)
    def chain_batches = COMBINE_CHAIN_METADATA(chain_jsons_done)

    def modelcif_inputs = model_ids.map { model_id -> tuple(model_id, db_file, chain_manifest_path) }
    def modelcif_inputs_result = EXPORT_MODELCIF_INPUT(modelcif_inputs)

    def pdb_channel = assets_for_structures.map { t -> tuple(t[0], t[1]) }
    def modelcif_ready = pdb_channel.join(modelcif_inputs_result.modelcif_input)

        def mmcif_files = RUN_MODELCIF_GEN(modelcif_ready).mmcif_file
        def dssp_files = RUN_DSSP(mmcif_files).dssp_mmcif

        // PDB header enrichment only needs the original PDB paired with its DSSP mmCIF
        def modelpdb_ready = dssp_files.join(pdb_channel)
        def modelpdb_results = RUN_MODELPDB(modelpdb_ready)
        RUN_CIF2BCIF(dssp_files)

        def convert_done = converted.plddt_json.map { path ->
            def base = path.getBaseName()
            tuple(base.replaceFirst(/-confidence_v1$/, ''), true)
        }
        def meta_cleanup_inputs = assets_for_meta_cleanup
            .map { t ->
                def id = t[0]
                def onDisk = shardPathForAf(id).resolve("${id}-meta_v1.json").toString()
                tuple(id, onDisk)
            }
            .join(convert_done)
        CLEANUP_META_JSON(meta_cleanup_inputs)

        def modelpdb_done = modelpdb_results.enriched_pdb.map { path ->
            def base = path.getBaseName()
            def suffix = "-model-${MODEL_VERSION}"
            def id = base.endsWith(suffix) ? base.substring(0, base.length() - suffix.length()) : base
            tuple(id, true)
        }
        def pdb_cleanup_inputs = assets_for_pdb_cleanup
            .map { t ->
                def id = t[0]
                def onDisk = shardPathForAf(id).resolve("${id}-model_v1.pdb").toString()
                tuple(id, onDisk)
            }
            .join(modelpdb_done)
        CLEANUP_PDB_FILE(pdb_cleanup_inputs)

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
