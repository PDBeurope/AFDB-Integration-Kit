nextflow.enable.dsl=2

def defaultParams = [
    db                 : 'uniprot/outputs/db/uniprot_2025_03.duckdb',
    config             : 'uniprot/config/dataset_config.json',
    mapping            : 'uniprot/config/uniprot_afid_mapping.csv',
    model_manifest     : 'uniprot/config/uniprot_model_metadata.csv',
    chunk_size         : 10_000,
    output_prefix      : 'AF-metadata',
    chain_output_prefix: 'AF-chain-metadata',
    python_cmd         : 'python3'
]

defaultParams.each { key, value ->
    if( !params.containsKey(key) || params.get(key) == null || params.get(key).toString().trim() == '' ) {
        params[key] = value
    }
}

def repoDir = projectDir.parent ?: projectDir
def mappingPath = file(params.mapping)
def configDir = mappingPath?.parent
def datasetDir = configDir?.parent ?: repoDir

if( !params.containsKey('per_model_outdir') || params.per_model_outdir == null || params.per_model_outdir.toString().trim() == '' ) {
    params.per_model_outdir = datasetDir.resolve('per_accession/models').toString()
}

if( !params.containsKey('per_chain_outdir') || params.per_chain_outdir == null || params.per_chain_outdir.toString().trim() == '' ) {
    params.per_chain_outdir = datasetDir.resolve('per_accession/chains').toString()
}

if( !params.containsKey('batch_model_outdir') || params.batch_model_outdir == null || params.batch_model_outdir.toString().trim() == '' ) {
    params.batch_model_outdir = datasetDir.resolve('batches/models').toString()
}

if( !params.containsKey('batch_chain_outdir') || params.batch_chain_outdir == null || params.batch_chain_outdir.toString().trim() == '' ) {
    params.batch_chain_outdir = datasetDir.resolve('batches/chains').toString()
}

def buildModelChannel() {
    Channel
        .fromPath(params.mapping)
        .splitCsv(header: true)
        .map { row ->
            def normalized = row.collectEntries { k, v ->
                [(k?.toString()?.toLowerCase()): v]
            }
            def modelId = normalized['model_entity_id'] ?: normalized['model_id']
            modelId?.toString()?.trim()
        }
        .filter { it }
        .unique()
}

process EXPORT_MODEL_METADATA {
    tag { model_id }

    publishDir params.per_model_outdir, mode: 'copy', overwrite: true

    input:
    val model_id

    output:
    path "${model_id}.json"

    script:
    """
    ${params.python_cmd} ${repoDir}/uniprot/scripts/export_model_metadata.py \
      --model-entity-id ${model_id} \
      --db ${file(params.db)} \
      --config ${file(params.config)} \
      --mapping ${file(params.mapping)} \
      --model-manifest ${file(params.model_manifest)} \
      --out ${model_id}.json
    """
}

process EXPORT_CHAIN_METADATA {
    tag { model_id }

    publishDir params.per_chain_outdir, mode: 'copy', overwrite: true

    input:
    val model_id

    output:
    path "${model_id}.json"

    script:
    """
    ${params.python_cmd} ${repoDir}/uniprot/scripts/export_chain_metadata.py \
      --model-entity-id ${model_id} \
      --db ${file(params.db)} \
      --config ${file(params.config)} \
      --mapping ${file(params.mapping)} \
      --model-manifest ${file(params.model_manifest)} \
      --out ${model_id}.json
    """
}

process COMBINE_MODEL_METADATA {
    publishDir params.batch_model_outdir, mode: 'copy', overwrite: true

    input:
    val metadata_files

    output:
    path "${params.output_prefix}-*.json"

    script:
    """
    ${params.python_cmd} ${repoDir}/uniprot/scripts/combine_metadata.py \
      --input-dir ${file(params.per_model_outdir)} \
      --output-dir ./ \
      --output-prefix ${params.output_prefix} \
      --chunk-size ${params.chunk_size}
    """
}

process COMBINE_CHAIN_METADATA {
    publishDir params.batch_chain_outdir, mode: 'copy', overwrite: true

    input:
    val metadata_files

    output:
    path "${params.chain_output_prefix}-*.json"

    script:
    """
    ${params.python_cmd} ${repoDir}/uniprot/scripts/combine_metadata.py \
      --input-dir ${file(params.per_chain_outdir)} \
      --output-dir ./ \
      --output-prefix ${params.chain_output_prefix} \
      --chunk-size ${params.chunk_size}
    """
}

workflow {
    main:
        def model_for_models = buildModelChannel()
        def model_for_chains = buildModelChannel()

        def model_exports = EXPORT_MODEL_METADATA(model_for_models)
        def chain_exports = EXPORT_CHAIN_METADATA(model_for_chains)

        def combine_model_ready = model_exports.collect()
        def combine_chain_ready = chain_exports.collect()

        COMBINE_MODEL_METADATA(combine_model_ready)
        COMBINE_CHAIN_METADATA(combine_chain_ready)
}
