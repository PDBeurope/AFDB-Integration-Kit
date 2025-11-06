nextflow.enable.dsl=2

params.db            = params.db            ?: 'uniprot/outputs/db/uniprot_2025_03.duckdb'
params.config        = params.config        ?: 'uniprot/config/dataset_config.json'
params.mapping       = params.mapping       ?: 'uniprot/config/uniprot_afid_mapping.csv'
params.per_outdir    = params.per_outdir    ?: 'uniprot/outputs/per_accession'
params.batch_outdir  = params.batch_outdir  ?: 'uniprot/outputs/batches'
params.chunk_size    = params.chunk_size    ?: 10_000
params.output_prefix = params.output_prefix ?: 'AF-metadata'
params.python_cmd    = params.python_cmd    ?: 'python3'

Channel
    .fromPath(params.mapping)
    .map { file -> file.newReader().withCloseable { reader ->
            reader.readLines()
                  .drop(1)
                  .collect { line ->
                      def cells = line.split(',')
                      cells.size() > 1 ? cells[1].trim() : null
                  }
        }
    }
    .flatten()
    .filter { it }
    .unique()
    .set { accession_ch }

process EXPORT_METADATA {
    tag { accession }

    publishDir params.per_outdir, mode: 'copy', overwrite: true

    input:
    val accession from accession_ch

    output:
    path "${accession}.json"

    script:
    """
    ${params.python_cmd} ${launchDir}/uniprot/scripts/export_metadata.py \
      --accession ${accession} \
      --db ${file(params.db)} \
      --config ${file(params.config)} \
      --mapping ${file(params.mapping)} \
      --out ${accession}.json
    """
}

process COMBINE_METADATA {
    publishDir params.batch_outdir, mode: 'copy', overwrite: true

    input:
    val ready from EXPORT_METADATA.out.collect()

    output:
    path "${params.output_prefix}-*.json"

    script:
    """
    ${params.python_cmd} ${launchDir}/uniprot/scripts/combine_metadata.py \
      --input-dir ${file(params.per_outdir)} \
      --output-dir ./ \
      --output-prefix ${params.output_prefix} \
      --chunk-size ${params.chunk_size}
    """
}

workflow {
    main:
        EXPORT_METADATA(accession_ch)
        COMBINE_METADATA()
}
