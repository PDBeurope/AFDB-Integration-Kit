nextflow.enable.dsl=2

params.db             = params.containsKey('db')             ? params.db             : 'uniprot/outputs/db/uniprot_2025_03.duckdb'
params.manifest       = params.containsKey('manifest')       ? params.manifest       : 'examples/complexes/config/uniprot_afid_mapping.csv'
params.template       = params.containsKey('template')       ? params.template       : 'uniprot/templates/modelcif_metadata.json'
params.output_dir     = params.containsKey('output_dir')     ? params.output_dir     : 'examples/complexes/modelcif_metadata'
params.python_cmd     = params.containsKey('python_cmd')     ? params.python_cmd     : 'python3'

Channel
    .fromPath(params.manifest)
    .flatMap { file ->
        def lines = file.readLines("UTF-8")
        if (!lines) return []
        lines.drop(1)
             .collect { line ->
                 def cells = line.split(",")
                 if (cells.size() < 4) {
                     return null
                 }
                 def model = cells[0].trim()
                 def entity = cells[1].trim()
                 def chain = cells[2].trim()
                 def uniprot = cells[3].trim()
                 if (!model || !entity || !chain || !uniprot) {
                     return null
                 }
                 tuple(model, [entity_id: entity, chain_id: chain, uniprot: uniprot])
             }
             .findAll { it != null }
    }
    .groupTuple()
    .map { model_id, rows -> tuple(model_id, rows.collect { it }) }
    .set { model_groups }

process EXPORT_MODEL_METADATA {
    tag { model_id }
    publishDir params.output_dir, mode: 'copy', overwrite: true

    input:
    tuple val(model_id), val(entries)

    output:
    path "${model_id}.json"

    script:
    def manifest_lines = entries.collect { rec ->
        "${model_id},${rec.entity_id},${rec.chain_id},${rec.uniprot}"
    }.join('\n')
    """
    cat <<'EOF' > manifest_${model_id}.csv
model_entity_id,entity_id,chain_id,uniprot_ac
${manifest_lines}
EOF

    ${params.python_cmd} ${launchDir}/uniprot/scripts/export_modelcif_input.py \
      --model-id ${model_id} \
      --manifest manifest_${model_id}.csv \
      --db ${file(params.db)} \
      --template ${file(params.template)} \
      --out ${model_id}.json
    """
}

workflow {
    main:
        EXPORT_MODEL_METADATA(model_groups)
}
