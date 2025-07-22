#!/usr/bin/env nextflow

// Include utility functions
include { getEntryDir } from './utils.nf'

/*
    Validate model metadata
*/
process validateModelMetadata {
    input:
        val entry

    script:
    """
    ${params.python_cmd} run-schema-validation \
        -i "${params.input_dir}/${getEntryDir(entry)}/${entry}-${params.version}.json" \
        -t model
    """
}

params.input_dir = "/input"
params.input_list = "${params.input_dir}/input.txt"
params.version = "v2"
params.python_cmd = "uv run /app/main.py"


/*
    Main workflow
*/
workflow {

    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs

    validateModelMetadata(input_channel)
}
