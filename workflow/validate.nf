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


/*
    Concatenate model metadata
    This process is designed to concatenate multiple model metadata files into a single JSON array.
    It uses jq to read multiple JSON files and combine them into a single array.
    The input is expected to be a tuple containing a counter and a list of entries.
*/
process concatModelMetadata {

    input:
        tuple val(counter), val(entries)


    script:
    """
    mkdir -p ${params.output_dir}/metadata
    jq -s '.' ${entries} > ${params.output_dir}/metadata/${params.metadata_prefix}-${counter}.json
    echo "Concatenated metadata for chunk ${counter} into ${params.metadata_prefix}-${counter}.json"
    """
}

params.input_dir = "/input"
params.output_dir = "/output"
params.input_list = "${params.input_dir}/input.txt"
params.version = "v1"
params.python_cmd = "uv run /app/main.py"
params.metadata_chunk_size = 100
params.metadata_prefix = "AF-metadata"

/*
    Main workflow
*/
workflow {

    def counter = 1
    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs

    validateModelMetadata(input_channel)

    chunked_entries = input_channel
        .map { entry -> "${params.input_dir}/${getEntryDir(entry)}/${entry}-${params.version}.json" }
        .collate(params.metadata_chunk_size)
        .map { chunk -> tuple(counter++, chunk.join(' ')) }

    concatModelMetadata(chunked_entries)

}
