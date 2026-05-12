#!/usr/bin/env nextflow


// Include utility functions
include { getEntryDir } from './utils.nf'


/*
  Create output directory for each entry
*/
process createEntryOutDir {

    publishDir "${params.output_dir}", mode: 'copy'

    input:
        val entry

    output:
        path "${getEntryDir(entry)}"

    script:
    """
    echo "Creating output directory for entry: ${entry}"
    mkdir -p "${getEntryDir(entry)}"
    """
}


/*
    Generate model CIF for each entry
*/
process runModelCifGenerator {

    publishDir "${params.output_dir}/${getEntryDir(entry)}", mode: 'copy'

    input:
        val entry

    output:
        path "${entry}-model-${params.version}.cif", emit: cif_file
        val entry, emit: entry

    script:
    """
    echo "Generating model CIF for entry: ${entry}"
    ${params.python_cmd} run-modelcif-gen \
        -p "${params.input_dir}/${getEntryDir(entry)}/${entry}-model-${params.version}.pdb" \
        -m "${params.input_dir}/${getEntryDir(entry)}/${entry}-${params.version}.cif.json" \
        -o "${entry}-model-${params.version}.cif"
    """
}

/*
    Convert CIF to BCIF for each entry
*/
process runCif2Bcif {

    publishDir "${params.output_dir}/${getEntryDir(entry)}", mode: 'copy'

    input:
        path cif_file
        val entry

    output:
        path "${entry}-model-${params.version}.bcif"

    script:
    """
    echo "Converting ${entry} CIF to BCIF"
    ${params.python_cmd} run-cif2bcif -i "${cif_file}" -o "${entry}-model-${params.version}.bcif"
    """
}


/*
    Run DSSP on CIF files for each entry
*/
process runDSSP {
    publishDir "${params.output_dir}/${getEntryDir(entry)}", mode: 'copy'

    input:
        path cif_file
        val entry

    output:
        path "${entry}-model-${params.version}.cif", emit: cif_file
        val entry, emit: entry

    script:
    """
    echo "Running DSSP on ${entry} CIF"
    ${params.python_cmd} run-dssp -i "${cif_file}" -o "${entry}-model-${params.version}.cif"
    """
}

/*
    Add AFDB-specific PDB headers from the DSSP-annotated CIF for each entry
*/
process runModelPdbGenerator {

    publishDir "${params.output_dir}/${getEntryDir(entry)}", mode: 'copy'

    input:
        path cif_file
        val entry

    output:
        path "${entry}-model-${params.version}.pdb"

    script:
    """
    echo "Generating enriched model PDB for entry: ${entry}"
    ${params.python_cmd} run-modelpdb-gen \
        -c "${cif_file}" \
        -p "${params.input_dir}/${getEntryDir(entry)}/${entry}-model-${params.version}.pdb" \
        -r "${params.provider_json}" \
        -o "${entry}-model-${params.version}.pdb"
    """
}

params.input_dir = "/input"
params.output_dir = "/output"
params.input_list = "${params.input_dir}/input.txt"
params.provider_json = "${params.input_dir}/provider.json"
params.version = "v1"
params.results_dir = "${params.output_dir}/results"
params.python_cmd = "uv run /app/main.py"


/*
    Main workflow
*/
workflow {

    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs

    createEntryOutDir(input_channel)

    model_cif = runModelCifGenerator(input_channel)

    dssp = runDSSP(model_cif.cif_file, model_cif.entry)

    runModelPdbGenerator(dssp.cif_file, dssp.entry)
    runCif2Bcif(dssp.cif_file, dssp.entry)
}
