

process runCif2Bcif {

    publishDir "${params.results_dir}", mode: 'copy'

    input:
        path cif_file
        val entry

    output:
        path "${entry}-model_${params.version}.bcif"

    script:
    """
    echo "Converting ${entry} CIF to BCIF"
    ${params.python_cmd} run-cif2bcif -i "${cif_file}" -o "${entry}-model_${params.version}.bcif"
    """
}

process runDSSP {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
        val entry

    output:
        path "${entry}-model_${params.version}.cif"
        val entry, emit: entry

    script:
    """
    echo "Running DSSP on ${entry} CIF"
    ${params.python_cmd} run-dssp -i "${params.input_dir}/${entry}-model_${params.version}.cif" -o "${entry}-model_${params.version}.cif"
    """
}

params.input_dir = "/input"
params.output_dir = "/output"
params.input_list = "${params.input_dir}/input.txt"
params.version = "v4"
params.results_dir = "${params.output_dir}/results"
params.python_cmd = "uv run /app/main.py"

workflow {

    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs

    runDSSP(input_channel)

    runCif2Bcif(runDSSP.out)
}
