

process runCif2Bcif {

    publishDir "${params.results_dir}", mode: 'copy'

    input:
        val entry

    output:
        path "${entry}-model_${params.version}.bcif", emit: bcif_file

    script:
    """
    echo "Converting ${entry} CIF to BCIF"
    cif2bcif "${params.input_dir}/${entry}-model_${params.version}.cif" "${entry}-model_${params.version}.bcif"
    """
}

process runDSSP {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
        val entry

    output:
        path "${entry}-model_${params.version}.cif", emit: cif_file

    script:
    """
    echo "Running DSSP on ${bcif_file}"
    dssp -i "${bcif_file}" -o "${bcif_file.baseName}.dssp"
    """
}

params.input_dir = "/input"
params.output_dir = "/output"
params.input_list = "${params.input_dir}/input.txt"
params.version = "v4"
params.results_dir = "${params.output_dir}/results"

workflow {

    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs



    runCif2Bcif(input_channel)
}
