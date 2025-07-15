
def getEntryDir(entry) {
    // Function to get the directory for a given entry
    // AF-1000000000000001 -> 1000/0000/0000/0001
    return entry.replaceFirst(/^AF-/, '').replaceAll(/(\d{4})(?=\d)/, '$1/')
}

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

params.input_dir = "/input"
params.output_dir = "/output"
params.input_list = "${params.input_dir}/input.txt"
params.version = "v1"
params.results_dir = "${params.output_dir}/results"
params.python_cmd = "uv run /app/main.py"

workflow {

    input_channel = Channel.fromPath(params.input_list)
        .splitCsv()
        .map { row -> row[0] } // Assuming the first column contains the entry IDs

    createEntryOutDir(input_channel)

    runModelCifGenerator(input_channel)

    runDSSP(runModelCifGenerator.out)

    runCif2Bcif(runDSSP.out)
}
