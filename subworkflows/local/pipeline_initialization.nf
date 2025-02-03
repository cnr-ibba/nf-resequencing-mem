//
// pipeline initialization checks
//


workflow PIPELINE_INITIALIZATION {
    take:
    input               // string: path to samplesheet
    multiqc_config      // file: multiqc config file
    genome_fasta        // file: genome fasta file
    genome_bwa_index    // file: genome bwa index file

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Check input path parameters to see if they exist
    def optionalFiles = [multiqc_config, genome_fasta, genome_bwa_index]
    optionalFiles.each { f ->
        if (f) {
            Channel.fromPath(f, checkIfExists: true)
        }
    }

    // this should be present
    Channel.fromPath(input, checkIfExists: true)
        .set { ch_input }

    emit:
    samplesheet = ch_input
}
