//
// pipeline initialization checks
//


workflow PIPELINE_INITIALIZATION {
    take:
    input               // string: path to samplesheet

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // this should be present
    Channel.fromPath(input, checkIfExists: true)
        .set { ch_input }

    // other input arguments are evaluated through nf-validation plugin and lib/WorkflowMain.groovy

    emit:
    samplesheet = ch_input
}

workflow NORMALIZATION_INITIALIZATION {
    take:
    input_vcf           // file: input vcf file
    input_tbi           // file: input tbi file
    genome_fasta        // file: genome fasta file

    main:
    // parameters are evaluated through nf-validation plugin and lib/WorkflowMain.groovy

    // Set channels for required files
    Channel.fromPath(input_vcf, checkIfExists: true)
        .map{ it -> [[id:"all-samples-normalized"], it] }
        .set { ch_input_vcf }

    Channel.fromPath(input_tbi, checkIfExists: true)
        .map{ it -> [[id:"all-samples-normalized"], it] }
        .set { ch_input_tbi }

    Channel.fromPath(genome_fasta, checkIfExists: true)
        .map{ it -> [[id:it[0].baseName], it] }
        .set { ch_genome_fasta }

    emit:
    vcf_ch = ch_input_vcf
    tbi_ch = ch_input_tbi
    fasta_ch = ch_genome_fasta
}
