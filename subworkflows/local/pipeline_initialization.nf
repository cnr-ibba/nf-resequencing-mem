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

workflow NORMALIZATION_INITIALIZATION {
    take:
    input_vcf           // file: input vcf file
    input_tbi           // file: input tbi file
    genome_fasta        // file: genome fasta file

    main:
    // Check input path parameters to see if they exist
    def require_parameter = [
        '--input_vcf': input_vcf,
        '--input_tbi': input_tbi,
        '--genome_fasta': genome_fasta
    ]
    require_parameter.each { key, value ->
        if (! value) {
            error "Required parameter '${key}' is missing"
        }
    }

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
