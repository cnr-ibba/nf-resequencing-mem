/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.genome_fasta, params.genome_bwa_index ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                         } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include {
            SEQKIT_RMDUP as SEQKIT_RMDUP_R1;
            SEQKIT_RMDUP as SEQKIT_RMDUP_R2;
                                            } from '../modules/cnr-ibba/seqkit/rmdup/main'
include { TRIMGALORE                        } from '../modules/nf-core/trimgalore/main'
include { BWA_MEM                           } from '../modules/nf-core/bwa/mem/main'
include { CRAM_FREEBAYES_PARALLEL           } from '../subworkflows/local/cram_freebayes_parallel/main'
include { BCFTOOLS_NORM                     } from '../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX                       } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS                    } from '../modules/nf-core/bcftools/stats/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CRAM_MARKDUPLICATES_PICARD        } from '../subworkflows/local/cram_markduplicates_picard/main'
include { SNPEFF_ANNOTATE                   } from '../subworkflows/local/snpeff_annotate'

// A workflow definition which does not declare any name is assumed to be the
// main workflow and it’s implicitly executed. Therefore it’s the entry point
// of the workflow application.
workflow RESEQUENCING_MEM {
  // collect software version
  ch_versions = Channel.empty()
    // collect software version
    ch_versions = Channel.empty()
    // multiqc stuff
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
        def meta_clone = meta.clone()
        tmp = meta_clone.id.split('_')
        if (tmp.size() > 1) {
            meta_clone.id = tmp[0..-2].join('_')
        }
        [ meta_clone, fastq ]
        }
        .groupTuple(by: [0])
        .branch {
        meta, fastq ->
            single  : fastq.size() == 1
            return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    genome_fasta = Channel.fromPath(params.genome_fasta).collect()
    genome_fasta_fai = params.genome_fasta_fai ? Channel.fromPath(params.genome_fasta_fai).collect() : Channel.empty()
    genome_bwa_index = params.genome_bwa_index ? Channel.fromPath(params.genome_bwa_index).collect() : Channel.empty()

    // Build indices if needed
    PREPARE_GENOME(
        genome_fasta,
        genome_fasta_fai,
        genome_bwa_index
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // call FASTQC from module
    FASTQC(INPUT_CHECK.out.reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // remove duplicates (if necessary)
    if (params.remove_fastq_duplicates) {
        // collect multiple files in one
        ch_cat_fastq
        .multiMap { meta, reads ->
            r1: [meta, reads[0]]
            r2: [meta, reads[1]]
        }.set{ ch_seqkit_input }

        SEQKIT_RMDUP_R1(ch_seqkit_input.r1)
        SEQKIT_RMDUP_R2(ch_seqkit_input.r2)
        ch_versions = ch_versions.mix(SEQKIT_RMDUP_R1.out.versions)

        ch_trimgalore_input = SEQKIT_RMDUP_R1.out.unique
        .join(SEQKIT_RMDUP_R2.out.unique)
        .map{ meta, r1, r2 -> [meta, [r1, r2]]}

        // Trimming reads
        TRIMGALORE(ch_trimgalore_input)
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    } else {
        // Trimming reads
        TRIMGALORE(ch_cat_fastq)
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    }

    // aligning with bwa: need reads in the same format used with FASTQC, a index file
    // which can be read from BWA_INDEX.out emit channel (https://www.nextflow.io/docs/edge/dsl2.html#process-named-output)
    // third params is if we want to sort data or not (true, so I don't need to SORT
    // with an additional step)
    BWA_MEM(TRIMGALORE.out.reads, PREPARE_GENOME.out.bwa_index, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // Perform Picard MarkDuplicates, index CRAM file and run samtools stats, flagstat and idxstats
    CRAM_MARKDUPLICATES_PICARD(BWA_MEM.out.cram, PREPARE_GENOME.out.genome_fasta, PREPARE_GENOME.out.genome_fasta_fai)
    ch_versions = ch_versions.mix(CRAM_MARKDUPLICATES_PICARD.out.versions)

    // prepare to call freebayes (multi) - get rid of meta.id
    freebayes_input_cram = CRAM_MARKDUPLICATES_PICARD.out.cram.map{ meta, cram -> [cram] }.collect().map{ it -> [[id: "all-samples"], it]}
    freebayes_input_crai = CRAM_MARKDUPLICATES_PICARD.out.crai.map{ meta, crai -> [crai] }.collect().map{ it -> [[id: "all-samples"], it]}

    // call freebayes paralle
    CRAM_FREEBAYES_PARALLEL(
        freebayes_input_cram,
        freebayes_input_crai,
        PREPARE_GENOME.out.genome_fasta,
        PREPARE_GENOME.out.genome_fasta_fai
    )
    ch_versions = ch_versions.mix(CRAM_FREEBAYES_PARALLEL.out.versions)

    // create bcftools channel. Freebayes multi will emit a single value for vcf and indexes.
    // join it and then change meta key to avoid file name collisions (meta is used to
    // determine output files)
    bcftools_ch = CRAM_FREEBAYES_PARALLEL.out.vcf
        .join(CRAM_FREEBAYES_PARALLEL.out.tbi)
        .map{ it -> [[id: "all-samples-normalized"], it[1], it[2]]}

    // normalize VCF (see https://github.com/freebayes/freebayes#normalizing-variant-representation)
    BCFTOOLS_NORM(
        bcftools_ch,
        PREPARE_GENOME.out.genome_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // index normalized VCF file
    TABIX_TABIX(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    // prepare input for bcftools stats
    bcftools_in_ch = BCFTOOLS_NORM.out.vcf
        .join(TABIX_TABIX.out.tbi)
        // .view()

    BCFTOOLS_STATS(
        bcftools_in_ch,
        [[], []],
        [[], []],
        [[], []],
        [[], []],
        [[], []]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    // for simplicity, I will collect SnpEff report in a new channel, in order
    // to use it as an empty chennel in MultiQC if SnpEff is not used
    snpeff_report = Channel.empty()

    // check for SnpEff annotation
    if (params.snpeff_database) {
        // annotate VCF with SnpEff
        SNPEFF_ANNOTATE(
        params.snpeff_database,
        BCFTOOLS_NORM.out.vcf,
        )

        // track version
        ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions)

        // update snpeff_report channel
        snpeff_report = SNPEFF_ANNOTATE.out.report
    }

    // get only the data I need for a MultiQC step
    ch_multiqc_files = FASTQC.out.html.map{it[1]}.ifEmpty([])
        .concat(FASTQC.out.zip.map{it[1]}.ifEmpty([]))
        .concat(TRIMGALORE.out.log.map{it[1]}.ifEmpty([]))
        .concat(CRAM_MARKDUPLICATES_PICARD.out.metrics.map{it[1]}.ifEmpty([]))
        .concat(CRAM_MARKDUPLICATES_PICARD.out.stats.map{it[1]}.ifEmpty([]))
        .concat(CRAM_MARKDUPLICATES_PICARD.out.idxstats.map{it[1]}.ifEmpty([]))
        .concat(CRAM_MARKDUPLICATES_PICARD.out.flagstat.map{it[1]}.ifEmpty([]))
        .concat(CRAM_MARKDUPLICATES_PICARD.out.coverage.map{it[1]}.ifEmpty([]))
        .concat(BCFTOOLS_STATS.out.stats.map{it[1]}.ifEmpty([]))
        .concat(snpeff_report.map{it[1]}.ifEmpty([]))
        // .view()

    // prepare multiqc_config file
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

    // TODO: requires nf-validation plugin and other stuff
    // summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    // calling MultiQC
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // return software version
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}
