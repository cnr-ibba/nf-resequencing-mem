#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { FASTQC } from './modules/nf-core/modules/fastqc/main' addParams( options: [:] )
include { MULTIQC } from './modules/nf-core/modules/multiqc/main' addParams( options: [:] )
// for trimgalore, publish only reports (txt - not trimmed files!)
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main' addParams( options: [publish_files: ['report.txt': '']] )
include { BWA_INDEX } from './modules/nf-core/modules/bwa/index/main' addParams( options: [publish_files: false] )
// override default publish_dir_mode: I don't want to copy a BAM file outside the "work" directory
include { BWA_MEM } from './modules/nf-core/modules/bwa/mem/main' addParams( options: [publish_files: false] )
include { SAMTOOLS_SORT } from './modules/nf-core/modules/samtools/sort/main' addParams( options: [publish_files: false, suffix: '.sort'] )
include { BAMADDRG } from './modules/cnr-ibba/nf-modules/bamaddrg/main' addParams( options: [:] )
include { PICARD_MARKDUPLICATES } from './modules/nf-core/modules/picard/markduplicates/main' addParams( options: [publish_files: false] )
include { SAMTOOLS_INDEX } from './modules/nf-core/modules/samtools/index/main' addParams( options: [publish_files: false] )
include { SAMTOOLS_FLAGSTAT } from './modules/nf-core/modules/samtools/flagstat/main' addParams( options: [:] )
include { FREEBAYES_SINGLE } from './modules/cnr-ibba/nf-modules/freebayes/single/main' addParams( options: [:] )
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'  addParams( options: [publish_files : ['_versions.yml':'']] )

// a function to read from reads file channel and convert this to the format used by
// imported workflows
def get_reads(reads_path) {
  // define a simple input channel derived form a glob path matching reads.
  reads = Channel
    .fromFilePairs( params.reads_path, checkIfExists: true )//.view()
    .ifEmpty { error "Cannot find any reads matching: '${ params.reads_path }'" }

  // create the meta key used in module. Is the sample name plus a boolean value
  // for single_end reads.
  // Functions implicitly return the result of the last evaluated statement. However:
  return reads.map( sample -> [ [id: sample[0], single_end: false], sample[1] ] )//.view()
}

// A workflow definition which does not declare any name is assumed to be the
// main workflow and it’s implicitly executed. Therefore it’s the entry point
// of the workflow application.
workflow {
  // collect software version
  ch_versions = Channel.empty()

  // get reads in the format required for fastqc module
  fastqc_input = get_reads(params.reads_path)

  // call FASTQC from module
  FASTQC(fastqc_input)

  ch_versions = ch_versions.mix(FASTQC.out.versions.first())

  // get only the data I need for a MultiQC step
  html_report = FASTQC.out.html.map( sample -> sample[1] )
  zip_report = FASTQC.out.zip.map( sample -> sample[1] )

  // combine two channel (mix) and the get only one emission
  multiqc_input = html_report.mix(zip_report).collect()//.view()

  // calling MultiQC
  MULTIQC(multiqc_input)

  // indexing genome
  BWA_INDEX(file(params.genome_path, checkIfExists: true))

  // trimming sequences: in DSL2 I can't use 'into' to duplicate the channel
  // like DSL1, but I can read the output from another workflow or open a new channel.
  // Call the same function I used before for FASTQC:
  trimgalore_input = get_reads(params.reads_path)

  // Trimming reads
  TRIMGALORE(trimgalore_input)

  // aligning with bwa: need reads in the same format used with FASTQC, a index file
  // which can be read from BWA_INDEX.out emit channel (https://www.nextflow.io/docs/edge/dsl2.html#process-named-output)
  BWA_MEM(TRIMGALORE.out.reads, BWA_INDEX.out.index)

  // BAM file need to be sorted in order to mark duplicates. Please note that the output
  // file has been renamed by modiying the nextflow core SAMTOOLS_SORT modules. This is
  // required since samtools can't sort using the same file names for input and output
  // as nextflow does in work directory
  SAMTOOLS_SORT(BWA_MEM.out.bam)

  // add sample names to BAM files. Required to name samples with freebayes
  BAMADDRG(SAMTOOLS_SORT.out.bam)

  // markduplicates step. It requires meta information + bam files, the same output of
  // SAMTOOLS_SORT step
  PICARD_MARKDUPLICATES(BAMADDRG.out.bam)

  // bam need to be indexed before doing flagstat
  SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)

  // now I can do the flagstat step. I need bam and bai files from markduplicates and
  // samtools index and the meta informations as three different input. From both channels,
  // I have an output like 'val(meta), path("*.bam")' and 'val(meta), path("*.bai")'
  // I can join two channels with the same key (https://www.nextflow.io/docs/latest/operator.html#join)
  // two options to check to have exactly the same keys with no duplications
  flagstat_input = PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true)//.view()

  // time to call flagstat
  SAMTOOLS_FLAGSTAT(flagstat_input)

  // prepare to call freebayes (single) - remove meta key
  FREEBAYES_SINGLE(PICARD_MARKDUPLICATES.out.bam, file(params.genome_path, checkIfExists: true))

  // return software version
  CUSTOM_DUMPSOFTWAREVERSIONS (
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
  )

}
