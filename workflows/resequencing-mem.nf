/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC } from './modules/nf-core/modules/fastqc/main'
include { MULTIQC } from './modules/nf-core/modules/multiqc/main'
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main'
include { BWA_INDEX } from './modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from './modules/nf-core/modules/bwa/mem/main'
include { BAMADDRG } from './modules/cnr-ibba/nf-modules/bamaddrg/main'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_FLAGSTAT } from './modules/nf-core/modules/samtools/flagstat/main'
include { FREEBAYES_SINGLE } from './modules/cnr-ibba/nf-modules/freebayes/single/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'

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
workflow RESEQUENCING_MEM {
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

  // prepare multiqc_config file
  multiqc_config = Channel.fromPath(params.multiqc_config)
  multiqc_config = multiqc_config.concat(Channel.fromPath(params.multiqc_logo))
  multiqc_config = multiqc_config.collect()

  // calling MultiQC
  MULTIQC(multiqc_input, multiqc_config)

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
  // third params is if we want to sort data or not (true, so I don't need to SORT
  // with an additional step)
  BWA_MEM(TRIMGALORE.out.reads, BWA_INDEX.out.index, true)

  // add sample names to BAM files. Required to name samples with freebayes
  BAMADDRG(BWA_MEM.out.bam)

  // BAM file need to be sorted in order to mark duplicates. This was don in BWA_MEM
  // step. Markduplicates requires meta information + bam files, the same output of
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
