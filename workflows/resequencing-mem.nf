/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.genome_path ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { FASTQC } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { BWA_INDEX } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from '../modules/nf-core/modules/bwa/mem/main'
include { BAMADDRG } from '../modules/cnr-ibba/nf-modules/bamaddrg/main'
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/modules/samtools/flagstat/main'
include { FREEBAYES_SINGLE } from '../modules/cnr-ibba/nf-modules/freebayes/single/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// A workflow definition which does not declare any name is assumed to be the
// main workflow and it’s implicitly executed. Therefore it’s the entry point
// of the workflow application.
workflow RESEQUENCING_MEM {
  // collect software version
  ch_versions = Channel.empty()

  //
  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
  //
  INPUT_CHECK (
    ch_input
  )
  ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

  // call FASTQC from module
  FASTQC(INPUT_CHECK.out.reads)
  ch_versions = ch_versions.mix(FASTQC.out.versions)

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
  ch_versions = ch_versions.mix(MULTIQC.out.versions)

  // indexing genome
  BWA_INDEX(file(params.genome_path, checkIfExists: true))
  ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

  // Trimming reads
  TRIMGALORE(INPUT_CHECK.out.reads)
  ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

  // aligning with bwa: need reads in the same format used with FASTQC, a index file
  // which can be read from BWA_INDEX.out emit channel (https://www.nextflow.io/docs/edge/dsl2.html#process-named-output)
  // third params is if we want to sort data or not (true, so I don't need to SORT
  // with an additional step)
  BWA_MEM(TRIMGALORE.out.reads, BWA_INDEX.out.index, true)
  ch_versions = ch_versions.mix(BWA_MEM.out.versions)

  // add sample names to BAM files. Required to name samples with freebayes
  BAMADDRG(BWA_MEM.out.bam)
  ch_versions = ch_versions.mix(BAMADDRG.out.versions)

  // BAM file need to be sorted in order to mark duplicates. This was don in BWA_MEM
  // step. Markduplicates requires meta information + bam files, the same output of
  // SAMTOOLS_SORT step
  PICARD_MARKDUPLICATES(BAMADDRG.out.bam)
  ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

  // bam need to be indexed before doing flagstat
  SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
  ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

  // now I can do the flagstat step. I need bam and bai files from markduplicates and
  // samtools index and the meta informations as three different input. From both channels,
  // I have an output like 'val(meta), path("*.bam")' and 'val(meta), path("*.bai")'
  // I can join two channels with the same key (https://www.nextflow.io/docs/latest/operator.html#join)
  // two options to check to have exactly the same keys with no duplications
  flagstat_input = PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true)//.view()

  // time to call flagstat
  SAMTOOLS_FLAGSTAT(flagstat_input)
  ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

  // prepare to call freebayes (single) - remove meta key
  FREEBAYES_SINGLE(PICARD_MARKDUPLICATES.out.bam, file(params.genome_path, checkIfExists: true))
  ch_versions = ch_versions.mix(FREEBAYES_SINGLE.out.versions)

  // return software version
  CUSTOM_DUMPSOFTWAREVERSIONS (
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
  )

}
