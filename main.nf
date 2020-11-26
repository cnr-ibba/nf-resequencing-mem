#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { BWA_INDEX } from './modules/nf-core/software/bwa/index/main' addParams( options: [:] )
// override default publish_dir_mode: I don't want to copy a BAM file outside the "work" directory
include { BWA_MEM } from './modules/nf-core/software/bwa/mem/main' addParams( options: [:], publish_dir_mode: "symlink" )

// a function to reads from reads file channel and convert this to the format used by
// imported workflows
def get_reads(reads_path) {
  // define a simple input channel derived form a glob path matching reads.
  reads = Channel
    .fromFilePairs( params.reads_path, checkIfExists: true )//.view()
    .ifEmpty { error "Cannot find any reads matching: '${ params.reads_path }'" }

  // Functions implicitly return the result of the last evaluated statement. However:
  return reads.map( sample -> [ [id: sample[0], single_end: false], sample[1] ] )//.view()
}

// A workflow definition which does not declare any name is assumed to be the
// main workflow and it’s implicitly executed. Therefore it’s the entry point
// of the workflow application.
workflow {
  // get reads in the format required for fastqc module
  fastqc_input = get_reads(params.reads_path)

  // call FASTQC from module
  FASTQC(fastqc_input)

  // reference to an emit channel to be used out FASTQC scope (https://www.nextflow.io/docs/edge/dsl2.html#process-named-output)
  // FASTQC.out.html.view()

  // TODO: add a multiqc step

  // indexing genome
  BWA_INDEX(file(params.genome_path, checkIfExists: true))

  // aligning reads to genome: in DSL2 I can't using into to duplicate the channel
  // like DSL1, but I can read the output from another workflow or open a new channel
  // call the same function as before for FASTQC:
  bwa_input = get_reads(params.reads_path)

  // aligning with bwa: need reads in the same format used with FASTQC, a index file
  // (which can be read from BWA_INDEX.out and the genome used)
  BWA_MEM(bwa_input, BWA_INDEX.out.index, file(params.genome_path, checkIfExists: true))
}
