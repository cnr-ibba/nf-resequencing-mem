#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2
version = '0.1.0'

include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )

// A workflow definition which does not declare any name is assumed to be the
// main workflow and it’s implicitly executed. Therefore it’s the entry point
// of the workflow application.
workflow {
  // Within this workfow, a simple input channel derived form a glob path matching reads is defined.
  Channel
    .fromFilePairs( params.reads_path, checkIfExists: true )//.view()
    .ifEmpty { error "Cannot find any reads matching: '${ params.reads_path }'" }
    .set { reads }

  input = reads.map( sample -> [ [id: sample[0], single_end: false], sample[1] ] )//.view()

  // call FASTQC from module
  FASTQC(input)

  // reference to an emit channel to be used out FASTQC scope (https://www.nextflow.io/docs/edge/dsl2.html#process-named-output)
  FASTQC.out.html.view()
}
