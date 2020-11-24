#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
