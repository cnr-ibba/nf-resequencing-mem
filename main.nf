#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cnr-ibba/nf-resequencing-mem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/cnr-ibba/nf-resequencing-mem
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp  } from 'plugin/nf-validation'
include { PIPELINE_INITIALIZATION         } from './subworkflows/local/pipeline_initialization'
include { NORMALIZATION_INITIALIZATION    } from './subworkflows/local/pipeline_initialization'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RESEQUENCING_MEM              } from './workflows/resequencing-mem'
include { NORMALIZE_VCF                 } from './subworkflows/local/normalize_vcf'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from './modules/nf-core/custom/dumpsoftwareversions/main'

//
// WORKFLOW: Run main cnr-ibba/nf-resequencing-mem analysis pipeline
//
workflow CNR_IBBA {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    RESEQUENCING_MEM (
        samplesheet
    )

    emit:
    multiqc_report = RESEQUENCING_MEM.out.multiqc_report // channel: /path/to/multiqc_report.html
}

workflow VCF_NORMALIZE {
    take:
    vcf_ch // channel: vcf file
    tbi_ch // channel: tbi file
    fasta_ch // channel: fasta file

    main:
    // collect software version
    ch_versions = Channel.empty()

    // calling the normalization workflow
    NORMALIZE_VCF(
        vcf_ch,
        tbi_ch,
        fasta_ch
    )
    ch_versions = ch_versions.mix(NORMALIZE_VCF.out.versions)

    // return software version
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN CNR_IBBA:RESEQUENCING_MEM WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    main:

    // Print help message if needed
    if (params.help) {
        def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
        def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
        def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome_fasta /path/to/genome.fasta -profile docker"
        log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
        System.exit(0)
    }

    // Validate input parameters
    if (params.validate_params) {
        validateParameters()
    }

    // Initialize the workflow and check specific parameters
    WorkflowMain.initialise(workflow, params, log)

    if (!params.normalization_only) {
        // doing the main analysis
        // Run initializations tasks
        PIPELINE_INITIALIZATION (
            params.input
        )

        // then run the main pipeline
        CNR_IBBA (
            PIPELINE_INITIALIZATION.out.samplesheet
        )
    } else {
        // doing only the normalization workflow
        // setting up
        NORMALIZATION_INITIALIZATION(
            params.input_vcf,
            params.input_tbi,
            params.genome_fasta
        )

        // run only the normalization workflow
        VCF_NORMALIZE (
            NORMALIZATION_INITIALIZATION.out.vcf_ch,
            NORMALIZATION_INITIALIZATION.out.tbi_ch,
            NORMALIZATION_INITIALIZATION.out.fasta_ch
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
