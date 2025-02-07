//
// This file holds several functions specific to the main.nf workflow in the nf-core/test pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check at least one input has been provided
        if (!params.normalization_only) {
            // check for mandatory input
            if (!params.input) {
                Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
            }

            // check for gvcf_chunk options and gvcf
            if ((params.gvcf_chunk || params.gvcf_dont_use_chunk) && !params.gvcf) {
                Nextflow.error("Please provide '--gvcf' option when providing '--gvcf_chunk' or '--gvcf_dont_use_chunk' parameters")
            } else if (params.gvcf_chunk && params.gvcf_dont_use_chunk) {
                Nextflow.error("Please provide only one of '--gvcf_chunk' or '--gvcf_dont_use_chunk' parameters")
            }
        }

        // doing the normalization workflow
        if (params.normalization_only) {
            if (!params.input_vcf || !params.input_tbi) {
                Nextflow.error("Please provide a VCF file and its index to the pipeline e.g. '--input_vcf input.vcf --input_tbi input.vcf.tbi' when using '--normalization_only'")
            }
            if (params.input) {
                log.warn("You choose to run the normalization workflow. The input samplesheet will be ignored.")
            }
        }
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }
}
