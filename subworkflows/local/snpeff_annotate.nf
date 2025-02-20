//
// Annotate VCF file with SnpEff
//

include { SNPEFF_DOWNLOAD } from '../../modules/nf-core/snpeff/download/main'
include { SNPEFF_SNPEFF   } from '../../modules/nf-core/snpeff/snpeff/main'


workflow SNPEFF_ANNOTATE {
    take:
    genome    // value: [mandatory] the SnpEff genome database to use
    vcf       // channel: [mandatory] the VCF file to annotate

    main:
    ch_versions = Channel.empty()
    cache = Channel.empty()

    // define configuration for SnpEff
    config = Channel.fromPath(params.snpeff_config, checkIfExists: true)
        .map { it -> [[id:genome], it] }
        .first()

    if (! params.snpeff_cachedir) {
        // create a input channel for SnpEff download
        snpeff_genome = [[id: genome], genome]
        SNPEFF_DOWNLOAD(snpeff_genome)

        // export snpeff cache
        cache = SNPEFF_DOWNLOAD.out.cache

        // track version
        ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)
    } else {
        // use the provided cache
        cache = Channel.fromPath(params.snpeff_cachedir, checkIfExists: true)
            .map { it -> [[id:genome], it] }
            .first()
    }

    // annotate the VCF file
    SNPEFF_SNPEFF(vcf, Channel.value(genome), cache, config)

    // track version
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

    emit:
    cache         = cache
    vcf           = SNPEFF_SNPEFF.out.vcf
    report        = SNPEFF_SNPEFF.out.report
    summary_html  = SNPEFF_SNPEFF.out.summary_html
    genes_txt     = SNPEFF_SNPEFF.out.genes_txt
    versions      = ch_versions
}
