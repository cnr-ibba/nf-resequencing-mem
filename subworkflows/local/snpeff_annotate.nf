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

    // create a input channel for SnpEff download
    snpeff_genome = [[id: genome], genome]
    SNPEFF_DOWNLOAD(snpeff_genome)

    // track version
    ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)

    // annotate the VCF file
    SNPEFF_SNPEFF(vcf, genome, SNPEFF_DOWNLOAD.out.cache)

    // track version
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

  emit:
    cache         = SNPEFF_DOWNLOAD.out.cache
    vcf           = SNPEFF_SNPEFF.out.vcf
    report        = SNPEFF_SNPEFF.out.report
    summary_html  = SNPEFF_SNPEFF.out.summary_html
    genes_txt     = SNPEFF_SNPEFF.out.genes_txt
    versions      = ch_versions
}
