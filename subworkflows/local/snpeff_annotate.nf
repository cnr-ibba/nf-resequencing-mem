//
// Annotate VCF file with snpEff
//

include { SNPEFF_DOWNLOAD } from '../../modules/nf-core/snpeff/download/main'
include { SNPEFF_SNPEFF   } from '../../modules/nf-core/snpeff/snpeff/main'


workflow SNPEFF_ANNOTATE {
  take:
    genome  // value: the SnpEff genome database to use

  main:
    ch_versions = Channel.empty()

    // create a input channel for SnpEff download
    snpeff_genome = [[id: genome], genome]
    SNPEFF_DOWNLOAD(snpeff_genome)

    // track version
    ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)

  emit:
    versions = ch_versions
}
