//
// Normalize VCF using freebayes and bcftools
//

include { VCFLIB_VCFWAVE                    } from '../../modules/local/vcflib_vcfwave'
include { BCFTOOLS_SORT                     } from '../../modules/nf-core/bcftools/sort/main'
include {
    TABIX_TABIX as BCFTOOLS_SORT_TABIX;
    TABIX_TABIX as BCFTOOLS_NORM_TABIX;
    TABIX_TABIX as BCFTOOLS_FILLTAGS_TABIX; } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM                     } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILLTAGS                 } from '../../modules/local/bcftools_filltags'


workflow NORMALIZE_VCF {
    take:
        vcf_ch    // channel: [mandatory] the VCF file to normalize
        tbi_ch    // channel: [mandatory] the index file for the VCF file
        ref_fasta // value: [mandatory] the reference fasta file coming from PREPARE_GENOME

    main:
        ch_versions = Channel.empty()

        // normalize input using vcflib/vcfwave (1st normalization)
        VCFLIB_VCFWAVE(vcf_ch.join(tbi_ch))
        ch_versions = ch_versions.mix(VCFLIB_VCFWAVE.out.versions)

        // sort VCF file
        BCFTOOLS_SORT(VCFLIB_VCFWAVE.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        // index sorted vcf
        BCFTOOLS_SORT_TABIX(BCFTOOLS_SORT.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT_TABIX.out.versions)

        // prepare to normalize with bcftools
        normalize_in_ch = BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT_TABIX.out.tbi)

        // normalize with bcftools (2nd normalization, after vcfwave)
        BCFTOOLS_NORM(
            normalize_in_ch,
            ref_fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

        // index normalized vcf
        BCFTOOLS_NORM_TABIX(BCFTOOLS_NORM.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_TABIX.out.versions)

        // add missing tags
        BCFTOOLS_FILLTAGS(BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM_TABIX.out.tbi))
        ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

        // index VCFs with tags
        BCFTOOLS_FILLTAGS_TABIX(BCFTOOLS_FILLTAGS.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS_TABIX.out.versions)

    emit:
        vcf           = BCFTOOLS_FILLTAGS.out.vcf
        tbi           = BCFTOOLS_FILLTAGS_TABIX.out.tbi
        versions      = ch_versions
}
