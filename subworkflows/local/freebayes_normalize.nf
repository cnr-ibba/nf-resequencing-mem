//
// Normalize VCF using freebayes and bcftools
//

include { FREEBAYES_NORM              } from '../../modules/local/freebayes_norm'
include { BCFTOOLS_SORT               } from '../../modules/nf-core/bcftools/sort/main'
include {
    TABIX_TABIX as BCFTOOLS_SORT_TBI;
    TABIX_TABIX as BCFTOOLS_NORM_TBI; } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM               } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILLTAGS           } from '../../modules/local/bcftools_filltags'


workflow FREEBAYES_NORMALIZE {
    take:
        vcf_ch    // channel: [mandatory] the VCF file to normalize
        ref_fasta // value: [mandatory] the reference fasta file coming from PREPARE_GENOME

    main:
        ch_versions = Channel.empty()

        // normalize input using vcfallelicprimitives
        FREEBAYES_NORM(vcf_ch)
        ch_versions = ch_versions.mix(FREEBAYES_NORM.out.versions)

        // sort VCF file
        BCFTOOLS_SORT(FREEBAYES_NORM.out.vcf.map{ meta, vcf -> [[id: "${meta.id}.sorted"], vcf] },)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        // index sorted vcf
        BCFTOOLS_SORT_TBI(BCFTOOLS_SORT.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT_TBI.out.versions)

        // prepare to normalize with bcftools
        normalize_in_ch = BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT_TBI.out.tbi)
            .map{ meta, vcf, tbi -> [[id: "${meta.id}.bcftools-normalized"], vcf, tbi]}

        // normalize with bcftools (2nd normalization, after vcfallelicprimitives)
        BCFTOOLS_NORM(
            normalize_in_ch,
            ref_fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

        // index normalized vcf
        BCFTOOLS_NORM_TBI(BCFTOOLS_NORM.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_TBI.out.versions)

        // add missing tags
        BCFTOOLS_FILLTAGS(BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM_TBI.out.tbi))
        ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

    emit:
        vcf           = BCFTOOLS_FILLTAGS.out.vcf
        tbi           = BCFTOOLS_FILLTAGS.out.index
        versions      = ch_versions
}
