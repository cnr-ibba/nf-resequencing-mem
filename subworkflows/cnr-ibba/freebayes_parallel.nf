//
// Prepare and run freebayes paralle
//

include { FREEBAYES_SPLITBAM } from '../../modules/cnr-ibba/freebayes/splitbam/main'
include { FREEBAYES_CHUNK } from '../../modules/cnr-ibba/freebayes/chunk/main'
include { FREEBAYES_MERGE } from '../../modules/cnr-ibba/freebayes/merge/main'

workflow FREEBAYES_PARALLEL {
    take:
    bam     // channel: [ val(meta), [ bam/cram ]]
    bai     // channel: [ val(meta), [ bai/crai ]]
    fasta   // channel: [ fasta ]
    fai     // channel: [ fai ]

    main:
    ch_versions = Channel.empty()

    // split fasta in chunks relying BAM size
    FREEBAYES_SPLITBAM ( bam, bai, fasta, fai )
    ch_versions = ch_versions.mix(FREEBAYES_SPLITBAM.out.versions)

    // create a channel from region list file
    regions_ch = FREEBAYES_SPLITBAM.out.regions
        .map{ it -> it[1]}
        .splitText()
        .map{ it -> [[id: it.trim()], it.trim()]}

    // call freebayes on each region
    FREEBAYES_CHUNK ( regions_ch, bam, bai, FREEBAYES_SPLITBAM.out.bam_list, fasta, fai )
    ch_versions = ch_versions.mix(FREEBAYES_CHUNK.out.versions)

    // merge freebayes chunks
    vcf_ch = FREEBAYES_CHUNK.out.vcf
        .collect{ it -> it[1]}
        .map{ it -> [[id: 'all-samples'], it]}

    FREEBAYES_MERGE ( vcf_ch )
    ch_versions = ch_versions.mix(FREEBAYES_MERGE.out.versions)

    emit:
    vcf      = FREEBAYES_MERGE.out.merged_vcf   // channel: [ val(meta), [ vcf ] ]
    index    = FREEBAYES_MERGE.out.merged_index // channel: [ val(meta), [ vcf ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
