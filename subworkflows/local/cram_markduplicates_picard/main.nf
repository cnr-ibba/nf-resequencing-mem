//
// Picard MarkDuplicates, index CRAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { CRAM_STATS_SAMTOOLS    } from '../cram_stats_samtools/main'

workflow CRAM_MARKDUPLICATES_PICARD {

    take:
    ch_cram  // channel: [ val(meta), path(cram) ]
    ch_fasta // channel: [ path(fasta) ]
    ch_fai   // channel: [ path(fai) ]

    main:

    ch_versions = Channel.empty()

    PICARD_MARKDUPLICATES ( ch_cram, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.cram )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_cram_crai = PICARD_MARKDUPLICATES.out.cram
        .join(SAMTOOLS_INDEX.out.crai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map{meta, cram, crai, csi ->
            if (crai) [ meta, cram, crai ]
            else [ meta, cram, csi ]
        }

    CRAM_STATS_SAMTOOLS ( ch_cram_crai, ch_fasta )
    ch_versions = ch_versions.mix(CRAM_STATS_SAMTOOLS.out.versions)

    emit:
    cram     = PICARD_MARKDUPLICATES.out.cram    // channel: [ val(meta), path(cram) ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_INDEX.out.crai          // channel: [ val(meta), path(crai) ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), path(csi) ]

    stats    = CRAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = CRAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = CRAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]
    coverage = CRAM_STATS_SAMTOOLS.out.coverage   // channel: [ val(meta), path(coverage) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
