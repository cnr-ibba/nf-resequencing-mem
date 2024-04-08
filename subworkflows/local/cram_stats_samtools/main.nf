//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_COVERAGE } from '../../../modules/nf-core/samtools/coverage/main'

workflow CRAM_STATS_SAMTOOLS {
    take:
    ch_cram_crai  // channel: [ val(meta), path(cram), path(crai) ]
    ch_fasta      // channel: [ val(meta), path(fasta) ]
    ch_faidx      // channel: [ val(meta), path(faidx) ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS ( ch_cram_crai, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    SAMTOOLS_FLAGSTAT ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    SAMTOOLS_IDXSTATS ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    SAMTOOLS_COVERAGE ( ch_cram_crai, ch_fasta, ch_faidx )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    coverage = SAMTOOLS_COVERAGE.out.coverage // channel: [ val(meta), path(coverage) ]

    versions = ch_versions                    // channel: [ path(versions.yml) ]
}
