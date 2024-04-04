//
// Sambamba MarkDuplicates, index CRAM file and run samtools stats, flagstat and idxstats
//

include { SAMBAMBA_MARKDUP      } from '../../../modules/nf-core/sambamba/markdup/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { CRAM_STATS_SAMTOOLS   } from '../cram_stats_samtools/main'

workflow CRAM_MARKDUPLICATES_SAMBAMBA {

    take:
    ch_cram  // channel: [ val(meta), path(cram) ]
    ch_fasta // channel: [ path(fasta) ]
    ch_fai   // channel: [ path(fai) ]

    main:

    ch_versions = Channel.empty()

    SAMBAMBA_MARKDUP ( ch_cram, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions.first())

    SAMTOOLS_INDEX ( SAMBAMBA_MARKDUP.out.cram )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_cram_crai = SAMBAMBA_MARKDUP.out.cram
        .join(SAMTOOLS_INDEX.out.crai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map{meta, cram, crai, csi ->
            if (crai) [ meta, cram, crai ]
            else [ meta, cram, csi ]
        }

    CRAM_STATS_SAMTOOLS ( ch_cram_crai, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix(CRAM_STATS_SAMTOOLS.out.versions)

    emit:
    cram     = SAMBAMBA_MARKDUP.out.cram          // channel: [ val(meta), path(cram) ]
    metrics  = SAMBAMBA_MARKDUP.out.metrics       // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_INDEX.out.crai            // channel: [ val(meta), path(crai) ]
    csi      = SAMTOOLS_INDEX.out.csi             // channel: [ val(meta), path(csi) ]

    stats    = CRAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = CRAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = CRAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]
    coverage = CRAM_STATS_SAMTOOLS.out.coverage   // channel: [ val(meta), path(coverage) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
