//
// Prepare and run freebayes parallel on CRAM files.
//

include { SAMTOOLS_DEPTH                        }from '../../../modules/nf-core/samtools/depth/main'
include { FREEBAYES_SPLITCRAM                   } from '../../../modules/local/freebayes_splitcram'
include { FREEBAYES_CHUNK                       } from '../../../modules/cnr-ibba/freebayes/chunk/main'
include { BCFTOOLS_CONCAT as FREEBAYES_CONCAT   } from '../../../modules/cnr-ibba/bcftools/concat/main'
include {
    TABIX_TABIX as FREEBAYES_CONCAT_TABIX;
    TABIX_TABIX as REMOVE_OVERLAP_TABIX         } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_OVERLAP       } from '../../../modules/nf-core/bcftools/norm/main'

workflow CRAM_FREEBAYES_PARALLEL {
    take:
    bam     // channel: [ val(meta), [ bam/cram ]]
    bai     // channel: [ val(meta), [ bai/crai ]]
    fasta   // channel: [ val(meta2), fasta ]
    fai     // channel: [ val(meta2), fai ]

    main:
    ch_versions = Channel.empty()

    // read chromosome list from fasta index
    chromosome_ch = fai.map{ it -> it[1] }
        .splitCsv(header: false, sep: '\t')
        .map{ row -> [[id: row[0], length: row[1]], row[0]] }
        // .view()

    // calculate total coverage depth for all samples by chromosome
    // bam is a value channel; chromosome_ch is a queue channel
    // for every chromosome
    SAMTOOLS_DEPTH( bam, bai, chromosome_ch )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    // split fasta in chunks relying BAM size
    FREEBAYES_SPLITCRAM ( SAMTOOLS_DEPTH.out.depth )
    ch_versions = ch_versions.mix(FREEBAYES_SPLITCRAM.out.versions)

    // create a channel from region list file
    regions_ch = FREEBAYES_SPLITCRAM.out.regions
        .map{ it -> it[1]}
        .splitText()
        .map{ it -> [[id: it.trim()], it.trim()]}

    // call freebayes on each region
    FREEBAYES_CHUNK (
        regions_ch,
        bam,
        bai,
        // converting to a value channel
        SAMTOOLS_DEPTH.out.bam_list.first(),
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(FREEBAYES_CHUNK.out.versions)

    // merge freebayes chunks
    vcf_ch = FREEBAYES_CHUNK.out.vcf
        .map{ meta, region -> {
            def chromosome = meta.id.tokenize(":")[0]
            [ [id: chromosome], region ]
        }}
        .groupTuple()
        // .view()
    tbi_ch = FREEBAYES_CHUNK.out.index
        .map{ meta, region -> {
            def chromosome = meta.id.tokenize(":")[0]
            [ [id: chromosome], region ]
        }}
        .groupTuple()
        // .view()

    FREEBAYES_CONCAT ( vcf_ch.join(tbi_ch) )
    ch_versions = ch_versions.mix(FREEBAYES_CONCAT.out.versions)

    // create index
    FREEBAYES_CONCAT_TABIX ( FREEBAYES_CONCAT.out.vcf )
    ch_versions = ch_versions.mix(FREEBAYES_CONCAT_TABIX.out.versions)

    // create bcftools channel. Freebayes multi will emit a single value for vcf and indexes.
    // join it and then change meta key to avoid file name collisions (meta is used to
    // determine output files)
    bcftools_in_ch = FREEBAYES_CONCAT.out.vcf
        .join(FREEBAYES_CONCAT_TABIX.out.tbi)

    // normalize VCF (see https://github.com/freebayes/freebayes#normalizing-variant-representation)
    // required to remove overlapping regions after concatenation
    REMOVE_OVERLAP(
        bcftools_in_ch,
        fasta
    )
    ch_versions = ch_versions.mix(REMOVE_OVERLAP.out.versions)

    // create index
    REMOVE_OVERLAP_TABIX ( REMOVE_OVERLAP.out.vcf )
    ch_versions = ch_versions.mix(REMOVE_OVERLAP_TABIX.out.versions)

    emit:
    vcf      = REMOVE_OVERLAP.out.vcf           // channel: [ val(meta), [ vcf ] ]
    tbi      = REMOVE_OVERLAP_TABIX.out.tbi     // channel: [ val(meta), [ tbi ] ]
    csi      = REMOVE_OVERLAP_TABIX.out.csi     // channel: [ val(meta), [ csi ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
