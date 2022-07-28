//
// PREPARE GENOME
//
// inspired from nf-core/sarek:subworkflows/local/prepare_genome.nf
//
// Initialize channels based on params or indices that were just built
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'
include { BWA_INDEX } from '../../modules/nf-core/modules/bwa/index/main'

workflow PREPARE_GENOME {
  take:
    genome_fasta      // channel: [mandatory] fasta
    genome_fasta_fai  // channel: [optional]  fasta_fai
    genome_bwa_index  // channel: [optional]  bwa index

  main:

    ch_versions = Channel.empty()

    // create fasta index if necessary
    if (! params.genome_fasta_fai) {
      SAMTOOLS_FAIDX(genome_fasta.map{ it -> [[id:it[0].getName()], it] })
      genome_fasta_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }

      // track version
      ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // indexing genome if necessary
    if (! params.genome_bwa_index) {
      BWA_INDEX(genome_fasta)
      genome_bwa_index = BWA_INDEX.out.index
      ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }


  emit:
    genome_fasta     = genome_fasta
    genome_fasta_fai = genome_fasta_fai
    bwa_index        = genome_bwa_index
    versions         = ch_versions
}
