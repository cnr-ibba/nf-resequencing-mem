
process FREEBAYES_CHUNK {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.8--h6a68c12_0':
        'biocontainers/freebayes:1.3.8--h6a68c12_0' }"

    input:
    tuple val(meta),  val(region)
    tuple val(meta2), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(genome_fasta)
    tuple val(meta3), path(genome_fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi") , emit: index
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy = params.ploidy == 2 ? "": "--ploidy ${params.ploidy}"
    def gvcf = params.gvcf ? "--gvcf" : ""
    def gvcf_chunk = params.gvcf_chunk ? "--gvcf-chunk ${params.gvcf_chunk}" : ""
    def gvcf_dont_use_chunk = params.gvcf_dont_use_chunk ? "--gvcf-dont-use-chunk true" : ""
    """
    ls $bam | xargs -n1 | sort > ${prefix}.list.txt

    freebayes \\
        $args \\
        $ploidy \\
        $gvcf \\
        $gvcf_chunk \\
        $gvcf_dont_use_chunk \\
        --bam-list ${prefix}.list.txt \\
        --standard-filters \\
        -f $genome_fasta \\
        --region $region \\
    | bgzip \\
        --threads $task.cpus \\
        --stdout > ${prefix}.vcf.gz

    tabix ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
