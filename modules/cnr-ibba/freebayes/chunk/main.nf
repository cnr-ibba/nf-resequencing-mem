
process FREEBAYES_CHUNK {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hb089aa1_0':
        'biocontainers/freebayes:1.3.6--hb089aa1_0' }"

    input:
    tuple val(meta),  val(region)
    tuple val(meta2), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta2), path(bam_list)
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
    """
    freebayes \\
        $args \\
        --bam-list $bam_list \\
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
