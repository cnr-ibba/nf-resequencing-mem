
process BAMADDRG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamaddrg:9baba65f88228e55639689a3cea38dd150e6284f--ha89c123_1':
        'biocontainers/bamaddrg:9baba65f88228e55639689a3cea38dd150e6284f--ha89c123_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '9baba65f88228e55639689a3cea38dd150e6284f' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    bamaddrg \\
        -b $bam \\
        -s $prefix \\
        > ${prefix}.rg.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamaddrg: $VERSION
    END_VERSIONS
    """
}
