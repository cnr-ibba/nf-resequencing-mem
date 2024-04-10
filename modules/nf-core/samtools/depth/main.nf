process SAMTOOLS_DEPTH {
    tag "$meta1.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta1), path(bam)
    tuple val(meta2), path(intervals)

    output:
    tuple val(meta1), path("*.tsv.gz"),     emit: depth
    tuple val(meta1), path("*.list.txt"),   emit: bam_list
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    def positions = intervals ? "-b ${intervals}" : ""
    """
    ls $bam | xargs -n1 > ${prefix}.list.txt

    samtools \\
        depth \\
        --threads ${task.cpus} \\
        $args \\
        $positions \\
        -H -a \\
        -f ${prefix}.list.txt \\
        | bgzip \\
        --threads ${task.cpus} \\
        --stdout > ${prefix}.depth.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
