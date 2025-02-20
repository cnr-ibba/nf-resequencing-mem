process SAMTOOLS_DEPTH {
    tag "$meta3.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta1), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), val(chromosome)

    output:
    tuple val(meta3), path("*.tsv.gz"),     emit: depth
    tuple val(meta1), path("*.list.txt"),   emit: bam_list
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    def positions = chromosome ? "-r ${chromosome}" : ""
    """
    ls $bam | xargs -n1 > ${prefix}.list.txt

    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        $positions \\
        -H \\
        -f ${prefix}.list.txt \\
        | gzip -c > ${prefix}.${meta3.id}.depth.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
