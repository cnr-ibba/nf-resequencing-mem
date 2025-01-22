
process FREEBAYES_NORM {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "bioconda::freebayes=1.3.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hb089aa1_0':
        'biocontainers/freebayes:1.3.6--hb089aa1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcfallelicprimitives \\
        $args \\
        -kg \\
        ${vcf} \\
    | bgzip \\
        $args2 \\
        --threads $task.cpus \\
        --stdout > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfallelicprimitives: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
