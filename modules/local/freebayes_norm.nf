
process FREEBAYES_NORM {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "bioconda::freebayes=1.3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.8--h6a68c12_0':
        'biocontainers/freebayes:1.3.8--h6a68c12_0' }"

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
        -k \\
        ${vcf} \\
    | bgzip \\
        $args2 \\
        --threads $task.cpus \\
        --stdout > ${prefix}.vcf.gz

    # collect exit statuses
    exit_status_vcfallelicprimitives=\${PIPESTATUS[0]:-0}
    exit_status_bgzip=\${PIPESTATUS[1]:-0}

    # Checking 1st exit status
    if [ "\$exit_status_vcfallelicprimitives" -ne 0 ]; then
        # a 132 status is a warning not an error
        if [ "\$exit_status_vcfallelicprimitives" -ne 132 ]; then
            echo "Critical error in vcfallelicprimitives! Exit status \$exit_status_vcfallelicprimitives" >&2
            exit 1
        fi
    fi

    # Checking 2nd exit status
    if [ "\$exit_status_bgzip" -ne 0 ]; then
        echo "Critical error in bgzip! Exit status \$exit_status_bgzip" >&2
        exit 1
    fi

    # turning on pipefail
    set -o pipefail

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfallelicprimitives: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
