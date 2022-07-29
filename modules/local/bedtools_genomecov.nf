
process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*cover*")       , emit: genomecov
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        $args \\
        -bga \\
        | gzip -c > ${prefix}.bga_coverage.gz

    bedtools \\
        genomecov \\
        -ibam $bam \\
        $args \\
        -d \\
        | gzip -c > ${prefix}.d_coverage.gz

    zcat ${prefix}.d_coverage.gz | awk '{total += \$3} END {print total}' > ${prefix}.seq_nr

    zcat ${prefix}.d_coverage.gz | awk 'END{print FNR}' > ${prefix}.positions

    zcat ${prefix}.d_coverage.gz | grep -w 0 | wc -l > ${prefix}.uncovered_positions

    awk 'BEGIN{getline to_add < "${prefix}.seq_nr"}{print \$0,to_add}' ${prefix}.positions > ${prefix}.cov_calculation

    awk '{cov = \$2/\$1} END {print cov}' ${prefix}.cov_calculation > ${prefix}.mean_coverage

    awk 'BEGIN{getline to_add < "${prefix}.uncovered_positions"}{print \$0,to_add}' ${prefix}.positions > ${prefix}.uncovered

    awk '{uncov = \$2*100/\$1} END {print uncov}' ${prefix}.uncovered > ${prefix}.uncovered_percentage

    rm ${prefix}.seq_nr
    rm ${prefix}.positions
    rm ${prefix}.uncovered_positions
    rm ${prefix}.cov_calculation
    rm ${prefix}.uncovered

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
