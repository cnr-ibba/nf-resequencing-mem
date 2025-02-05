
process VCFLIB_VCFWAVE {
    tag "$meta.id"
    label 'process_low'

    // TODO: conda version is 1.0.10 while container version is v1.0.12
    conda "bioconda::vcflib=1.0.10"
    container "docker.io/bunop/vcflib:0.1"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.10' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vcfwave \\
        $args \\
        $vcf \\
        | bgzip --threads $task.cpus -c $args2 > ${prefix}.vcfwave.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.10' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo | gzip > ${prefix}.vcfwave.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
