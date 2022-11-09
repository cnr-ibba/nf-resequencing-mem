
process FREEBAYES_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hb089aa1_0':
        'quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.merged.vcf.gz")     , emit: merged_vcf
    tuple val(meta), path("*.merged.vcf.gz.tbi") , emit: merged_index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $vcf \\
        | vcffirstheader \\
        | vcfstreamsort -w 1000 \\
        | vcfuniq \\
        | bgzip \\
        --threads $task.cpus \\
        --stdout > ${prefix}.merged.vcf.gz

    tabix ${prefix}.merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
