process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-faf8276e4ce95c7ed33ebd0a33c53621fc720774:fe25dd03855090ac1ab4a90a97f2fff19641c7da-0' :
        'biocontainers/mulled-v2-faf8276e4ce95c7ed33ebd0a33c53621fc720774:fe25dd03855090ac1ab4a90a97f2fff19641c7da-0' }"

    input:
    tuple val(meta),    path(cram)
    tuple val(meta2),   path(fasta)
    tuple val(meta3),   path(fai)

    output:
    tuple val(meta), path("*.cram")       , emit: cram
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$cram" == "${prefix}.cram") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    mkfifo picard_output

    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $cram \\
        --OUTPUT picard_output \\
        --REFERENCE_SEQUENCE $fasta \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt &
    samtools view $args2 --threads $task.cpus -O cram -o ${prefix}.cram picard_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$cram" == "${prefix}.cram") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.cram
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
