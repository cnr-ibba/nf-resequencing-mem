process SNPEFF_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f:c51b2e46bf63786b2d9a7a7d23680791163ab39a-0' :
        'biocontainers/mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f:c51b2e46bf63786b2d9a7a7d23680791163ab39a-0' }"

    input:
    tuple val(meta), path(vcf)
    val   db
    tuple val(meta2), path(cache)
    tuple val(meta2), path(config)

    output:
    tuple val(meta), path("*.ann.vcf.gz"),      emit: vcf
    tuple val(meta), path("*.tbi"),             optional:true, emit: tbi
    tuple val(meta), path("*.csi"),             optional:true, emit: csi
    tuple val(meta), path("*.csv"),             emit: report
    tuple val(meta), path("*.html"),            emit: summary_html
    tuple val(meta), path("*.genes.txt"),       emit: genes_txt
    path "versions.yml"                 ,       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    def config_command = config ? "-c ${config}" : ""
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        $db \\
        $args \\
        -csvStats ${prefix}.csv \\
        $cache_command \\
        $config_command \\
        $vcf \\
        | bgzip -c $args2 -@${task.cpus} > ${prefix}.ann.vcf.gz

    tabix $args3 ${prefix}.ann.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf.gz
    touch ${prefix}.ann.vcf.gz.tbi
    touch ${prefix}.ann.csv
    touch ${prefix}.ann.genes.txt
    touch snpEff_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
