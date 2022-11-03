
process FREEBAYES_MULTI {
    tag "freebayes.multi"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hb089aa1_0' :
        'quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path(genome_fasta)
    path(genome_fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls $bam | xargs -n1 > bam_list.txt

    freebayes-parallel \\
        <(fasta_generate_regions.py $genome_fasta_fai 100000) $task.cpus \\
        $args \\
        --bam-list bam_list.txt \\
        --standard-filters \\
        -f $genome_fasta \\
    | bgzip \\
        --threads $task.cpus \\
        --stdout > ${prefix}.vcf.gz

    tabix ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
