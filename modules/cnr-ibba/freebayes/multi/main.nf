
process FREEBAYES_MULTI {
    tag "freebayes.multi"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6 main::numpy main::scipy" : null)
    container "bunop/freebayes:v0.1"

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
    def args2 = task.ext.args2 ?: ''
    def prefix = prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls $bam | xargs -n1 > bam_list.txt

    freebayes-parallel \\
        <(split_ref_by_bai_datasize.py --bam-list bam_list.txt --reference-fai $genome_fasta_fai $args2 -v | awk -F '[[:space:]]+' '{printf("%s:%d-%d\\n", \$1, \$2, \$3)}') \\
        $task.cpus \\
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
