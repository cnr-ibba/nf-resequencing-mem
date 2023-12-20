
process FREEBAYES_SPLITBAM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/bunop/freebayes:v0.1"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta2), path(genome_fasta)
    tuple val(meta2), path(genome_fasta_fai)

    output:
    tuple val(meta), path("*.list.txt"), emit: bam_list
    tuple val(meta), path("*.regions.txt"), emit: regions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ls $bam | xargs -n1 > ${prefix}.list.txt

    split_ref_by_bai_datasize.py \\
        --bam-list ${prefix}.list.txt \\
        --reference-fai $genome_fasta_fai \\
        $args \\
        -v | awk -F '[[:space:]]+' '{printf("%s:%d-%d\\n", \$1, \$2, \$3)}' > ${prefix}.regions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
