// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FREEBAYES_SINGLE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3"
    } else {
        container "quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3"
    }

    input:
    tuple val(meta), path(bam)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi") , emit: index
    path  "versions.yml"                  , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    freebayes $options.args -b $bam --standard-filters -f $genome_fasta | bgzip --threads $task.cpus --stdout > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(freebayes --version 2>&1 | sed 's/^version: * v//')
    END_VERSIONS
    """
}
