// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bwa=0.7.17=hed695b0_7 bioconda::samtools=1.10=h9402c20_2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  fasta

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    bwa mem \\
        $options.args \\
        $read_group \\
        -t $task.cpus \\
        $fasta \\
        $reads \\
        | samtools view $options.args2 -@ $task.cpus -bS -o ${prefix}.bam -

    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
