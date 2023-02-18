process combine_lanes {
    tag "$meta.library"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:
        tuple val(meta), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("${meta.library}.1.fastq.gz"), path("${meta.library}.2.fastq.gz")    , emit: fastq

    script:
    """
    cat ${fastq1} > ${meta.library}.1.fastq.gz
    cat ${fastq2} > ${meta.library}.2.fastq.gz
    """
}