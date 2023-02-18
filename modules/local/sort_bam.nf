process sort_bam {
    tag "$meta.library"
    label 'process_high'
    label 'per_library'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (bam)
        path (reference_fasta) // genome fasta

    output:
        tuple val(meta), path ("${meta.library}.cram") , path ("${meta.library}.cram.crai"), emit: cram_crai
        path ("versions.yml")                          , emit: versions

    script:
        """

        samtools sort -@ ${task.cpus} --reference ${reference_fasta} -o ${meta.library}.cram ${bam} && samtools index -@ ${task.cpus} ${meta.library}.cram

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}