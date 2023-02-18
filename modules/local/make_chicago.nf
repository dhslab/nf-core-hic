process make_chicago {
    tag "$meta.sample"
    label 'process_medium'
    label 'per_sample'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (cram), path(crai)
        path (reference_fasta) // genome fasta

    output:
        tuple val(meta), path ("${meta.sample}.chicago.bam") , emit: chicago_bam
        path ("versions.yml")                                , emit: versions

    script:
        """
        if (( ${task.cpus} > 1 )); then MAXTHREADS=\$(( ${task.cpus} - 1 )) ; else MAXTHREADS=1 ; fi

        samtools view -@ 1 -T ${reference_fasta} -Shu -F 2048 ${cram} | samtools sort -n -T tmp/temp.bam --threads \$MAXTHREADS -o ${meta.sample}.chicago.bam -
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}