process merge_crams {
    tag "$meta.sample"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (crams)
        path (reference_fasta) // genome fasta

    output:
        tuple val(meta), path ("${meta.sample}.pairs.cram") , path ("${meta.sample}.pairs.cram.crai"), emit: cram_crai
        path ("versions.yml")                         , emit: versions

    script:
        """
        samtools merge -@ ${task.cpus} --reference ${reference_fasta} --write-index -o ${meta.sample}.pairs.cram##idx##${meta.sample}.pairs.cram.crai ${crams}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}
