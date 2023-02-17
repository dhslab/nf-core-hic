process MAKEBIGWIG {
    tag "$meta.sample"
    label 'process_high'
    label 'per_sample'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0' :
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0' }"

    input:

        tuple val(meta), path (cram), path(crai)
        path (reference_fasta) // genome fasta

    output:
        tuple val(meta), path ("${meta.sample}.bw") , emit: bigwig
        path ("versions.yml")                       , emit: versions

    script:
        """
        samtools view -b -T ${reference_fasta} ${cram} | \\
        bamCoverage -b - -o ${meta.sample}.bw

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
        END_VERSIONS
        """
}