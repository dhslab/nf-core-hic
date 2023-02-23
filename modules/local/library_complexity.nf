process library_complexity {
    tag "$meta.library"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (cram), path(crai)
        path (reference_fasta) // genome fasta

    output:
        tuple val(meta), path ("${meta.library}.preseq.txt") , emit: qc
        path ("versions.yml")                                , emit: versions

    script:
        """
        samtools view -b -T ${reference_fasta} ${cram} | \\
        bamToBed -i - | preseq lc_extrap /dev/stdin -pe -extrap 2100000000 -step 100000000 -seg_len 1000000000 -output ${meta.library}.preseq.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
}
