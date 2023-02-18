process merge_pairs {
    tag "$meta.sample"
    label 'process_high'
    label 'per_sample'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (pairs)

    output:
        tuple val(meta), path ("${meta.sample}.pairs.gz")  , emit: pairs
        path ("versions.yml")                              , emit: versions

    script:
        """
        pairtools merge --nproc ${task.cpus} --tmpdir tmp -o ${meta.sample}.pairs.gz ${pairs}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
        END_VERSIONS
        """
}