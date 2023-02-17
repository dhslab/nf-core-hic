process PAIRS2MCOOL {
    tag "$meta.sample"
    label 'process_high'
    label 'per_sample'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (pairs)
        path (chromsizes)

    output:
        tuple val(meta), path ("*.mcool")  , emit: mcool
        path ("versions.yml")              , emit: versions

    script:
        """
        MAPQ_LIST="1 30"
        for MAPQ in \$MAPQ_LIST; do
            bgzip -cd -@ ${task.cpus} ${pairs} | pairtools select "(mapq1>=\${MAPQ}) and (mapq2>=\${MAPQ})" | \\
                cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly ${params.genome} ${chromsizes}:${params.min_res} - ${meta.sample}.mapq_\${MAPQ}.cool && \\
                cooler zoomify --nproc ${task.cpus} --out ${meta.sample}.mapq_\${MAPQ}.mcool --resolutions ${params.resolutions} --balance ${meta.sample}.mapq_\${MAPQ}.cool
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
            cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
        END_VERSIONS
        """
}