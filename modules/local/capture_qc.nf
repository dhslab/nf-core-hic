process capture_qc {
    tag "$meta.library"
    label 'process_medium'
    label 'per_library'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (cram), path(crai), path(dedup_stats)
        path (reference_fasta) // genome fasta
        path (bed) // capture baits bed

    output:
        tuple val(meta), path ("${meta.library}.capture_qc.txt") , emit: capture_qc
        path ("versions.yml")                                    , emit: versions

    script:
        """
        get_qc.py -p ${dedup_stats} > ${meta.library}".capture_qc.txt"
        mosdepth -f ${reference_fasta} -t ${task.cpus} -b ${bed} -x ${meta.library} -n ${cram}
        echo -e 'Overall coverage depth\t'\$(tail -n 2 ${meta.library}.mosdepth.summary.txt | head -n 1) >> ${meta.library}.capture_qc.txt
        echo -e 'Target coverage depth\t'\$(tail -n 1 ${meta.library}.mosdepth.summary.txt) >> ${meta.library}.capture_qc.txt


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
        END_VERSIONS
        """
}