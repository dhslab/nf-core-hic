process FASTQ2PAIRS {
    tag "$meta.library"
    label 'process_high'
    label 'per_library'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-hic' :
        'ghcr.io/dhslab/docker-hic' }"

    input:

        tuple val(meta), path (fastq1), path (fastq2)
        path (bwa_index_files) // staged index files
        path (reference_fasta) // genome fasta
        path (chromsizes)

    output:
        tuple val(meta), path ("${meta.library}.fastp.json")                                            , emit: json
        tuple val(meta), path ("${meta.library}.fastp.html")                                            , emit: html
        tuple val(meta), path ("${meta.library}.bam")                                                   , emit: bam
        tuple val(meta), path ("${meta.library}.pairs.gz")                                              , emit: pairs
        tuple val(meta), path ("${meta.library}.pairs.cram") , path ("${meta.library}.pairs.cram.crai") , emit: cram_crai
        tuple val(meta), path ("${meta.library}.dedup.stats.txt")                                       , emit: stats
        path ("versions.yml")                                                                           , emit: versions

    script:
        """
        if (( ${task.cpus} > 3 )); then MAXTHREADS=\$(( ${task.cpus} - 3 )) ; else MAXTHREADS=1 ; fi

        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        fastp -q ${params.trim_qual} --json ${meta.library}.fastp.json --html ${meta.library}.fastp.html -i ${fastq1} -I ${fastq2} --stdout | \\
            bwa mem -p -t \$MAXTHREADS -5 -T0 -SP \$INDEX - | tee >(samtools view -bS -o ${meta.library}.bam) | \\
            pairtools parse --min-mapq ${params.parsemq} --walks-policy ${params.parse_walks_policy} --max-inter-align-gap ${params.parse_max_gap} --add-columns mapq --nproc-in 1 --nproc-out 1 -c ${chromsizes} | \\
            pairtools sort --nproc 1 --tmpdir tmp -o ${meta.library}.sorted.pairs.gz &&

        pairtools dedup --max-mismatch ${params.max_mismatch} --mark-dups --output-stats ${meta.library}.dedup.stats.txt ${meta.library}.sorted.pairs.gz | \\
            pairtools split --nproc-in 1 --nproc-out 1 --output-pairs ${meta.library}.pairs.gz --output-sam - | \\
            samtools view -bS -@ 1 | \\
            samtools sort -@ \$MAXTHREADS --reference ${reference_fasta} --write-index -o ${meta.library}.pairs.cram##idx##${meta.library}.pairs.cram.crai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
        END_VERSIONS
        """
}