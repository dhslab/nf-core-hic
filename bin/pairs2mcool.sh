#!/bin/bash

GENOME=${params.genome}
CHROMSIZES=${chromsizes}

# this could be variable, but we should go down to at least 5000, maybe 1000.
RESOLUTIONS=${params.resolutions}
# this is the minimum of the above
MINRES=${params.min_res}

SAMPLE=${meta.sample}
PAIRS=${pairs}
MAPQ=1

# could be increased to 3
THREADS=${task.cpus}

bgzip -cd -@ $THREADS $PAIRS | pairtools select "(mapq1>=$MAPQ) and (mapq2>=$MAPQ)" | \
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly $GENOME $CHROMSIZES:$MINRES - $SAMPLE.mapq_$MAPQ.cool && \
    cooler zoomify --nproc $THREADS --out $SAMPLE.mapq_$MAPQ.mcool --resolutions $RESOLUTIONS --balance $SAMPLE.mapq_$MAPQ.cool


        # for MAPQ in "1 30"
        #     do
        #         bgzip -cd -@ ${task.cpus} ${pairs} | pairtools select "(mapq1>=\${MAPQ}) and (mapq2>=\${MAPQ})" | \\
        #             cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly ${params.genome} ${chromsizes}:${params.min_res} - ${meta.sample}.mapq_\${MAPQ}.cool && \\
        #             cooler zoomify --nproc ${task.cpus} --out ${meta.sample}.mapq_\${MAPQ}.mcool --resolutions ${params.resolutions} --balance ${meta.sample}.mapq_\${MAPQ}.cool
        #     done