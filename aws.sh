# Run default Hic in RIS
NXF_HOME=/scratch1/fs1/dspencer/mohamed/.nextflow \
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME" \
bsub -g /dspencer/nextflow -G compute-dspencer -q dspencer \
-e nextflow_launcher.err -o nextflow_launcher.log -We 2:00 -n 2 -M 12GB \
-R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" \
-a "docker(ghcr.io/dhslab/docker-nextflow)" \
nextflow run dhslab/nf-core-hic \
-r dev -latest \
-profile dhslab,ris \
-c /storage1/fs1/dspencer/Active/spencerlab/mohamed/github/nf-core-hic/conf/test_dhs.config \
--outdir results


# Run default Hic in AWS (using igenomes hosted on aws)
nextflow run dhslab/nf-core-hic \
-r dev -latest \
-profile dhslab,aws \
-c /storage1/fs1/dspencer/Active/spencerlab/mohamed/github/nf-core-hic/conf/test.config \
--outdir results \
-bucket-dir s3://dhs-lab-data/nxf_tmp/

# Run default Hic in hybrid mode (using igenomes hosted on aws)
nextflow run dhslab/nf-core-hic \
-r dev -latest \
-profile dhslab,hybrid \
-c /storage1/fs1/dspencer/Active/spencerlab/mohamed/github/nf-core-hic/conf/test.config \
--outdir results \
-bucket-dir s3://dhs-lab-data/nxf_tmp/