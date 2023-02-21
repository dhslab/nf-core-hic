# hic
## Pipeline for Hi-C/Capture-C data analysis


[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core-hic** is a bioinformatics best-practice analysis pipeline for Hi-C/Capture-C data analysis. This pipeline is optimal for large scale analysis in High Performance Computing Clusters (HPCs) and Cloud Computing environments (eg. AWS)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

## Workflow Summary
### 1) Hi-C Workflow (default)

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. **fastq2pair (per library)**:
    1. **Preprocessing** ([`fastp`](https://github.com/OpenGene/fastp)) >> `library.html` & `library.json`
    1. **Alignment** ([`bwa`](https://bio-bwa.sourceforge.net/bwa.shtml)) >> `library.cram`/`library.cram.crai`
    1. **Extract ligation junctions** ([`pairtools`](https://pairtools.readthedocs.io/en/latest/))
    1. **Remove PCR/optical duplicates** ([`pairtools`](https://pairtools.readthedocs.io/en/latest/)) >> `library.pairs.gz` & `library.dedup.stats.txt`
    1. **Make pairs cram file** ([`pairtools`](https://pairtools.readthedocs.io/en/latest/) & [`samtools`](http://www.htslib.org/doc/samtools.html))  >> `library.pairs.cram`/`library.pairs.cram.crai`
    - These steps are based on this [Dovetail tutorial](https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html). Check the link for more details
1. **Merge all** `library.pairs.gz` & `library.pairs.cram` for libraries per individual sample >> `sample.pairs.gz` & `sample.pairs.cram`
1. **Make `.mcool` file** ([`cooler`](https://cooler.readthedocs.io/en/latest/quickstart.html)) >> `sample.mcool`

### 2) Capture-C Workflow
**- Initial steps Similar to Hi-C Workflow (steps 1-3)** \
4\. **QC for Capture** (Baits regions coverage) \
5\. **Make bam** file compatible with [CHiCAGO algorithm](https://bioconductor.org/packages/release/bioc/html/Chicago.html) ([`samtools`](http://www.htslib.org/doc/samtools.html))

### 3) QC Workflow
- This workflow is intended to check **library Complexity** from shallow-depth sequencing for QC before doing deep sequencing. it is based on this [Dovetail tutorial](https://micro-c.readthedocs.io/en/latest/library_qc.html#library-complexity).
1. **fastq2pair (per library)**: Same steps as in HiC and Capture-C workflows.
2. **Estimate library complexity** ([`preseq`](https://github.com/smithlabcode/preseq)) >> `sample.preseq.txt`. For interpretation of this results refer to [Dovetail tutorial](https://micro-c.readthedocs.io/en/latest/library_qc.html#library-complexity)


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(**this pipeline can NOT be run with conda**))_. This requirement is not needed for running the pipeline in WashU RIS cluster. This pipeline is also successfully tested using Amazon Cloud Computing (AWS). For details on how to run nextflow pipelines in AWS refer to [nextflow documentation](https://www.nextflow.io/docs/latest/aws.html) and to [this excellent tutorial](https://staphb.org/resources/2020-04-29-nextflow_batch.html).


3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run dhslab/nf-core-hic -profile test,YOURPROFILE(S) --outdir <OUTDIR>
   ```


4. Start running your own analysis!
    1. Hi-C workflow (default)
   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run dhslab/nf-core-hic -r dev \
   -profile YOURPROFILE(S) \
   --input <SAMPLESHEET> \
   --fasta <FASTA> \
   --bwa_index <INDEX_PREFIX> \
   --chromsizes <CHROMSIZES> \
   --genome <GENOME_NAME> \
   --outdir <OUTDIR> 
   ```

    2. Capture-C workflow
   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run dhslab/nf-core-hic -r dev \
   -entry capture \
   -profile YOURPROFILE(S) \
   --input <SAMPLESHEET> \
   --fasta <FASTA> \
   --bwa_index <INDEX_PREFIX> \
   --chromsizes <CHROMSIZES> \
   --genome <GENOME_NAME> \
   --baits_bed <BAITS_BED> \
   --outdir <OUTDIR>
   ```
    3. QC workflow
   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run dhslab/nf-core-hic -r dev \
   -entry qc
   -profile YOURPROFILE(S) \
   --input <SAMPLESHEET> \
   --fasta <FASTA> \
   --bwa_index <INDEX_PREFIX> \
   --chromsizes <CHROMSIZES> \
   --genome <GENOME_NAME> \
   --outdir <OUTDIR> 
   ```

## Usage
### Required parameters:
1. **Input** `samplesheet.cvs` which provides directory paths for `fastq1`, `fastq2` raw reads and their metadata (`id`, `sample`, `library`, `flowcell`). this can be provided either in a configuration file or as `--input path/to/samplesheet.cvs` command line parameter. Example sheet located in `assets/samplesheet.csv`.
2. **Genome fasta**, either in a configuration file or as `--fasta path/to/genome.fasta` command line parameter.
3. **BWA index**, either in a configuration file or as `--bwa_index path/to/bwa_index/with_prefix` command line parameter. It is important to provide the full path including **index prefix**.
4. **Chromosome sizes file**, either in a configuration file or as `--chromsizes path/to/chromsizes` command line parameter.
5. **Genome name (eg. hg38)**, either in a configuration file or as `--fasta path/to/genome.fasta` command line parameter.
6. **Capture-C Baits bed file (Only in Capture-C workflow)** either in a configuration file or as `--baits_bed path/to/baits_bed` command line parameter.

### Tools specific parameters:
The following parameters are set to the shown default values, but should be modified when required in command line, or in user-provided config files:
#### Preprocessing options for fastp


| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `trim_qual` | fastp `-q` option for the quality value that a base is qualified | `integer` | 15 |

#### pairtools options

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `parsemq` | pairtools parse `--min-mapq` option for the minimal MAPQ score to consider a read as uniquely mapped | `integer` | 1 |
| `parse_walks_policy` | pairtools parse `--walks-policy` option. See pairtools documentation for details | `string` | 5unique |
| `parse_max_gap` | pairtools parse `--max-inter-align-gap` option. See pairtools documentation for details | `integer` | 30 |
| `max_mismatch` | pairtools dedup `--max-mismatch` option. Pairs with both sides mapped within this distance (bp) from each other are considered duplicates | `integer` | 1 |

#### mcool file options

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `resolutions` | cooler zoomify `--resolutions` option: Comma-separated list of target resolutions | `string` | 1000000,500000,250000,100000,50000,20000,10000,5000 |
| `min_res` | Minimum resolution for the mcool file (from the resolutions list provided) | `integer` | 5000 |
| `mcool_mapq_threshold` | Minimum resolution for the mcool file (from the resolutions list provided) | `string` | 1 30 |

#### Capture-C options

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `baits_bed` | Bed file for regions targeted  by Capture baits. Required only for Capture-C workflow | `string` | None |


## Example for running a test on WashU RIS HPC:

```bash
# Run default Hi-C workflow
NXF_HOME=/scratch1/fs1/dspencer/.nextflow \
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME" \
bsub -g /dspencer/nextflow -G compute-dspencer -q dspencer \
-e nextflow_launcher.err -o nextflow_launcher.log -We 2:00 -n 2 -M 12GB \
-R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" \
-a "docker(ghcr.io/dhslab/docker-nextflow)" \
nextflow run dhslab/nf-core-hic \
-r dev \
-profile dhslab,ris \
-c /storage1/fs1/dspencer/Active/spencerlab/mohamed/github/nf-core-hic/conf/test_dhs.config \
--outdir results
```
**Notice the profiles which are used here:**
1. **`dhslab`** -> to set **lab-specific** cluster/cloud configuration (by importing `conf/dhslab.config`, see `nextflow.config`)
2. **`ris`** -> to set **general** configuration for RIS LSF cluster
- any number of profiles/config-files can be used. Just consider how configuration priorities are set in nextflow as documented [here](https://www.nextflow.io/docs/latest/config.html)

### **Directory tree for test run output (default workflow):**

```
.
├── pipeline_info
│   ├── execution_report_2023-02-17_23-57-34.html
│   ├── execution_timeline_2023-02-17_23-57-34.html
│   ├── execution_trace_2023-02-17_23-57-34.txt
│   ├── pipeline_dag_2023-02-17_23-57-34.html
│   ├── samplesheet.valid.csv
│   └── software_versions.yml
└── samples
    └── TEST
        ├── fastq2pairs
        │   ├── TESTA
        │   │   ├── TESTA.cram
        │   │   ├── TESTA.cram.crai
        │   │   ├── TESTA.dedup.stats.txt
        │   │   ├── TESTA.fastp.html
        │   │   ├── TESTA.fastp.json
        │   │   ├── TESTA.pairs.cram
        │   │   ├── TESTA.pairs.cram.crai
        │   │   └── TESTA.pairs.gz
        │   ├── TESTB
        │   │   ├── TESTB.cram
        │   │   ├── TESTB.cram.crai
        │   │   ├── TESTB.dedup.stats.txt
        │   │   ├── TESTB.fastp.html
        │   │   ├── TESTB.fastp.json
        │   │   ├── TESTB.pairs.cram
        │   │   ├── TESTB.pairs.cram.crai
        │   │   └── TESTB.pairs.gz
        │   ├── TESTC
        │   │   ├── TESTC.cram
        │   │   ├── TESTC.cram.crai
        │   │   ├── TESTC.dedup.stats.txt
        │   │   ├── TESTC.fastp.html
        │   │   ├── TESTC.fastp.json
        │   │   ├── TESTC.pairs.cram
        │   │   ├── TESTC.pairs.cram.crai
        │   │   └── TESTC.pairs.gz
        │   └── merged
        │       ├── TEST.pairs.cram
        │       ├── TEST.pairs.cram.crai
        │       └── TEST.pairs.gz
        └── mcool
            ├── TEST.mapq_1.mcool
            └── TEST.mapq_30.mcool


```


### Notes:
- The pipeline is developed and optimized to be run in WashU RIS (LSF) HPC, but could be deployed in any [`HPC environment supported by Nextflow`](https://www.nextflow.io/docs/latest/executor.html).
- The pipeline does NOT support conda.
- The Test workflow can be run on personal computer, but is not advised. It is recommended to do the testing in environment with at least 16 GB memory. If the test workflow failed (especially at **fastq2pair** step ), try re-run with more allocated resources. Such errors are likely because of broken pipes due to maxed-out memory. The pipeline is designed with several pipe steps to avoid making large intermediate files.
