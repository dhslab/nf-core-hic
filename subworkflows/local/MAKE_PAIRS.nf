/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHic.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.bwa_index ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from './INPUT_CHECK'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
// include { FASTQC                      } from '../../modules/nf-core/fastqc/main'
// include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { combine_lanes                   } from '../../modules/local/combine_lanes.nf'
include { fastq2pairs                    } from '../../modules/local/fastq2pairs.nf'
include { sort_bam                        } from '../../modules/local/sort_bam.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAKE_PAIRS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // CHANNEL: Channel operation to merge reads of similar libraries and different lanes/flowcells together
    //
    // Group Read 1 in one channel
    INPUT_CHECK.out.reads
    .map { meta, fastq -> [[sample: meta.sample, library: meta.library] , fastq[0]] }
    .groupTuple(by: 0)
    .set { ch_read_1 }
    // Group Read 2 in one channel
    INPUT_CHECK.out.reads
    .map { meta, fastq -> [[sample: meta.sample, library: meta.library] , fastq[1]] }
    .groupTuple(by: 0)
    .set { ch_read_2 }
    // Merge them
    ch_read_1
    .join(ch_read_2)
    .set { ch_fastq }

    //
    // MODULE: Combine fastqs for same library and different lanes
    //
    combine_lanes (
        ch_fastq
    )

    //
    // MODULE: Convert Fastq to .pairs files
    //
    fastq2pairs (
        combine_lanes.out.fastq,
        file(params.bwa_index + '.*'), // channel to stage all bwa index files
        file(params.fasta), // for genome fasta
        file(params.chromsizes)

    )
    ch_versions = ch_versions.mix(fastq2pairs.out.versions)

    //
    // MODULE: Sort mapping bam and convert it to indexed cram
    //
    sort_bam (
        fastq2pairs.out.bam,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(sort_bam.out.versions)

    emit:
        aligned_bam     = fastq2pairs.out.bam
        pairs           = fastq2pairs.out.pairs
        pairs_cram      = fastq2pairs.out.cram_crai
        dedup_stats      = fastq2pairs.out.stats
        ch_versions

}