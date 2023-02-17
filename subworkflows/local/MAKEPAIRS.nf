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
include { INPUT_CHECK } from './input_check'

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
include { COMBINELANES                   } from '../../modules/local/COMBINELANES.nf'
include { FASTQ2PAIRS                    } from '../../modules/local/FASTQ2PAIRS.nf'
include { SORTBAM                        } from '../../modules/local/SORTBAM.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MAKEPAIRS {

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
    COMBINELANES (
        ch_fastq
    )

    //
    // MODULE: Convert Fastq to .pairs files
    //
    FASTQ2PAIRS (
        COMBINELANES.out.fastq,
        file(params.bwa_index + '.*'), // channel to stage all bwa index files
        file(params.fasta), // for genome fasta
        file(params.chromsizes)

    )
    ch_versions = ch_versions.mix(FASTQ2PAIRS.out.versions)

    //
    // MODULE: Sort mapping bam and convert it to indexed cram
    //
    SORTBAM (
        FASTQ2PAIRS.out.bam,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(SORTBAM.out.versions)

    emit:
        aligned_bam     = FASTQ2PAIRS.out.bam
        aligned_cram    = SORTBAM.out.cram_crai
        pairs           = FASTQ2PAIRS.out.pairs
        pairs_cram      = FASTQ2PAIRS.out.cram_crai
        ch_versions

}