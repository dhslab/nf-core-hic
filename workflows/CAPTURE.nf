/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE Baits Bed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Check bait bed parameter is provided
def checkPathParamList = [ params.baits_bed ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check bait bed path if it exist
if (params.baits_bed) { ch_baits_bed = file(params.baits_bed) } else { exit 1, 'Baits Bed file is not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAKE_MCOOL                   } from '../subworkflows/local/MAKE_MCOOL.nf'
include { capture_qc                   } from '../modules/local/capture_qc.nf'
include { make_chicago                 } from '../modules/local/make_chicago.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAPTURE {
    MAKE_MCOOL()
    
    ch_versions = MAKE_MCOOL.out.ch_versions

    //
    // CHANNEL: Channel operation to mix cram files and deduplications stats per individual library to use it as input for capture qc
    //
    MAKE_MCOOL.out.pairs_cram
    .join(MAKE_MCOOL.out.dedup_stats)
    .set { ch_cram_stats }
    
    //
    // MODULE: Do QC for capture quality
    //
    capture_qc (
        ch_cram_stats,
        file(params.fasta),
        ch_baits_bed
    )
    ch_versions = ch_versions.mix(capture_qc.out.versions)

    //
    // MODULE: make chicago-compatible bam file
    //
    make_chicago (
        MAKE_MCOOL.out.ch_merged_cram,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(make_chicago.out.versions)   

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    ) 
}