
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKEPAIRS            } from '../subworkflows/local/MAKEPAIRS.nf'
include { MERGECRAMS           } from '../modules/local/MERGECRAMS.nf'
// include { MAKEBIGWIG           } from '../modules/local/MAKEBIGWIG.nf'
include { MERGEPAIRS           } from '../modules/local/MERGEPAIRS.nf'
include { PAIRS2MCOOL          } from '../modules/local/PAIRS2MCOOL.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow HIC {
    MAKEPAIRS()

    ch_versions = MAKEPAIRS.out.ch_versions

    //
    // CHANNEL: Channel operation to mix all crams per sample in one channel
    //
    MAKEPAIRS.out.pairs_cram
    .map { meta, cram, crai -> [[sample: meta.sample] , cram] }
    .groupTuple(by: 0)
    .set { ch_cram_per_sample }

    //
    // MODULE: Merge all crams per sample
    //
    MERGECRAMS (
        ch_cram_per_sample,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(MERGECRAMS.out.versions)

    //
    // MODULE: Make BigWig file from merged cram file
    //
    // MAKEBIGWIG (
    //     MERGECRAMS.out.cram_crai,
    //     file(params.fasta)
    // )
    // ch_versions = ch_versions.mix(MAKEBIGWIG.out.versions)

    //
    // CHANNEL: Channel operation to mix all .pairs files per sample in one channel
    //
    MAKEPAIRS.out.pairs
    .map { meta, pairs -> [[sample: meta.sample] , pairs] }
    .groupTuple(by: 0)
    .set { ch_pairs_per_sample }


    //
    // MODULE: Merge all .pairs files per sample
    //
    MERGEPAIRS (
        ch_pairs_per_sample
    )
    ch_versions = ch_versions.mix(MERGEPAIRS.out.versions)

    //
    // MODULE: make mcool file from pairs file
    //
    PAIRS2MCOOL (
        MERGEPAIRS.out.pairs,
        file(params.chromsizes)
    )
    ch_versions = ch_versions.mix(PAIRS2MCOOL.out.versions)

    emit: ch_pairs_per_sample
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
