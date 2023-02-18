
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKE_PAIRS             } from './MAKE_PAIRS.nf'
include { merge_crams           } from '../../modules/local/merge_crams.nf'
include { merge_pairs           } from '../../modules/local/merge_pairs.nf'
include { pairs2mcool           } from '../../modules/local/pairs2mcool.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow MAKE_MCOOL {
    MAKE_PAIRS()

    ch_versions = MAKE_PAIRS.out.ch_versions

    //
    // CHANNEL: Channel operation to mix all crams per sample in one channel
    //
    MAKE_PAIRS.out.pairs_cram
    .map { meta, cram, crai -> [[sample: meta.sample] , cram] }
    .groupTuple(by: 0)
    .set { ch_cram_per_sample }

    //
    // CHANNEL: Channel operation to mix all .pairs files per sample in one channel
    //
    MAKE_PAIRS.out.pairs
    .map { meta, pairs -> [[sample: meta.sample] , pairs] }
    .groupTuple(by: 0)
    .set { ch_pairs_per_sample }

    //
    // MODULE: Merge all crams per sample
    //
    merge_crams (
        ch_cram_per_sample,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(merge_crams.out.versions)

    //
    // MODULE: Merge all .pairs files per sample
    //
    merge_pairs (
        ch_pairs_per_sample
    )
    ch_versions = ch_versions.mix(merge_pairs.out.versions)

    //
    // MODULE: make mcool file from pairs file
    //
    pairs2mcool (
        merge_pairs.out.pairs,
        file(params.chromsizes)
    )
    ch_versions = ch_versions.mix(pairs2mcool.out.versions)

    emit:
        pairs_cram       = MAKE_PAIRS.out.pairs_cram
        dedup_stats      = MAKE_PAIRS.out.dedup_stats
        ch_merged_cram   = merge_crams.out.cram_crai
        ch_versions


}