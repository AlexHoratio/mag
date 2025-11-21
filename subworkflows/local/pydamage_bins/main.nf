/*
 * SHORTREAD_PREPROCESSING: Preprocessing and QC for short reads
 */

include { SUMMARISEPYDAMAGE } from '../../../modules/local/summarisepydamage/main'

workflow PYDAMAGE_BINS {
    take:
    ch_contig_pydamage_results
    ch_input_for_postbinning

    main:
    ch_versions = channel.empty()

    // Separate bins per sample into distinct entries
    ch_binned_contig_assignments = ch_input_for_postbinning.transpose()

    // Save the assembly meta for each bin in a separate channel for later re-joining
    ch_input_for_bins_metaonly = ch_binned_contig_assignments.map { meta, binfile -> [[bin_id: binfile.name], meta] }

    // Extract contig ID per bin
    ch_binned_contig_assignments_meta = ch_binned_contig_assignments
        .map { meta, binfile -> [meta + [bin_id: binfile.name], binfile] }
        .splitFasta(record: [header: true])
        .map { meta, record -> [[id: meta.id, contig_id: record.header.split(' ')[0]], [bin_id: meta.bin_id]] }
    // CONCOCT appears to be the only binner that does not remove everything after a space in the contig header
    // So we split it ourselves here to allow re-joining as with all others.

    // Split up the pydamage results on a per-contig to allow re-grouping
    // corresponding to the contigs in each bin
    ch_contig_pydamage_forreordering_results = ch_contig_pydamage_results
        .splitCsv(header: true)
        .map { meta, pydamage_stats -> [[id: meta.id, contig_id: pydamage_stats.reference], meta, pydamage_stats] }

    // Combine the split pydamage results with the contig IDs plus bin meta,
    // so we can group the pyDamage contig rows into their respective bins
    ch_reordered_pydamage_stats = ch_contig_pydamage_forreordering_results
        .combine(ch_binned_contig_assignments_meta, by: 0)
        .map { _coremeta, full_meta, pydamage_stats, bin_meta -> [bin_meta, full_meta, pydamage_stats] }
        .groupTuple(by: [0, 1])

    // Convert contents of the reordered contigs to a CSV file
    // and re-attach the original sample/assembly meta based on bin_id
    // This also saves the 'reordered' pydamage results to allow inspection
    // with the final summarised stats
    ch_pydamage_to_bins = ch_reordered_pydamage_stats
        .map { bin_id, _full_meta, data ->
            def header = data[0].keySet().join(',')
            def rows = data.collect { row -> row.values().join(',') }.join('\n')
            def content = header + '\n' + rows
            [bin_id, content]
        }
        .collectFile(
            [storeDir: "${params.outdir}/GenomeBinning/QC/pydamage/analyze_bins/"]
        ) { meta, content ->
            ["${meta.bin_id}_pydamagebins.csv", content]
        }
        .map { file ->
            def meta = [:]
            meta.bin_id = file.getName() - '_pydamagebins.csv'
            [meta, file]
        }
        .join(ch_input_for_bins_metaonly)
        .map { bin_id, file, meta ->
            def meta_new = meta + bin_id
            [meta_new, file]
        }

    // Generate the per-bin median stats summary
    SUMMARISEPYDAMAGE(ch_pydamage_to_bins)
    ch_versions = ch_versions.mix(SUMMARISEPYDAMAGE.out.versions)

    // Stick per bin-summary into a single summary CSV
    // that can be used to bind to the bin table later on
    ch_aggregate_summaries = SUMMARISEPYDAMAGE.out.summary_tsv
        .map { _meta, file -> [file] }
        .splitCsv(header: true, sep: '\t')
        .map { content ->
            content[0].id = content[0].id - '_pydamagebins'
            [[id: 'all'], content.flatten()]
        }
        .groupTuple()
        .map { _meta, contents ->
            // Get keys of first row to act as header
            def header = contents[0][0].keySet().join('\t')
            // Get values of remaining rows for cells
            def rows = contents.collect { row -> row[0].values().join('\t') }.join('\n')
            header + '\n' + rows
        }
        .collectFile(
            name: "pydamage_bin_summary.tsv",
            newLine: true,
            sort: false,
            storeDir: "${params.outdir}/GenomeBinning/QC/",
        )

    emit:
    tsv      = ch_aggregate_summaries
    versions = ch_versions
}
