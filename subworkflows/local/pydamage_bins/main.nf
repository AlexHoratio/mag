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

    // Get sample ID for each contig (record header)
    ch_binned_contig_assignments = ch_input_for_postbinning.transpose()

    // Get meta only version of input bins for later re-binding
    ch_input_for_bins_metaonly = ch_binned_contig_assignments.map { meta, binfile -> [[bin_id: binfile.name], meta] }

    // Extract contig ID per bin
    // TODO: CONCOCT adds additional information in header that needs to be cleaned to merge
    ch_binned_contig_assignments_meta = ch_binned_contig_assignments
        .map { meta, binfile -> [meta + [bin_id: binfile.name], binfile] }
        .splitFasta(record: [header: true])
        .map { meta, record -> [[id: meta.id, contig_id: record.header], [bin_id: meta.bin_id]] }

    // Make useable sample id/contig name meta map for merging
    ch_contig_pydamage_forreordering_results = ch_contig_pydamage_results
        .splitCsv(header: true)
        .map { meta, pydamage_stats -> [[id: meta.id, contig_id: pydamage_stats.reference], meta, pydamage_stats] }

    //Merge sample ID to pydamage results so we can group by contig AND node name
    // We need to `combine` as a contig from one a assembly maybe present in multiple bins
    ch_reordered_pydamage_stats = ch_contig_pydamage_forreordering_results
        .dump(tag: 'pre-combine')
        .combine(ch_binned_contig_assignments_meta, by: 0)
        .dump(tag: 'post-combine')
        .map { _coremeta, bin_meta, _full_meta, pydamage_stats -> [bin_meta, pydamage_stats] }
        .dump(tag: 'pregroup')
        .groupTuple(by: 0)
        .dump(tag: 'ch_reordered_pydamage_stats')

    // Convert contents of the reordered contigs to a CSV file, and re-attach meta based on bin_id (i.e., from bin file name)
    ch_pydamage_to_bins = ch_reordered_pydamage_stats
        .map { bin_id, data ->
            // Process your data to create CSV content TODO THIS IS BROKEN IT ONLY HAS A SINGLE COLUMN
            def header = data[0].keySet().join(',')
            def rows = data.collect { row -> row.values().join(',') }.join('\n')
            def content = header + '\n' + rows
            [bin_id, content]
        }
        .dump(tag: 'precollect')
        .collectFile(
            [storeDir: "${params.outdir}/GenomeBinning/QC/pydamage/analyze_bins/"]
        ) { meta, content ->
            ["${meta.bin_id}_pydamagebins.csv", content]
        }
        .dump(tag: 'postcollect')
        .map { file ->
            def meta = [:]
            meta.bin_id = file.getName() - '_pydamagebins.csv'
            [meta, file]
        }
        .dump(tag: 'prejoin')
        .join(ch_input_for_bins_metaonly)
        .dump(tag: 'postjoin')
        .map { bin_id, file, meta ->
            def meta_new = meta + bin_id
            [meta_new, file]
        }
        .dump(tag: 'ch_pydamage_to_bins')

    // Generate the per bin-summary
    SUMMARISEPYDAMAGE(ch_pydamage_to_bins)
    ch_versions = ch_versions.mix(SUMMARISEPYDAMAGE.out.versions)

    // Stick per bin-summary into a single summary CSV
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
