// --- classify reads by Kraken2 taxid, up to (not including) consensus ----------
// Split out of subworkflows/mapping.nf: this is the Kraken2/Kraken2Ref half of what
// mapping.nf used to do end to end. Everything from GENERATE_CONSENSUS onward now lives
// in mapping.nf, which is shared with subworkflows/classifying_index.nf -- so a species
// found by either classifier goes through one consensus/Nextclade/subtyping/report path
// instead of two parallel ones.
//
// The interface mapping.nf consumes is deliberately identical for both classifiers:
//   sample_taxid_ch                tuple(meta, [read_1, read_2], reference_fasta)
//   sample_report_with_join_key_ch [join_key, report_meta]  (join_key == meta.id)
include {SORT_READS_BY_REF} from '../workflows/SORT_READS_BY_REF.nf'

workflow CLASSIFYING_KRAKEN2 {
    /*
    -----------------------------------------------------------------
    Classifies preprocessed reads against a Kraken2 database, sorts
    them per selected reference taxid (Kraken2Ref), and resolves each
    taxid's reference FASTA -- leaving reads ready for consensus
    generation in mapping.nf.
    -----------------------------------------------------------------
    */

    take:
        preprocessed_3tuple_ch // tuple (meta, read1, read2)

    main:
        // reconstruct the tuple(meta, [read1, read2]) shape SORT_READS_BY_REF expects
        preprocessed_3tuple_ch
            .map { meta, read1, read2 -> [meta, [read1, read2]] }
            .set { sort_reads_in_ch }

        SORT_READS_BY_REF(sort_reads_in_ch)

        // One row per (sample, selected reference taxid), keyed so mapping.nf can join it
        // against GENERATE_CONSENSUS's own per-consensus output. `id` is set on the row
        // itself as well as used as the key: mapping.nf's Nextclade input carries this
        // report row forward as the base of its meta, and downstream steps
        // (RUN_NEXTCLADE's own per-consensus JSON keying) expect meta.id to be present.
        SORT_READS_BY_REF.out.sample_pre_report_ch
            .map { row ->
                def join_key = "${row.sample_id}.${row.selected_taxid}".toString()
                [join_key, row + [id: join_key]]
            }
            .set { sample_report_with_join_key_ch }

        // --- rvi_integration_1: per-sample "species already identified by Kraken2" ---
        // Consumed by CLASSIFYING_INDEX to decide which of its own species calls are
        // genuinely new. Built straight off SORT_READS_BY_REF's raw per-sample
        // pre-report FILE (one element per sample, available as soon as THAT sample's
        // Kraken2/k2r pass finishes) rather than the exploded sample_pre_report_ch +
        // groupTuple() -- groupTuple() can't emit a group until its whole upstream
        // channel closes, which would mean waiting for every sample in the run, not
        // just this one.
        //
        // Match on virus_name, NOT ref_selected. Both are free-text names lifted from the
        // Kraken2 report by bin/k2r_report.py, but they sit at different ranks and only
        // virus_name shares a vocabulary with the sequence-index methods'
        // species_label/species fields:
        //   virus_name   <- kraken2ref's source_taxid, rank S      e.g. "Betacoronavirus pandemicum"
        //   ref_selected <- the selected reference, rank S1/S2/S3  e.g. "Severe acute
        //                                                          respiratory syndrome coronavirus 2"
        // Themisto2/Metagraph label the RVDB index with ICTV species binomials, so comparing
        // against ref_selected alone never matched: every species Kraken2 had already
        // found looked "new", and the feature called a redundant consensus for it rather
        // than a new-species one. Proven on the farm -- see INSTRUCT.md's
        // "new-species consensus: matched on the wrong column" section.
        //
        // virus_name is reliably species rank: kraken2ref decomposes species -> below-species
        // by construction, and its decomposed JSON carries the rank code in its own `source`
        // field ('S' for all 68 reference selections across the 10 samples run so far).
        //
        // ref_selected is still folded into the same set, as a widening rather than a
        // replacement: a candidate whose label happens to match a strain-level reference
        // Kraken2 already selected is also genuinely already covered, so suppressing it is
        // correct too. Adding it can only ever suppress a candidate, never invent one.
        SORT_READS_BY_REF.out.raw_sample_pre_report_ch
            .filter { it -> it.size() > 1 } // mirror SORT_READS_BY_REF's own empty-file guard
            .map { report_file ->
                def lines = report_file.readLines()
                def header = lines[0].split('\t')
                def sample_id_idx    = header.findIndexOf { String col -> col == 'sample_id' }
                def virus_name_idx   = header.findIndexOf { String col -> col == 'virus_name' }
                def ref_selected_idx = header.findIndexOf { String col -> col == 'ref_selected' }
                // Fail loudly rather than silently mis-matching: Groovy's row[-1] returns the
                // LAST field, so a renamed/removed column would quietly compare against
                // report_name instead of erroring.
                if (sample_id_idx < 0 || virus_name_idx < 0 || ref_selected_idx < 0) {
                    error("pre-report ${report_file} lacks one of the sample_id/virus_name/" +
                          "ref_selected columns (header: ${header}). bin/k2r_report.py's output " +
                          "format has changed -- update identified_species_ch in " +
                          "subworkflows/classifying_kraken2.nf.")
                }
                // split('\t') drops trailing empty fields, so a row is only usable if it
                // actually reaches the columns being read.
                def max_idx = [sample_id_idx, virus_name_idx, ref_selected_idx].max()
                def rows = lines[1..-1]
                    .collect { line -> line.split('\t') }
                    .findAll { row -> row.size() > max_idx }
                if (!rows) {
                    // Not reachable via bin/k2r_report.py today: a sample with nothing
                    // selected yields a column-less 1-byte file, already dropped by the
                    // filter above, and any real row carries all 11 columns. Fail loudly
                    // instead of letting rows[0] throw an opaque IndexOutOfBounds.
                    error("pre-report ${report_file} has a header but no parseable data rows " +
                          "(need >${max_idx} tab-separated fields). Check bin/k2r_report.py's output.")
                }
                def sample_id = rows[0][sample_id_idx]
                def species = rows
                    .collectMany { row -> [row[virus_name_idx], row[ref_selected_idx]] }
                    .collect { name -> name.trim().toLowerCase() }
                    .findAll { name -> name }
                    .unique()
                [sample_id, species]
            }
            .set { identified_species_ch }

    emit:
        sample_taxid_ch = SORT_READS_BY_REF.out.sample_taxid_ch // tuple (meta, [read_1, read_2], reference_fasta)
        sample_report_with_join_key_ch                          // [join_key, report_meta]
        identified_species_ch                                   // [sample_id, [normalized_species_name, ...]]
}
