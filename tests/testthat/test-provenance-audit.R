toy_provenance_raw_sequences <- function() {
	tibble::tibble(
		strain_type = c("A", "A", "A", "A"),
		strain_subtype = c("H1N1", "H1N1", "H3N2", "H3N2"),
		strain_name = c(
			"A/H1N1/Full/1/2001",
			"A/H1N1/Short/2/2002",
			"A/H3N2/Unused/3/2003",
			"A/H3N2/Included/4/2004"
		),
		short_name = c("one", "two", "unused", "four"),
		protein_sequence_source = c("GenBank", "Gisaid", "GenBank", ""),
		protein_sequence_accession_number = c("AAA1", "EPI2", "BBB3", ""),
		full_length = c(TRUE, FALSE, TRUE, TRUE),
		protein_sequence = c("AAAA", "AA", "CCCC", "DDDD")
	)
}

toy_provenance_cartography_antigens <- function() {
	list(
		h1 = c("H1N1-Full-2001", "H1N1-Short-2002"),
		h3 = c("H3N2-Included-2004")
	)
}

test_that("strain provenance records label inclusion and source completeness", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	alignment_completeness <- tibble::tibble(
		subtype = c("h1", "h1", "h3"),
		short_name = c("one", "two", "four"),
		protein_non_gap_length = c(4L, 2L, 4L),
		alignment_full_length = c(TRUE, FALSE, TRUE)
	)

	out <- make_strain_provenance_records(
		toy_provenance_raw_sequences(),
		toy_provenance_cartography_antigens(),
		alignment_completeness
	)

	expect_equal(nrow(out), 4)
	expect_true(all(c(
		"inclusion_status",
		"exclusion_reason",
		"source_status",
		"accession_status",
		"protein_sequence_status",
		"alignment_full_length"
	) %in% names(out)))
	expect_equal(
		out$inclusion_status,
		c("included_in_analysis", "included_in_analysis", "excluded_from_analysis", "included_in_analysis")
	)
	expect_equal(out$exclusion_reason[out$short_name == "unused"], "absent_from_cartography_maps")
	expect_equal(out$source_status[out$short_name == "four"], "missing_source")
	expect_equal(out$accession_status[out$short_name == "four"], "missing_accession")
	expect_equal(out$alignment_full_length[out$short_name == "two"], FALSE)
})

test_that("pair count audit checks unique off-diagonal pair counts by subtype and method", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	clean_sequences <- tibble::tibble(
		subtype = c("h1", "h1", "h1", "h3", "h3"),
		short_name = c("a", "b", "c", "d", "e")
	)
	distance_table <- tibble::tibble(
		subtype = c(rep("h1", 6), "h3"),
		method = c(rep(c("year", "cart"), each = 3), "year"),
		Var1 = c("b", "c", "c", "b", "c", "c", "e"),
		Var2 = c("a", "a", "b", "a", "a", "b", "d"),
		d = seq_len(7)
	)

	out <- make_pair_count_audit(distance_table, clean_sequences)

	expect_equal(
		names(out),
		c("subtype", "method", "included_strains", "expected_unique_pairs", "observed_unique_pairs", "pair_count_complete")
	)
	expect_equal(out$expected_unique_pairs[out$subtype == "h1"], c(3, 3))
	expect_true(all(out$pair_count_complete[out$subtype == "h1"]))
	expect_equal(out$expected_unique_pairs[out$subtype == "h3"], 1)
	expect_true(out$pair_count_complete[out$subtype == "h3"])
})

test_that("strain flow summary reports raw, included, excluded, incomplete, and pair counts", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	provenance <- make_strain_provenance_records(
		toy_provenance_raw_sequences(),
		toy_provenance_cartography_antigens(),
		tibble::tibble(
			subtype = c("h1", "h1", "h3"),
			short_name = c("one", "two", "four"),
			protein_non_gap_length = c(4L, 2L, 4L),
			alignment_full_length = c(TRUE, FALSE, TRUE)
		)
	)
	pair_counts <- tibble::tibble(
		subtype = c("h1", "h3"),
		method = c("year", "year"),
		included_strains = c(2L, 1L),
		expected_unique_pairs = c(1L, 0L),
		observed_unique_pairs = c(1L, 0L),
		pair_count_complete = c(TRUE, TRUE)
	)

	out <- make_strain_flow_summary(provenance, pair_counts)

	expect_true(all(c("stage", "category", "subtype", "n", "detail") %in% names(out)))
	expect_equal(out$n[out$category == "raw_sequence_rows" & is.na(out$subtype)], 4)
	expect_equal(out$n[out$category == "analysis_sequence_rows" & is.na(out$subtype)], 3)
	expect_equal(out$n[out$category == "excluded_absent_from_cartography_maps" & is.na(out$subtype)], 1)
	expect_match(out$detail[out$category == "analysis_non_full_length_rows"], "two", fixed = TRUE)
	expect_equal(out$n[out$category == "unique_pairwise_comparisons" & out$subtype == "h1"], 1)
})
