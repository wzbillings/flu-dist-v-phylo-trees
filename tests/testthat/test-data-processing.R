toy_raw_sequences <- function() {
	tibble::tibble(
		strain_name = c(" Full One ", "full two", "full three"),
		protein_sequence = c("AA AA", "BBBB", "CCCC"),
		nucleotide_sequence_type = c("rna", "rna", "rna"),
		nucleotide_sequence_source = c("source", "source", "source"),
		nucleotide_sequence = c("UU UU", "TTTT", "CCCC")
	)
}

toy_virus_metadata <- function() {
	tibble::tibble(
		subtype = c("h3", "h1", "h1", "h1"),
		analysis_name = c("A/Test/3/2003", "A/Test/1/2001", "A/Test/2/2002", "A/Test/4/2004"),
		genbank_strain_name = c("full three", "full one", "full two", "unused"),
		short_name = c("three", "one", "MI/15", "unused"),
		factor_order = c(3L, 1L, 2L, 4L),
		vaccine_strain = c(FALSE, TRUE, FALSE, FALSE)
	) |>
		dplyr::mutate(metadata_join_key = normalize_join_key(.data$genbank_strain_name))
}

test_that("clean_sequence_data joins metadata, orders records, and normalizes sequences", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	clean <- clean_sequence_data(toy_raw_sequences(), toy_virus_metadata())

	expect_equal(clean$short_name, c("one", "MI/15", "three"))
	expect_equal(clean$subtype, c("h1", "h1", "h3"))
	expect_equal(clean$protein_sequence, c("aaaa", "bbbb", "cccc"))
	expect_equal(clean$nucleotide_sequence, c("uuuu", "tttt", "cccc"))
	expect_equal(clean$protein_length, c(4L, 4L, 4L))
	expect_equal(clean$nucleotide_length, c(4L, 4L, 4L))
	expect_equal(clean$full_length, c(TRUE, FALSE, TRUE))
})

test_that("clean_sequence_data fails on missing metadata and duplicate join keys", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	raw <- toy_raw_sequences()
	metadata <- toy_virus_metadata()

	expect_error(
		clean_sequence_data(raw[-1, ], metadata[-3, ]),
		"Sequence workbook strain\\(s\\) missing from virus metadata: full two"
	)

	duplicate_metadata <- dplyr::bind_rows(metadata, metadata[1, ])
	expect_error(
		clean_sequence_data(raw, duplicate_metadata),
		"virus metadata strain names contains duplicate value\\(s\\)"
	)
})

test_that("sequence-cleaning audit reports row counts and excluded metadata", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	raw <- toy_raw_sequences()
	metadata <- toy_virus_metadata()
	clean <- clean_sequence_data(raw, metadata)
	audit <- sequence_cleaning_audit(raw, metadata, clean)

	expect_true(all(c("check", "value", "detail") %in% names(audit)))
	expect_equal(audit$value[audit$check == "raw_sequence_rows"], 3)
	expect_equal(audit$value[audit$check == "clean_sequence_rows"], 3)
	expect_equal(audit$value[audit$check == "flagged_not_full_length"], 1)
	expect_match(audit$detail[audit$check == "flagged_not_full_length"], "MI/15", fixed = TRUE)
	expect_match(audit$detail[audit$check == "metadata_rows_not_in_sequence_workbook"], "unused", fixed = TRUE)
})

test_that("split_sequences_by_subtype returns subtype groups", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	clean <- clean_sequence_data(toy_raw_sequences(), toy_virus_metadata())
	split <- split_sequences_by_subtype(clean)

	expect_equal(sort(names(split)), c("h1", "h3"))
	expect_equal(nrow(split$h1), 2)
	expect_equal(nrow(split$h3), 1)
})
