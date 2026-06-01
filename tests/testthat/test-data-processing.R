toy_raw_sequences <- function() {
	tibble::tibble(
		strain_type = c("A", "A", "A", "A"),
		strain_subtype = c("H1N1", "H1N1", "H3N2", "H1N1"),
		strain_name = c(
			"A/H1N1/Full/1/2001",
			"A/H1N1/Full/2/2002",
			"A/H3N2/Full/3/2003",
			"A/H1N1/Unused/4/2004"
		),
		short_name = c("one", "two", "three", "unused"),
		protein_sequence_source = c("GenBank", "Gisaid", "GenBank", "GenBank"),
		protein_sequence_accession_number = c("AAA1", "EPI2", "BBB3", "CCC4"),
		full_length = c(TRUE, FALSE, TRUE, TRUE),
		protein_sequence = c("AA AA", "BBBB", "CCCC", "DDDD")
	)
}

toy_cartography_antigens <- function() {
	list(
		h1 = c("H1N1-Full-2001", "H1N1-Full-2002"),
		h3 = c("H3N2-Full-2003")
	)
}

test_that("clean_sequence_data keeps protein-only sequence records in cartography maps", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	clean <- clean_sequence_data(toy_raw_sequences(), toy_cartography_antigens())

	expect_equal(clean$short_name, c("one", "two", "three"))
	expect_equal(clean$subtype, c("h1", "h1", "h3"))
	expect_equal(clean$protein_sequence, c("aaaa", "bbbb", "cccc"))
	expect_equal(clean$protein_length, c(4L, 4L, 4L))
	expect_equal(clean$full_length, c(TRUE, FALSE, TRUE))
	expect_equal(clean$protein_sequence_source, c("GenBank", "Gisaid", "GenBank"))
	expect_equal(clean$protein_sequence_accession_number, c("AAA1", "EPI2", "BBB3"))
	expect_equal(clean$analysis_name, c("H1N1-Full-2001", "H1N1-Full-2002", "H3N2-Full-2003"))
})

test_that("clean_sequence_data fails on missing columns and duplicate short names", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	raw <- toy_raw_sequences()

	expect_error(
		clean_sequence_data(dplyr::select(raw, -"protein_sequence"), toy_cartography_antigens()),
		"raw sequence CSV is missing required column\\(s\\): protein_sequence"
	)

	duplicate_raw <- dplyr::bind_rows(raw, raw[1, ])
	expect_error(
		clean_sequence_data(duplicate_raw, toy_cartography_antigens()),
		"raw sequence short names contains duplicate value\\(s\\)"
	)
})

test_that("cartography name resolution handles small stored-map spelling differences", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	raw <- tibble::tibble(
		strain_type = "A",
		strain_subtype = "H3N2",
		strain_name = "A/H3N2/Shandong/9/1993",
		short_name = "Shan/93",
		protein_sequence_source = "GenBank",
		protein_sequence_accession_number = "AAA1",
		full_length = TRUE,
		protein_sequence = "AAAA"
	)
	cartography_antigens <- list(h1 = character(), h3 = "H3N2-Shangdong-1993")

	clean <- clean_sequence_data(raw, cartography_antigens)

	expect_equal(clean$analysis_name, "H3N2-Shangdong-1993")
	expect_equal(clean$cartography_match_status, "fuzzy_name_resolution")
})

test_that("sequence-cleaning audit reports row counts, cartography exclusions, and full-length flags", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	raw <- toy_raw_sequences()
	cartography_antigens <- toy_cartography_antigens()
	clean <- clean_sequence_data(raw, cartography_antigens)
	audit <- sequence_cleaning_audit(raw, cartography_antigens, clean)

	expect_true(all(c("check", "value", "detail") %in% names(audit)))
	expect_equal(audit$value[audit$check == "raw_sequence_rows"], 4)
	expect_equal(audit$value[audit$check == "analysis_sequence_rows"], 3)
	expect_equal(audit$value[audit$check == "analysis_sequence_rows_h1"], 2)
	expect_equal(audit$value[audit$check == "source_non_full_length_rows"], 1)
	expect_match(audit$detail[audit$check == "source_non_full_length_rows"], "two", fixed = TRUE)
	expect_equal(audit$value[audit$check == "sequence_rows_not_in_cartography_maps"], 1)
	expect_match(audit$detail[audit$check == "sequence_rows_not_in_cartography_maps"], "unused", fixed = TRUE)
})

test_that("split_sequences_by_subtype returns subtype groups", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	clean <- clean_sequence_data(toy_raw_sequences(), toy_cartography_antigens())
	split <- split_sequences_by_subtype(clean)

	expect_equal(sort(names(split)), c("h1", "h3"))
	expect_equal(nrow(split$h1), 2)
	expect_equal(nrow(split$h3), 1)
})
