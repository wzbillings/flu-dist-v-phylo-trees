test_that("alignment completeness is computed from aligned protein lengths without hard-coded strain names", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	alignment_result <- list(
		subtype = "h3",
		aligned_sequences = tibble::tibble(
			short_name = c("full one", "short one", "full two"),
			pro_aligned = c("AAAA--", "AA----", "AAAA-A"),
			source_full_length = c(TRUE, FALSE, TRUE)
		)
	)

	out <- alignment_completeness_audit(alignment_result)

	expect_equal(out$short_name, c("full one", "short one", "full two"))
	expect_equal(out$protein_non_gap_length, c(4L, 2L, 5L))
	expect_equal(out$reference_full_length_min, c(4L, 4L, 4L))
	expect_equal(out$alignment_full_length, c(TRUE, FALSE, TRUE))
	expect_equal(out$full_length_consistent, c(TRUE, TRUE, TRUE))
})

test_that("alignment_audit reports protein-only widths and non-full-length details", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	alignment_result <- list(
		subtype = "h3",
		aligned_sequences = tibble::tibble(
			short_name = c("MI/85", "HK/68"),
			pro_aligned = c("AA----", "AAAAAA"),
			source_full_length = c(FALSE, TRUE)
		)
	)

	audit <- alignment_audit(alignment_result)

	expect_equal(names(audit), c(
		"subtype",
		"n_sequences",
		"protein_alignment_width",
		"source_non_full_length",
		"alignment_non_full_length",
		"non_full_length_detail"
	))
	expect_equal(audit$n_sequences, 2)
	expect_equal(audit$protein_alignment_width, 6)
	expect_equal(audit$source_non_full_length, 1)
	expect_equal(audit$alignment_non_full_length, 1)
	expect_match(audit$non_full_length_detail, "MI/85", fixed = TRUE)
})
