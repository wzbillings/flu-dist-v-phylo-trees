toy_sensitivity_matrix <- function(values, names = letters[1:4]) {
	mat <- matrix(0, nrow = length(names), ncol = length(names), dimnames = list(names, names))
	mat[lower.tri(mat)] <- values
	mat + t(mat)
}

test_that("sequence character audit classifies residues, gaps, ambiguity, and unexpected codes", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	alignment_result <- list(
		subtype = "h3",
		aligned_sequences = tibble::tibble(
			short_name = c("one", "two"),
			pro_aligned = c("ACD-X?", "BJZUO@"),
			source_full_length = c(TRUE, TRUE)
		)
	)

	out <- make_sequence_character_audit(list(h3 = alignment_result))

	expect_equal(out$n[out$residue_class == "standard_residue"], 3)
	expect_equal(out$n[out$residue_class == "gap"], 1)
	expect_equal(out$n[out$residue_class == "ambiguous_or_missing"], 5)
	expect_equal(out$n[out$residue_class == "known_nonstandard_residue"], 2)
	expect_equal(out$n[out$residue_class == "unexpected"], 1)
	expect_equal(out$characters[out$residue_class == "unexpected"], "@")
})

test_that("complete-sequence sensitivity removes non-full-length included strains from matrices", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	cophenetic <- toy_sensitivity_matrix(c(1, 2, 3, 4, 5, 6))
	year <- toy_sensitivity_matrix(c(1, 2, 3, 4, 5, 20))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = year)
	)
	provenance <- tibble::tibble(
		subtype = "h1",
		short_name = letters[1:4],
		inclusion_status = "included_in_analysis",
		source_full_length = c(TRUE, FALSE, TRUE, TRUE),
		alignment_full_length = c(TRUE, FALSE, TRUE, TRUE)
	)

	out <- calculate_complete_sequence_matrix_sensitivity(
		distances,
		provenance,
		methods = c(year = "Temporal distance"),
		threshold = 0.05
	)

	expected_full <- stats::cor(lower_triangle_values(cophenetic), lower_triangle_values(year))
	expected_complete <- stats::cor(
		lower_triangle_values(cophenetic[c("a", "c", "d"), c("a", "c", "d")]),
		lower_triangle_values(year[c("a", "c", "d"), c("a", "c", "d")])
	)

	expect_equal(nrow(out), 1)
	expect_equal(out$`Removed non-full-length strains`, "b")
	expect_equal(out$`Full pairwise comparisons`, 6)
	expect_equal(out$`Complete-sequence pairwise comparisons`, 3)
	expect_equal(out$`Full estimate`, expected_full)
	expect_equal(out$`Complete-sequence estimate`, expected_complete)
	expect_equal(out$`Estimate change`, expected_complete - expected_full)
	expect_equal(out$`Sensitivity flag`, "yes")
})

test_that("sequence missing-data sensitivity compares pairwise and complete deletion distances", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	alignment_result <- list(
		subtype = "h1",
		aligned_sequences = tibble::tibble(
			short_name = c("one", "two", "three", "four"),
			pro_aligned = c("AAAA", "AVXA", "AV-A", "VVVA"),
			source_full_length = c(TRUE, TRUE, TRUE, TRUE)
		)
	)
	cophenetic <- toy_sensitivity_matrix(c(1, 2, 3, 4, 5, 6), names = c("one", "two", "three", "four"))
	primary <- list(
		h1 = list(
			cophenetic = cophenetic,
			grantham = dist_grantham(extract_aligned_proteins(alignment_result))
		)
	)

	out <- calculate_sequence_missing_data_sensitivity(
		list(h1 = alignment_result),
		primary,
		methods = c(grantham = "Grantham distance")
	)

	expect_equal(nrow(out), 1)
	expect_equal(out$Comparison, "Grantham distance")
	expect_equal(out$`Primary deletion rule`, "pairwise")
	expect_equal(out$`Sensitivity deletion rule`, "complete")
	expect_equal(out$`Total candidate sites`, 4)
	expect_equal(out$`Complete-deletion retained sites`, 3)
	expect_equal(out$`Sensitivity flag`, "no")
	expect_true(out$`Maximum absolute distance change` > 0)
})
