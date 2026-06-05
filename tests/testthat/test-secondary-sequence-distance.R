test_that("amino-acid hamming distance uses pairwise deletion and preserves names", {
	seqs <- c(one = "AA-X", two = "AVCX", three = "VVC?")

	distances <- dist_amino_acid_hamming(seqs)

	expect_equal(dim(distances), c(3L, 3L))
	expect_equal(rownames(distances), names(seqs))
	expect_equal(colnames(distances), names(seqs))
	expect_equal(unname(diag(distances)), c(0, 0, 0))
	expect_equal(distances, t(distances))
	expect_equal(distances["one", "two"], 1 / 2)
	expect_equal(distances["one", "three"], 1)
})

test_that("optimal string alignment distance is normalized after pairwise deletion", {
	expect_equal(optimal_string_alignment_pair_distance("AR", "RA"), 1 / 2)
	expect_equal(optimal_string_alignment_pair_distance("AA-X", "AVCX"), 1 / 2)
	expect_true(is.na(optimal_string_alignment_pair_distance("----", "XXXX")))

	seqs <- c(one = "ARCD", two = "RACD", three = "AVCD")
	distances <- dist_optimal_string_alignment(seqs)

	expect_equal(distances["one", "two"], 1 / 4)
	expect_equal(distances["one", "three"], 1 / 4)
	expect_equal(distances, t(distances))
})

test_that("p-all-epitope averages epitope-specific hamming distances", {
	seq_one <- paste(rep("A", 326), collapse = "")
	seq_two_chars <- rep("A", 326)
	seq_two_chars[118] <- "C"
	seq_two_chars[124] <- "C"
	seq_two <- paste(seq_two_chars, collapse = "")
	seqs <- c(one = seq_one, two = seq_two)

	expect_equal(p_all_epitope(seq_one, seq_two, "h1n1"), mean(c(1 / 24, 1 / 22, 0, 0, 0)))

	distances <- dist_p_all_epitope(seqs, "h1n1")

	expect_equal(distances["one", "two"], mean(c(1 / 24, 1 / 22, 0, 0, 0)))
	expect_equal(distances, t(distances))
})

test_that("BLOSUM62 distance is zero for identities and lower for conservative changes", {
	expect_equal(blosum62_pair_distance("AAAA", "AAAA"), 0)
	expect_lt(blosum62_pair_distance("I", "L"), blosum62_pair_distance("C", "W"))
	expect_true(is.na(blosum62_pair_distance("----", "XXXX")))

	seqs <- c(one = "II", two = "IL", three = "CW")
	distances <- dist_blosum62(seqs)

	expect_equal(distances["one", "two"], blosum62_pair_distance("II", "IL"))
	expect_equal(distances, t(distances))
})

test_that("secondary sequence distances are supplementary and use expected comparators", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("purrr")
	skip_if_not_installed("flextable")

	alignment_result <- list(
		subtype = "h1",
		aligned_sequences = tibble::tibble(
			short_name = c("one", "two", "three", "four"),
			pro_aligned = c(
				paste(rep("A", 326), collapse = ""),
				paste(rep("V", 326), collapse = ""),
				paste(rep("S", 326), collapse = ""),
				paste(rep("T", 326), collapse = "")
			),
			source_full_length = c(TRUE, TRUE, TRUE, TRUE)
		)
	)
	secondary <- calculate_secondary_sequence_distances(list(h1 = alignment_result))

	expect_equal(names(secondary$h1), c("hamming", "osa", "p_all_epitope", "blosum62"))

	make_mat <- function(values) {
		mat <- matrix(0, nrow = 4, ncol = 4, dimnames = list(c("one", "two", "three", "four"), c("one", "two", "three", "four")))
		mat[lower.tri(mat)] <- values
		mat + t(mat)
	}
	primary <- list(
		h1 = list(
			cophenetic = make_mat(c(1, 2, 3, 4, 5, 6)),
			grantham = make_mat(c(1, 2, 3, 4, 5, 6)),
			pepi = make_mat(c(1, 1, 2, 3, 5, 8))
		)
	)
	secondary_for_sensitivity <- list(
		h1 = list(
			hamming = make_mat(c(1, 2, 3, 4, 5, 6)),
			osa = make_mat(c(6, 5, 4, 3, 2, 1)),
			p_all_epitope = make_mat(c(1, 1, 2, 3, 5, 8)),
			blosum62 = make_mat(c(1, 2, 3, 4, 5, 6))
		)
	)

	out <- calculate_secondary_sequence_distance_sensitivity(
		primary,
		secondary_for_sensitivity,
		settings = make_analysis_settings("test"),
		threshold = 0.10
	)

	expect_equal(nrow(out), 4)
	expect_equal(out$Scope, rep("H1N1", 4))
	expect_equal(out$`Primary comparator`, c(
		"Grantham distance",
		"Grantham distance",
		"p-Epitope distance",
		"Grantham distance"
	))
	expect_equal(out$`Sensitivity flag`[out$method == "osa"], "yes")
	expect_false(any(names(out) %in% c("Permutation p-value", "SH p-value", "RF distance")))
	expect_s3_class(make_secondary_sequence_distance_sensitivity_table(out), "flextable")
})
