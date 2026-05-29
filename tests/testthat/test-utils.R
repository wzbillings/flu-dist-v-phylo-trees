test_that("required column checks report missing columns", {
	skip_if_not_installed("tibble")
	data <- tibble::tibble(a = 1, b = 2)

	expect_invisible(check_required_columns(data, c("a", "b"), "toy data"))
	expect_error(
		check_required_columns(data, c("a", "c"), "toy data"),
		"toy data is missing required column\\(s\\): c",
		fixed = FALSE
	)
})

test_that("join keys and duplicate validation are deterministic", {
	skip_if_not_installed("stringr")

	expect_identical(normalize_join_key(c(" A / Test  ", "MiXeD  Case")), c("a / test", "mixed case"))
	expect_invisible(validate_unique_values(c("a", "b", "c"), "ids"))
	expect_error(validate_unique_values(c("a", "b", "a"), "ids"), "ids contains duplicate value\\(s\\): a")
})

test_that("distance matrix validation catches structural problems", {
	good <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
	bad_order <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("a", "b"), c("b", "a")))

	expect_invisible(validate_distance_matrix(good, "good"))
	expect_error(validate_distance_matrix(1:3, "bad"), "bad must be a matrix.")
	expect_error(validate_distance_matrix(matrix(1:6, nrow = 2), "bad"), "bad must be square.")
	expect_error(validate_distance_matrix(matrix(1:4, nrow = 2), "bad"), "bad must have row and column names.")
	expect_error(validate_distance_matrix(bad_order, "bad"), "bad row and column names must be identical and ordered.")
})

test_that("distance matrix tidying supports full and unique-pair outputs", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("tidyr")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("forcats")

	mat <- matrix(
		c(0, 1, 2, 1, 0, 3, 2, 3, 0),
		nrow = 3,
		dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
	)

	full <- tidy_dist_mat(mat)
	unique_pairs <- tidy_unique_dist_mat(mat)

	expect_equal(nrow(full), 9)
	expect_equal(nrow(unique_pairs), 3)
	expect_equal(as.character(unique_pairs$Var1), c("b", "c", "c"))
	expect_equal(as.character(unique_pairs$Var2), c("a", "a", "b"))
	expect_equal(unique_pairs$d, c(1, 2, 3))
})

test_that("normalization handles varying and constant values", {
	mat <- matrix(c(2, 4, 6, 8), nrow = 2)

	expect_equal(normalize_matrix(mat), matrix(c(0, 1 / 3, 2 / 3, 1), nrow = 2))
	expect_equal(normalize_matrix(matrix(5, nrow = 2, ncol = 2)), matrix(0, nrow = 2, ncol = 2))
	expect_equal(normalize_vector(c(2, 4, 6)), c(0, 0.5, 1))
	expect_equal(normalize_vector(c(7, 7, 7)), c(0, 0, 0))
})

test_that("strain-name replacement and year distances use explicit metadata", {
	skip_if_not_installed("tibble")

	virus_info <- tibble::tibble(
		analysis_name = c("A/Test/1/2001", "A/Test/2/2004"),
		genbank_strain_name = c("full one", "full two"),
		short_name = c("one", "two"),
		subtype = c("h1", "h1")
	)

	expect_equal(replace_strain_names(c("one", "two"), virus_info, from = "short", to = "analysis"), virus_info$analysis_name)
	expect_error(replace_strain_names("missing", virus_info, from = "short", to = "analysis"), "absent from virus metadata")

	year_dist <- dist_year(c("one", "two"), virus_info, format = "short")
	expect_equal(rownames(year_dist), c("one", "two"))
	expect_equal(year_dist["one", "two"], 3)
	expect_error(dist.year(c("one", "two")), "`virus_info` is required")
})

test_that("seeded perturbation is reproducible and preserves distance-matrix structure", {
	mat <- matrix(
		c(0, 1, 2, 1, 0, 3, 2, 3, 0),
		nrow = 3,
		dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
	)

	p1 <- perturb_matrix(mat, mag = 0.01, seed = 123)
	p2 <- perturb_matrix(mat, mag = 0.01, seed = 123)

	expect_equal(p1, p2)
	expect_equal(unname(diag(p1)), c(0, 0, 0))
	expect_equal(p1, t(p1))
	expect_equal(dimnames(p1), dimnames(mat))
})

test_that("analysis settings expose fast test and full modes", {
	test_settings <- make_analysis_settings("test")
	full_settings <- make_analysis_settings("FULL")

	expect_equal(test_settings$mode, "test")
	expect_equal(full_settings$mode, "full")
	expect_true(test_settings$sh_bootstrap < full_settings$sh_bootstrap)
	expect_true(test_settings$mantel_permutations < full_settings$mantel_permutations)
	expect_true(test_settings$mantel_bootstrap_reps < full_settings$mantel_bootstrap_reps)
	expect_error(make_analysis_settings("other"), "FLU_TARGETS_MODE must be either 'test' or 'full'.")
})
