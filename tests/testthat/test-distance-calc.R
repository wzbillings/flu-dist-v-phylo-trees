toy_distance_matrix <- function(order = c("a", "b", "c")) {
	base <- matrix(
		c(0, 1, 2, 1, 0, 3, 2, 3, 0),
		nrow = 3,
		dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
	)
	base[order, order, drop = FALSE]
}

test_that("distance sets validate common names and ordering", {
	skip_if_not_installed("purrr")

	distances <- list(first = toy_distance_matrix(), second = toy_distance_matrix())
	expect_invisible(validate_distance_set(distances, subtype = "toy"))

	mismatched <- distances
	mismatched$second <- toy_distance_matrix(c("b", "a", "c"))
	expect_error(validate_distance_set(mismatched, subtype = "toy"), "do not share the same strain names/order: second")
})

test_that("standardize_distance_order reorders compatible matrices", {
	skip_if_not_installed("purrr")

	distances <- list(shuffled = toy_distance_matrix(c("c", "a", "b")))
	out <- standardize_distance_order(distances, reference_names = c("a", "b", "c"), subtype = "toy")

	expect_equal(rownames(out$shuffled), c("a", "b", "c"))
	expect_equal(colnames(out$shuffled), c("a", "b", "c"))
	expect_equal(out$shuffled["b", "c"], 3)

	bad <- list(extra = toy_distance_matrix())
	rownames(bad$extra)[1] <- "missing"
	colnames(bad$extra)[1] <- "missing"
	expect_error(standardize_distance_order(bad, c("a", "b", "c"), subtype = "toy"), "incompatible strain names")
})

test_that("combined distance tables preserve subtype and method labels", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("tidyr")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("forcats")
	skip_if_not_installed("purrr")

	distances <- list(
		h1 = list(year = toy_distance_matrix(), grantham = toy_distance_matrix()),
		h3 = list(year = toy_distance_matrix())
	)

	out <- combine_distance_tables(distances, unique_pairs = TRUE)

	expect_equal(nrow(out), 9)
	expect_equal(names(out), c("subtype", "method", "Var1", "Var2", "d"))
	expect_equal(sort(unique(out$subtype)), c("h1", "h3"))
	expect_equal(sort(unique(out$method)), c("grantham", "year"))
})

test_that("distance matrix alignment uses shared strains and requires enough overlap", {
	mat1 <- toy_distance_matrix(c("a", "b", "c"))
	mat2 <- toy_distance_matrix(c("c", "b", "a"))
	aligned <- align_distance_matrices(mat1, mat2)

	expect_equal(rownames(aligned$mat1), c("a", "b", "c"))
	expect_equal(rownames(aligned$mat2), c("a", "b", "c"))
	expect_equal(aligned$mat2["a", "c"], 2)

	small <- mat2[c("a", "b"), c("a", "b")]
	expect_error(align_distance_matrices(mat1, small), "Need at least three shared strains")
})

test_that("paired lower-triangle distances exclude diagonals and duplicate pairs", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	mat1 <- toy_distance_matrix()
	mat2 <- toy_distance_matrix()
	pairs <- paired_lower_triangle_distances(
		mat1,
		mat2,
		subtype = "toy",
		x_name = "cophenetic",
		y_name = "distance"
	)

	expect_equal(nrow(pairs), 3)
	expect_equal(as.character(pairs$Var1), c("b", "c", "c"))
	expect_equal(as.character(pairs$Var2), c("a", "a", "b"))
	expect_false(any(as.character(pairs$Var1) == as.character(pairs$Var2)))
	expect_false(any(paste(pairs$Var1, pairs$Var2) %in% paste(pairs$Var2, pairs$Var1)))
	expect_equal(pairs$cophenetic, c(1, 2, 3))
	expect_equal(pairs$distance, c(1, 2, 3))
})

test_that("mantel_permutation_test returns deterministic result structure", {
	skip_if_not_installed("tibble")

	mat1 <- toy_distance_matrix()
	mat2 <- toy_distance_matrix(c("b", "c", "a"))
	result <- mantel_permutation_test(mat1, mat2, permutations = 9, seed = 123)

	expect_equal(names(result), c("estimate", "p_value", "permutations", "n_pairs", "method"))
	expect_equal(result$permutations, 9)
	expect_equal(result$n_pairs, 3)
	expect_equal(result$method, "mantel_pearson")
	expect_true(result$p_value >= 0 && result$p_value <= 1)
	expect_error(mantel_permutation_test(mat1, mat2, permutations = 0), "positive numeric scalar")
})

test_that("stratified Mantel permutation combines subtype-specific unique pairs", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("purrr")
	skip_if_not_installed("dplyr")

	distances <- list(
		h1 = list(cophenetic = toy_distance_matrix(), year = toy_distance_matrix()),
		h3 = list(cophenetic = toy_distance_matrix(), year = toy_distance_matrix(c("b", "c", "a")))
	)

	result <- stratified_mantel_permutation_test(
		distances,
		"cophenetic",
		"year",
		permutations = 9,
		seed = 123
	)

	expect_equal(names(result), c("estimate", "p_value", "permutations", "n_pairs", "method"))
	expect_equal(result$n_pairs, 6)
	expect_equal(result$method, "stratified_mantel_pearson")
	expect_true(result$p_value >= 0 && result$p_value <= 1)
})
