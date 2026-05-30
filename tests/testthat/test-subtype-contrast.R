toy_contrast_matrix <- function(values, names = letters[1:4]) {
	mat <- matrix(0, nrow = length(names), ncol = length(names), dimnames = list(names, names))
	mat[lower.tri(mat)] <- values
	mat + t(mat)
}

test_that("cophenetic subtype contrast summary estimates H3N2 minus H1N1", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("bayesboot")

	cophenetic <- toy_contrast_matrix(c(1, 2, 3, 4, 5, 6))
	h1_year <- toy_contrast_matrix(c(1, 3, 2, 5, 4, 6))
	h3_year <- toy_contrast_matrix(c(1, 2, 3, 4, 5, 6))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = h1_year),
		h3 = list(cophenetic = cophenetic, year = h3_year)
	)
	settings <- make_analysis_settings("test")
	settings$mantel_permutations <- 9L
	settings$subtype_contrast_bootstrap_reps <- 7L

	out <- calculate_cophenetic_subtype_contrast_summary(
		distances,
		settings = settings,
		methods = c(year = "Temporal distance")
	)

	expect_equal(nrow(out), 1)
	expect_equal(out$Comparison, "Temporal distance")
	expect_equal(out$`H1N1 pairwise comparisons`, 6)
	expect_equal(out$`H3N2 pairwise comparisons`, 6)
	expect_equal(out$`H1N1 Mantel r`, stats::cor(lower_triangle_values(cophenetic), lower_triangle_values(h1_year)))
	expect_equal(out$`H3N2 Mantel r`, 1)
	expect_equal(out$`Difference (H3N2 - H1N1)`, out$`H3N2 Mantel r` - out$`H1N1 Mantel r`)
	expect_equal(out$`Bootstrap draws`, 7)
	expect_equal(out$`CI method`, "Bayesian bootstrap contrast over strains")
	expect_true(out$`CI lower` <= out$`CI upper`)
})

test_that("subtype contrast sensitivity summaries report Fisher and permutation alternatives", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	cophenetic <- toy_contrast_matrix(c(1, 2, 3, 4, 5, 6))
	h1_year <- toy_contrast_matrix(c(1, 3, 2, 5, 4, 6))
	h3_year <- toy_contrast_matrix(c(1, 2, 3, 4, 5, 6))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = h1_year),
		h3 = list(cophenetic = cophenetic, year = h3_year)
	)
	settings <- make_analysis_settings("test")
	settings$subtype_contrast_permutations <- 9L

	out <- calculate_cophenetic_subtype_contrast_sensitivity(
		distances,
		settings = settings,
		methods = c(year = "Temporal distance")
	)

	expect_equal(nrow(out), 2)
	expect_equal(sort(unique(out$`Sensitivity method`)), c("Fisher z approximation", "Independent Mantel permutation contrast"))
	expect_true(all(out$`Sensitivity p-value` >= 0 & out$`Sensitivity p-value` <= 1))
	expect_equal(unique(out$`Difference (H3N2 - H1N1)`), stats::cor(lower_triangle_values(cophenetic), lower_triangle_values(h3_year)) - stats::cor(lower_triangle_values(cophenetic), lower_triangle_values(h1_year)))
	expect_true(all(nzchar(out$Interpretation)))
})

test_that("permutation contrast sensitivity returns NA when the observed contrast is undefined", {
	skip_if_not_installed("tibble")

	cophenetic <- toy_contrast_matrix(c(1, 2, 3, 4, 5, 6))
	constant <- toy_contrast_matrix(rep(1, 6))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = constant),
		h3 = list(cophenetic = cophenetic, year = cophenetic)
	)

	out <- permutation_subtype_contrast_sensitivity(
		distances,
		comparison_method = "year",
		observed_difference = NA_real_,
		permutations = 9L,
		seed = 370L
	)

	expect_equal(out$statistic, NA_real_)
	expect_equal(out$p_value, NA_real_)
})

test_that("subtype contrast table and plot constructors return manuscript-ready objects", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")
	skip_if_not_installed("ggplot2")

	summary <- tibble::tibble(
		method = c("year", "cart"),
		Comparison = c("Temporal distance", "Cartographic distance"),
		`H1N1 pairwise comparisons` = c(153L, 153L),
		`H3N2 pairwise comparisons` = c(210L, 210L),
		`H1N1 Mantel r` = c(0.25, 0.32),
		`H3N2 Mantel r` = c(0.99, 0.84),
		`Difference (H3N2 - H1N1)` = c(0.74, 0.52),
		`CI lower` = c(0.50, 0.20),
		`CI upper` = c(0.90, 0.75),
		`Bootstrap draws` = c(199L, 199L),
		`CI method` = "Bayesian bootstrap contrast over strains",
		`Bootstrap seed` = c(3701L, 3704L)
	)

	expect_s3_class(make_subtype_contrast_table(summary), "flextable")
	expect_s3_class(plot_subtype_contrasts(summary), "ggplot")
})
