alignment_sensitivity_toy_matrix <- function(values, names = letters[1:4]) {
	mat <- matrix(0, nrow = length(names), ncol = length(names), dimnames = list(names, names))
	mat[lower.tri(mat)] <- values
	mat + t(mat)
}

test_that("alignment distance sensitivity compares primary and alternative alignment estimates", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("purrr")
	skip_if_not_installed("flextable")

	cophenetic <- alignment_sensitivity_toy_matrix(c(1, 2, 4, 3, 5, 6))
	primary_grantham <- alignment_sensitivity_toy_matrix(c(1, 2, 4, 3, 5, 6))
	alternative_grantham <- alignment_sensitivity_toy_matrix(c(6, 5, 3, 4, 2, 1))
	pepi <- alignment_sensitivity_toy_matrix(c(1, 1, 2, 3, 5, 8))

	primary <- list(
		h1 = list(cophenetic = cophenetic, grantham = primary_grantham, pepi = pepi),
		h3 = list(cophenetic = cophenetic, grantham = primary_grantham, pepi = pepi)
	)
	alternative <- list(
		ClustalW = list(
			h1 = list(grantham = alternative_grantham, pepi = pepi),
			h3 = list(grantham = primary_grantham, pepi = pepi)
		)
	)

	out <- calculate_alignment_distance_sensitivity(
		primary,
		alternative,
		settings = make_analysis_settings("test"),
		methods = c(grantham = "Grantham distance"),
		threshold = 0.10
	)

	expect_equal(nrow(out), 3)
	expect_equal(unique(out$Comparison), "Grantham distance")
	expect_equal(unique(out$`Alternative alignment`), "ClustalW")
	expect_equal(out$Scope, c("Overall", "H1N1", "H3N2"))
	expect_equal(out$`Primary pairwise comparisons`, c(12, 6, 6))
	expect_equal(out$`Alternative pairwise comparisons`, c(12, 6, 6))
	expect_equal(out$`Sensitivity flag`[out$Scope == "H1N1"], "yes")
	expect_equal(out$`Sensitivity flag`[out$Scope == "H3N2"], "no")
	expect_s3_class(make_alignment_distance_sensitivity_table(out), "flextable")
})

test_that("alignment model sensitivity flags primary model only when tolerance is exceeded", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("purrr")
	skip_if_not_installed("flextable")

	primary_tests <- list(
		h1 = tibble::tibble(Model = c("FLU", "FLU+G"), AICc = c(100, 110), logLik = c(-100, -105)),
		h3 = tibble::tibble(Model = c("FLU", "FLU+G"), AICc = c(100, 110), logLik = c(-100, -105))
	)
	alternative_tests <- list(
		ClustalW = list(
			h1 = tibble::tibble(Model = c("FLU", "FLU+G"), AICc = c(104, 100), logLik = c(-105, -100)),
			h3 = tibble::tibble(Model = c("FLU", "FLU+G"), AICc = c(125, 100), logLik = c(-120, -100))
		)
	)
	primary_choice <- list(
		selected_models = c(h1 = "FLU", h3 = "FLU")
	)

	out <- calculate_alignment_model_sensitivity(
		primary_tests,
		alternative_tests,
		primary_choice,
		settings = make_analysis_settings("test"),
		performance_tolerance = 0.10
	)

	expect_equal(nrow(out), 2)
	expect_equal(out$`Selected model changed`, c("yes", "yes"))
	expect_equal(out$`Model sensitivity flag`, c("no", "yes"))
	expect_equal(out$`Primary selected model loss fraction`, c(0.05, 0.20))
	expect_s3_class(make_alignment_model_sensitivity_table(out), "flextable")
})

#### END OF FILE ####
