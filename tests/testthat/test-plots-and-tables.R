test_that("display helpers expose stable labels and palette names", {
	skip_if_not_installed("dplyr")

	expect_equal(names(subtype_color_palette()), c("H1N1", "H3N2"))
	expect_equal(names(distance_method_labels()), c("year", "grantham", "pepi", "cart"))
	expect_equal(subtype_display_label(c("h1", "h3", "other")), c("H1N1", "H3N2", "other"))
})

test_that("plotting distance normalization groups by method", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	full_distance_table <- tibble::tibble(
		subtype = c("h1", "h1", "h1", "h1"),
		method = c("year", "year", "grantham", "grantham"),
		Var1 = c("b", "c", "b", "c"),
		Var2 = c("a", "a", "a", "a"),
		d = c(2, 4, 10, 20)
	)

	out <- normalize_distances_for_plot(full_distance_table)

	expect_equal(out$d_normalized, c(0, 1, 0, 1))
	expect_equal(names(out), c(names(full_distance_table), "d_normalized"))
})

test_that("weighted correlations match unweighted correlations for equal weights", {
	x <- c(1, 2, 3, NA)
	y <- c(1, 4, 9, 16)
	weights <- c(1, 1, 1, 1)

	expect_equal(weighted_pearson_correlation(x, y), stats::cor(x[1:3], y[1:3]))
	expect_equal(weighted_pearson_correlation(x, y, weights), stats::cor(x[1:3], y[1:3]))
	expect_true(is.na(weighted_pearson_correlation(c(1, 2, 3), c(3, 2, 1), c(0, 0, 0))))
	expect_true(is.na(weighted_pearson_correlation(c(1, 1, 1), c(3, 2, 1), c(1, 1, 1))))
})

test_that("wald correlation intervals return guardrail missing values", {
	valid <- wald_correlation_ci(estimate = 0.5, n = 10, conf_level = 0.95)

	expect_equal(names(valid), c("lower", "upper"))
	expect_true(valid[["lower"]] < 0.5)
	expect_true(valid[["upper"]] > 0.5)
	expect_true(all(is.na(wald_correlation_ci(estimate = NA_real_, n = 10))))
	expect_true(all(is.na(wald_correlation_ci(estimate = 0.5, n = 3))))
	expect_true(all(is.na(wald_correlation_ci(estimate = 1, n = 10))))
})

test_that("strain bootstrap units and pair ids preserve subtype strata", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	data <- tibble::tibble(
		subtype = c("h1", "h1", "h3"),
		Var1 = c("a", "b", "a"),
		Var2 = c("b", "c", "b"),
		cophenetic = c(1, 2, 3),
		distance = c(2, 3, 4)
	)

	ids <- pair_strain_ids(data)
	units <- make_strain_bootstrap_units(data)

	expect_equal(ids$strain1_id, c("h1::a", "h1::b", "h3::a"))
	expect_equal(ids$strain2_id, c("h1::b", "h1::c", "h3::b"))
	expect_equal(units$strain_id, c("h1::a", "h1::b", "h1::c", "h3::a", "h3::b"))
	expect_equal(units$stratum, c("h1", "h1", "h1", "h3", "h3"))
	expect_error(
		make_strain_bootstrap_units(tibble::tibble(Var1 = "a")),
		"correlation pair data is missing required column\\(s\\): Var2"
	)
})

test_that("strain weights normalize within subtype and map to pair weights", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	data <- tibble::tibble(
		subtype = c("h1", "h1", "h3"),
		Var1 = c("a", "b", "a"),
		Var2 = c("b", "c", "b"),
		cophenetic = c(1, 2, 3),
		distance = c(2, 3, 4)
	)
	units <- make_strain_bootstrap_units(data)
	weights <- normalize_strain_weights_by_stratum(units, c(1, 1, 2, 3, 1))

	expect_equal(sum(weights[startsWith(names(weights), "h1::")]), 1)
	expect_equal(sum(weights[startsWith(names(weights), "h3::")]), 1)
	expect_equal(unname(weights["h1::c"]), 0.5)
	expect_error(
		normalize_strain_weights_by_stratum(units, c(1, 2)),
		"Strain weights must match the number of strain bootstrap units."
	)

	pair_weights <- strain_pair_weights(data, weights)
	expect_equal(pair_weights, unname(c(
		weights["h1::a"] * weights["h1::b"],
		weights["h1::b"] * weights["h1::c"],
		weights["h3::a"] * weights["h3::b"]
	)))

	expect_error(
		strain_pair_weights(data, c("h1::a" = 1)),
		"Strain weights are missing endpoint\\(s\\)"
	)
})

test_that("correlation_ci falls back to Wald intervals without pair columns", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	data <- tibble::tibble(
		cophenetic = c(1, 2, 3, 4, 5),
		distance = c(1.1, 1.9, 3.2, 3.8, 5.1)
	)
	estimate <- stats::cor(data$cophenetic, data$distance)

	ci <- correlation_ci(data, estimate, reps = 5, conf_level = 0.95, seed = 1)

	expect_equal(names(ci), c("ci_lower", "ci_upper", "ci_method"))
	expect_equal(ci$ci_method, "Fisher z / Wald")
	expect_true(ci$ci_lower < estimate)
	expect_true(ci$ci_upper > estimate)
})

test_that("cophenetic Mantel summary reports overall and subtype-specific matrix results", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("purrr")
	skip_if_not_installed("flextable")

	make_mat <- function(values) {
		mat <- matrix(0, nrow = 4, ncol = 4, dimnames = list(letters[1:4], letters[1:4]))
		mat[lower.tri(mat)] <- values
		mat <- mat + t(mat)
		mat
	}
	cophenetic <- make_mat(c(1, 2, 4, 3, 5, 6))
	year <- make_mat(c(1, 3, 5, 2, 4, 6))
	grantham <- make_mat(c(2, 1, 4, 3, 6, 5))
	pepi <- make_mat(c(1, 1, 2, 3, 5, 8))
	cart <- make_mat(c(6, 5, 3, 4, 2, 1))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = year, grantham = grantham, pepi = pepi, cart = cart),
		h3 = list(cophenetic = cophenetic, year = cart, grantham = pepi, pepi = grantham, cart = year)
	)
	settings <- make_analysis_settings("test")
	settings$mantel_permutations <- 9L
	settings$mantel_bootstrap_reps <- 5L

	pair_data <- make_cophenetic_pair_table(distances["h1"], "year")
	summary <- calculate_cophenetic_mantel_summary(distances, settings)

	expect_equal(nrow(pair_data), 6)
	expect_false(any(pair_data$Var1 == pair_data$Var2))
	expect_equal(nrow(summary), 12)
	expect_equal(unique(summary$Scope), c("Overall", "H1N1", "H3N2"))
	expect_true(all(summary$Permutations == 9))
	expect_true(all(summary$`Bootstrap draws` == 5))
	expect_false("Pearson r" %in% names(summary))
	expect_s3_class(make_cophenetic_mantel_table(summary), "flextable")
})

test_that("descriptive cophenetic correlations use raw distance values", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("tidyr")

	pair_values <- tibble::tibble(
		subtype = c("h1", "h1", "h3", "h3"),
		Var1 = c("b", "c", "b", "c"),
		Var2 = "a",
		cophenetic = c(1, 2, 100, 110),
		year = c(10, 20, 30, 40),
		grantham = c(1, 2, 3, 4),
		pepi = c(2, 3, 4, 5),
		cart = c(5, 6, 7, 8)
	)
	full_distance_table <- pair_values |>
		tidyr::pivot_longer(
			cols = c("cophenetic", "year", "grantham", "pepi", "cart"),
			names_to = "method",
			values_to = "d"
		)
	settings <- make_analysis_settings("test")
	settings$correlation_bootstrap_reps <- 5L

	summary <- calculate_cophenetic_correlation_summary(full_distance_table, settings)
	overall <- summary |>
		dplyr::filter(.data$Comparison == "Temporal distance", .data$Scope == "Overall")

	expect_equal(overall$`Pearson r`, stats::cor(c(1, 2, 100, 110), c(10, 20, 30, 40)))
})

test_that("distance scale audit is supplement-ready and labels normalized outputs", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")

	audit <- make_distance_scale_audit()

	expect_true(all(c(
		"Output",
		"Reported statistic",
		"Distance scale",
		"Normalization scope",
		"Reason",
		"Supplement placement"
	) %in% names(audit)))
	expect_equal(
		audit$`Distance scale`[audit$Output == "Neighbor-joining distance trees"],
		"Raw distances"
	)
	expect_equal(
		audit$`Distance scale`[audit$Output == "Correlation plot"],
		"Min-max normalized distances"
	)
	expect_equal(
		audit$`Normalization scope`[audit$Output == "Correlation plot"],
		"By metric across all displayed unique off-diagonal pairs"
	)
	expect_true(all(audit$`Supplement placement` %in% c("Main text", "Supplement-ready audit")))
	expect_s3_class(make_distance_scale_audit_table(audit), "flextable")
})

make_toy_tree_analysis <- function(subtype, sh_diffs, rf_distances, model_test = NULL) {
	method_names <- c("ml", names(distance_method_labels()))
	sh_test <- cbind(
		"Trees" = seq_along(method_names),
		"ln L" = -1000 - sh_diffs,
		"Diff ln L" = sh_diffs,
		"p-value" = c(0.95, 0.12, 0.34, 0.001, 0.56)
	)
	rf_matrix <- matrix(
		0,
		nrow = length(method_names),
		ncol = length(method_names),
		dimnames = list(method_names, method_names)
	)
	rf_matrix[, "ml"] <- rf_distances
	rf_matrix["ml", ] <- rf_distances

	list(
		subtype = subtype,
		model_test = model_test,
		selected_model = NULL,
		sh_test = sh_test,
		tree_distance_metrics = list(
			SPR = rf_matrix + 1,
			RF = rf_matrix,
			wRF = rf_matrix / 10,
			KF = rf_matrix / 100,
			path = rf_matrix * 2
		)
	)
}

test_that("tree comparison summary centers delta log likelihood and keeps SH/RF secondary", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")

	h1 <- make_toy_tree_analysis(
		"h1",
		sh_diffs = c(0, 40, 2, 8, 55),
		rf_distances = c(0, 12, 2, 4, 14)
	)
	h3 <- make_toy_tree_analysis(
		"h3",
		sh_diffs = c(0, 5, 4, 12, 60),
		rf_distances = c(0, 6, 8, 10, 16)
	)

	summary <- make_tree_comparison_summary(h1, h3)

	expect_equal(nrow(summary), 8)
	expect_equal(
		names(summary),
		c(
			"Subtype",
			"Distance tree",
			"Delta log likelihood",
			"Tree log likelihood",
			"ML log likelihood",
			"SH p-value",
			"RF distance"
		)
	)
	expect_false("Maximum likelihood baseline" %in% summary$`Distance tree`)
	expect_equal(summary$`Distance tree`[1:4], unname(distance_method_labels()))
	expect_equal(summary$`Delta log likelihood`[summary$Subtype == "H1N1"], c(40, 2, 8, 55))
	expect_equal(summary$`RF distance`[summary$Subtype == "H3N2"], c(6, 8, 10, 16))
	expect_s3_class(make_tree_comparison_table(summary), "flextable")
})

test_that("topology distance summary reports supplemental metrics against the ML tree", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")

	h1 <- make_toy_tree_analysis(
		"h1",
		sh_diffs = c(0, 40, 2, 8, 55),
		rf_distances = c(0, 12, 2, 4, 14)
	)
	h3 <- make_toy_tree_analysis(
		"h3",
		sh_diffs = c(0, 5, 4, 12, 60),
		rf_distances = c(0, 6, 8, 10, 16)
	)

	summary <- make_topology_distance_summary(h1, h3)

	expect_equal(nrow(summary), 8)
	expect_equal(
		names(summary),
		c(
			"Subtype",
			"Distance tree",
			"SPR distance",
			"RF distance",
			"Weighted RF distance",
			"Branch-score distance",
			"Path distance"
		)
	)
	expect_equal(summary$`SPR distance`[1], 13)
	expect_equal(summary$`Weighted RF distance`[1], 1.2)
	expect_equal(summary$`Branch-score distance`[1], 0.12)
	expect_equal(summary$`Path distance`[1], 24)
	expect_s3_class(make_topology_distance_table(summary), "flextable")
})

test_that("model selection prefers a common model within likelihood-loss tolerance", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	h1_model_test <- tibble::tibble(
		Model = c("FLU+G(4)", "FLU+G(4)+I", "FLU"),
		df = c(34, 35, 33),
		logLik = c(-1000, -1005, -1300),
		AICc = c(2068, 2079, 2667),
		AICcw = c(0.99, 0.01, 0),
		BIC = c(2100, 2110, 2700)
	)
	h3_model_test <- tibble::tibble(
		Model = c("FLU+G(4)", "FLU+G(4)+I", "FLU"),
		df = c(40, 41, 39),
		logLik = c(-910, -900, -1200),
		AICc = c(1900, 1882, 2478),
		AICcw = c(0.01, 0.99, 0),
		BIC = c(1950, 1932, 2520)
	)

	choice <- choose_tree_model(
		list(h1 = h1_model_test, h3 = h3_model_test),
		performance_tolerance = 0.10
	)
	summary <- make_model_selection_summary(list(h1 = h1_model_test, h3 = h3_model_test), choice)

	expect_equal(choice$strategy, "common")
	expect_equal(choice$selected_model, "FLU+G(4)+I")
	expect_equal(choice$selected_models[["h1"]], "FLU+G(4)+I")
	expect_equal(choice$selected_models[["h3"]], "FLU+G(4)+I")
	expect_equal(unique(summary$`Selected for ML tree`[summary$Model == "FLU+G(4)+I"]), TRUE)
	expect_true(all(summary$`Log-likelihood loss fraction`[summary$`Selected for ML tree`] <= 0.10))
})

test_that("model selection falls back to subtype-specific models when common loss is too large", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	h1_model_test <- tibble::tibble(
		Model = c("FLU+G(4)", "FLU+G(4)+I"),
		df = c(34, 35),
		logLik = c(-1000, -1400),
		AICc = c(2068, 2870),
		AICcw = c(1, 0),
		BIC = c(2100, 2910)
	)
	h3_model_test <- tibble::tibble(
		Model = c("FLU+G(4)", "FLU+G(4)+I"),
		df = c(40, 41),
		logLik = c(-1400, -900),
		AICc = c(2880, 1882),
		AICcw = c(0, 1),
		BIC = c(2920, 1932)
	)

	choice <- choose_tree_model(
		list(h1 = h1_model_test, h3 = h3_model_test),
		performance_tolerance = 0.10
	)

	expect_equal(choice$strategy, "subtype-specific")
	expect_true(is.na(choice$selected_model))
	expect_equal(choice$selected_models[["h1"]], "FLU+G(4)")
	expect_equal(choice$selected_models[["h3"]], "FLU+G(4)+I")
})

make_toy_support_result <- function(subtype) {
	list(
		branch_support_summary = tibble::tibble(
			subtype = subtype,
			bootstrap_replicates = 10L,
			internal_branches = 3L,
			mean_support_percent = 82,
			median_support_percent = 80,
			min_support_percent = 60,
			branches_ge_70_percent = 2L,
			branches_ge_90_percent = 1L,
			branches_ge_95_percent = 0L
		),
		topology_stability_summary = tibble::tibble(
			subtype = subtype,
			bootstrap_replicates = 10L,
			usable_topology_replicates = 10L,
			distance_failures = 0L,
			identical_topology_fraction = 0.3,
			median_normalized_rf_distance = 0.2,
			mean_normalized_rf_distance = 0.25,
			max_normalized_rf_distance = 0.6,
			median_weighted_rf_distance = 0.4,
			median_branch_score_distance = 0.5,
			median_path_distance = 2
		),
		branch_support_detail = tibble::tibble(
			subtype = subtype,
			node = c(6L, 7L, 8L),
			bootstrap_replicates = 10L,
			support_percent = c(100, 80, 60),
			support_category = c(">=95%", "70-89%", "<70%")
		)
	)
}

test_that("ML support tables combine branch support and topology stability", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")

	h1 <- make_toy_support_result("h1")
	h3 <- make_toy_support_result("h3")

	summary <- make_ml_tree_support_summary(h1, h3)
	detail <- make_ml_branch_support_detail(h1, h3)

	expect_equal(nrow(summary), 2)
	expect_equal(
		names(summary),
		c(
			"Subtype",
			"Bootstrap replicates",
			"Internal branches",
			"Mean branch support (%)",
			"Median branch support (%)",
			"Minimum branch support (%)",
			"Branches >=70%",
			"Branches >=90%",
			"Branches >=95%",
			"Usable topology replicates",
			"Topology distance failures",
			"Identical topology fraction",
			"Median normalized RF distance",
			"Mean normalized RF distance",
			"Maximum normalized RF distance",
			"Median weighted RF distance",
			"Median branch-score distance",
			"Median path distance"
		)
	)
	expect_equal(nrow(detail), 6)
	expect_equal(unique(detail$Subtype), c("H1N1", "H3N2"))
	expect_s3_class(make_ml_tree_support_table(summary), "flextable")
	expect_s3_class(make_ml_branch_support_detail_table(detail), "flextable")
})
