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
