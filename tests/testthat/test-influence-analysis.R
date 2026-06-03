toy_influence_matrix <- function(values, names = letters[1:4]) {
	mat <- matrix(0, nrow = length(names), ncol = length(names), dimnames = list(names, names))
	mat[lower.tri(mat)] <- values
	mat + t(mat)
}

test_that("cophenetic influence summary recomputes estimates after removing each strain", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	cophenetic <- toy_influence_matrix(c(1, 2, 3, 4, 5, 6))
	year <- toy_influence_matrix(c(1, 2, 3, 4, 5, 20))
	distances <- list(
		h1 = list(cophenetic = cophenetic, year = year)
	)

	out <- calculate_cophenetic_influence_summary(
		distances,
		methods = c(year = "Temporal distance"),
		estimate_methods = c("mantel", "pearson"),
		threshold = 0.05
	)

	expect_equal(nrow(out), 8)
	expect_equal(
		names(out),
		c(
			"method",
			"Comparison",
			"subtype",
			"Subtype",
			"Removed strain",
			"Estimate method",
			"Full pairwise comparisons",
			"Leave-one-out pairwise comparisons",
			"Full estimate",
			"Leave-one-out estimate",
			"Estimate change",
			"Absolute change",
			"Flag threshold",
			"Influence flag"
		)
	)
	expect_equal(unique(out$Comparison), "Temporal distance")
	expect_equal(unique(out$Subtype), "H1N1")
	expect_equal(sort(unique(out$`Estimate method`)), c("Descriptive Pearson r", "Mantel r"))
	expect_equal(unique(out$`Full pairwise comparisons`), 6)
	expect_equal(unique(out$`Leave-one-out pairwise comparisons`), 3)

	expected_full <- stats::cor(lower_triangle_values(cophenetic), lower_triangle_values(year))
	expected_without_d <- stats::cor(
		lower_triangle_values(cophenetic[c("a", "b", "c"), c("a", "b", "c")]),
		lower_triangle_values(year[c("a", "b", "c"), c("a", "b", "c")])
	)
	without_d <- out |>
		dplyr::filter(.data$`Removed strain` == "d", .data$`Estimate method` == "Mantel r")

	expect_equal(without_d$`Full estimate`, expected_full)
	expect_equal(without_d$`Leave-one-out estimate`, expected_without_d)
	expect_equal(without_d$`Estimate change`, expected_without_d - expected_full)
	expect_true(any(out$`Influence flag` == "yes"))
})

test_that("influence analysis validates threshold and required matrices", {
	cophenetic <- toy_influence_matrix(c(1, 2, 3, 4, 5, 6))
	year <- toy_influence_matrix(c(1, 2, 3, 4, 5, 6))

	expect_error(
		calculate_cophenetic_influence_summary(
			list(h1 = list(cophenetic = cophenetic, year = year)),
			methods = c(year = "Temporal distance"),
			threshold = -0.1
		),
		"Influence threshold must be a non-negative numeric scalar"
	)
	expect_error(
		calculate_cophenetic_influence_summary(
			list(h1 = list(cophenetic = cophenetic)),
			methods = c(year = "Temporal distance")
		),
		"must include cophenetic and year matrices"
	)
})

test_that("influence table and plot constructors return supplement-ready objects", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("flextable")
	skip_if_not_installed("ggplot2")

	summary <- tibble::tibble(
		method = c("year", "cart"),
		Comparison = c("Temporal distance", "Cartographic distance"),
		subtype = c("h1", "h3"),
		Subtype = c("H1N1", "H3N2"),
		`Removed strain` = c("a", "b"),
		`Estimate method` = c("Mantel r", "Mantel r"),
		`Full pairwise comparisons` = c(153L, 210L),
		`Leave-one-out pairwise comparisons` = c(136L, 190L),
		`Full estimate` = c(0.25, 0.80),
		`Leave-one-out estimate` = c(0.40, 0.74),
		`Estimate change` = c(0.15, -0.06),
		`Absolute change` = c(0.15, 0.06),
		`Flag threshold` = c(0.10, 0.10),
		`Influence flag` = c("yes", "no")
	)

	expect_s3_class(make_cophenetic_influence_table(summary), "flextable")
	expect_s3_class(plot_cophenetic_influence(summary), "ggplot")
})
