toy_panel_vaccine_source <- function() {
	tibble::tibble(
		season = c(
			"2018-2019", "2018-2019", "2018-2019", "2018-2019",
			"2020-2021", "2020-2021", "2020-2021", "2020-2021"
		),
		vaccine = "Fluzone",
		subtype = c("H1N1", "H3N2", "B-Vic", "B-Yam", "H1N1", "H3N2", "B-Vic", "B-Yam"),
		strain = c("MI/15", "Sing/16", "CO/17", "PH/13", "GD/19", "HK/19", "WA/19", "PH/13"),
		excluded_from_high_dose = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
	)
}

toy_panel_strain_provenance <- function() {
	tibble::tibble(
		subtype = c("h1", "h1", "h3", "h3", "h3"),
		short_name = c("MI/15", "Bris/18", "Sing/16", "HK/19", "old/68"),
		strain_name = c(
			"A/H1N1/Michigan/45/2015",
			"A/H1N1/Brisbane/02/2018",
			"A/H3N2/Singapore/infimh-16-0019/2016",
			"A/H3N2/Hong Kong/2671/2019",
			"A/H3N2/Hong Kong/8/1968"
		),
		inclusion_status = "included_in_analysis"
	)
}

test_that("vaccine source validation requires canonical columns and logical high-dose exclusions", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	source <- toy_panel_vaccine_source()
	validated <- validate_vaccine_strain_source(source)

	expect_equal(names(validated), names(source))
	expect_error(
		validate_vaccine_strain_source(dplyr::select(source, -"strain")),
		"Fluzone vaccine strain source is missing required column\\(s\\): strain"
	)
	expect_error(
		validate_vaccine_strain_source(dplyr::mutate(source, excluded_from_high_dose = "False")),
		"excluded_from_high_dose must be logical"
	)
	source_with_missing <- source
	source_with_missing$excluded_from_high_dose[1] <- NA
	expect_error(
		validate_vaccine_strain_source(source_with_missing),
		"excluded_from_high_dose must not contain missing values"
	)
})

test_that("vaccine source records classify high-dose formulation and match standardized panel names", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	source <- toy_panel_vaccine_source()
	panel <- toy_panel_strain_provenance()
	records <- make_vaccine_strain_records(source, panel)

	expect_equal(unique(records$subtype), c("h1", "h3"))
	expect_equal(
		records$high_dose_formulation[records$season == "2018-2019"],
		rep("trivalent", 2)
	)
	expect_equal(
		records$high_dose_formulation[records$season == "2020-2021"],
		rep("quadrivalent", 2)
	)
	expect_equal(
		records$match_status[records$strain == "Sing/16"],
		"matched_analysis_panel"
	)
	expect_equal(
		records$match_status[records$strain == "GD/19"],
		"absent_from_analysis_panel"
	)
	expect_false(any(records$strain %in% c("SG/16", "BR/18")))
})

test_that("panel strain status summarizes vaccine seasons for included strains", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	vaccine_records <- make_vaccine_strain_records(
		toy_panel_vaccine_source(),
		toy_panel_strain_provenance()
	)
	status <- make_panel_strain_status(toy_panel_strain_provenance(), vaccine_records)

	expect_equal(status$is_vaccine_component, c(TRUE, FALSE, TRUE, TRUE, FALSE))
	expect_equal(status$vaccine_seasons[status$short_name == "MI/15"], "2018-2019")
	expect_equal(status$vaccine_seasons[status$short_name == "HK/19"], "2020-2021")
	expect_equal(status$vaccine_seasons[status$short_name == "old/68"], "not a cohort-season Fluzone component")
	expect_equal(status$isolation_year, c(2015L, 2018L, 2016L, 2019L, 1968L))
})

test_that("panel summary reports subtype counts, vaccine counts, and pair completeness", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	vaccine_records <- make_vaccine_strain_records(
		toy_panel_vaccine_source(),
		toy_panel_strain_provenance()
	)
	status <- make_panel_strain_status(toy_panel_strain_provenance(), vaccine_records)
	pair_counts <- tibble::tibble(
		subtype = c("h1", "h1", "h3", "h3"),
		method = c("year", "cart", "year", "cart"),
		included_strains = c(2L, 2L, 3L, 3L),
		expected_unique_pairs = c(1L, 1L, 3L, 3L),
		observed_unique_pairs = c(1L, 1L, 3L, 2L),
		pair_count_complete = c(TRUE, TRUE, TRUE, FALSE)
	)

	summary <- make_panel_summary(status, pair_counts)

	expect_equal(summary$included_strains, c(2L, 3L))
	expect_equal(summary$vaccine_component_strains, c(1L, 2L))
	expect_equal(summary$earliest_year, c(2015L, 1968L))
	expect_equal(summary$latest_year, c(2018L, 2019L))
	expect_equal(summary$expected_unique_pairs, c(1L, 3L))
	expect_equal(summary$all_distance_methods_complete, c(TRUE, FALSE))
})

test_that("panel space coverage summarizes cartographic and phylogenetic ranges", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	full_distance_table <- tibble::tibble(
		subtype = rep(c("h1", "h3"), each = 6),
		method = rep(c("cophenetic", "cophenetic", "cophenetic", "cart", "cart", "cart"), 2),
		Var1 = rep(c("b", "c", "c"), 4),
		Var2 = rep(c("a", "a", "b"), 4),
		d = c(0.1, 0.3, 0.2, 1, 3, 2, 0.4, 0.8, 0.6, 5, 9, 7)
	)
	cartography <- tibble::tibble(
		subtype = c("h1", "h3"),
		map_dimensions = c(2L, 2L),
		antigen_count = c(2L, 3L),
		serum_count = c(10L, 12L),
		missing_titer_percent = c(20, 30)
	)

	coverage <- make_panel_space_coverage(full_distance_table, cartography)

	expect_equal(coverage$cartographic_distance_min, c(1, 5))
	expect_equal(coverage$cartographic_distance_max, c(3, 9))
	expect_equal(coverage$cophenetic_distance_median, c(0.2, 0.6))
	expect_equal(coverage$serum_count, c(10L, 12L))
})

test_that("temporal coverage plot data uses stable vaccine-component labels", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("stringr")

	vaccine_records <- make_vaccine_strain_records(
		toy_panel_vaccine_source(),
		toy_panel_strain_provenance()
	)
	status <- make_panel_strain_status(toy_panel_strain_provenance(), vaccine_records)
	plot_data <- make_panel_temporal_coverage_plot_data(status)

	expect_equal(
		unique(plot_data$vaccine_component_label),
		c("Cohort-season Fluzone component", "Not a cohort-season Fluzone component")
	)
	expect_equal(plot_data$subtype_label, c("H1N1", "H1N1", "H3N2", "H3N2", "H3N2"))
	expect_equal(plot_data$isolation_year, c(2015L, 2018L, 2016L, 2019L, 1968L))
	expect_error(
		make_panel_temporal_coverage_plot_data(dplyr::select(status, -"vaccine_seasons")),
		"panel strain status is missing required column\\(s\\): vaccine_seasons"
	)
})
