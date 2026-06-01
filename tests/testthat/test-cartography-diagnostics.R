test_that("titer diagnostics count missing and lower-bound titer cells", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	titers <- matrix(
		c("20", "<10", "*", "80", "", NA_character_),
		nrow = 2,
		dimnames = list(c("ag1", "ag2"), c("sr1", "sr2", "sr3"))
	)

	summary <- summarise_titer_table(titers)

	expect_equal(summary$titer_cells, 6)
	expect_equal(summary$measured_titers, 3)
	expect_equal(summary$missing_titers, 3)
	expect_equal(summary$lower_bound_titers, 1)
	expect_equal(summary$minimum_titer, 10)
	expect_equal(summary$maximum_titer, 80)
	expect_equal(summary$missing_titer_percent, 50)
	expect_equal(summary$lower_bound_titer_percent, 100 / 3)
})

test_that("cartography diagnostics summarize map dimensions and unavailable optimizer history", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	ag_coords <- matrix(
		c(0, 0, 1, 1, 2, 1),
		ncol = 2,
		byrow = TRUE,
		dimnames = list(c("ag1", "ag2", "ag3"), c("x", "y"))
	)
	sr_coords <- matrix(
		c(0, 2, 2, 0),
		ncol = 2,
		byrow = TRUE,
		dimnames = list(c("sr1", "sr2"), c("x", "y"))
	)
	titers <- matrix(
		c("20", "<10", "*", "40", "80", "160"),
		nrow = 3,
		dimnames = list(rownames(ag_coords), rownames(sr_coords))
	)

	diagnostics <- summarise_cartography_diagnostics(
		subtype = "h1",
		antigen_coords = ag_coords,
		serum_coords = sr_coords,
		titer_table = titers,
		map_stress = 12.5,
		projection_count = 1L,
		optimization_metadata = NA_character_
	)

	expect_equal(nrow(diagnostics), 1)
	expect_equal(diagnostics$subtype, "h1")
	expect_equal(diagnostics$subtype_label, "H1N1")
	expect_equal(diagnostics$map_dimensions, 2)
	expect_equal(diagnostics$antigen_count, 3)
	expect_equal(diagnostics$serum_count, 2)
	expect_equal(diagnostics$stored_map_stress, 12.5)
	expect_equal(diagnostics$stored_projection_count, 1L)
	expect_equal(diagnostics$optimization_metadata, "not available in stored .ace input")
	expect_equal(diagnostics$measured_titers, 5)
	expect_equal(diagnostics$missing_titers, 1)
})

test_that("cartography summary records source files supplied by pipeline targets", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	make_map_summary <- function(subtype, source_file = NA_character_) {
		tibble::tibble(
			subtype = subtype,
			subtype_label = subtype_display_label(subtype),
			source_file = source_file,
			map_dimensions = 2L,
			antigen_count = 1L,
			serum_count = 1L,
			stored_projection_count = 1L,
			stored_map_stress = 1,
			optimization_metadata = "not available in stored .ace input",
			titer_cells = 1L,
			measured_titers = 1L,
			missing_titers = 0L,
			missing_titer_percent = 0,
			lower_bound_titers = 0L,
			lower_bound_titer_percent = 0,
			minimum_titer = 10,
			maximum_titer = 10
		)
	}

	helper_env <- environment(make_cartography_diagnostics_summary)
	old_extract <- get("extract_cartography_diagnostics", envir = helper_env)
	on.exit(assign("extract_cartography_diagnostics", old_extract, envir = helper_env), add = TRUE)
	assign(
		"extract_cartography_diagnostics",
		function(cartography_map, subtype, source_file = NA_character_) {
			make_map_summary(subtype, source_file)
		},
		envir = helper_env
	)

	summary <- make_cartography_diagnostics_summary(
		h1_cartography_map = list(),
		h3_cartography_map = list(),
		h1_cartography_file = "custom/h1.ace",
		h3_cartography_file = "custom/h3.ace"
	)

	expect_equal(summary$source_file, c("custom/h1.ace", "custom/h3.ace"))
})

test_that("cartography diagnostics table uses audit-oriented column labels", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")
	skip_if_not_installed("flextable")

	diagnostics <- tibble::tibble(
		subtype = c("h1", "h3"),
		subtype_label = c("H1N1", "H3N2"),
		map_dimensions = c(2L, 2L),
		antigen_count = c(18L, 21L),
		serum_count = c(1524L, 1909L),
		stored_projection_count = c(1L, 1L),
		stored_map_stress = c(30449.36, 35301.42),
		titer_cells = c(27432L, 40089L),
		measured_titers = c(20858L, 27071L),
		missing_titers = c(6574L, 13018L),
		missing_titer_percent = c(23.96544, 32.47275),
		lower_bound_titers = c(5895L, 4180L),
		lower_bound_titer_percent = c(28.26254, 15.44014),
		minimum_titer = c(10, 10),
		maximum_titer = c(5120, 10240),
		optimization_metadata = rep("not available in stored .ace input", 2)
	)

	table_data <- format_cartography_diagnostics_summary(diagnostics)

	expect_equal(
		names(table_data),
		c(
			"Subtype",
			"Dimensions",
			"Strains",
			"Sera",
			"Stored projections",
			"Stored stress",
			"Titer cells",
			"Measured titers",
			"Missing titers",
			"Lower-bound titers",
			"Titer range",
			"Optimization metadata"
		)
	)
	expect_equal(table_data$`Missing titers`, c("6,574 (24.0%)", "13,018 (32.5%)"))
	expect_equal(table_data$`Lower-bound titers`, c("5,895 (28.3% of measured)", "4,180 (15.4% of measured)"))
	expect_s3_class(make_cartography_diagnostics_table(diagnostics), "flextable")
})
