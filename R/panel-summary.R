###
# Strain panel coverage and vaccine-status summaries
# Zane Billings
###

validate_vaccine_strain_source <- function(vaccine_source) {
	check_required_columns(
		vaccine_source,
		c("season", "vaccine", "subtype", "strain", "excluded_from_high_dose"),
		"Fluzone vaccine strain source"
	)
	if (!is.logical(vaccine_source$excluded_from_high_dose)) {
		stop("excluded_from_high_dose must be logical.", call. = FALSE)
	}
	invisible(vaccine_source)
}

load_vaccine_strain_source <- function(path) {
	validate_file_exists(path, "Fluzone vaccine strain source")
	readr::read_csv(
		path,
		show_col_types = FALSE,
		col_types = readr::cols(
			season = readr::col_character(),
			vaccine = readr::col_character(),
			subtype = readr::col_character(),
			strain = readr::col_character(),
			excluded_from_high_dose = readr::col_logical()
		)
	) |>
		validate_vaccine_strain_source()
}

normalize_vaccine_subtype <- function(subtype) {
	dplyr::case_when(
		stringr::str_to_lower(subtype) == "h1n1" ~ "h1",
		stringr::str_to_lower(subtype) == "h3n2" ~ "h3",
		TRUE ~ NA_character_
	)
}

make_high_dose_formulation_summary <- function(vaccine_source) {
	validate_vaccine_strain_source(vaccine_source)
	vaccine_source |>
		dplyr::group_by(.data$season) |>
		dplyr::summarise(
			high_dose_formulation = dplyr::if_else(
				any(.data$excluded_from_high_dose),
				"trivalent",
				"quadrivalent"
			),
			.groups = "drop"
		)
}

make_vaccine_strain_records <- function(vaccine_source, panel_strains) {
	validate_vaccine_strain_source(vaccine_source)
	check_required_columns(panel_strains, c("subtype", "short_name"), "panel strains")
	
	panel_lookup_base <- if ("inclusion_status" %in% names(panel_strains)) {
		dplyr::filter(panel_strains, .data$inclusion_status == "included_in_analysis")
	} else {
		panel_strains
	}
	panel_lookup <- panel_lookup_base |>
		dplyr::distinct(.data$subtype, .data$short_name) |>
		dplyr::mutate(in_analysis_panel = TRUE)
	
	formulations <- make_high_dose_formulation_summary(vaccine_source)
	
	vaccine_source |>
		dplyr::mutate(
			dplyr::across(dplyr::where(is.character), stringr::str_squish),
			subtype = normalize_vaccine_subtype(.data$subtype),
			strain = stringr::str_squish(.data$strain)
		) |>
		dplyr::filter(!is.na(.data$subtype)) |>
		dplyr::left_join(formulations, by = "season") |>
		dplyr::left_join(
			panel_lookup,
			by = c("subtype", "strain" = "short_name")
		) |>
		dplyr::mutate(
			subtype_label = subtype_display_label(.data$subtype),
			analysis_short_name = .data$strain,
			analysis_short_name = dplyr::if_else(
				is.na(.data$in_analysis_panel),
				NA_character_,
				.data$analysis_short_name
			),
			match_status = dplyr::if_else(
				is.na(.data$in_analysis_panel),
				"absent_from_analysis_panel",
				"matched_analysis_panel"
			)
		) |>
		dplyr::select(
			"season",
			"vaccine",
			"subtype",
			"subtype_label",
			"strain",
			"analysis_short_name",
			"excluded_from_high_dose",
			"high_dose_formulation",
			"match_status"
		) |>
		dplyr::arrange(.data$season, .data$subtype, .data$strain)
}

collapse_panel_values <- function(values, empty = NA_character_) {
	values <- unique(as.character(values[!is.na(values) & values != ""]))
	if (length(values) == 0) {
		return(empty)
	}
	paste(values, collapse = ", ")
}

make_panel_strain_status <- function(strain_provenance, vaccine_strain_records) {
	check_required_columns(
		strain_provenance,
		c("subtype", "short_name", "strain_name", "inclusion_status"),
		"strain provenance records"
	)
	check_required_columns(
		vaccine_strain_records,
		c(
			"season", "subtype", "analysis_short_name",
			"excluded_from_high_dose", "high_dose_formulation", "match_status"
		),
		"vaccine strain records"
	)
	
	vaccine_matches <- vaccine_strain_records |>
		dplyr::filter(.data$match_status == "matched_analysis_panel") |>
		dplyr::group_by(.data$subtype, .data$analysis_short_name) |>
		dplyr::summarise(
			vaccine_seasons = collapse_panel_values(.data$season),
			high_dose_seasons = collapse_panel_values(
				.data$season[!.data$excluded_from_high_dose],
				empty = "excluded from high-dose formulations"
			),
			vaccine_formulations = collapse_panel_values(.data$high_dose_formulation),
			cohort_season_count = dplyr::n_distinct(.data$season),
			.groups = "drop"
		)
	
	strain_provenance |>
		dplyr::mutate(.panel_order = dplyr::row_number()) |>
		dplyr::filter(.data$inclusion_status == "included_in_analysis") |>
		dplyr::mutate(
			isolation_year = if ("strain_year" %in% names(strain_provenance)) {
				as.integer(.data$strain_year)
			} else {
				extract_strain_year(.data$strain_name)
			}
		) |>
		dplyr::left_join(
			vaccine_matches,
			by = c("subtype", "short_name" = "analysis_short_name")
		) |>
		dplyr::mutate(
			subtype_label = subtype_display_label(.data$subtype),
			is_vaccine_component = !is.na(.data$cohort_season_count),
			vaccine_seasons = dplyr::if_else(
				.data$is_vaccine_component,
				.data$vaccine_seasons,
				"not a cohort-season Fluzone component"
			),
			high_dose_seasons = dplyr::if_else(
				.data$is_vaccine_component,
				.data$high_dose_seasons,
				"not a cohort-season Fluzone component"
			),
			vaccine_formulations = dplyr::if_else(
				.data$is_vaccine_component,
				.data$vaccine_formulations,
				"not a cohort-season Fluzone component"
			),
			cohort_season_count = dplyr::coalesce(.data$cohort_season_count, 0L),
			vaccine_match_status = dplyr::if_else(
				.data$is_vaccine_component,
				"matched_analysis_panel",
				"not_matched_vaccine_source"
			)
		) |>
		dplyr::arrange(.data$.panel_order) |>
		dplyr::select(
			"subtype",
			"subtype_label",
			"short_name",
			"strain_name",
			"isolation_year",
			"is_vaccine_component",
			"vaccine_seasons",
			"high_dose_seasons",
			"vaccine_formulations",
			"cohort_season_count",
			"vaccine_match_status"
		)
}

make_panel_summary <- function(panel_strain_status, pair_counts) {
	check_required_columns(
		panel_strain_status,
		c("subtype", "subtype_label", "short_name", "isolation_year", "is_vaccine_component"),
		"panel strain status"
	)
	check_required_columns(
		pair_counts,
		c("subtype", "expected_unique_pairs", "observed_unique_pairs", "pair_count_complete"),
		"pair counts"
	)
	
	strain_summary <- panel_strain_status |>
		dplyr::group_by(.data$subtype, .data$subtype_label) |>
		dplyr::summarise(
			included_strains = dplyr::n_distinct(.data$short_name),
			unique_isolation_years = dplyr::n_distinct(.data$isolation_year),
			earliest_year = min(.data$isolation_year, na.rm = TRUE),
			latest_year = max(.data$isolation_year, na.rm = TRUE),
			year_span = .data$latest_year - .data$earliest_year,
			vaccine_component_strains = sum(.data$is_vaccine_component),
			vaccine_component_percent = 100 * .data$vaccine_component_strains / .data$included_strains,
			.groups = "drop"
		)
	
	pair_summary <- pair_counts |>
		dplyr::group_by(.data$subtype) |>
		dplyr::summarise(
			expected_unique_pairs = max(.data$expected_unique_pairs, na.rm = TRUE),
			observed_unique_pairs_min = min(.data$observed_unique_pairs, na.rm = TRUE),
			observed_unique_pairs_max = max(.data$observed_unique_pairs, na.rm = TRUE),
			all_distance_methods_complete = all(.data$pair_count_complete),
			.groups = "drop"
		)
	
	strain_summary |>
		dplyr::left_join(pair_summary, by = "subtype") |>
		dplyr::arrange(.data$subtype)
}

summarise_panel_distance_range <- function(full_distance_table, method, prefix) {
	full_distance_table |>
		dplyr::filter(.data$method == !!method) |>
		dplyr::group_by(.data$subtype) |>
		dplyr::summarise(
			distance_min = min(.data$d, na.rm = TRUE),
			distance_median = stats::median(.data$d, na.rm = TRUE),
			distance_max = max(.data$d, na.rm = TRUE),
			.groups = "drop"
		) |>
		dplyr::rename_with(
			\(name) paste0(prefix, "_", name),
			dplyr::all_of(c("distance_min", "distance_median", "distance_max"))
		)
}

make_panel_space_coverage <- function(full_distance_table, cartography_diagnostics_summary) {
	check_required_columns(
		full_distance_table,
		c("subtype", "method", "d"),
		"full distance table"
	)
	check_required_columns(
		cartography_diagnostics_summary,
		c("subtype", "map_dimensions", "antigen_count", "serum_count", "missing_titer_percent"),
		"cartography diagnostics summary"
	)
	
	subtypes <- tibble::tibble(subtype = sort(unique(full_distance_table$subtype)))
	cartographic <- summarise_panel_distance_range(full_distance_table, "cart", "cartographic")
	cophenetic <- summarise_panel_distance_range(full_distance_table, "cophenetic", "cophenetic")
	
	subtypes |>
		dplyr::left_join(cartographic, by = "subtype") |>
		dplyr::left_join(cophenetic, by = "subtype") |>
		dplyr::left_join(
			cartography_diagnostics_summary |>
				dplyr::select(
					"subtype",
					"map_dimensions",
					"antigen_count",
					"serum_count",
					"missing_titer_percent"
				),
			by = "subtype"
		) |>
		dplyr::mutate(subtype_label = subtype_display_label(.data$subtype), .after = "subtype")
}

make_panel_temporal_coverage_plot_data <- function(panel_strain_status) {
	check_required_columns(
		panel_strain_status,
		c("subtype_label", "short_name", "isolation_year", "is_vaccine_component"),
		"panel strain status"
	)
	panel_strain_status |>
		dplyr::mutate(
			vaccine_component_label = dplyr::if_else(
				.data$is_vaccine_component,
				"Cohort-season Fluzone component",
				"Not a cohort-season Fluzone component"
			)
		) |>
		dplyr::select(
			"subtype",
			"subtype_label",
			"short_name",
			"isolation_year",
			"is_vaccine_component",
			"vaccine_component_label",
			"vaccine_seasons"
		)
}

plot_panel_temporal_coverage <- function(panel_strain_status) {
	plot_data <- make_panel_temporal_coverage_plot_data(panel_strain_status)
	ggplot2::ggplot(
		plot_data,
		ggplot2::aes(
			x = .data$isolation_year,
			y = .data$subtype_label,
			color = .data$vaccine_component_label,
			shape = .data$vaccine_component_label
		)
	) +
		ggplot2::geom_point(size = 2.8, alpha = 0.85) +
		ggplot2::scale_color_manual(
			values = c(
				"Cohort-season Fluzone component" = "#D55E00",
				"Not a cohort-season Fluzone component" = "#0072B2"
			),
			name = "Panel strain status"
		) +
		ggplot2::scale_shape_manual(
			values = c(
				"Cohort-season Fluzone component" = 17,
				"Not a cohort-season Fluzone component" = 16
			),
			name = "Panel strain status"
		) +
		ggplot2::labs(
			x = "Isolation year",
			y = "Subtype"
		) +
		ggplot2::theme_bw(base_size = 11) +
		ggplot2::theme(
			legend.position = "bottom",
			panel.grid.minor = ggplot2::element_blank()
		)
}

write_panel_temporal_coverage_plot <- function(panel_strain_status, path) {
	ensure_dir(path)
	plot <- plot_panel_temporal_coverage(panel_strain_status)
	ggplot2::ggsave(path, plot = plot, width = 7.5, height = 4.25, dpi = 300)
	path
}

make_panel_summary_table <- function(panel_summary) {
	panel_summary |>
		dplyr::mutate(
			`Vaccine component strains` = paste0(
				.data$vaccine_component_strains,
				" (",
				sprintf("%.1f", .data$vaccine_component_percent),
				"%)"
			),
			`Observed unique pairs` = dplyr::if_else(
				.data$observed_unique_pairs_min == .data$observed_unique_pairs_max,
				as.character(.data$observed_unique_pairs_min),
				paste0(.data$observed_unique_pairs_min, "-", .data$observed_unique_pairs_max)
			),
			`All methods complete` = dplyr::if_else(.data$all_distance_methods_complete, "yes", "no")
		) |>
		dplyr::select(
			"Subtype" = "subtype_label",
			"Included strains" = "included_strains",
			"Unique isolation years" = "unique_isolation_years",
			"Earliest year" = "earliest_year",
			"Latest year" = "latest_year",
			"Year span" = "year_span",
			"Vaccine component strains",
			"Expected unique pairs" = "expected_unique_pairs",
			"Observed unique pairs",
			"All methods complete"
		) |>
		flextable::flextable() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Subtype-level summary of the expert-selected H1N1/H3N2 analysis panel."
		))
}

make_panel_strain_status_table <- function(panel_strain_status) {
	panel_strain_status |>
		dplyr::mutate(
			`Vaccine component` = dplyr::if_else(.data$is_vaccine_component, "yes", "no")
		) |>
		dplyr::select(
			"Subtype" = "subtype_label",
			"Short name" = "short_name",
			"Strain name" = "strain_name",
			"Isolation year" = "isolation_year",
			"Vaccine component",
			"Vaccine seasons" = "vaccine_seasons",
			"High-dose seasons" = "high_dose_seasons",
			"High-dose formulation(s)" = "vaccine_formulations"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Included panel strains and cohort-season Fluzone vaccine-component status.")
}

make_fluzone_vaccine_strain_match_table <- function(vaccine_strain_records) {
	vaccine_strain_records |>
		dplyr::mutate(
			`Excluded from high dose` = dplyr::if_else(.data$excluded_from_high_dose, "yes", "no"),
			`Analysis panel status` = dplyr::case_when(
				.data$match_status == "matched_analysis_panel" ~ "matched",
				TRUE ~ "absent from analysis panel"
			)
		) |>
		dplyr::select(
			"Season" = "season",
			"Vaccine" = "vaccine",
			"Subtype" = "subtype_label",
			"Vaccine strain" = "strain",
			"High-dose formulation" = "high_dose_formulation",
			"Excluded from high dose",
			"Analysis panel status"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Season", "High-dose formulation")) |>
		flextable::valign(j = c("Season", "High-dose formulation"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Fluzone H1N1/H3N2 vaccine components and analysis-panel matching status.")
}

make_panel_space_coverage_table <- function(panel_space_coverage) {
	panel_space_coverage |>
		dplyr::mutate(
			`Cartographic distance range` = paste0(
				sprintf("%.2f", .data$cartographic_distance_min),
				"-",
				sprintf("%.2f", .data$cartographic_distance_max)
			),
			`ML-tree cophenetic range` = paste0(
				sprintf("%.3f", .data$cophenetic_distance_min),
				"-",
				sprintf("%.3f", .data$cophenetic_distance_max)
			),
			`Median cartographic distance` = sprintf("%.2f", .data$cartographic_distance_median),
			`Median ML-tree cophenetic distance` = sprintf("%.3f", .data$cophenetic_distance_median),
			`Missing titers (%)` = sprintf("%.1f", .data$missing_titer_percent)
		) |>
		dplyr::select(
			"Subtype" = "subtype_label",
			"Map dimensions" = "map_dimensions",
			"Map antigens" = "antigen_count",
			"Sera" = "serum_count",
			"Missing titers (%)",
			"Cartographic distance range",
			"Median cartographic distance",
			"ML-tree cophenetic range",
			"Median ML-tree cophenetic distance"
		) |>
		flextable::flextable() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Available cartographic and phylogenetic distance coverage for the analysis panel."
		))
}

write_panel_summary_table <- function(panel_summary, path) {
	make_panel_summary_table(panel_summary) |>
		write_rds_target(path)
}

write_panel_strain_status_table <- function(panel_strain_status, path) {
	make_panel_strain_status_table(panel_strain_status) |>
		write_rds_target(path)
}

write_fluzone_vaccine_strain_match_table <- function(vaccine_strain_records, path) {
	make_fluzone_vaccine_strain_match_table(vaccine_strain_records) |>
		write_rds_target(path)
}

write_panel_space_coverage_table <- function(panel_space_coverage, path) {
	make_panel_space_coverage_table(panel_space_coverage) |>
		write_rds_target(path)
}

#### END OF FILE ####
