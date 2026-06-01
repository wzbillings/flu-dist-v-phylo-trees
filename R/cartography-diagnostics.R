###
# Antigenic cartography diagnostics
# Zane Billings
###

missing_titer_values <- function(titers) {
	titer_text <- trimws(as.character(titers))
	is.na(titers) | titer_text %in% c("", "*", "NA")
}

numeric_titer_values <- function(titers) {
	titer_text <- trimws(as.character(titers))
	suppressWarnings(as.numeric(gsub("[^0-9.]+", "", titer_text)))
}

summarise_titer_table <- function(titer_table) {
	if (!is.matrix(titer_table)) {
		stop("Titer table must be a matrix.", call. = FALSE)
	}
	titers <- as.vector(titer_table)
	missing <- missing_titer_values(titers)
	measured <- !missing
	lower_bound <- measured & startsWith(trimws(as.character(titers)), "<")
	numeric_titers <- numeric_titer_values(titers)
	measured_numeric <- numeric_titers[measured & !is.na(numeric_titers)]
	titer_cells <- length(titers)
	measured_titers <- sum(measured)

	tibble::tibble(
		titer_cells = titer_cells,
		measured_titers = measured_titers,
		missing_titers = sum(missing),
		missing_titer_percent = if (titer_cells > 0) 100 * sum(missing) / titer_cells else NA_real_,
		lower_bound_titers = sum(lower_bound),
		lower_bound_titer_percent = if (measured_titers > 0) {
			100 * sum(lower_bound) / measured_titers
		} else {
			NA_real_
		},
		minimum_titer = if (length(measured_numeric) > 0) min(measured_numeric) else NA_real_,
		maximum_titer = if (length(measured_numeric) > 0) max(measured_numeric) else NA_real_
	)
}

normalize_unavailable_metadata <- function(value) {
	if (length(value) == 0 || all(is.na(value)) || identical(trimws(as.character(value[[1]])), "")) {
		return("not available in stored .ace input")
	}
	as.character(value[[1]])
}

summarise_cartography_diagnostics <- function(
		subtype,
		antigen_coords,
		serum_coords,
		titer_table,
		map_stress = NA_real_,
		projection_count = NA_integer_,
		optimization_metadata = NA_character_,
		source_file = NA_character_
	) {
	if (!is.matrix(antigen_coords) || !is.matrix(serum_coords)) {
		stop("Antigen and serum coordinates must be matrices.", call. = FALSE)
	}
	if (ncol(antigen_coords) != ncol(serum_coords)) {
		stop("Antigen and serum coordinates must have the same dimensionality.", call. = FALSE)
	}

	titer_summary <- summarise_titer_table(titer_table)
	tibble::tibble(
		subtype = subtype,
		subtype_label = subtype_display_label(subtype),
		source_file = source_file,
		map_dimensions = ncol(antigen_coords),
		antigen_count = nrow(antigen_coords),
		serum_count = nrow(serum_coords),
		stored_projection_count = as.integer(projection_count[[1]]),
		stored_map_stress = as.numeric(map_stress[[1]]),
		optimization_metadata = normalize_unavailable_metadata(optimization_metadata)
	) |>
		dplyr::bind_cols(titer_summary)
}

extract_cartography_diagnostics <- function(cartography_map, subtype, source_file = NA_character_) {
	all_stresses <- tryCatch(
		Racmacs::allMapStresses(cartography_map),
		error = function(e) numeric()
	)
	map_stress <- tryCatch(
		Racmacs::mapStress(cartography_map),
		error = function(e) {
			if (length(all_stresses) > 0) {
				all_stresses[[1]]
			} else {
				NA_real_
			}
		}
	)
	projection_count <- if (length(all_stresses) > 0) length(all_stresses) else NA_integer_

	summarise_cartography_diagnostics(
		subtype = subtype,
		antigen_coords = Racmacs::agCoords(cartography_map),
		serum_coords = Racmacs::srCoords(cartography_map),
		titer_table = Racmacs::titerTable(cartography_map),
		map_stress = map_stress,
		projection_count = projection_count,
		optimization_metadata = NA_character_,
		source_file = source_file
	)
}

make_cartography_diagnostics_summary <- function(
		h1_cartography_map,
		h3_cartography_map,
		h1_cartography_file,
		h3_cartography_file
	) {
	dplyr::bind_rows(
		extract_cartography_diagnostics(h1_cartography_map, "h1", h1_cartography_file),
		extract_cartography_diagnostics(h3_cartography_map, "h3", h3_cartography_file)
	)
}

format_cartography_count <- function(x) {
	dplyr::if_else(
		is.na(x),
		"NA",
		formatC(as.integer(x), format = "d", big.mark = ",")
	)
}

format_cartography_number <- function(x, digits = 1) {
	dplyr::if_else(
		is.na(x),
		"NA",
		formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
	)
}

format_cartography_count_percent <- function(count, percent, suffix = NULL) {
	out <- paste0(format_cartography_count(count), " (", format_cartography_number(percent, digits = 1), "%")
	if (!is.null(suffix)) {
		out <- paste0(out, " ", suffix)
	}
	paste0(out, ")")
}

format_cartography_diagnostics_summary <- function(cartography_diagnostics_summary) {
	cartography_diagnostics_summary |>
		dplyr::transmute(
			Subtype = .data$subtype_label,
			Dimensions = format_cartography_count(.data$map_dimensions),
			Strains = format_cartography_count(.data$antigen_count),
			Sera = format_cartography_count(.data$serum_count),
			`Stored projections` = format_cartography_count(.data$stored_projection_count),
			`Stored stress` = format_cartography_number(.data$stored_map_stress, digits = 2),
			`Titer cells` = format_cartography_count(.data$titer_cells),
			`Measured titers` = format_cartography_count(.data$measured_titers),
			`Missing titers` = format_cartography_count_percent(
				.data$missing_titers,
				.data$missing_titer_percent
			),
			`Lower-bound titers` = format_cartography_count_percent(
				.data$lower_bound_titers,
				.data$lower_bound_titer_percent,
				suffix = "of measured"
			),
			`Titer range` = paste0(
				format_cartography_number(.data$minimum_titer, digits = 0),
				"-",
				format_cartography_number(.data$maximum_titer, digits = 0)
			),
			`Optimization metadata` = .data$optimization_metadata
		)
}

make_cartography_diagnostics_table <- function(cartography_diagnostics_summary) {
	format_cartography_diagnostics_summary(cartography_diagnostics_summary) |>
		flextable::flextable() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Antigenic-cartography diagnostics extracted from the stored .ace ",
			"inputs used as first-pass cartography maps. Missing titers are ",
			"unmeasured cells in the titer table; lower-bound titers are censored ",
			"values such as <10."
		))
}

write_cartography_diagnostics_table <- function(cartography_diagnostics_summary, path) {
	table <- make_cartography_diagnostics_table(cartography_diagnostics_summary)
	write_rds_target(table, path)
}

#### END OF FILE ####
