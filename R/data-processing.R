###
# Sequence data import, validation, and cleaning
# Zane Billings
###

load_raw_sequences <- function(path) {
	validate_file_exists(path, "sequence CSV")
	readr::read_csv(
		path,
		show_col_types = FALSE,
		col_types = readr::cols(
			`Strain Type` = readr::col_character(),
			`Strain Subtype` = readr::col_character(),
			`Strain Name` = readr::col_character(),
			`Short Name` = readr::col_character(),
			`Protein Sequence Source` = readr::col_character(),
			`Protein Sequence Accession Number` = readr::col_character(),
			`Full length?` = readr::col_logical(),
			`Protein Sequence` = readr::col_character()
		)
	) |>
		janitor::clean_names()
}

validate_source_inputs <- function(sequence_file, ace_files) {
	files <- c(
		sequence_file = sequence_file,
		ace_files
	)
	
	tibble::tibble(
		input = names(files),
		path = unname(files)
	) |>
		dplyr::mutate(
			absolute_path = purrr::map_chr(.data$path, validate_file_exists),
			size_bytes = file.info(.data$absolute_path)$size
		)
}

normalize_sequence_subtype <- function(subtype) {
	dplyr::case_when(
		stringr::str_to_lower(subtype) == "h1n1" ~ "h1",
		stringr::str_to_lower(subtype) == "h3n2" ~ "h3",
		TRUE ~ stringr::str_to_lower(subtype)
	)
}

extract_strain_year <- function(strain_name) {
	stringr::str_extract(strain_name, "[0-9]{4}$") |>
		as.integer()
}

extract_strain_location <- function(strain_name) {
	purrr::map_chr(
		stringr::str_split(strain_name, "/", simplify = FALSE),
		\(parts) {
			if (length(parts) < 4) {
				return(NA_character_)
			}
			location <- parts[[3]]
			stringr::str_remove(location, "\\s+[0-9]+$")
		}
	)
}

make_cartography_analysis_name <- function(strain_name, subtype) {
	location <- extract_strain_location(strain_name)
	year <- extract_strain_year(strain_name)
	subtype_label <- dplyr::case_when(
		subtype == "h1" ~ "H1N1",
		subtype == "h3" ~ "H3N2",
		TRUE ~ subtype
	)
	dplyr::if_else(
		is.na(location) | is.na(year),
		NA_character_,
		paste(subtype_label, location, year, sep = "-")
	)
}

extract_cartography_antigen_names <- function(h1_cartography_map, h3_cartography_map) {
	list(
		h1 = Racmacs::agNames(h1_cartography_map),
		h3 = Racmacs::agNames(h3_cartography_map)
	)
}

resolve_cartography_name <- function(candidate, subtype, cartography_antigens, max_distance = 2L) {
	antigens <- cartography_antigens[[subtype]]
	if (is.null(antigens) || is.na(candidate)) {
		return(NA_character_)
	}
	if (candidate %in% antigens) {
		return(candidate)
	}
	candidate_year <- extract_year(candidate)
	antigen_years <- extract_year(antigens)
	eligible <- antigens[antigen_years == candidate_year]
	if (length(eligible) == 0) {
		return(NA_character_)
	}
	distances <- utils::adist(stringr::str_to_lower(candidate), stringr::str_to_lower(eligible))
	best <- which.min(distances)
	if (length(best) == 0 || distances[[best]] > max_distance) {
		return(NA_character_)
	}
	eligible[[best]]
}

cartography_match_status <- function(candidate, resolved) {
	dplyr::case_when(
		is.na(resolved) ~ "not_in_cartography_map",
		identical(candidate, resolved) ~ "exact",
		TRUE ~ "fuzzy_name_resolution"
	)
}

prepare_sequence_records <- function(raw_sequences, cartography_antigens) {
	check_required_columns(
		raw_sequences,
		c(
			"strain_type", "strain_subtype", "strain_name", "short_name",
			"protein_sequence_source", "protein_sequence_accession_number",
			"full_length", "protein_sequence"
		),
		"raw sequence CSV"
	)
	
	out <- raw_sequences |>
		dplyr::mutate(
			dplyr::across(dplyr::where(is.character), stringr::str_squish),
			subtype = normalize_sequence_subtype(.data$strain_subtype),
			short_name = stringr::str_squish(.data$short_name),
			sequence_join_key = normalize_join_key(.data$short_name),
			strain_year = extract_strain_year(.data$strain_name),
			candidate_analysis_name = make_cartography_analysis_name(
				.data$strain_name,
				.data$subtype
			),
			analysis_name = purrr::pmap_chr(
				list(.data$candidate_analysis_name, .data$subtype),
				\(candidate, subtype) resolve_cartography_name(
					candidate,
					subtype,
					cartography_antigens
				)
			),
			cartography_match_status = purrr::map2_chr(
				.data$candidate_analysis_name,
				.data$analysis_name,
				cartography_match_status
			),
			protein_sequence = .data$protein_sequence |>
				stringr::str_to_lower() |>
				stringr::str_remove_all("\\s"),
			protein_length = nchar(.data$protein_sequence),
			factor_order = dplyr::row_number(),
			genbank_strain_name = .data$strain_name,
			vaccine_strain = FALSE
		)
	
	validate_unique_values(out$sequence_join_key, "raw sequence short names")
	validate_unique_values(
		paste(out$protein_sequence_source, out$protein_sequence_accession_number),
		"raw sequence source accession values"
	)
	out
}

clean_sequence_data <- function(raw_sequences, cartography_antigens) {
	prepare_sequence_records(raw_sequences, cartography_antigens) |>
		dplyr::filter(!is.na(.data$analysis_name)) |>
		dplyr::arrange(.data$subtype, .data$factor_order) |>
		dplyr::select(
			dplyr::all_of(c(
				"strain_type", "strain_subtype", "subtype", "strain_name",
				"genbank_strain_name", "analysis_name", "candidate_analysis_name",
				"short_name", "factor_order", "vaccine_strain",
				"protein_sequence_source", "protein_sequence_accession_number",
				"protein_sequence", "protein_length", "full_length",
				"strain_year", "cartography_match_status"
			))
		)
}

sequence_cleaning_audit <- function(raw_sequences, cartography_antigens, clean_sequences) {
	check_required_columns(
		clean_sequences,
		c("subtype", "short_name", "full_length"),
		"clean sequences"
	)
	
	source_records <- prepare_sequence_records(raw_sequences, cartography_antigens)
	excluded_from_analysis <- source_records |>
		dplyr::filter(is.na(.data$analysis_name))
	non_full_source <- source_records |>
		dplyr::filter(!.data$full_length)
	non_full_analysis <- clean_sequences |>
		dplyr::filter(!.data$full_length)
	
	dplyr::bind_rows(
		tibble::tibble(
			check = "raw_sequence_rows",
			value = nrow(raw_sequences),
			detail = NA_character_
		),
		tibble::tibble(
			check = "analysis_sequence_rows",
			value = nrow(clean_sequences),
			detail = NA_character_
		),
		clean_sequences |>
			dplyr::count(.data$subtype, name = "value") |>
			dplyr::transmute(
				check = paste0("analysis_sequence_rows_", .data$subtype),
				value = .data$value,
				detail = NA_character_
			),
		tibble::tibble(
			check = "source_non_full_length_rows",
			value = nrow(non_full_source),
			detail = paste(non_full_source$short_name, collapse = ", ")
		),
		tibble::tibble(
			check = "analysis_non_full_length_rows",
			value = nrow(non_full_analysis),
			detail = paste(non_full_analysis$short_name, collapse = ", ")
		),
		tibble::tibble(
			check = "sequence_rows_not_in_cartography_maps",
			value = nrow(excluded_from_analysis),
			detail = paste(excluded_from_analysis$short_name, collapse = ", ")
		)
	)
}

split_sequences_by_subtype <- function(clean_sequences) {
	check_required_columns(clean_sequences, c("subtype", "short_name"), "clean sequences")
	split(clean_sequences, clean_sequences$subtype)
}

#### END OF FILE ####
