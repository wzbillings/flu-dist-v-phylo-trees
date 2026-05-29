###
# Sequence data import, validation, and cleaning
# Zane Billings
###

load_virus_metadata <- function(path) {
	validate_file_exists(path, "virus metadata CSV")
	
	virus_info <- readr::read_csv(
		path,
		col_types = readr::cols(
			subtype = readr::col_character(),
			analysis_name = readr::col_character(),
			genbank_strain_name = readr::col_character(),
			short_name = readr::col_character(),
			factor_order = readr::col_integer(),
			vaccine_strain = readr::col_logical()
		)
	)
	
	check_required_columns(
		virus_info,
		c(
			"subtype", "analysis_name", "genbank_strain_name", "short_name",
			"factor_order", "vaccine_strain"
		),
		"virus metadata"
	)
	
	virus_info |>
		dplyr::mutate(
			dplyr::across(
				c("subtype", "analysis_name", "genbank_strain_name", "short_name"),
				stringr::str_squish
			),
			subtype = stringr::str_to_lower(.data$subtype),
			metadata_join_key = normalize_join_key(.data$genbank_strain_name)
		)
}

load_raw_sequences <- function(path) {
	validate_file_exists(path, "sequence workbook")
	readxl::read_xlsx(path) |>
		janitor::clean_names()
}

validate_source_inputs <- function(sequence_file, virus_name_file, ace_files) {
	files <- c(
		sequence_file = sequence_file,
		virus_name_file = virus_name_file,
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

clean_sequence_data <- function(raw_sequences, virus_metadata) {
	check_required_columns(
		raw_sequences,
		c(
			"strain_name", "protein_sequence", "nucleotide_sequence_type",
			"nucleotide_sequence_source", "nucleotide_sequence"
		),
		"raw sequence workbook"
	)
	check_required_columns(
		virus_metadata,
		c(
			"subtype", "analysis_name", "genbank_strain_name", "short_name",
			"factor_order", "vaccine_strain", "metadata_join_key"
		),
		"virus metadata"
	)
	
	validate_unique_values(virus_metadata$metadata_join_key, "virus metadata strain names")
	
	raw_prepped <- raw_sequences |>
		dplyr::mutate(
			strain_name = stringr::str_squish(.data$strain_name),
			sequence_join_key = normalize_join_key(.data$strain_name)
		)
	validate_unique_values(raw_prepped$sequence_join_key, "raw sequence strain names")
	
	missing_metadata <- setdiff(raw_prepped$sequence_join_key, virus_metadata$metadata_join_key)
	if (length(missing_metadata) > 0) {
		missing_names <- raw_prepped$strain_name[match(missing_metadata, raw_prepped$sequence_join_key)]
		stop(
			"Sequence workbook strain(s) missing from virus metadata: ",
			paste(missing_names, collapse = ", "),
			call. = FALSE
		)
	}
	
	metadata_for_join <- virus_metadata |>
		dplyr::rename(sequence_join_key = metadata_join_key)
	
	raw_prepped |>
		dplyr::left_join(
			metadata_for_join,
			by = "sequence_join_key",
			relationship = "many-to-one"
		) |>
		dplyr::arrange(.data$subtype, .data$factor_order) |>
		dplyr::mutate(
			dplyr::across(
				c("protein_sequence", "nucleotide_sequence"),
				\(x) x |>
					stringr::str_to_lower() |>
					stringr::str_remove_all("\\s")
			),
			protein_length = nchar(.data$protein_sequence),
			nucleotide_length = nchar(.data$nucleotide_sequence),
			full_length = .data$short_name != "MI/15"
		) |>
		dplyr::select(
			dplyr::all_of(c(
				"strain_name", "protein_sequence", "nucleotide_sequence_type",
				"nucleotide_sequence_source", "nucleotide_sequence", "subtype",
				"analysis_name", "genbank_strain_name", "short_name", "factor_order",
				"vaccine_strain", "protein_length", "nucleotide_length", "full_length"
			))
		)
}

sequence_cleaning_audit <- function(raw_sequences, virus_metadata, clean_sequences) {
	check_required_columns(clean_sequences, c("subtype", "short_name", "full_length"), "clean sequences")
	
	used_metadata_keys <- normalize_join_key(clean_sequences$genbank_strain_name)
	unused_metadata <- virus_metadata |>
		dplyr::filter(!.data$metadata_join_key %in% used_metadata_keys)
	
	dplyr::bind_rows(
		tibble::tibble(
			check = "raw_sequence_rows",
			value = nrow(raw_sequences),
			detail = NA_character_
		),
		tibble::tibble(
			check = "clean_sequence_rows",
			value = nrow(clean_sequences),
			detail = NA_character_
		),
		clean_sequences |>
			dplyr::count(.data$subtype, name = "value") |>
			dplyr::transmute(
				check = paste0("clean_sequence_rows_", .data$subtype),
				value = .data$value,
				detail = NA_character_
			),
		tibble::tibble(
			check = "flagged_not_full_length",
			value = sum(!clean_sequences$full_length),
			detail = paste(clean_sequences$short_name[!clean_sequences$full_length], collapse = ", ")
		),
		tibble::tibble(
			check = "metadata_rows_not_in_sequence_workbook",
			value = nrow(unused_metadata),
			detail = paste(unused_metadata$short_name, collapse = ", ")
		)
	)
}

split_sequences_by_subtype <- function(clean_sequences) {
	check_required_columns(clean_sequences, c("subtype", "short_name"), "clean sequences")
	split(clean_sequences, clean_sequences$subtype)
}

#### END OF FILE ####
