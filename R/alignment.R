###
# Multiple sequence alignment helpers shared across influenza subtypes
# Zane Billings
###

align_subtype_sequences <- function(clean_sequences, subtype, settings = make_analysis_settings()) {
	check_required_columns(
		clean_sequences,
		c("subtype", "short_name", "protein_sequence", "full_length"),
		"clean sequences"
	)
	subtype_requested <- stringr::str_to_lower(subtype)
	subtype_data <- clean_sequences |>
		dplyr::filter(.data$subtype == subtype_requested) |>
		dplyr::arrange(.data$factor_order)

	if (nrow(subtype_data) == 0) {
		stop("No sequences found for subtype: ", subtype_requested, call. = FALSE)
	}

	protein_sequences <- with(
		subtype_data,
		rlang::set_names(protein_sequence, short_name)
	)

	protein_msa <- msa::msa(
		protein_sequences,
		method = settings$alignment_method,
		type = "protein",
		order = "input",
		verbose = FALSE
	)

	aligned_sequences <- tibble::tibble(
		short_name = subtype_data$short_name,
		pro_aligned = alignment_to_character(protein_msa),
		source_full_length = subtype_data$full_length
	)

	list(
		subtype = subtype_requested,
		protein_msa = protein_msa,
		aligned_sequences = aligned_sequences
	)
}

alignment_completeness_audit <- function(alignment_result) {
	aligned_sequences <- alignment_result$aligned_sequences
	check_required_columns(
		aligned_sequences,
		c("short_name", "pro_aligned", "source_full_length"),
		"aligned sequences"
	)

	audit <- aligned_sequences |>
		dplyr::mutate(
			protein_non_gap_length = nchar(stringr::str_remove_all(.data$pro_aligned, "-"))
		)
	reference_full_length_min <- audit |>
		dplyr::filter(.data$source_full_length) |>
		dplyr::summarise(reference = min(.data$protein_non_gap_length, na.rm = TRUE)) |>
		dplyr::pull("reference")
	if (length(reference_full_length_min) == 0 || is.infinite(reference_full_length_min)) {
		reference_full_length_min <- NA_integer_
	}

	audit |>
		dplyr::mutate(
			reference_full_length_min = as.integer(reference_full_length_min),
			alignment_full_length = dplyr::if_else(
				is.na(.data$reference_full_length_min),
				NA,
				.data$protein_non_gap_length >= .data$reference_full_length_min
			),
			full_length_consistent = .data$source_full_length == .data$alignment_full_length
		)
}

alignment_audit <- function(alignment_result) {
	aligned_sequences <- alignment_result$aligned_sequences
	completeness <- alignment_completeness_audit(alignment_result)
	protein_width <- unique(nchar(aligned_sequences$pro_aligned))
	if (length(protein_width) != 1) {
		stop("Protein alignment widths should be identical within an alignment.", call. = FALSE)
	}
	non_full_length <- completeness |>
		dplyr::filter(!.data$source_full_length | !.data$alignment_full_length)

	tibble::tibble(
		subtype = alignment_result$subtype,
		n_sequences = nrow(aligned_sequences),
		protein_alignment_width = protein_width[[1]],
		source_non_full_length = sum(!completeness$source_full_length),
		alignment_non_full_length = sum(!completeness$alignment_full_length, na.rm = TRUE),
		non_full_length_detail = paste(non_full_length$short_name, collapse = ", ")
	)
}

#### END OF FILE ####
