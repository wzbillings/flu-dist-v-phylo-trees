###
# Multiple sequence alignment helpers shared across influenza subtypes
# Zane Billings
###

align_subtype_sequences <- function(clean_sequences, subtype, settings = make_analysis_settings()) {
	check_required_columns(
		clean_sequences,
		c("subtype", "short_name", "nucleotide_sequence", "protein_sequence"),
		"clean sequences"
	)
	subtype <- stringr::str_to_lower(subtype)
	subtype_data <- clean_sequences |>
		dplyr::filter(.data$subtype == subtype) |>
		dplyr::arrange(.data$factor_order)
	
	if (nrow(subtype_data) == 0) {
		stop("No sequences found for subtype: ", subtype, call. = FALSE)
	}
	
	nucleotide_sequences <- with(
		subtype_data,
		rlang::set_names(nucleotide_sequence, short_name)
	) |>
		gsub(pattern = "t", replacement = "u")
	
	nucleotide_msa <- msa::msa(
		nucleotide_sequences,
		method = settings$alignment_method,
		type = "rna",
		order = "input",
		verbose = FALSE
	)
	
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
		nuc_aligned = alignment_to_character(nucleotide_msa),
		pro_aligned = alignment_to_character(protein_msa)
	)
	
	list(
		subtype = subtype,
		nucleotide_msa = nucleotide_msa,
		protein_msa = protein_msa,
		aligned_sequences = aligned_sequences
	)
}

alignment_audit <- function(alignment_result) {
	aligned_sequences <- alignment_result$aligned_sequences
	tibble::tibble(
		subtype = alignment_result$subtype,
		n_sequences = nrow(aligned_sequences),
		nucleotide_alignment_width = unique(nchar(aligned_sequences$nuc_aligned)),
		protein_alignment_width = unique(nchar(aligned_sequences$pro_aligned))
	)
}

#### END OF FILE ####
