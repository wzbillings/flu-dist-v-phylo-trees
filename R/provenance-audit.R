###
# Strain provenance, inclusion, and cleaning audit outputs
# Zane Billings
###

is_missing_text <- function(x) {
	is.na(x) | stringr::str_squish(as.character(x)) == ""
}

make_alignment_completeness_records <- function(alignments_by_subtype) {
	purrr::imap_dfr(
		alignments_by_subtype,
		\(alignment_result, subtype) alignment_completeness_audit(alignment_result) |>
			dplyr::mutate(subtype = subtype, .before = 1)
	)
}

make_strain_provenance_records <- function(
		raw_sequences,
		cartography_antigens,
		alignment_completeness = NULL
	) {
	source_records <- prepare_sequence_records(raw_sequences, cartography_antigens)
	
	if (is.null(alignment_completeness)) {
		alignment_completeness <- tibble::tibble(
			subtype = character(),
			short_name = character(),
			protein_non_gap_length = integer(),
			alignment_full_length = logical()
		)
	}
	check_required_columns(
		alignment_completeness,
		c("subtype", "short_name", "protein_non_gap_length", "alignment_full_length"),
		"alignment completeness"
	)
	
	source_records |>
		dplyr::mutate(
			in_analysis = !is.na(.data$analysis_name),
			inclusion_status = dplyr::if_else(
				.data$in_analysis,
				"included_in_analysis",
				"excluded_from_analysis"
			),
			exclusion_reason = dplyr::if_else(
				.data$in_analysis,
				NA_character_,
				"absent_from_cartography_maps"
			),
			source_status = dplyr::if_else(
				is_missing_text(.data$protein_sequence_source),
				"missing_source",
				"source_recorded"
			),
			accession_status = dplyr::if_else(
				is_missing_text(.data$protein_sequence_accession_number),
				"missing_accession",
				"accession_recorded"
			),
			protein_sequence_status = dplyr::if_else(
				is_missing_text(.data$protein_sequence),
				"missing_protein_sequence",
				"protein_sequence_recorded"
			)
		) |>
		dplyr::left_join(
			alignment_completeness |>
				dplyr::select(
					"subtype",
					"short_name",
					"protein_non_gap_length",
					"alignment_full_length"
				),
			by = c("subtype", "short_name")
		) |>
		dplyr::transmute(
			.data$strain_type,
			.data$strain_subtype,
			.data$subtype,
			.data$strain_name,
			.data$short_name,
			.data$candidate_analysis_name,
			.data$analysis_name,
			.data$cartography_match_status,
			.data$inclusion_status,
			.data$exclusion_reason,
			.data$protein_sequence_source,
			.data$protein_sequence_accession_number,
			.data$source_status,
			.data$accession_status,
			.data$protein_sequence_status,
			source_full_length = .data$full_length,
			.data$protein_length,
			.data$protein_non_gap_length,
			.data$alignment_full_length
		)
}

make_pair_count_audit <- function(distance_table, clean_sequences) {
	check_required_columns(
		distance_table,
		c("subtype", "method", "Var1", "Var2"),
		"distance table"
	)
	check_required_columns(clean_sequences, c("subtype", "short_name"), "clean sequences")
	
	included_counts <- clean_sequences |>
		dplyr::count(.data$subtype, name = "included_strains") |>
		dplyr::mutate(
			expected_unique_pairs = as.integer(.data$included_strains * (.data$included_strains - 1L) / 2L)
		)
	
	distance_table |>
		dplyr::distinct(.data$subtype, .data$method, .data$Var1, .data$Var2) |>
		dplyr::count(.data$subtype, .data$method, name = "observed_unique_pairs") |>
		dplyr::left_join(included_counts, by = "subtype") |>
		dplyr::mutate(
			pair_count_complete = .data$observed_unique_pairs == .data$expected_unique_pairs
		) |>
		dplyr::select(
			"subtype",
			"method",
			"included_strains",
			"expected_unique_pairs",
			"observed_unique_pairs",
			"pair_count_complete"
		) |>
		dplyr::arrange(.data$subtype, .data$method)
}

collapse_short_names <- function(data) {
	values <- data$short_name
	if (length(values) == 0) {
		return(NA_character_)
	}
	paste(values, collapse = ", ")
}

make_count_row <- function(stage, category, subtype, n, detail = NA_character_) {
	tibble::tibble(
		stage = stage,
		category = category,
		subtype = subtype,
		n = as.integer(n),
		detail = detail
	)
}

make_strain_flow_summary <- function(strain_provenance, pair_counts) {
	check_required_columns(
		strain_provenance,
		c(
			"subtype", "short_name", "inclusion_status", "exclusion_reason",
			"source_status", "accession_status", "protein_sequence_status",
			"source_full_length", "alignment_full_length"
		),
		"strain provenance records"
	)
	check_required_columns(
		pair_counts,
		c("subtype", "method", "observed_unique_pairs"),
		"pair counts"
	)
	
	included <- strain_provenance |>
		dplyr::filter(.data$inclusion_status == "included_in_analysis")
	excluded_cartography <- strain_provenance |>
		dplyr::filter(.data$exclusion_reason == "absent_from_cartography_maps")
	source_non_full <- strain_provenance |>
		dplyr::filter(!.data$source_full_length)
	analysis_non_full <- included |>
		dplyr::filter(!.data$source_full_length)
	alignment_non_full <- included |>
		dplyr::filter(.data$alignment_full_length %in% FALSE)
	
	dplyr::bind_rows(
		make_count_row("Source sequence CSV", "raw_sequence_rows", NA_character_, nrow(strain_provenance)),
		strain_provenance |>
			dplyr::count(.data$subtype, name = "n") |>
			dplyr::transmute(
				stage = "Source sequence CSV",
				category = "raw_sequence_rows",
				subtype = .data$subtype,
				n = .data$n,
				detail = NA_character_
			),
		make_count_row(
			"Source sequence CSV",
			"missing_source_rows",
			NA_character_,
			sum(strain_provenance$source_status == "missing_source"),
			collapse_short_names(dplyr::filter(strain_provenance, .data$source_status == "missing_source"))
		),
		make_count_row(
			"Source sequence CSV",
			"missing_accession_rows",
			NA_character_,
			sum(strain_provenance$accession_status == "missing_accession"),
			collapse_short_names(dplyr::filter(strain_provenance, .data$accession_status == "missing_accession"))
		),
		make_count_row(
			"Source sequence CSV",
			"missing_protein_sequence_rows",
			NA_character_,
			sum(strain_provenance$protein_sequence_status == "missing_protein_sequence"),
			collapse_short_names(dplyr::filter(strain_provenance, .data$protein_sequence_status == "missing_protein_sequence"))
		),
		make_count_row("Analysis inclusion", "analysis_sequence_rows", NA_character_, nrow(included)),
		included |>
			dplyr::count(.data$subtype, name = "n") |>
			dplyr::transmute(
				stage = "Analysis inclusion",
				category = "analysis_sequence_rows",
				subtype = .data$subtype,
				n = .data$n,
				detail = NA_character_
			),
		make_count_row(
			"Analysis inclusion",
			"excluded_absent_from_cartography_maps",
			NA_character_,
			nrow(excluded_cartography),
			collapse_short_names(excluded_cartography)
		),
		make_count_row(
			"Sequence completeness",
			"source_non_full_length_rows",
			NA_character_,
			nrow(source_non_full),
			collapse_short_names(source_non_full)
		),
		make_count_row(
			"Sequence completeness",
			"analysis_non_full_length_rows",
			NA_character_,
			nrow(analysis_non_full),
			collapse_short_names(analysis_non_full)
		),
		make_count_row(
			"Sequence completeness",
			"alignment_non_full_length_rows",
			NA_character_,
			nrow(alignment_non_full),
			collapse_short_names(alignment_non_full)
		),
		pair_counts |>
			dplyr::transmute(
				stage = "Pairwise distance matrices",
				category = "unique_pairwise_comparisons",
				subtype = .data$subtype,
				n = .data$observed_unique_pairs,
				detail = .data$method
			)
	)
}

make_sequence_source_summary <- function(strain_provenance) {
	check_required_columns(
		strain_provenance,
		c("subtype", "protein_sequence_source", "inclusion_status"),
		"strain provenance records"
	)
	
	strain_provenance |>
		dplyr::mutate(
			protein_sequence_source = dplyr::if_else(
				is_missing_text(.data$protein_sequence_source),
				"Missing",
				.data$protein_sequence_source
			)
		) |>
		dplyr::count(
			.data$subtype,
			.data$inclusion_status,
			.data$protein_sequence_source,
			name = "n"
		) |>
		dplyr::arrange(.data$subtype, .data$inclusion_status, .data$protein_sequence_source)
}

make_strain_flow_table <- function(strain_flow_summary) {
	strain_flow_summary |>
		dplyr::mutate(
			subtype = dplyr::if_else(is.na(.data$subtype), "all", .data$subtype),
			detail = dplyr::if_else(is.na(.data$detail), "", .data$detail)
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "stage") |>
		flextable::valign(j = "stage", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Strain inclusion and cleaning audit summary.")
}

make_strain_accession_table <- function(strain_provenance) {
	strain_provenance |>
		dplyr::transmute(
			Subtype = subtype_display_label(.data$subtype),
			`Short name` = .data$short_name,
			`Strain name` = .data$strain_name,
			`Analysis name` = .data$analysis_name,
			`Sequence source` = .data$protein_sequence_source,
			Accession = .data$protein_sequence_accession_number,
			`Source full length` = .data$source_full_length,
			`Protein length` = .data$protein_length,
			`Aligned non-gap length` = .data$protein_non_gap_length,
			`Alignment full length` = .data$alignment_full_length,
			`Inclusion status` = .data$inclusion_status,
			`Exclusion reason` = .data$exclusion_reason
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Sequence source, accession, and analysis-inclusion provenance.")
}

make_pair_count_table <- function(pair_counts) {
	pair_counts |>
		dplyr::mutate(
			Subtype = subtype_display_label(.data$subtype),
			Method = unname(distance_method_labels()[.data$method]),
			`Pair count complete` = dplyr::if_else(.data$pair_count_complete, "yes", "no")
		) |>
		dplyr::select(
			"Subtype",
			"Method",
			"Included strains" = "included_strains",
			"Expected unique pairs" = "expected_unique_pairs",
			"Observed unique pairs" = "observed_unique_pairs",
			"Pair count complete"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Unique off-diagonal pair counts by subtype and distance metric.")
}

make_sequence_source_summary_table <- function(sequence_source_summary) {
	sequence_source_summary |>
		dplyr::mutate(Subtype = subtype_display_label(.data$subtype)) |>
		dplyr::select(
			"Subtype",
			"Inclusion status" = "inclusion_status",
			"Sequence source" = "protein_sequence_source",
			"Rows" = "n"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption("Sequence source counts by subtype and inclusion status.")
}

write_strain_flow_table <- function(strain_flow_summary, path) {
	make_strain_flow_table(strain_flow_summary) |>
		write_rds_target(path)
}

write_strain_accession_table <- function(strain_provenance, path) {
	make_strain_accession_table(strain_provenance) |>
		write_rds_target(path)
}

write_pair_count_table <- function(pair_counts, path) {
	make_pair_count_table(pair_counts) |>
		write_rds_target(path)
}

write_sequence_source_summary_table <- function(sequence_source_summary, path) {
	make_sequence_source_summary_table(sequence_source_summary) |>
		write_rds_target(path)
}

#### END OF FILE ####
