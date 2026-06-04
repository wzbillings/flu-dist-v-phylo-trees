###
# Sequence completeness and missing-data sensitivity audit
# Zane Billings
###

collapse_unique_values <- function(x) {
	x <- sort(unique(x[!is.na(x) & x != ""]))
	if (length(x) == 0) {
		return("")
	}
	paste(x, collapse = ", ")
}

make_sequence_character_audit_one <- function(alignment_result, subtype = alignment_result$subtype) {
	aligned_sequences <- alignment_result$aligned_sequences
	check_required_columns(
		aligned_sequences,
		c("short_name", "pro_aligned"),
		"aligned sequences"
	)

	purrr::map_dfr(
		seq_len(nrow(aligned_sequences)),
		\(row_index) {
			chars <- split_sequence_characters(aligned_sequences$pro_aligned[[row_index]])
			tibble::tibble(
				subtype = subtype,
				short_name = aligned_sequences$short_name[[row_index]],
				character = chars,
				character_upper = normalize_amino_acid_code(chars),
				residue_class = classify_amino_acid_code(chars)
			)
		}
	) |>
		dplyr::group_by(.data$subtype, .data$residue_class) |>
		dplyr::summarise(
			n = dplyr::n(),
			characters = collapse_unique_values(.data$character_upper),
			strains = collapse_unique_values(.data$short_name),
			.groups = "drop"
		) |>
		dplyr::arrange(.data$subtype, .data$residue_class)
}

make_sequence_character_audit <- function(alignments_by_subtype) {
	purrr::imap_dfr(
		alignments_by_subtype,
		\(alignment_result, subtype) make_sequence_character_audit_one(alignment_result, subtype)
	)
}

make_sequence_character_audit_table <- function(sequence_character_audit) {
	sequence_character_audit |>
		dplyr::mutate(
			Subtype = subtype_display_label(.data$subtype),
			`Residue class` = dplyr::recode(
				.data$residue_class,
				standard_residue = "Standard amino-acid residue",
				gap = "Alignment gap",
				ambiguous_or_missing = "Ambiguous or missing residue code",
				known_nonstandard_residue = "Known nonstandard residue code",
				unexpected = "Unexpected character"
			)
		) |>
		dplyr::select(
			"Subtype",
			"Residue class",
			"Characters" = "characters",
			"Aligned character count" = "n",
			"Strains" = "strains"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Audit of aligned sequence characters by residue class. Standard ",
			"residues are used in sequence-distance denominators; gaps, ambiguous ",
			"or missing codes, and known nonstandard residue codes are excluded. ",
			"Unexpected characters require source-data review."
		))
}

write_sequence_character_audit_table <- function(sequence_character_audit, path) {
	make_sequence_character_audit_table(sequence_character_audit) |>
		write_rds_target(path)
}

remove_strains_from_matrix <- function(mat, removed_strains) {
	validate_distance_matrix(mat)
	removed_strains <- intersect(removed_strains, rownames(mat))
	if (length(removed_strains) == 0) {
		return(mat)
	}
	keep <- setdiff(rownames(mat), removed_strains)
	mat[keep, keep, drop = FALSE]
}

safe_matrix_association_estimate <- function(mat1, mat2) {
	aligned <- align_distance_matrices(mat1, mat2)
	x <- lower_triangle_values(aligned$mat1)
	y <- lower_triangle_values(aligned$mat2)
	complete <- stats::complete.cases(x, y)
	n_pairs <- sum(complete)
	estimate <- if (
		n_pairs < 2 ||
			stats::sd(x[complete]) == 0 ||
			stats::sd(y[complete]) == 0
	) {
		NA_real_
	} else {
		stats::cor(x[complete], y[complete])
	}
	tibble::tibble(estimate = estimate, n_pairs = n_pairs)
}

non_full_length_strains <- function(strain_provenance, subtype) {
	check_required_columns(
		strain_provenance,
		c("subtype", "short_name", "inclusion_status", "source_full_length", "alignment_full_length"),
		"strain provenance records"
	)
	strain_provenance |>
		dplyr::filter(
			.data$subtype == !!subtype,
			.data$inclusion_status == "included_in_analysis",
			!(.data$source_full_length %in% TRUE & .data$alignment_full_length %in% TRUE)
		) |>
		dplyr::pull("short_name")
}

calculate_complete_sequence_matrix_sensitivity <- function(
		distances_with_tree_by_subtype,
		strain_provenance,
		methods = distance_method_labels(),
		threshold = make_analysis_settings()$influence_threshold
	) {
	validate_influence_threshold(threshold)
	purrr::imap_dfr(
		distances_with_tree_by_subtype,
		\(distances, subtype) {
			removed <- non_full_length_strains(strain_provenance, subtype)
			purrr::imap_dfr(
				methods,
				\(comparison_label, comparison_method) {
					if (!all(c("cophenetic", comparison_method) %in% names(distances))) {
						stop(
							"Distance set for ",
							subtype,
							" must include cophenetic and ",
							comparison_method,
							" matrices.",
							call. = FALSE
						)
					}
					full <- safe_matrix_association_estimate(
						distances$cophenetic,
						distances[[comparison_method]]
					)
					complete <- safe_matrix_association_estimate(
						remove_strains_from_matrix(distances$cophenetic, removed),
						remove_strains_from_matrix(distances[[comparison_method]], removed)
					)
					change <- complete$estimate - full$estimate
					absolute_change <- abs(change)
					tibble::tibble(
						method = comparison_method,
						Comparison = comparison_label,
						subtype = subtype,
						Subtype = subtype_display_label(subtype),
						`Removed non-full-length strains` = collapse_unique_values(removed),
						`Full pairwise comparisons` = full$n_pairs,
						`Complete-sequence pairwise comparisons` = complete$n_pairs,
						`Full estimate` = full$estimate,
						`Complete-sequence estimate` = complete$estimate,
						`Estimate change` = change,
						`Absolute change` = absolute_change,
						`Flag threshold` = threshold,
						`Sensitivity flag` = dplyr::case_when(
							is.na(absolute_change) ~ "not estimable",
							absolute_change >= threshold ~ "yes",
							TRUE ~ "no"
						)
					)
				}
			)
		}
	)
}

make_complete_sequence_sensitivity_table <- function(complete_sequence_sensitivity) {
	complete_sequence_sensitivity |>
		dplyr::mutate(
			`Removed non-full-length strains` = dplyr::if_else(
				.data$`Removed non-full-length strains` == "",
				"none",
				.data$`Removed non-full-length strains`
			),
			`Full estimate` = sprintf("%.2f", .data$`Full estimate`),
			`Complete-sequence estimate` = sprintf("%.2f", .data$`Complete-sequence estimate`),
			`Estimate change` = sprintf("%+.2f", .data$`Estimate change`),
			`Absolute change` = sprintf("%.2f", .data$`Absolute change`),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`)
		) |>
		dplyr::select(
			"Comparison",
			"Subtype",
			"Removed non-full-length strains",
			"Full pairwise comparisons",
			"Complete-sequence pairwise comparisons",
			"Full estimate",
			"Complete-sequence estimate",
			"Estimate change",
			"Absolute change",
			"Flag threshold",
			"Sensitivity flag"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Comparison", "Subtype")) |>
		flextable::valign(j = c("Comparison", "Subtype"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Complete-sequence-only matrix sensitivity. Estimates compare the ",
			"full current distance matrices with matrices after excluding included ",
			"strains that are not full-length by source or alignment audit. This ",
			"sensitivity does not refit alignments, ML trees, or cartography maps."
		))
}

write_complete_sequence_sensitivity_table <- function(complete_sequence_sensitivity, path) {
	make_complete_sequence_sensitivity_table(complete_sequence_sensitivity) |>
		write_rds_target(path)
}

candidate_sites_for_method <- function(aligned_proteins, method, subtype) {
	width <- unique(nchar(aligned_proteins))
	if (length(width) != 1) {
		stop("Aligned proteins must have a common width.", call. = FALSE)
	}
	if (identical(method, "grantham")) {
		return(seq_len(width[[1]]))
	}
	if (identical(method, "pepi")) {
		sites <- get_pepitope_sites(subtype, sites = c("a", "b", "c", "d", "e")) |>
			unlist(use.names = FALSE) |>
			unique() |>
			sort()
		return(sites[sites <= width[[1]]])
	}
	stop("Missing-data sensitivity is only defined for Grantham and p-epitope distances.", call. = FALSE)
}

complete_deletion_retained_sites <- function(aligned_proteins, sites) {
	residue_matrix <- sequence_residue_matrix(aligned_proteins)
	site_matrix <- residue_matrix[, sites, drop = FALSE]
	apply(site_matrix, 2, \(site) all(is_standard_amino_acid(site)))
}

complete_deletion_distance <- function(aligned_proteins, method, subtype) {
	switch(
		method,
		grantham = dist_grantham(aligned_proteins, deletion = "complete"),
		pepi = dist.pepi(aligned_proteins, subtype = subtype, deletion = "complete"),
		stop("Missing-data sensitivity is only defined for Grantham and p-epitope distances.", call. = FALSE)
	)
}

distance_matrix_difference_summary <- function(primary, sensitivity) {
	aligned <- align_distance_matrices(primary, sensitivity)
	x <- lower_triangle_values(aligned$mat1)
	y <- lower_triangle_values(aligned$mat2)
	complete <- stats::complete.cases(x, y)
	diff <- y[complete] - x[complete]
	distance_correlation <- if (
		sum(complete) < 2 ||
			stats::sd(x[complete]) == 0 ||
			stats::sd(y[complete]) == 0
	) {
		NA_real_
	} else {
		stats::cor(x[complete], y[complete])
	}
	tibble::tibble(
		`Compared pairwise distances` = sum(complete),
		`Distance-matrix correlation` = distance_correlation,
		`Mean absolute distance change` = if (length(diff) == 0) NA_real_ else mean(abs(diff)),
		`Maximum absolute distance change` = if (length(diff) == 0) NA_real_ else max(abs(diff))
	)
}

calculate_sequence_missing_data_sensitivity <- function(
		alignments_by_subtype,
		distances_with_tree_by_subtype,
		methods = c(grantham = "Grantham distance", pepi = "p-Epitope distance"),
		threshold = make_analysis_settings()$influence_threshold
	) {
	validate_influence_threshold(threshold)
	purrr::imap_dfr(
		methods,
		\(comparison_label, method) purrr::imap_dfr(
			alignments_by_subtype,
			\(alignment_result, subtype) {
				distances <- distances_with_tree_by_subtype[[subtype]]
				if (is.null(distances) || !all(c("cophenetic", method) %in% names(distances))) {
					stop(
						"Distance set for ",
						subtype,
						" must include cophenetic and ",
						method,
						" matrices.",
						call. = FALSE
					)
				}
				aligned_proteins <- extract_aligned_proteins(alignment_result)
				candidate_sites <- candidate_sites_for_method(aligned_proteins, method, subtype)
				retained_sites <- complete_deletion_retained_sites(aligned_proteins, candidate_sites)
				sensitivity <- complete_deletion_distance(aligned_proteins, method, subtype)
				primary_association <- safe_matrix_association_estimate(
					distances$cophenetic,
					distances[[method]]
				)
				sensitivity_association <- safe_matrix_association_estimate(
					distances$cophenetic,
					sensitivity
				)
				change <- sensitivity_association$estimate - primary_association$estimate
				absolute_change <- abs(change)

				tibble::tibble(
					method = method,
					Comparison = comparison_label,
					subtype = subtype,
					Subtype = subtype_display_label(subtype),
					`Primary deletion rule` = "pairwise",
					`Sensitivity deletion rule` = "complete",
					`Total candidate sites` = length(candidate_sites),
					`Complete-deletion retained sites` = sum(retained_sites),
					`Complete-deletion excluded sites` = length(candidate_sites) - sum(retained_sites),
					`Primary pairwise comparisons` = primary_association$n_pairs,
					`Sensitivity pairwise comparisons` = sensitivity_association$n_pairs,
					`Primary estimate` = primary_association$estimate,
					`Sensitivity estimate` = sensitivity_association$estimate,
					`Estimate change` = change,
					`Absolute estimate change` = absolute_change,
					`Flag threshold` = threshold,
					`Sensitivity flag` = dplyr::case_when(
						is.na(absolute_change) ~ "not estimable",
						absolute_change >= threshold ~ "yes",
						TRUE ~ "no"
					)
				) |>
					dplyr::bind_cols(distance_matrix_difference_summary(distances[[method]], sensitivity))
			}
		)
	)
}

make_sequence_missing_data_sensitivity_table <- function(sequence_missing_data_sensitivity) {
	sequence_missing_data_sensitivity |>
		dplyr::mutate(
			`Primary estimate` = sprintf("%.2f", .data$`Primary estimate`),
			`Sensitivity estimate` = sprintf("%.2f", .data$`Sensitivity estimate`),
			`Estimate change` = sprintf("%+.2f", .data$`Estimate change`),
			`Absolute estimate change` = sprintf("%.2f", .data$`Absolute estimate change`),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`),
			`Distance-matrix correlation` = sprintf("%.2f", .data$`Distance-matrix correlation`),
			`Mean absolute distance change` = sprintf("%.3f", .data$`Mean absolute distance change`),
			`Maximum absolute distance change` = sprintf("%.3f", .data$`Maximum absolute distance change`)
		) |>
		dplyr::select(
			"Comparison",
			"Subtype",
			"Primary deletion rule",
			"Sensitivity deletion rule",
			"Total candidate sites",
			"Complete-deletion retained sites",
			"Complete-deletion excluded sites",
			"Primary estimate",
			"Sensitivity estimate",
			"Estimate change",
			"Absolute estimate change",
			"Flag threshold",
			"Sensitivity flag",
			"Distance-matrix correlation",
			"Mean absolute distance change",
			"Maximum absolute distance change"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Comparison") |>
		flextable::valign(j = "Comparison", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Cheap missing-data sensitivity for sequence-derived distances. The ",
			"primary rule uses pairwise deletion of gaps, ambiguous or missing ",
			"codes, and known nonstandard residue codes. The sensitivity rule uses ",
			"complete deletion of non-comparable sites across all included strains ",
			"within subtype."
		))
}

write_sequence_missing_data_sensitivity_table <- function(sequence_missing_data_sensitivity, path) {
	make_sequence_missing_data_sensitivity_table(sequence_missing_data_sensitivity) |>
		write_rds_target(path)
}

#### END OF FILE ####
