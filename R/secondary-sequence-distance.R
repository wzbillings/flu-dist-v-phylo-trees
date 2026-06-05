###
# Supplementary sequence-distance sensitivity metrics
# Zane Billings
###

secondary_sequence_distance_labels <- function() {
	c(
		hamming = "Hamming distance",
		osa = "Optimal string alignment distance",
		p_all_epitope = "p-all-epitope distance",
		blosum62 = "BLOSUM62 distance"
	)
}

secondary_sequence_distance_comparators <- function() {
	c(
		hamming = "grantham",
		osa = "grantham",
		p_all_epitope = "pepi",
		blosum62 = "grantham"
	)
}

comparable_residue_pairs <- function(seq_1, seq_2, site_mask = NULL) {
	validate_aligned_sequence_pair(seq_1, seq_2)
	residues_1 <- split_residues(seq_1)
	residues_2 <- split_residues(seq_2)
	comparable <- is_standard_amino_acid(residues_1) & is_standard_amino_acid(residues_2)
	if (!is.null(site_mask)) {
		if (length(site_mask) != length(residues_1)) {
			stop("Site mask must match the aligned sequence length.", call. = FALSE)
		}
		comparable <- comparable & site_mask
	}
	list(
		residues_1 = residues_1[comparable],
		residues_2 = residues_2[comparable],
		n_comparable = sum(comparable)
	)
}

validate_secondary_sequence_vector <- function(seqs, label = "secondary sequence") {
	if (!is.character(seqs)) {
		stop("seqs must be a character vector of aligned amino-acid sequences.", call. = FALSE)
	}
	if (is.null(names(seqs)) || any(is.na(names(seqs))) || any(names(seqs) == "")) {
		stop("seqs must have non-missing strain names.", call. = FALSE)
	}
	validate_unique_values(names(seqs), paste(label, "names"))
	invisible(seqs)
}

secondary_sequence_distance_matrix <- function(seqs, pair_distance, label, ...) {
	validate_secondary_sequence_vector(seqs, label)
	res <- matrix(
		NA_real_,
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	diag(res) <- 0
	
	if (length(seqs) > 1) {
		for (i in 2:nrow(res)) {
			for (j in 1:(i - 1)) {
				distance <- pair_distance(seqs[[i]], seqs[[j]], ...)
				res[[i, j]] <- distance
				res[[j, i]] <- distance
			}
		}
	}
	
	res
}

amino_acid_hamming_pair_distance <- function(seq_1, seq_2, site_mask = NULL) {
	pairs <- comparable_residue_pairs(seq_1, seq_2, site_mask = site_mask)
	if (pairs$n_comparable == 0) {
		return(NA_real_)
	}
	mean(pairs$residues_1 != pairs$residues_2)
}

dist_amino_acid_hamming <- function(seqs) {
	secondary_sequence_distance_matrix(
		seqs,
		amino_acid_hamming_pair_distance,
		label = "Hamming sequence"
	)
}

optimal_string_alignment_pair_distance <- function(seq_1, seq_2, site_mask = NULL) {
	pairs <- comparable_residue_pairs(seq_1, seq_2, site_mask = site_mask)
	if (pairs$n_comparable == 0) {
		return(NA_real_)
	}
	filtered_1 <- paste(pairs$residues_1, collapse = "")
	filtered_2 <- paste(pairs$residues_2, collapse = "")
	stringdist::stringdist(filtered_1, filtered_2, method = "osa") / pairs$n_comparable
}

dist_optimal_string_alignment <- function(seqs) {
	secondary_sequence_distance_matrix(
		seqs,
		optimal_string_alignment_pair_distance,
		label = "optimal string alignment sequence"
	)
}

p_all_epitope <- function(seq_1, seq_2, subtype, site_masks = NULL) {
	p_epi_sites <- get_pepitope_sites(subtype, sites = c("a", "b", "c", "d", "e"))
	epi_dists <- purrr::imap_dbl(
		p_epi_sites,
		\(current_sites, site_name) pepitope_epitope_distance(
			seq_1,
			seq_2,
			current_sites,
			site_mask = site_masks[[site_name]]
		)
	)
	if (all(is.na(epi_dists))) {
		return(NA_real_)
	}
	mean(epi_dists, na.rm = TRUE)
}

dist_p_all_epitope <- function(seqs, subtype, deletion = c("pairwise", "complete")) {
	deletion <- match.arg(deletion)
	validate_secondary_sequence_vector(seqs, label = "p-all-epitope sequence")
	site_masks <- if (identical(deletion, "complete")) {
		pepitope_complete_site_masks(seqs, subtype)
	} else {
		NULL
	}
	
	res <- matrix(
		NA_real_,
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	diag(res) <- 0
	
	if (length(seqs) > 1) {
		for (i in 2:nrow(res)) {
			for (j in 1:(i - 1)) {
				distance <- p_all_epitope(seqs[[i]], seqs[[j]], subtype, site_masks = site_masks)
				res[[i, j]] <- distance
				res[[j, i]] <- distance
			}
		}
	}
	
	res
}

blosum62_substitution_matrix <- function() {
	amino_acids <- c(
		"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
		"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
	)
	values <- c(
		4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0,
		-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3,
		-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
		-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3,
		0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
		-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
		-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
		0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3,
		-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3,
		-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3,
		-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1,
		-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2,
		-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1,
		-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1,
		-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2,
		1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
		0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0,
		-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3,
		-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1,
		0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4
	)
	matrix(values, nrow = length(amino_acids), dimnames = list(amino_acids, amino_acids), byrow = TRUE)
}

blosum62_pair_distance <- function(
		seq_1,
		seq_2,
		matrix = blosum62_substitution_matrix(),
		diag_scores = diag(matrix),
		site_mask = NULL
	) {
	pairs <- comparable_residue_pairs(seq_1, seq_2, site_mask = site_mask)
	if (pairs$n_comparable == 0) {
		return(NA_real_)
	}
	scores <- matrix[cbind(pairs$residues_1, pairs$residues_2)]
	distances <- diag_scores[pairs$residues_1] + diag_scores[pairs$residues_2] - 2 * scores
	mean(distances)
}

dist_blosum62 <- function(seqs) {
	substitution_matrix <- blosum62_substitution_matrix()
	diag_scores <- diag(substitution_matrix)
	secondary_sequence_distance_matrix(
		seqs,
		\(seq_1, seq_2) blosum62_pair_distance(
			seq_1,
			seq_2,
			matrix = substitution_matrix,
			diag_scores = diag_scores
		),
		label = "BLOSUM62 sequence"
	)
}

calculate_secondary_sequence_distances_one <- function(alignment_result, subtype = alignment_result$subtype) {
	aligned_proteins <- extract_aligned_proteins(alignment_result)
	distances <- list(
		hamming = dist_amino_acid_hamming(aligned_proteins),
		osa = dist_optimal_string_alignment(aligned_proteins),
		p_all_epitope = dist_p_all_epitope(aligned_proteins, subtype = subtype),
		blosum62 = dist_blosum62(aligned_proteins)
	)
	distances <- standardize_distance_order(
		distances,
		reference_names = names(aligned_proteins),
		subtype = subtype
	)
	validate_distance_set(distances, subtype = subtype)
	distances
}

calculate_secondary_sequence_distances <- function(alignments_by_subtype) {
	purrr::imap(
		alignments_by_subtype,
		\(alignment_result, subtype) calculate_secondary_sequence_distances_one(alignment_result, subtype)
	)
}

secondary_sequence_sensitivity_scopes <- function(distances_by_subtype) {
	subtype_names <- names(distances_by_subtype)
	subtype_rows <- tibble::tibble(
		Scope = subtype_display_label(subtype_names),
		subtype = subtype_names
	)
	if (all(c("h1", "h3") %in% subtype_names)) {
		return(dplyr::bind_rows(
			tibble::tibble(Scope = "Overall", subtype = NA_character_),
			subtype_rows
		))
	}
	subtype_rows
}

calculate_secondary_sequence_distance_sensitivity <- function(
		primary_distances_with_tree_by_subtype,
		secondary_distances_by_subtype,
		settings = make_analysis_settings(),
		methods = secondary_sequence_distance_labels(),
		comparators = secondary_sequence_distance_comparators(),
		threshold = settings$influence_threshold
	) {
	validate_influence_threshold(threshold)
	comparators <- comparators[names(methods)]
	if (any(is.na(comparators))) {
		stop("Each secondary sequence metric must have a primary comparator.", call. = FALSE)
	}
	primary_labels <- distance_method_labels()
	secondary_with_tree <- add_primary_tree_distances(
		secondary_distances_by_subtype,
		primary_distances_with_tree_by_subtype
	)
	scopes <- secondary_sequence_sensitivity_scopes(secondary_distances_by_subtype)
	
	purrr::imap_dfr(
		methods,
		\(comparison_label, comparison_method) {
			primary_method <- comparators[[comparison_method]]
			primary_label <- unname(primary_labels[[primary_method]])
			purrr::pmap_dfr(
				scopes,
				\(Scope, subtype) {
					subtype_arg <- if (is.na(subtype)) NULL else subtype
					primary <- safe_cophenetic_association_estimate(
						primary_distances_with_tree_by_subtype,
						primary_method,
						subtype = subtype_arg
					)
					secondary <- safe_cophenetic_association_estimate(
						secondary_with_tree,
						comparison_method,
						subtype = subtype_arg
					)
					change <- secondary$estimate - primary$estimate
					absolute_change <- abs(change)
					tibble::tibble(
						method = comparison_method,
						Comparison = comparison_label,
						Scope = Scope,
						`Primary comparator` = primary_label,
						`Primary pairwise comparisons` = primary$n_pairs,
						`Secondary pairwise comparisons` = secondary$n_pairs,
						`Primary estimate` = primary$estimate,
						`Secondary estimate` = secondary$estimate,
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

make_secondary_sequence_distance_sensitivity_table <- function(secondary_sequence_distance_sensitivity) {
	secondary_sequence_distance_sensitivity |>
		dplyr::mutate(
			`Primary estimate` = sprintf("%.2f", .data$`Primary estimate`),
			`Secondary estimate` = sprintf("%.2f", .data$`Secondary estimate`),
			`Estimate change` = sprintf("%+.2f", .data$`Estimate change`),
			`Absolute change` = sprintf("%.2f", .data$`Absolute change`),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`)
		) |>
		dplyr::select(
			"Comparison",
			"Scope",
			"Primary comparator",
			"Primary pairwise comparisons",
			"Secondary pairwise comparisons",
			"Primary estimate",
			"Secondary estimate",
			"Estimate change",
			"Absolute change",
			"Flag threshold",
			"Sensitivity flag"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Comparison") |>
		flextable::valign(j = "Comparison", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Supplementary sequence-distance sensitivity metrics. Secondary ",
			"metrics are compared with the primary ML-tree cophenetic distances ",
			"and with their primary sequence-distance comparator; they are not ",
			"used for primary interpretation or distance-tree construction."
		))
}

write_secondary_sequence_distance_sensitivity_table <- function(secondary_sequence_distance_sensitivity, path) {
	make_secondary_sequence_distance_sensitivity_table(secondary_sequence_distance_sensitivity) |>
		write_rds_target(path)
}

#### END OF FILE ####
