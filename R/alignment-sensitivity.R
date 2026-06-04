###
# Alignment-method sensitivity helpers
# Zane Billings
###

alignment_sensitive_distance_labels <- function() {
	c(
		grantham = "Grantham distance",
		pepi = "p-Epitope distance"
	)
}

alignment_sensitivity_scopes <- function() {
	tibble::tibble(
		Scope = c("Overall", "H1N1", "H3N2"),
		subtype = c(NA_character_, "h1", "h3")
	)
}

calculate_alignment_sensitivity_alignments <- function(
		clean_sequences,
		settings = make_analysis_settings()
	) {
	methods <- alignment_sensitivity_methods(settings)
	stats::setNames(
		purrr::map(
			methods,
			\(method) {
				method_settings <- settings_with_alignment_method(settings, method)
				list(
					h1 = align_h1_sequences(clean_sequences, method_settings),
					h3 = align_h3_sequences(clean_sequences, method_settings)
				)
			}
		),
		methods
	)
}

make_alignment_sensitivity_alignment_summary <- function(alignments_by_method) {
	purrr::imap_dfr(
		alignments_by_method,
		\(alignments_by_subtype, alignment_method) purrr::imap_dfr(
			alignments_by_subtype,
			\(alignment_result, subtype) alignment_audit(alignment_result) |>
				dplyr::mutate(
					alignment_method = alignment_method,
					alignment_method_key = alignment_method_key(alignment_method),
					.before = 1
				)
		)
	)
}

calculate_alignment_sensitivity_distances <- function(
		alignments_by_method,
		h1_cartography_map,
		h3_cartography_map,
		virus_metadata
	) {
	cartography_maps <- list(h1 = h1_cartography_map, h3 = h3_cartography_map)
	purrr::imap(
		alignments_by_method,
		\(alignments_by_subtype, alignment_method) purrr::imap(
			alignments_by_subtype,
			\(alignment_result, subtype) calculate_subtype_distances(
				alignment_result,
				cartography_maps[[subtype]],
				virus_metadata
			)
		)
	)
}

calculate_alignment_sensitivity_model_tests <- function(
		alignments_by_method,
		settings = make_analysis_settings()
	) {
	purrr::map(
		alignments_by_method,
		\(alignments_by_subtype) purrr::map(
			alignments_by_subtype,
			\(alignment_result) calculate_model_test(alignment_result, settings)
		)
	)
}

add_primary_tree_distances <- function(alternative_distances_by_subtype, primary_distances_with_tree_by_subtype) {
	purrr::imap(
		alternative_distances_by_subtype,
		\(distances, subtype) {
			if (!"cophenetic" %in% names(primary_distances_with_tree_by_subtype[[subtype]])) {
				stop("Primary distance set for ", subtype, " must include cophenetic distances.", call. = FALSE)
			}
			c(
				list(cophenetic = primary_distances_with_tree_by_subtype[[subtype]]$cophenetic),
				distances
			)
		}
	)
}

safe_cophenetic_association_estimate <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		subtype = NULL
	) {
	selected <- if (is.null(subtype)) {
		distances_with_tree_by_subtype
	} else {
		distances_with_tree_by_subtype[subtype]
	}
	pair_data <- purrr::imap_dfr(
		selected,
		\(distances, current_subtype) {
			if (!all(c("cophenetic", comparison_method) %in% names(distances))) {
				stop(
					"Distance set for ",
					current_subtype,
					" must include cophenetic and ",
					comparison_method,
					" matrices.",
					call. = FALSE
				)
			}
			paired_lower_triangle_distances(
				distances$cophenetic,
				distances[[comparison_method]],
				subtype = current_subtype,
				x_name = "cophenetic",
				y_name = "distance"
			)
		}
	)
	complete <- stats::complete.cases(pair_data$cophenetic, pair_data$distance)
	n_pairs <- sum(complete)
	estimate <- if (
		n_pairs < 2 ||
			stats::sd(pair_data$cophenetic[complete]) == 0 ||
			stats::sd(pair_data$distance[complete]) == 0
	) {
		NA_real_
	} else {
		stats::cor(pair_data$cophenetic[complete], pair_data$distance[complete])
	}
	tibble::tibble(estimate = estimate, n_pairs = n_pairs)
}

calculate_alignment_distance_sensitivity <- function(
		primary_distances_with_tree_by_subtype,
		alternative_distances_by_method,
		settings = make_analysis_settings(),
		methods = alignment_sensitive_distance_labels(),
		threshold = settings$influence_threshold
	) {
	validate_influence_threshold(threshold)
	primary_alignment <- canonical_alignment_method(settings$alignment_method)
	scopes <- alignment_sensitivity_scopes()
	purrr::imap_dfr(
		alternative_distances_by_method,
		\(alternative_distances_by_subtype, alternative_alignment) {
			alternative_with_tree <- add_primary_tree_distances(
				alternative_distances_by_subtype,
				primary_distances_with_tree_by_subtype
			)
			purrr::imap_dfr(
				methods,
				\(comparison_label, comparison_method) purrr::pmap_dfr(
					scopes,
					\(Scope, subtype) {
						subtype_arg <- if (is.na(subtype)) NULL else subtype
						primary <- safe_cophenetic_association_estimate(
							primary_distances_with_tree_by_subtype,
							comparison_method,
							subtype = subtype_arg
						)
						alternative <- safe_cophenetic_association_estimate(
							alternative_with_tree,
							comparison_method,
							subtype = subtype_arg
						)
						change <- alternative$estimate - primary$estimate
						absolute_change <- abs(change)
						tibble::tibble(
							method = comparison_method,
							Comparison = comparison_label,
							Scope = Scope,
							`Primary alignment` = primary_alignment,
							`Alternative alignment` = alternative_alignment,
							`Primary pairwise comparisons` = primary$n_pairs,
							`Alternative pairwise comparisons` = alternative$n_pairs,
							`Primary estimate` = primary$estimate,
							`Alternative estimate` = alternative$estimate,
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
			)
		}
	)
}

make_alignment_distance_sensitivity_table <- function(alignment_distance_sensitivity) {
	alignment_distance_sensitivity |>
		dplyr::mutate(
			`Primary estimate` = sprintf("%.2f", .data$`Primary estimate`),
			`Alternative estimate` = sprintf("%.2f", .data$`Alternative estimate`),
			`Estimate change` = sprintf("%+.2f", .data$`Estimate change`),
			`Absolute change` = sprintf("%.2f", .data$`Absolute change`),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`)
		) |>
		dplyr::select(
			"Comparison",
			"Scope",
			"Primary alignment",
			"Alternative alignment",
			"Primary pairwise comparisons",
			"Alternative pairwise comparisons",
			"Primary estimate",
			"Alternative estimate",
			"Estimate change",
			"Absolute change",
			"Flag threshold",
			"Sensitivity flag"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Comparison", "Primary alignment", "Alternative alignment")) |>
		flextable::valign(j = c("Comparison", "Primary alignment", "Alternative alignment"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Matrix-only alignment sensitivity for sequence-derived distances. ",
			"Alternative alignments are compared with the primary MUSCLE workflow ",
			"using the primary ML-tree cophenetic distances; ML trees, distance ",
			"trees, and cartography maps are not refit."
		))
}

write_alignment_distance_sensitivity_table <- function(alignment_distance_sensitivity, path) {
	make_alignment_distance_sensitivity_table(alignment_distance_sensitivity) |>
		write_rds_target(path)
}

calculate_alignment_model_sensitivity <- function(
		primary_model_tests_by_subtype,
		alternative_model_tests_by_method,
		primary_model_choice,
		settings = make_analysis_settings(),
		performance_tolerance = settings$model_performance_tolerance
	) {
	purrr::imap_dfr(
		alternative_model_tests_by_method,
		\(alternative_model_tests_by_subtype, alternative_alignment) {
			alternative_model_choice <- choose_tree_model(
				alternative_model_tests_by_subtype,
				performance_tolerance = performance_tolerance
			)
			purrr::imap_dfr(
				alternative_model_tests_by_subtype,
				\(alternative_model_test, subtype) {
					primary_selected_model <- primary_model_choice$selected_models[[subtype]]
					alternative_selected_model <- alternative_model_choice$selected_models[[subtype]]
					loss_fraction <- model_log_likelihood_loss_fraction(
						alternative_model_test,
						primary_selected_model
					)
					tibble::tibble(
						subtype = subtype,
						Subtype = subtype_display_label(subtype),
						`Primary alignment` = canonical_alignment_method(settings$alignment_method),
						`Alternative alignment` = alternative_alignment,
						`Primary selected model` = primary_selected_model,
						`Primary AICc-best model` = best_model_by_aicc(primary_model_tests_by_subtype[[subtype]]),
						`Alternative selected model` = alternative_selected_model,
						`Alternative AICc-best model` = best_model_by_aicc(alternative_model_test),
						`Alternative selection strategy` = alternative_model_choice$strategy,
						`Selected model changed` = dplyr::if_else(
							identical(primary_selected_model, alternative_selected_model),
							"no",
							"yes"
						),
						`Primary selected model loss fraction` = loss_fraction,
						`Flag threshold` = performance_tolerance,
						`Model sensitivity flag` = dplyr::case_when(
							is.na(loss_fraction) ~ "not estimable",
							loss_fraction > performance_tolerance ~ "yes",
							TRUE ~ "no"
						)
					)
				}
			)
		}
	)
}

make_alignment_model_sensitivity_table <- function(alignment_model_sensitivity) {
	alignment_model_sensitivity |>
		dplyr::mutate(
			`Primary selected model loss fraction` = sprintf(
				"%.3f",
				.data$`Primary selected model loss fraction`
			),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`)
		) |>
		dplyr::select(
			"Subtype",
			"Primary alignment",
			"Alternative alignment",
			"Primary selected model",
			"Alternative selected model",
			"Alternative AICc-best model",
			"Alternative selection strategy",
			"Selected model changed",
			"Primary selected model loss fraction",
			"Flag threshold",
			"Model sensitivity flag"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Primary alignment", "Alternative alignment")) |>
		flextable::valign(j = c("Primary alignment", "Alternative alignment"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Alignment sensitivity for FLU-family model selection. The primary ",
			"selected model is flagged only when its log-likelihood loss under an ",
			"alternative alignment exceeds the model-performance tolerance."
		))
}

write_alignment_model_sensitivity_table <- function(alignment_model_sensitivity, path) {
	make_alignment_model_sensitivity_table(alignment_model_sensitivity) |>
		write_rds_target(path)
}

#### END OF FILE ####
