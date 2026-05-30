###
# Subtype contrasts for distance-matrix association estimates
# Zane Billings
###

subtype_contrast_labels <- function() {
	c(h1 = "H1N1", h3 = "H3N2")
}

subtype_contrast_setting <- function(settings, name, fallback_name) {
	if (!is.null(settings[[name]])) {
		return(settings[[name]])
	}
	settings[[fallback_name]]
}

complete_pair_correlation <- function(data) {
	complete <- data |>
		dplyr::filter(stats::complete.cases(.data$cophenetic, .data$distance))
	tibble::tibble(
		n_pairs = nrow(complete),
		estimate = stats::cor(complete$cophenetic, complete$distance)
	)
}

subtype_pair_correlation <- function(pair_data, subtype) {
	pair_data |>
		dplyr::filter(.data$subtype == !!subtype) |>
		complete_pair_correlation()
}

weighted_subtype_correlation <- function(pair_data, strain_weights, subtype) {
	subtype_data <- pair_data |>
		dplyr::filter(.data$subtype == !!subtype)
	pair_weights <- strain_pair_weights(subtype_data, strain_weights)
	weighted_pearson_correlation(
		subtype_data$cophenetic,
		subtype_data$distance,
		pair_weights
	)
}

weighted_subtype_contrast <- function(pair_data, strain_weights) {
	weighted_subtype_correlation(pair_data, strain_weights, "h3") -
		weighted_subtype_correlation(pair_data, strain_weights, "h1")
}

bayesboot_subtype_contrast_ci <- function(pair_data, reps, conf_level, seed) {
	if (!requireNamespace("bayesboot", quietly = TRUE)) {
		stop(
			"The bayesboot package is required for subtype contrast intervals.",
			call. = FALSE
		)
	}
	alpha <- 1 - conf_level
	strain_units <- make_strain_bootstrap_units(pair_data)
	draws <- with_seed(seed, {
		bayesboot::bayesboot(
			strain_units,
			statistic = function(strains, weights) {
				strain_weights <- normalize_strain_weights_by_stratum(strains, weights)
				weighted_subtype_contrast(pair_data, strain_weights)
			},
			R = reps,
			R2 = reps,
			use.weights = TRUE,
			.progress = "none"
		)
	})
	values <- draws[[1]]
	quantiles <- stats::quantile(
		values,
		probs = c(alpha / 2, 1 - alpha / 2),
		na.rm = TRUE,
		names = FALSE
	)
	tibble::tibble(
		`CI lower` = as.numeric(quantiles[[1]]),
		`CI upper` = as.numeric(quantiles[[2]]),
		`CI method` = "Bayesian bootstrap contrast over strains"
	)
}

summarise_one_subtype_contrast <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		comparison_label,
		settings,
		seed_offset = 0L
	) {
	pair_data <- make_cophenetic_pair_table(
		distances_with_tree_by_subtype,
		comparison_method,
		subtypes = names(subtype_contrast_labels())
	)
	h1 <- subtype_pair_correlation(pair_data, "h1")
	h3 <- subtype_pair_correlation(pair_data, "h3")
	bootstrap_reps <- subtype_contrast_setting(
		settings,
		"subtype_contrast_bootstrap_reps",
		"mantel_bootstrap_reps"
	)
	bootstrap_seed <- settings$seed + 20000L + seed_offset
	ci <- bayesboot_subtype_contrast_ci(
		pair_data,
		reps = bootstrap_reps,
		conf_level = settings$correlation_ci_level,
		seed = bootstrap_seed
	)
	tibble::tibble(
		method = comparison_method,
		Comparison = comparison_label,
		`H1N1 pairwise comparisons` = h1$n_pairs,
		`H3N2 pairwise comparisons` = h3$n_pairs,
		`H1N1 Mantel r` = h1$estimate,
		`H3N2 Mantel r` = h3$estimate,
		`Difference (H3N2 - H1N1)` = h3$estimate - h1$estimate,
		`Bootstrap draws` = bootstrap_reps,
		`Bootstrap seed` = bootstrap_seed
	) |>
		dplyr::bind_cols(ci)
}

calculate_cophenetic_subtype_contrast_summary <- function(
		distances_with_tree_by_subtype,
		settings = make_analysis_settings(),
		methods = distance_method_labels()
	) {
	required <- names(subtype_contrast_labels())
	if (!all(required %in% names(distances_with_tree_by_subtype))) {
		stop(
			"Subtype contrast requires distance sets named: ",
			paste(required, collapse = ", "),
			call. = FALSE
		)
	}
	purrr::imap_dfr(
		methods,
		\(comparison_label, comparison_method) summarise_one_subtype_contrast(
			distances_with_tree_by_subtype,
			comparison_method = comparison_method,
			comparison_label = comparison_label,
			settings = settings,
			seed_offset = match(comparison_method, names(methods))
		)
	)
}

clamp_correlation_for_fisher_z <- function(x) {
	epsilon <- sqrt(.Machine$double.eps)
	pmin(pmax(x, -1 + epsilon), 1 - epsilon)
}

fisher_z_subtype_contrast <- function(h1_estimate, h3_estimate, n_h1, n_h3) {
	if (n_h1 <= 3 || n_h3 <= 3 || any(is.na(c(h1_estimate, h3_estimate)))) {
		return(tibble::tibble(statistic = NA_real_, p_value = NA_real_))
	}
	z_h1 <- atanh(clamp_correlation_for_fisher_z(h1_estimate))
	z_h3 <- atanh(clamp_correlation_for_fisher_z(h3_estimate))
	se <- sqrt(1 / (n_h1 - 3) + 1 / (n_h3 - 3))
	statistic <- (z_h3 - z_h1) / se
	tibble::tibble(
		statistic = statistic,
		p_value = 2 * stats::pnorm(-abs(statistic))
	)
}

permuted_subtype_mantel_estimate <- function(distances, comparison_method) {
	aligned <- align_distance_matrices(distances$cophenetic, distances[[comparison_method]])
	x <- lower_triangle_values(aligned$mat1)
	perm <- sample(seq_len(nrow(aligned$mat2)))
	permuted <- aligned$mat2[perm, perm, drop = FALSE]
	y_perm <- lower_triangle_values(permuted)
	complete <- stats::complete.cases(x, y_perm)
	stats::cor(x[complete], y_perm[complete])
}

permutation_subtype_contrast_sensitivity <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		observed_difference,
		permutations,
		seed
	) {
	null <- with_seed(seed, {
		replicate(permutations, {
			permuted_subtype_mantel_estimate(
				distances_with_tree_by_subtype$h3,
				comparison_method
			) -
				permuted_subtype_mantel_estimate(
					distances_with_tree_by_subtype$h1,
					comparison_method
				)
		})
	})
	tibble::tibble(
		statistic = observed_difference,
		p_value = (sum(abs(null) >= abs(observed_difference), na.rm = TRUE) + 1) /
			(permutations + 1)
	)
}

summarise_one_subtype_contrast_sensitivity <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		comparison_label,
		settings,
		seed_offset = 0L
	) {
	pair_data <- make_cophenetic_pair_table(
		distances_with_tree_by_subtype,
		comparison_method,
		subtypes = names(subtype_contrast_labels())
	)
	h1 <- subtype_pair_correlation(pair_data, "h1")
	h3 <- subtype_pair_correlation(pair_data, "h3")
	observed_difference <- h3$estimate - h1$estimate
	fisher <- fisher_z_subtype_contrast(
		h1$estimate,
		h3$estimate,
		h1$n_pairs,
		h3$n_pairs
	)
	permutations <- subtype_contrast_setting(
		settings,
		"subtype_contrast_permutations",
		"mantel_permutations"
	)
	permutation_seed <- settings$seed + 30000L + seed_offset
	permutation <- permutation_subtype_contrast_sensitivity(
		distances_with_tree_by_subtype,
		comparison_method = comparison_method,
		observed_difference = observed_difference,
		permutations = permutations,
		seed = permutation_seed
	)
	dplyr::bind_rows(
		tibble::tibble(
			method = comparison_method,
			Comparison = comparison_label,
			`Sensitivity method` = "Fisher z approximation",
			`Difference (H3N2 - H1N1)` = observed_difference,
			`Sensitivity statistic` = fisher$statistic,
			`Sensitivity p-value` = fisher$p_value,
			Permutations = NA_integer_,
			`Sensitivity seed` = NA_integer_,
			Interpretation = paste0(
				"Naive independent-correlation approximation using unique pair counts; ",
				"included only as a sensitivity check because pairwise distances share strains."
			)
		),
		tibble::tibble(
			method = comparison_method,
			Comparison = comparison_label,
			`Sensitivity method` = "Independent Mantel permutation contrast",
			`Difference (H3N2 - H1N1)` = observed_difference,
			`Sensitivity statistic` = permutation$statistic,
			`Sensitivity p-value` = permutation$p_value,
			Permutations = permutations,
			`Sensitivity seed` = permutation_seed,
			Interpretation = paste0(
				"Independent within-subtype matrix-label permutations; this tests the ",
				"contrast against a joint no-association null rather than a direct ",
				"equality-of-subtype-associations null."
			)
		)
	)
}

calculate_cophenetic_subtype_contrast_sensitivity <- function(
		distances_with_tree_by_subtype,
		settings = make_analysis_settings(),
		methods = distance_method_labels()
	) {
	purrr::imap_dfr(
		methods,
		\(comparison_label, comparison_method) summarise_one_subtype_contrast_sensitivity(
			distances_with_tree_by_subtype,
			comparison_method = comparison_method,
			comparison_label = comparison_label,
			settings = settings,
			seed_offset = match(comparison_method, names(methods))
		)
	)
}

make_subtype_contrast_table <- function(subtype_contrast_summary) {
	subtype_contrast_summary |>
		dplyr::mutate(
			`H1N1 Mantel r` = sprintf("%.2f", .data$`H1N1 Mantel r`),
			`H3N2 Mantel r` = sprintf("%.2f", .data$`H3N2 Mantel r`),
			`Difference (H3N2 - H1N1)` = sprintf("%.2f", .data$`Difference (H3N2 - H1N1)`),
			`95% CI` = paste0(
				"(",
				sprintf("%.2f", .data$`CI lower`),
				", ",
				sprintf("%.2f", .data$`CI upper`),
				")"
			)
		) |>
		dplyr::select(
			"Comparison",
			"H1N1 pairwise comparisons",
			"H3N2 pairwise comparisons",
			"H1N1 Mantel r",
			"H3N2 Mantel r",
			"Difference (H3N2 - H1N1)",
			"95% CI",
			"Bootstrap draws"
		) |>
		flextable::flextable() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Subtype contrasts in Mantel correlation between each distance metric ",
			"and ML-tree cophenetic distance. Differences are H3N2 minus H1N1; ",
			"intervals are Bayesian bootstrap intervals over strain units."
		))
}

format_optional_p_value <- function(p_value) {
	dplyr::case_when(
		is.na(p_value) ~ "NA",
		p_value < 0.001 ~ "< 0.001",
		TRUE ~ sprintf("%.3f", p_value)
	)
}

make_subtype_contrast_sensitivity_table <- function(subtype_contrast_sensitivity) {
	subtype_contrast_sensitivity |>
		dplyr::mutate(
			`Difference (H3N2 - H1N1)` = sprintf("%.2f", .data$`Difference (H3N2 - H1N1)`),
			`Sensitivity statistic` = sprintf("%.2f", .data$`Sensitivity statistic`),
			`Sensitivity p-value` = format_optional_p_value(.data$`Sensitivity p-value`)
		) |>
		dplyr::select(
			"Comparison",
			"Sensitivity method",
			"Difference (H3N2 - H1N1)",
			"Sensitivity statistic",
			"Sensitivity p-value",
			"Permutations"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Comparison") |>
		flextable::valign(j = "Comparison", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Sensitivity summaries for subtype contrasts in Mantel correlation. ",
			"These alternatives are not the primary contrast estimand."
		))
}

plot_subtype_contrasts <- function(subtype_contrast_summary) {
	plot_data <- subtype_contrast_summary |>
		dplyr::mutate(
			Comparison = factor(.data$Comparison, levels = rev(distance_method_labels()))
		)
	ggplot2::ggplot(plot_data) +
		ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
		ggplot2::geom_segment(
			ggplot2::aes(
				x = .data$`CI lower`,
				xend = .data$`CI upper`,
				y = .data$Comparison,
				yend = .data$Comparison
			),
			linewidth = 0.8,
			color = "#4D4D4D"
		) +
		ggplot2::geom_point(
			ggplot2::aes(
				x = .data$`Difference (H3N2 - H1N1)`,
				y = .data$Comparison
			),
			size = 3,
			color = "#0072B2"
		) +
		ggplot2::labs(
			x = "Mantel r difference (H3N2 - H1N1)",
			y = NULL
		) +
		ggplot2::theme_bw(base_size = 11) +
		ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

write_subtype_contrast_table <- function(subtype_contrast_summary, path) {
	table <- make_subtype_contrast_table(subtype_contrast_summary)
	write_rds_target(table, path)
}

write_subtype_contrast_sensitivity_table <- function(subtype_contrast_sensitivity, path) {
	table <- make_subtype_contrast_sensitivity_table(subtype_contrast_sensitivity)
	write_rds_target(table, path)
}

write_subtype_contrast_plot <- function(subtype_contrast_summary, path) {
	ensure_dir(path)
	plot <- plot_subtype_contrasts(subtype_contrast_summary)
	ggplot2::ggsave(path, plot = plot, width = 7, height = 4, dpi = 300)
	path
}

#### END OF FILE ####
