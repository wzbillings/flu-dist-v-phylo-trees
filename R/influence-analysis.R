###
# Leave-one-strain-out influence analysis for matrix associations
# Zane Billings
###

influence_estimate_labels <- function() {
	c(
		mantel = "Mantel r",
		pearson = "Descriptive Pearson r"
	)
}

validate_influence_threshold <- function(threshold) {
	if (
		!is.numeric(threshold) ||
			length(threshold) != 1 ||
			is.na(threshold) ||
			threshold < 0
	) {
		stop("Influence threshold must be a non-negative numeric scalar.", call. = FALSE)
	}
	invisible(threshold)
}

validate_influence_estimate_methods <- function(estimate_methods) {
	labels <- influence_estimate_labels()
	bad <- setdiff(estimate_methods, names(labels))
	if (length(bad) > 0) {
		stop(
			"Unknown influence estimate method(s): ",
			paste(bad, collapse = ", "),
			call. = FALSE
		)
	}
	invisible(estimate_methods)
}

influence_threshold_from_settings <- function(settings, threshold = NULL) {
	if (!is.null(threshold)) {
		validate_influence_threshold(threshold)
		return(threshold)
	}
	if (!is.null(settings$influence_threshold)) {
		validate_influence_threshold(settings$influence_threshold)
		return(settings$influence_threshold)
	}
	0.10
}

matrix_association_estimate <- function(mat1, mat2) {
	aligned <- align_distance_matrices(mat1, mat2)
	x <- lower_triangle_values(aligned$mat1)
	y <- lower_triangle_values(aligned$mat2)
	complete <- stats::complete.cases(x, y)
	n_pairs <- sum(complete)
	estimate <- if (n_pairs < 2) {
		NA_real_
	} else {
		stats::cor(x[complete], y[complete])
	}
	tibble::tibble(
		estimate = estimate,
		n_pairs = n_pairs
	)
}

remove_strain_from_matrix <- function(mat, removed_strain) {
	validate_distance_matrix(mat)
	if (!removed_strain %in% rownames(mat)) {
		stop("Removed strain is absent from distance matrix: ", removed_strain, call. = FALSE)
	}
	keep <- setdiff(rownames(mat), removed_strain)
	mat[keep, keep, drop = FALSE]
}

summarise_one_influence_estimate <- function(
		distances,
		subtype,
		comparison_method,
		comparison_label,
		estimate_method,
		threshold
	) {
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
	aligned <- align_distance_matrices(distances$cophenetic, distances[[comparison_method]])
	if (nrow(aligned$mat1) < 4) {
		stop("Leave-one-strain-out influence requires at least four shared strains.", call. = FALSE)
	}
	full <- matrix_association_estimate(aligned$mat1, aligned$mat2)
	estimate_label <- unname(influence_estimate_labels()[[estimate_method]])

	purrr::map_dfr(
		rownames(aligned$mat1),
		\(removed_strain) {
			loo <- matrix_association_estimate(
				remove_strain_from_matrix(aligned$mat1, removed_strain),
				remove_strain_from_matrix(aligned$mat2, removed_strain)
			)
			change <- loo$estimate - full$estimate
			absolute_change <- abs(change)
			tibble::tibble(
				method = comparison_method,
				Comparison = comparison_label,
				subtype = subtype,
				Subtype = subtype_display_label(subtype),
				`Removed strain` = removed_strain,
				`Estimate method` = estimate_label,
				`Full pairwise comparisons` = full$n_pairs,
				`Leave-one-out pairwise comparisons` = loo$n_pairs,
				`Full estimate` = full$estimate,
				`Leave-one-out estimate` = loo$estimate,
				`Estimate change` = change,
				`Absolute change` = absolute_change,
				`Flag threshold` = threshold,
				`Influence flag` = dplyr::case_when(
					is.na(absolute_change) ~ "not estimable",
					absolute_change >= threshold ~ "yes",
					TRUE ~ "no"
				)
			)
		}
	)
}

calculate_cophenetic_influence_summary <- function(
		distances_with_tree_by_subtype,
		settings = make_analysis_settings(),
		methods = distance_method_labels(),
		estimate_methods = c("mantel", "pearson"),
		threshold = NULL
	) {
	threshold <- influence_threshold_from_settings(settings, threshold)
	validate_influence_estimate_methods(estimate_methods)

	purrr::imap_dfr(
		distances_with_tree_by_subtype,
		\(distances, subtype) purrr::imap_dfr(
			methods,
			\(comparison_label, comparison_method) purrr::map_dfr(
				estimate_methods,
				\(estimate_method) summarise_one_influence_estimate(
					distances,
					subtype = subtype,
					comparison_method = comparison_method,
					comparison_label = comparison_label,
					estimate_method = estimate_method,
					threshold = threshold
				)
			)
		)
	) |>
		dplyr::mutate(
			Comparison = factor(.data$Comparison, levels = unname(methods)),
			Subtype = factor(.data$Subtype, levels = c("H1N1", "H3N2")),
			`Estimate method` = factor(.data$`Estimate method`, levels = unname(influence_estimate_labels()))
		) |>
		dplyr::arrange(
			.data$Comparison,
			.data$Subtype,
			.data$`Removed strain`,
			.data$`Estimate method`
		) |>
		dplyr::mutate(
			Comparison = as.character(.data$Comparison),
			Subtype = as.character(.data$Subtype),
			`Estimate method` = as.character(.data$`Estimate method`)
		)
}

make_cophenetic_influence_table <- function(
		cophenetic_influence_summary,
		estimate_method = "Mantel r"
	) {
	cophenetic_influence_summary |>
		dplyr::filter(.data$`Estimate method` == !!estimate_method) |>
		dplyr::mutate(
			`Full estimate` = sprintf("%.2f", .data$`Full estimate`),
			`Leave-one-out estimate` = sprintf("%.2f", .data$`Leave-one-out estimate`),
			`Estimate change` = sprintf("%+.2f", .data$`Estimate change`),
			`Absolute change` = sprintf("%.2f", .data$`Absolute change`),
			`Flag threshold` = sprintf("%.2f", .data$`Flag threshold`)
		) |>
		dplyr::select(
			"Comparison",
			"Subtype",
			"Removed strain",
			"Full pairwise comparisons",
			"Leave-one-out pairwise comparisons",
			"Full estimate",
			"Leave-one-out estimate",
			"Estimate change",
			"Absolute change",
			"Flag threshold",
			"Influence flag"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Comparison", "Subtype")) |>
		flextable::valign(j = c("Comparison", "Subtype"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Leave-one-strain-out influence on Mantel correlations between each ",
			"distance metric and ML-tree cophenetic distance. Changes are ",
			"leave-one-out estimates minus full-panel estimates; flagged rows meet ",
			"or exceed the prespecified absolute-change threshold."
		))
}

influence_plot_y_separator <- function() {
	"___facet___"
}

influence_plot_y_label <- function(x) {
	sub(paste0(influence_plot_y_separator(), ".*$"), "", x)
}

prepare_cophenetic_influence_plot_data <- function(
		cophenetic_influence_summary,
		estimate_method = "Mantel r",
		top_n = 6L
	) {
	top_n <- as.integer(top_n)
	if (length(top_n) != 1 || is.na(top_n) || top_n < 1) {
		stop("`top_n` must be a positive integer.", call. = FALSE)
	}
	plot_data <- cophenetic_influence_summary |>
		dplyr::filter(.data$`Estimate method` == !!estimate_method) |>
		dplyr::group_by(.data$Comparison, .data$Subtype) |>
		dplyr::slice_max(.data$`Absolute change`, n = top_n, with_ties = FALSE) |>
		dplyr::ungroup() |>
		dplyr::mutate(
			Comparison = factor(.data$Comparison, levels = rev(distance_method_labels())),
			`Influence flag` = factor(.data$`Influence flag`, levels = c("no", "yes", "not estimable")),
			removed_strain_facet = paste(
				.data$`Removed strain`,
				.data$Comparison,
				.data$Subtype,
				sep = influence_plot_y_separator()
			)
		)

	facet_levels <- plot_data |>
		dplyr::arrange(.data$Comparison, .data$Subtype, .data$`Absolute change`) |>
		dplyr::pull("removed_strain_facet")

	plot_data |>
		dplyr::mutate(
			removed_strain_facet = factor(.data$removed_strain_facet, levels = unique(facet_levels))
		)
}

plot_cophenetic_influence <- function(
		cophenetic_influence_summary,
		estimate_method = "Mantel r",
		top_n = 6L
	) {
	plot_data <- prepare_cophenetic_influence_plot_data(
		cophenetic_influence_summary,
		estimate_method = estimate_method,
		top_n = top_n
	)

	ggplot2::ggplot(
		plot_data,
		ggplot2::aes(
			x = .data$`Estimate change`,
			y = .data$removed_strain_facet,
			color = .data$`Influence flag`
		)
	) +
		ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
		ggplot2::geom_point(size = 1.8, alpha = 0.85) +
		ggplot2::facet_grid(
			rows = ggplot2::vars(.data$Comparison),
			cols = ggplot2::vars(.data$Subtype),
			scales = "free_y",
			space = "free_y"
		) +
		ggplot2::scale_color_manual(
			values = c(no = "#4D4D4D", yes = "#D55E00", `not estimable` = "#999999"),
			name = "Flagged"
		) +
		ggplot2::scale_y_discrete(labels = influence_plot_y_label) +
		ggplot2::labs(
			x = "Leave-one-out estimate change",
			y = "Removed strain"
		) +
		ggplot2::theme_bw(base_size = 10) +
		ggplot2::theme(
			legend.position = "bottom",
			panel.grid.minor = ggplot2::element_blank()
		)
}

write_cophenetic_influence_table <- function(cophenetic_influence_summary, path) {
	table <- make_cophenetic_influence_table(cophenetic_influence_summary)
	write_rds_target(table, path)
}

write_cophenetic_influence_plot <- function(cophenetic_influence_summary, path) {
	ensure_dir(path)
	plot <- plot_cophenetic_influence(cophenetic_influence_summary)
	ggplot2::ggsave(path, plot = plot, width = 9, height = 8, dpi = 300)
	path
}

#### END OF FILE ####
