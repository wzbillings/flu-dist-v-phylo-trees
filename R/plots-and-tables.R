###
# Manuscript figures and tables
# Zane Billings
###

make_full_distance_table <- function(distances_by_subtype, tree_analyses) {
	distances_with_tree <- purrr::imap(
		distances_by_subtype,
		\(distances, subtype) c(
			list(cophenetic = tree_analyses[[subtype]]$ml_tree_distances),
			distances
		)
	)
	combine_distance_tables(distances_with_tree, unique_pairs = TRUE)
}

subtype_color_palette <- function() {
	c(
		H1N1 = "#0072B2",
		H3N2 = "#D55E00"
	)
}

distance_method_labels <- function() {
	c(
		year = "Temporal distance",
		grantham = "Grantham distance",
		pepi = "p-Epitope distance",
		cart = "Cartographic distance"
	)
}

subtype_display_label <- function(subtype) {
	dplyr::if_else(
		subtype == "h1",
		"H1N1",
		dplyr::if_else(subtype == "h3", "H3N2", subtype)
	)
}

normalize_distances_for_plot <- function(full_distance_table) {
	full_distance_table |>
		dplyr::group_by(.data$method) |>
		dplyr::mutate(d_normalized = normalize_vector(.data$d)) |>
		dplyr::ungroup()
}

plot_tree_distance_correlations <- function(
		full_distance_table,
		settings = make_analysis_settings()
	) {
	method_labels <- distance_method_labels()
	
	plot_df <- full_distance_table |>
		normalize_distances_for_plot() |>
		dplyr::mutate(method = as.character(.data$method)) |>
		tidyr::pivot_wider(
			id_cols = c("subtype", "Var1", "Var2"),
			names_from = "method",
			values_from = "d_normalized"
		) |>
		tidyr::pivot_longer(
			cols = dplyr::all_of(names(method_labels)),
			names_to = "method",
			values_to = "distance"
		) |>
		dplyr::mutate(
			method_label = factor(.data$method, levels = names(method_labels), labels = method_labels),
			subtype_label = subtype_display_label(.data$subtype)
		)
	
	ggplot2::ggplot(
		plot_df,
		ggplot2::aes(
			x = .data$cophenetic,
			y = .data$distance,
			color = .data$subtype_label
		)
	) +
		ggplot2::geom_point(
			size = 2,
			alpha = 0.65,
			position = ggplot2::position_jitter(
				width = 0.01,
				height = 0.01,
				seed = settings$jitter_seed
			)
		) +
		ggplot2::facet_wrap(ggplot2::vars(.data$method_label), nrow = 2) +
		ggplot2::scale_color_manual(values = subtype_color_palette(), name = "Subtype") +
		ggplot2::labs(
			x = "Normalized ML-tree cophenetic distance",
			y = "Normalized distance"
		) +
		ggplot2::theme_bw(base_size = 11) +
		ggplot2::theme(
			legend.position = "bottom",
			strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
			panel.grid.minor = ggplot2::element_blank()
		)
}

write_correlation_plot <- function(
		full_distance_table,
		path,
		settings = make_analysis_settings()
	) {
	ensure_dir(path)
	plot <- plot_tree_distance_correlations(full_distance_table, settings)
	ggplot2::ggsave(path, plot = plot, width = 9, height = 7.5, dpi = 300)
	path
}

make_cophenetic_pair_table <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		subtypes = names(distances_with_tree_by_subtype)
	) {
	selected <- distances_with_tree_by_subtype[subtypes]
	purrr::imap_dfr(
		selected,
		\(distances, subtype) {
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
			paired_lower_triangle_distances(
				distances$cophenetic,
				distances[[comparison_method]],
				subtype = subtype,
				x_name = "cophenetic",
				y_name = "distance"
			)
		}
	)
}

summarise_cophenetic_mantel <- function(
		distances_with_tree_by_subtype,
		comparison_method,
		comparison_label,
		scope,
		settings,
		seed_offset = 0L,
		subtype = NULL
	) {
	perm_seed <- settings$seed + seed_offset
	ci_seed <- settings$seed + 10000L + seed_offset

	if (is.null(subtype)) {
		selected <- distances_with_tree_by_subtype
		mantel <- stratified_mantel_permutation_test(
			selected,
			"cophenetic",
			comparison_method,
			permutations = settings$mantel_permutations,
			seed = perm_seed,
			method = "pearson"
		)
		pair_data <- make_cophenetic_pair_table(selected, comparison_method)
	} else {
		selected <- distances_with_tree_by_subtype[subtype]
		mantel <- mantel_permutation_test(
			selected[[subtype]]$cophenetic,
			selected[[subtype]][[comparison_method]],
			permutations = settings$mantel_permutations,
			seed = perm_seed,
			method = "pearson"
		)
		pair_data <- make_cophenetic_pair_table(selected, comparison_method)
	}

	ci <- correlation_ci(
		pair_data,
		estimate = mantel$estimate,
		reps = settings$mantel_bootstrap_reps,
		conf_level = settings$correlation_ci_level,
		seed = ci_seed
	)

	tibble::tibble(
		method = comparison_method,
		Comparison = comparison_label,
		Scope = scope,
		`Pairwise comparisons` = mantel$n_pairs,
		`Mantel r` = mantel$estimate,
		`CI lower` = ci$ci_lower,
		`CI upper` = ci$ci_upper,
		`Permutation p-value` = mantel$p_value,
		Permutations = mantel$permutations,
		`Bootstrap draws` = settings$mantel_bootstrap_reps,
		`CI method` = ci$ci_method,
		`Permutation seed` = perm_seed,
		`Bootstrap seed` = ci_seed
	)
}

calculate_cophenetic_mantel_summary <- function(
		distances_with_tree_by_subtype,
		settings = make_analysis_settings()
	) {
	method_labels <- distance_method_labels()
	scopes <- tibble::tibble(
		Scope = c("Overall", "H1N1", "H3N2"),
		subtype = c(NA_character_, "h1", "h3"),
		scope_index = seq_len(3)
	)

	purrr::imap_dfr(
		method_labels,
		\(comparison_label, comparison_method) {
			method_index <- match(comparison_method, names(method_labels))
			purrr::pmap_dfr(
				scopes,
				\(Scope, subtype, scope_index) summarise_cophenetic_mantel(
					distances_with_tree_by_subtype,
					comparison_method = comparison_method,
					comparison_label = comparison_label,
					scope = Scope,
					settings = settings,
					seed_offset = 1000L * scope_index + method_index,
					subtype = if (is.na(subtype)) NULL else subtype
				)
			)
		}
	) |>
		dplyr::mutate(
			Comparison = factor(.data$Comparison, levels = unname(method_labels)),
			Scope = factor(.data$Scope, levels = c("Overall", "H1N1", "H3N2"))
		) |>
		dplyr::arrange(.data$Comparison, .data$Scope) |>
		dplyr::mutate(
			Comparison = as.character(.data$Comparison),
			Scope = as.character(.data$Scope)
		)
}

format_p_value <- function(p_value) {
	dplyr::case_when(
		is.na(p_value) ~ "NA",
		p_value < 0.001 ~ "< 0.001",
		TRUE ~ sprintf("%.3f", p_value)
	)
}

make_cophenetic_mantel_table <- function(cophenetic_mantel_summary) {
	cophenetic_mantel_summary |>
		dplyr::mutate(
			`Mantel r` = sprintf("%.2f", .data$`Mantel r`),
			`95% CI` = paste0(
				"(",
				sprintf("%.2f", .data$`CI lower`),
				", ",
				sprintf("%.2f", .data$`CI upper`),
				")"
			),
			`Permutation p-value` = format_p_value(.data$`Permutation p-value`)
		) |>
		dplyr::select(
			"Comparison",
			"Scope",
			"Pairwise comparisons",
			"Mantel r",
			"95% CI",
			"Permutation p-value",
			"Permutations"
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Comparison") |>
		flextable::valign(j = "Comparison", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Mantel correlations between each distance metric and ML-tree ",
			"cophenetic distance, using unique off-diagonal strain pairs. ",
			"Intervals are 95% Bayesian bootstrap intervals over sampled strains; ",
			"p-values are two-sided matrix-permutation p-values."
		))
}

write_cophenetic_mantel_table <- function(cophenetic_mantel_summary, path) {
	table <- make_cophenetic_mantel_table(cophenetic_mantel_summary)
	write_rds_target(table, path)
}

weighted_pearson_correlation <- function(x, y, weights = NULL) {
	complete <- if (is.null(weights)) {
		stats::complete.cases(x, y)
	} else {
		stats::complete.cases(x, y, weights)
	}
	x <- x[complete]
	y <- y[complete]
	if (is.null(weights)) {
		return(stats::cor(x, y))
	}
	weights <- weights[complete]
	if (sum(weights) <= 0) {
		return(NA_real_)
	}
	weights <- weights / sum(weights)
	x_centered <- x - sum(weights * x)
	y_centered <- y - sum(weights * y)
	cov_xy <- sum(weights * x_centered * y_centered)
	var_x <- sum(weights * x_centered^2)
	var_y <- sum(weights * y_centered^2)
	if (var_x <= 0 || var_y <= 0) {
		return(NA_real_)
	}
	cov_xy / sqrt(var_x * var_y)
}

wald_correlation_ci <- function(estimate, n, conf_level = 0.95) {
	alpha <- 1 - conf_level
	if (is.na(estimate) || n <= 3 || abs(estimate) >= 1) {
		return(c(lower = NA_real_, upper = NA_real_))
	}
	z <- stats::qnorm(1 - alpha / 2)
	se <- 1 / sqrt(n - 3)
	bounds <- tanh(atanh(estimate) + c(-1, 1) * z * se)
	stats::setNames(bounds, c("lower", "upper"))
}

pair_strain_ids <- function(data, strata_col = "subtype") {
	check_required_columns(data, c("Var1", "Var2"), "correlation pair data")
	stratum <- if (strata_col %in% names(data)) {
		as.character(data[[strata_col]])
	} else {
		rep("all", nrow(data))
	}
	tibble::tibble(
		strain1_id = paste(stratum, as.character(data$Var1), sep = "::"),
		strain2_id = paste(stratum, as.character(data$Var2), sep = "::")
	)
}

make_strain_bootstrap_units <- function(data, strata_col = "subtype") {
	check_required_columns(data, c("Var1", "Var2"), "correlation pair data")
	stratum <- if (strata_col %in% names(data)) {
		as.character(data[[strata_col]])
	} else {
		rep("all", nrow(data))
	}
	units <- dplyr::bind_rows(
		tibble::tibble(stratum = stratum, strain = as.character(data$Var1)),
		tibble::tibble(stratum = stratum, strain = as.character(data$Var2))
	) |>
		dplyr::distinct() |>
		dplyr::arrange(.data$stratum, .data$strain) |>
		dplyr::mutate(strain_id = paste(.data$stratum, .data$strain, sep = "::")) |>
		dplyr::select("strain_id", "stratum", "strain")
	
	validate_unique_values(units$strain_id, "strain bootstrap units")
	units
}

normalize_strain_weights_by_stratum <- function(strain_units, weights) {
	if (length(weights) != nrow(strain_units)) {
		stop("Strain weights must match the number of strain bootstrap units.", call. = FALSE)
	}
	out <- strain_units |>
		dplyr::mutate(weight = as.numeric(weights)) |>
		dplyr::group_by(.data$stratum) |>
		dplyr::mutate(
			stratum_total = sum(.data$weight),
			weight = dplyr::if_else(
				.data$stratum_total > 0,
				.data$weight / .data$stratum_total,
				0
			)
		) |>
		dplyr::ungroup()
	stats::setNames(out$weight, out$strain_id)
}

strain_pair_weights <- function(
		data,
		strain_weights,
		strata_col = "subtype",
		preserve_stratum_pair_weight = FALSE
	) {
	ids <- pair_strain_ids(data, strata_col)
	endpoint1 <- ids$strain1_id
	endpoint2 <- ids$strain2_id
	
	if (!all(c(endpoint1, endpoint2) %in% names(strain_weights))) {
		raw_endpoint1 <- as.character(data$Var1)
		raw_endpoint2 <- as.character(data$Var2)
		if (!all(c(raw_endpoint1, raw_endpoint2) %in% names(strain_weights))) {
			missing <- setdiff(unique(c(endpoint1, endpoint2)), names(strain_weights))
			stop(
				"Strain weights are missing endpoint(s): ",
				paste(missing, collapse = ", "),
				call. = FALSE
			)
		}
		endpoint1 <- raw_endpoint1
		endpoint2 <- raw_endpoint2
	}
	
	pair_weights <- unname(strain_weights[endpoint1] * strain_weights[endpoint2])
	
	if (isTRUE(preserve_stratum_pair_weight) && strata_col %in% names(data)) {
		pair_weights <- tibble::tibble(
			stratum = as.character(data[[strata_col]]),
			weight = pair_weights
		) |>
			dplyr::group_by(.data$stratum) |>
			dplyr::mutate(
				stratum_total = sum(.data$weight),
				weight = dplyr::if_else(
					.data$stratum_total > 0,
					.data$weight / .data$stratum_total * dplyr::n(),
					0
				)
			) |>
			dplyr::ungroup() |>
			dplyr::pull("weight")
	}
	
	pair_weights
}

strain_weighted_correlation <- function(data, strain_weights) {
	pair_weights <- strain_pair_weights(
		data,
		strain_weights,
		preserve_stratum_pair_weight = TRUE
	)
	weighted_pearson_correlation(data$cophenetic, data$distance, pair_weights)
}

bayesboot_correlation_ci <- function(data, reps, conf_level, seed) {
	if (!requireNamespace("bayesboot", quietly = TRUE)) {
		return(NULL)
	}
	if (!all(c("Var1", "Var2") %in% names(data))) {
		return(NULL)
	}
	alpha <- 1 - conf_level
	strain_units <- make_strain_bootstrap_units(data)
	draws <- with_seed(seed, {
		bayesboot::bayesboot(
			strain_units,
			statistic = function(strains, weights) {
				strain_weights <- normalize_strain_weights_by_stratum(strains, weights)
				strain_weighted_correlation(data, strain_weights)
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
	c(
		lower = quantiles[[1]],
		upper = quantiles[[2]],
		method = "Bayesian bootstrap over strains"
	)
}

boot_bca_correlation_ci <- function(data, reps, conf_level, seed) {
	if (!requireNamespace("boot", quietly = TRUE)) {
		return(NULL)
	}
	if (!all(c("Var1", "Var2") %in% names(data))) {
		return(NULL)
	}
	strain_units <- make_strain_bootstrap_units(data)
	boot_result <- with_seed(seed, {
		boot::boot(
			strain_units,
			statistic = function(strains, indices) {
				counts <- tabulate(indices, nbins = nrow(strains))
				strain_weights <- normalize_strain_weights_by_stratum(strains, counts)
				strain_weighted_correlation(data, strain_weights)
			},
			R = reps,
			strata = strain_units$stratum
		)
	})
	ci <- tryCatch(
		boot::boot.ci(
			boot_result,
			conf = conf_level,
			type = "bca"
		),
		error = function(e) NULL
	)
	if (is.null(ci) || is.null(ci$bca)) {
		return(NULL)
	}
	c(lower = ci$bca[4], upper = ci$bca[5], method = "BCa strain bootstrap")
}

correlation_ci <- function(data, estimate, reps, conf_level, seed) {
	data <- data |>
		dplyr::select(
			dplyr::any_of(c("subtype", "Var1", "Var2")),
			dplyr::all_of(c("cophenetic", "distance"))
		) |>
		dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::everything())))
	
	bayesboot_ci <- bayesboot_correlation_ci(data, reps, conf_level, seed)
	if (!is.null(bayesboot_ci)) {
		return(tibble::tibble(
			ci_lower = as.numeric(bayesboot_ci[["lower"]]),
			ci_upper = as.numeric(bayesboot_ci[["upper"]]),
			ci_method = unname(bayesboot_ci[["method"]])
		))
	}
	
	bca_ci <- boot_bca_correlation_ci(data, reps, conf_level, seed)
	if (!is.null(bca_ci)) {
		return(tibble::tibble(
			ci_lower = as.numeric(bca_ci[["lower"]]),
			ci_upper = as.numeric(bca_ci[["upper"]]),
			ci_method = unname(bca_ci[["method"]])
		))
	}
	
	wald_ci <- wald_correlation_ci(estimate, nrow(data), conf_level)
	tibble::tibble(
		ci_lower = as.numeric(wald_ci[["lower"]]),
		ci_upper = as.numeric(wald_ci[["upper"]]),
		ci_method = "Fisher z / Wald"
	)
}

summarise_cophenetic_correlation <- function(data, scope, settings, seed_offset = 0L) {
	complete <- data |>
		dplyr::filter(stats::complete.cases(.data$cophenetic, .data$distance))
	estimate <- stats::cor(complete$cophenetic, complete$distance)
	ci <- correlation_ci(
		complete,
		estimate = estimate,
		reps = settings$correlation_bootstrap_reps,
		conf_level = settings$correlation_ci_level,
		seed = settings$seed + seed_offset
	)
	tibble::tibble(
		scope = scope,
		n_pairs = nrow(complete),
		pearson_r = estimate
	) |>
		dplyr::bind_cols(ci)
}

calculate_cophenetic_correlation_summary <- function(
		full_distance_table,
		settings = make_analysis_settings()
	) {
	method_labels <- distance_method_labels()
	
	long_table <- full_distance_table |>
		dplyr::mutate(method = as.character(.data$method)) |>
		tidyr::pivot_wider(
			id_cols = c("subtype", "Var1", "Var2"),
			names_from = "method",
			values_from = "d"
		) |>
		tidyr::pivot_longer(
			cols = dplyr::all_of(names(method_labels)),
			names_to = "method",
			values_to = "distance"
		) |>
		dplyr::mutate(
			comparison = unname(method_labels[.data$method]),
			subtype_label = subtype_display_label(.data$subtype)
		)
	
	overall <- long_table |>
		dplyr::group_by(.data$method, .data$comparison) |>
		dplyr::group_modify(
			\(data, key) summarise_cophenetic_correlation(
				data,
				scope = "Overall",
				settings = settings,
				seed_offset = match(key$method, names(method_labels))
			)
		) |>
		dplyr::ungroup()
	
	by_subtype <- long_table |>
		dplyr::group_by(.data$method, .data$comparison, .data$subtype_label) |>
		dplyr::group_modify(
			\(data, key) summarise_cophenetic_correlation(
				data,
				scope = key$subtype_label,
				settings = settings,
				seed_offset = 100L * match(key$subtype_label, c("H1N1", "H3N2")) +
					match(key$method, names(method_labels))
			)
		) |>
		dplyr::ungroup() |>
		dplyr::select(-dplyr::all_of("subtype_label"))
	
	dplyr::bind_rows(overall, by_subtype) |>
		dplyr::ungroup() |>
		dplyr::mutate(
			comparison = factor(.data$comparison, levels = unname(method_labels)),
			scope = factor(.data$scope, levels = c("Overall", "H1N1", "H3N2"))
		) |>
		dplyr::arrange(.data$comparison, .data$scope) |>
		dplyr::transmute(
			Comparison = as.character(.data$comparison),
			Scope = as.character(.data$scope),
			`Pairwise comparisons` = .data$n_pairs,
			`Pearson r` = .data$pearson_r,
			`CI lower` = .data$ci_lower,
			`CI upper` = .data$ci_upper,
			`CI method` = .data$ci_method
		)
}

make_cophenetic_correlation_table <- function(cophenetic_correlation_summary) {
	cophenetic_correlation_summary |>
		dplyr::mutate(
			`Pearson r` = sprintf("%.2f", .data$`Pearson r`),
			`95% CI` = paste0(
				"(",
				sprintf("%.2f", .data$`CI lower`),
				", ",
				sprintf("%.2f", .data$`CI upper`),
				")"
			)
		) |>
		dplyr::select(-dplyr::all_of(c("CI lower", "CI upper", "CI method"))) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Comparison") |>
		flextable::valign(j = "Comparison", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Descriptive Pearson correlations between each distance metric and ",
			"ML-tree cophenetic distance, using unique off-diagonal strain pairs. ",
			"Intervals are 95% Bayesian bootstrap intervals over sampled strains ",
			"when the bayesboot package is available."
		))
}

write_cophenetic_correlation_table <- function(cophenetic_correlation_summary, path) {
	table <- make_cophenetic_correlation_table(cophenetic_correlation_summary)
	write_rds_target(table, path)
}

make_ml_tree_plot <- function(tree_analysis, title) {
	phangorn::midpoint(tree_analysis$ml_tree$tree) |>
		ggtree::ggtree(ladderize = FALSE, layout = "circular") +
		ggtree::geom_tiplab() +
		ggplot2::labs(title = title) +
		ggplot2::theme(
			plot.margin = ggplot2::margin(-0.5, -0.5, -0.5, -0.5, "cm"),
			plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", vjust = -15, size = 20)
		)
}

write_ml_tree_plot <- function(tree_analysis, path, title) {
	ensure_dir(path)
	plot <- make_ml_tree_plot(tree_analysis, title)
	if (identical(tree_analysis$subtype, "h3")) {
		plot <- ggtree::rotate_tree(plot, 90)
	}
	ggplot2::ggsave(path, plot = plot, width = 6.5, height = 6.5, dpi = 300)
	path
}

write_combined_ml_tree_plot <- function(h1_tree_analysis, h3_tree_analysis, path) {
	ensure_dir(path)
	h1_plot <- make_ml_tree_plot(h1_tree_analysis, "H1N1")
	h3_plot <- make_ml_tree_plot(h3_tree_analysis, "H3N2") |>
		ggtree::rotate_tree(90)
	
	combined <- cowplot::plot_grid(h1_plot, h3_plot, nrow = 1) +
		ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"))
	
	ggplot2::ggsave(path, plot = combined, width = 13, height = 6.5, dpi = 300)
	
	if (requireNamespace("magick", quietly = TRUE)) {
		magick::image_read(path) |>
			magick::image_crop(geometry = "3861x1930", gravity = "center") |>
			magick::image_write(path)
	}
	
	path
}

tree_comparison_tree_labels <- function() {
	c("ml" = "Maximum likelihood baseline", distance_method_labels())
}

tree_analysis_method_order <- function() {
	names(tree_comparison_tree_labels())
}

extract_tree_metric_to_ml <- function(tree_analysis, metric) {
	metric_matrix <- tree_analysis$tree_distance_metrics[[metric]]
	if (is.null(metric_matrix)) {
		return(rep(NA_real_, length(tree_analysis_method_order())))
	}
	if (!"ml" %in% colnames(metric_matrix)) {
		stop(metric, " tree-distance matrix must include an 'ml' column.", call. = FALSE)
	}
	method_order <- tree_analysis_method_order()
	if (!all(method_order %in% rownames(metric_matrix))) {
		missing <- setdiff(method_order, rownames(metric_matrix))
		stop(
			metric,
			" tree-distance matrix is missing row(s): ",
			paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
	as.numeric(metric_matrix[method_order, "ml"])
}

make_tree_comparison_summary_one <- function(tree_analysis, subtype_label = NULL) {
	sh_table <- tibble::as_tibble(tree_analysis$sh_test)
	tree_labels <- tree_comparison_tree_labels()
	method_order <- tree_analysis_method_order()
	
	if (nrow(sh_table) != length(method_order)) {
		stop(
			"SH-test table must have ",
			length(method_order),
			" rows in ML/tree-method order.",
			call. = FALSE
		)
	}
	if (is.null(subtype_label)) {
		subtype_label <- subtype_display_label(tree_analysis$subtype)
	}
	
	sh_table |>
		dplyr::mutate(
			Subtype = subtype_label,
			tree_key = method_order,
			`Distance tree` = unname(tree_labels[.data$tree_key]),
			log_likelihood = .data[["ln L"]],
			delta_log_likelihood = .data[["Diff ln L"]],
			sh_p_value = .data[["p-value"]],
			rf_distance = extract_tree_metric_to_ml(tree_analysis, "RF"),
			ml_log_likelihood = .data$log_likelihood[.data$tree_key == "ml"],
			.before = 1
		) |>
		dplyr::filter(.data$tree_key != "ml") |>
		dplyr::transmute(
			.data$Subtype,
			.data$`Distance tree`,
			`Delta log likelihood` = .data$delta_log_likelihood,
			`Tree log likelihood` = .data$log_likelihood,
			`ML log likelihood` = .data$ml_log_likelihood,
			`SH p-value` = .data$sh_p_value,
			`RF distance` = .data$rf_distance
		)
}

make_tree_comparison_summary <- function(h1_tree_analysis, h3_tree_analysis) {
	dplyr::bind_rows(
		make_tree_comparison_summary_one(h1_tree_analysis, "H1N1"),
		make_tree_comparison_summary_one(h3_tree_analysis, "H3N2")
	)
}

format_tree_comparison_summary <- function(tree_comparison_summary) {
	tree_comparison_summary |>
		dplyr::mutate(
			`Delta log likelihood` = sprintf("%.1f", .data$`Delta log likelihood`),
			`Tree log likelihood` = sprintf("%.1f", .data$`Tree log likelihood`),
			`ML log likelihood` = sprintf("%.1f", .data$`ML log likelihood`),
			`SH p-value` = format_p_value(.data$`SH p-value`),
			`RF distance` = as.character(.data$`RF distance`)
		)
}

make_tree_comparison_table <- function(tree_comparison_summary) {
	table_data <- format_tree_comparison_summary(tree_comparison_summary)
	
	table_data |>
		flextable::flextable() |>
		flextable::merge_v(j = c("Subtype", "ML log likelihood")) |>
		flextable::valign(j = c("Subtype", "ML log likelihood"), valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Distance-tree comparison against the optimized maximum-likelihood ",
			"tree. Delta log likelihood is the primary comparison metric; ",
			"Shimodaira-Hasegawa p-values and Robinson-Foulds distances are ",
			"secondary checks."
		))
}

write_tree_comparison_table <- function(h1_tree_analysis, h3_tree_analysis, path) {
	table <- make_tree_comparison_summary(h1_tree_analysis, h3_tree_analysis) |>
		make_tree_comparison_table()
	write_rds_target(table, path)
}

topology_metric_labels <- function() {
	c(
		SPR = "SPR distance",
		RF = "RF distance",
		wRF = "Weighted RF distance",
		KF = "Branch-score distance",
		path = "Path distance"
	)
}

make_topology_distance_summary_one <- function(tree_analysis, subtype_label = NULL) {
	method_order <- tree_analysis_method_order()
	if (is.null(subtype_label)) {
		subtype_label <- subtype_display_label(tree_analysis$subtype)
	}

	base <- tibble::tibble(
		tree_key = method_order,
		Subtype = subtype_label,
		`Distance tree` = unname(tree_comparison_tree_labels()[method_order])
	) |>
		dplyr::filter(.data$tree_key != "ml")

	metric_values <- purrr::imap(
		topology_metric_labels(),
		\(metric_label, metric_name) {
			stats::setNames(
				list(extract_tree_metric_to_ml(tree_analysis, metric_name)[method_order != "ml"]),
				metric_label
			)
		}
	) |>
		purrr::list_c()

	dplyr::bind_cols(
		base |> dplyr::select(-dplyr::all_of("tree_key")),
		tibble::as_tibble(metric_values)
	)
}

make_topology_distance_summary <- function(h1_tree_analysis, h3_tree_analysis) {
	dplyr::bind_rows(
		make_topology_distance_summary_one(h1_tree_analysis, "H1N1"),
		make_topology_distance_summary_one(h3_tree_analysis, "H3N2")
	)
}

make_topology_distance_table <- function(topology_distance_summary) {
	topology_distance_summary |>
		dplyr::mutate(
			dplyr::across(
				dplyr::all_of(unname(topology_metric_labels())),
				\(x) sprintf("%.2f", x)
			)
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Supplemental topology-distance metrics between each distance tree ",
			"and the optimized maximum-likelihood tree."
		))
}

write_topology_distance_table <- function(h1_tree_analysis, h3_tree_analysis, path) {
	table <- make_topology_distance_summary(h1_tree_analysis, h3_tree_analysis) |>
		make_topology_distance_table()
	write_rds_target(table, path)
}

make_model_selection_summary <- function(model_tests_by_subtype, model_choice) {
	purrr::imap_dfr(
		model_tests_by_subtype,
		\(model_test, subtype) {
			best_model <- best_model_by_aicc(model_test)
			selected_model <- model_choice$selected_models[[subtype]]
			if (!"model_test_gamma_categories" %in% names(model_test)) {
				model_test$model_test_gamma_categories <- NA_integer_
			}
			model_test |>
				dplyr::mutate(
					Subtype = subtype_display_label(subtype),
					`Delta AICc` = .data$AICc - min(.data$AICc, na.rm = TRUE),
					`Best by AICc` = .data$Model == best_model,
					`Selected for ML tree` = .data$Model == selected_model,
					`Log-likelihood loss fraction` = purrr::map_dbl(
						.data$Model,
						\(model) model_log_likelihood_loss_fraction(model_test, model)
					),
					.before = 1
				)
		}
	) |>
		dplyr::select(
			"Subtype",
			"Model",
			"model_test_gamma_categories",
			"df",
			"logLik",
			"AICc",
			"Delta AICc",
			"AICcw",
			"BIC",
			"Log-likelihood loss fraction",
			"Best by AICc",
			"Selected for ML tree"
		)
}

make_model_selection_table <- function(model_selection_summary, model_choice) {
	strategy_note <- if (identical(model_choice$strategy, "common")) {
		paste0(
			"Common model selected for both subtypes: ",
			model_choice$selected_model,
			"."
		)
	} else {
		"Subtype-specific models selected because no common model met the likelihood-loss tolerance."
	}

	model_selection_summary |>
		dplyr::mutate(
			logLik = sprintf("%.1f", .data$logLik),
			AICc = sprintf("%.1f", .data$AICc),
			`Delta AICc` = sprintf("%.1f", .data$`Delta AICc`),
			AICcw = sprintf("%.3f", .data$AICcw),
			BIC = sprintf("%.1f", .data$BIC),
			`Log-likelihood loss fraction` = sprintf("%.4f", .data$`Log-likelihood loss fraction`),
			`Best by AICc` = dplyr::if_else(.data$`Best by AICc`, "yes", ""),
			`Selected for ML tree` = dplyr::if_else(.data$`Selected for ML tree`, "yes", "")
		) |>
		flextable::flextable() |>
		flextable::merge_v(j = "Subtype") |>
		flextable::valign(j = "Subtype", valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"FLU-family amino-acid substitution model comparison by subtype. ",
			"Models are selected by AICc within subtype, then a common model is ",
			"preferred when its log-likelihood loss is within ",
			sprintf("%.0f%%", 100 * model_choice$performance_tolerance),
			" for each subtype. ",
			strategy_note
		))
}

write_model_selection_table <- function(
		model_tests_by_subtype,
		model_choice,
		path
	) {
	table <- make_model_selection_summary(model_tests_by_subtype, model_choice) |>
		make_model_selection_table(model_choice)
	write_rds_target(table, path)
}

#### END OF FILE ####
