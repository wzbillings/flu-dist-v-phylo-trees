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

plot_tree_distance_correlations <- function(
		full_distance_table,
		distance_comparisons,
		settings = make_analysis_settings()
	) {
	method_labels <- c(
		year = "Temporal distance",
		hamming = "Hamming distance",
		pepi = "p-Epitope distance",
		cart = "Cartographic distance"
	)
	
	plot_df <- full_distance_table |>
		dplyr::mutate(method = as.character(.data$method)) |>
		tidyr::pivot_wider(
			id_cols = c(.data$subtype, .data$Var1, .data$Var2),
			names_from = .data$method,
			values_from = .data$d
		) |>
		tidyr::pivot_longer(
			cols = dplyr::all_of(names(method_labels)),
			names_to = "method",
			values_to = "distance"
		) |>
		dplyr::mutate(
			method_label = factor(.data$method, levels = names(method_labels), labels = method_labels),
			subtype_label = dplyr::case_match(
				.data$subtype,
				"h1" ~ "H1N1",
				"h3" ~ "H3N2",
				.default = .data$subtype
			)
		)
	
	label_df <- distance_comparisons |>
		dplyr::filter(.data$method == "descriptive_pearson") |>
		dplyr::filter(.data$m1 == "cophenetic" | .data$m2 == "cophenetic") |>
		dplyr::mutate(
			method = dplyr::if_else(.data$m1 == "cophenetic", .data$m2, .data$m1),
			method_label = factor(.data$method, levels = names(method_labels), labels = method_labels),
			subtype_label = dplyr::case_match(
				.data$subtype,
				"h1" ~ "H1N1",
				"h3" ~ "H3N2",
				.default = .data$subtype
			),
			label = paste0("R = ", sprintf("%.2f", .data$estimate))
		)
	
	ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$cophenetic, y = .data$distance)) +
		ggplot2::geom_point(
			size = 2,
			alpha = 0.55,
			position = ggplot2::position_jitter(
				width = 0.01,
				height = 0.01,
				seed = settings$jitter_seed
			)
		) +
		ggplot2::geom_label(
			data = label_df,
			ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
			inherit.aes = FALSE,
			hjust = -0.05,
			vjust = 1.1,
			label.size = 0,
			fill = grDevices::adjustcolor("white", alpha.f = 0.75)
		) +
		ggplot2::facet_grid(.data$subtype_label ~ .data$method_label, scales = "free_y") +
		ggplot2::labs(x = "ML-tree cophenetic distance", y = NULL) +
		ggplot2::theme_bw(base_size = 11) +
		ggplot2::theme(
			strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
			panel.grid.minor = ggplot2::element_blank()
		)
}

write_correlation_plot <- function(
		full_distance_table,
		distance_comparisons,
		path,
		settings = make_analysis_settings()
	) {
	ensure_dir(path)
	plot <- plot_tree_distance_correlations(full_distance_table, distance_comparisons, settings)
	ggplot2::ggsave(path, plot = plot, width = 13, height = 7.5, dpi = 300)
	path
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

format_sh_table <- function(tree_analysis, subtype_label) {
	sh_table <- tibble::as_tibble(tree_analysis$sh_test)
	tree_labels <- c(
		"Maximum likelihood baseline",
		"Temporal distance",
		"Hamming distance",
		"p-Epitope distance",
		"Cartographic distance"
	)
	
	rf_distance <- rep(NA_real_, nrow(sh_table))
	if (!is.null(tree_analysis$tree_distance_metrics$RF)) {
		rf_matrix <- tree_analysis$tree_distance_metrics$RF
		if ("ml" %in% colnames(rf_matrix)) {
			rf_distance <- as.numeric(rf_matrix[, "ml"])
		}
	}
	
	sh_table |>
		dplyr::mutate(
			Subtype = subtype_label,
			Tree = tree_labels[seq_len(dplyr::n())],
			log_likelihood = .data[["ln L"]],
			delta_log_likelihood = .data[["Diff ln L"]],
			sh_p_value = .data[["p-value"]],
			rf_distance = rf_distance[seq_len(dplyr::n())],
			.before = 1
		) |>
		dplyr::transmute(
			.data$Subtype,
			.data$Tree,
			`log likelihood` = sprintf("%.1f", .data$log_likelihood),
			`Delta log likelihood` = sprintf("%.1f", .data$delta_log_likelihood),
			`SH p-value` = dplyr::if_else(
				.data$sh_p_value < 0.001,
				"< 0.001",
				sprintf("%.3f", .data$sh_p_value)
			),
			`RF distance` = dplyr::if_else(
				.data$Tree == "Maximum likelihood baseline",
				"NA",
				as.character(.data$rf_distance)
			)
		)
}

make_stat_table <- function(h1_tree_analysis, h3_tree_analysis) {
	table_data <- dplyr::bind_rows(
		format_sh_table(h1_tree_analysis, "H1N1"),
		format_sh_table(h3_tree_analysis, "H3N2")
	)
	
	table_data |>
		flextable::flextable() |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit() |>
		flextable::set_caption(paste0(
			"Log likelihood of constructed trees, change in log likelihood from ",
			"the ML tree, Shimodaira-Hasegawa test p-value, and ",
			"Robinson-Foulds distance from the ML tree."
		))
}

write_stat_table <- function(h1_tree_analysis, h3_tree_analysis, path) {
	table <- make_stat_table(h1_tree_analysis, h3_tree_analysis)
	write_rds_target(table, path)
}

#### END OF FILE ####
