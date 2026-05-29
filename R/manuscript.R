###
# Manuscript rendering helpers
# Zane Billings
###

render_quarto_file <- function(qmd_path, output_file) {
	quarto <- Sys.which("quarto")
	if (identical(unname(quarto), "")) {
		stop("Quarto CLI is required to render the manuscript.", call. = FALSE)
	}

	ensure_dir(output_file)
	render_env <- list(
		R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep),
		RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE"
	)
	old_env <- Sys.getenv(names(render_env), unset = NA_character_)
	on.exit({
		unset_names <- names(old_env)[is.na(old_env)]
		if (length(unset_names) > 0) {
			Sys.unsetenv(unset_names)
		}
		restore_env <- as.list(old_env[!is.na(old_env)])
		if (length(restore_env) > 0) {
			do.call(Sys.setenv, restore_env)
		}
	}, add = TRUE)
	do.call(Sys.setenv, render_env)

	status <- system2(
		quarto,
		args = c("render", qmd_path, "--to", "docx"),
		stdout = TRUE,
		stderr = TRUE
	)
	if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
		stop(paste(status, collapse = "\n"), call. = FALSE)
	}
	validate_file_exists(output_file, "rendered Quarto document")
	output_file
}

render_quarto_manuscript <- function(
		qmd_path,
		correlation_plot_file,
		ml_tree_plot_file,
		cophenetic_correlation_table_file,
		stat_table_file,
		bibliography_file,
		csl_file,
		output_file = "products/manuscript.docx"
	) {
	validate_file_exists(qmd_path, "manuscript qmd")
	validate_file_exists(correlation_plot_file, "correlation figure")
	validate_file_exists(ml_tree_plot_file, "ML tree figure")
	validate_file_exists(cophenetic_correlation_table_file, "distance correlation table")
	validate_file_exists(stat_table_file, "statistical table")
	validate_file_exists(bibliography_file, "bibliography")
	validate_file_exists(csl_file, "CSL file")

	render_quarto_file(qmd_path, output_file)
}

render_quarto_supplement <- function(
		qmd_path,
		mantel_summary_file,
		mantel_correlation_table_file,
		pearson_correlation_table_file,
		bibliography_file,
		csl_file,
		output_file = "products/supplement.docx"
	) {
	validate_file_exists(qmd_path, "supplement qmd")
	validate_file_exists(mantel_summary_file, "Mantel distance correlation summary")
	validate_file_exists(mantel_correlation_table_file, "Mantel distance correlation table")
	validate_file_exists(pearson_correlation_table_file, "descriptive Pearson distance correlation table")
	validate_file_exists(bibliography_file, "bibliography")
	validate_file_exists(csl_file, "CSL file")

	render_quarto_file(qmd_path, output_file)
}

#### END OF FILE ####
