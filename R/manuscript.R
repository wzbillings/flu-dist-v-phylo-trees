###
# Manuscript rendering helpers
# Zane Billings
###

render_quarto_manuscript <- function(
		qmd_path,
		correlation_plot_file,
		ml_tree_plot_file,
		stat_table_file,
		bibliography_file,
		csl_file,
		output_file = "products/manuscript.docx"
	) {
	validate_file_exists(qmd_path, "manuscript qmd")
	validate_file_exists(correlation_plot_file, "correlation figure")
	validate_file_exists(ml_tree_plot_file, "ML tree figure")
	validate_file_exists(stat_table_file, "statistical table")
	validate_file_exists(bibliography_file, "bibliography")
	validate_file_exists(csl_file, "CSL file")
	
	quarto <- Sys.which("quarto")
	if (identical(unname(quarto), "")) {
		stop("Quarto CLI is required to render the manuscript.", call. = FALSE)
	}
	
	ensure_dir(output_file)
	status <- system2(
		quarto,
		args = c("render", qmd_path, "--to", "docx"),
		stdout = TRUE,
		stderr = TRUE
	)
	if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
		stop(paste(status, collapse = "\n"), call. = FALSE)
	}
	validate_file_exists(output_file, "rendered manuscript")
	output_file
}

#### END OF FILE ####
