###
# Reusable pipeline utilities
# Zane Billings
###

ensure_dir <- function(path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	path
}

validate_file_exists <- function(path, label = path) {
	if (!file.exists(path)) {
		stop(label, " does not exist at: ", path, call. = FALSE)
	}
	normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_rds_target <- function(object, path) {
	ensure_dir(path)
	readr::write_rds(object, path)
	path
}

check_required_columns <- function(data, required, data_name = "data") {
	missing <- setdiff(required, names(data))
	if (length(missing) > 0) {
		stop(
			data_name,
			" is missing required column(s): ",
			paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
	invisible(data)
}

normalize_join_key <- function(x) {
	x |>
		stringr::str_squish() |>
		stringr::str_to_lower()
}

validate_unique_values <- function(x, label) {
	duplicated_values <- unique(x[duplicated(x)])
	if (length(duplicated_values) > 0) {
		stop(
			label,
			" contains duplicate value(s): ",
			paste(duplicated_values, collapse = ", "),
			call. = FALSE
		)
	}
	invisible(x)
}

# Convert an MSA object to a named character vector.
alignment_to_character <- function(msa_alignment) {
	as.character(msa_alignment@unmasked)
}

normalize_matrix <- function(mat) {
	mmax <- max(mat, na.rm = TRUE)
	mmin <- min(mat, na.rm = TRUE)
	if (isTRUE(all.equal(mmax, mmin))) {
		return(mat * 0)
	}
	(mat - mmin) / (mmax - mmin)
}

validate_distance_matrix <- function(mat, label = "distance matrix") {
	if (!is.matrix(mat)) {
		stop(label, " must be a matrix.", call. = FALSE)
	}
	if (nrow(mat) != ncol(mat)) {
		stop(label, " must be square.", call. = FALSE)
	}
	if (is.null(rownames(mat)) || is.null(colnames(mat))) {
		stop(label, " must have row and column names.", call. = FALSE)
	}
	if (!identical(rownames(mat), colnames(mat))) {
		stop(label, " row and column names must be identical and ordered.", call. = FALSE)
	}
	invisible(mat)
}

tidy_dist_mat <- function(d, unique_pairs = FALSE, include_diagonal = TRUE) {
	validate_distance_matrix(d)
	row_names <- rownames(d)
	col_names <- colnames(d)
	out <- d |>
		tibble::as_tibble(rownames = "Var1") |>
		tidyr::pivot_longer(
			cols = -Var1,
			names_to = "Var2",
			values_to = "d"
		) |>
		dplyr::mutate(
			Var1 = as.character(.data$Var1),
			Var2 = as.character(.data$Var2),
			.row_index = match(.data$Var1, row_names),
			.col_index = match(.data$Var2, col_names)
		)
	
	if (isTRUE(unique_pairs)) {
		out <- dplyr::filter(out, .data$.row_index > .data$.col_index)
	} else if (!isTRUE(include_diagonal)) {
		out <- dplyr::filter(out, .data$.row_index != .data$.col_index)
	}
	
	out |>
		dplyr::select(-dplyr::all_of(c(".row_index", ".col_index"))) |>
		dplyr::mutate(
			Var1 = forcats::fct_inorder(.data$Var1),
			Var2 = forcats::fct_inorder(.data$Var2) |> forcats::fct_rev()
		)
}

tidy_unique_dist_mat <- function(d) {
	tidy_dist_mat(d, unique_pairs = TRUE, include_diagonal = FALSE)
}

racmacs_map_to_distances <- function(map) {
	coords <- Racmacs::agCoords(map)
	strains <- rownames(coords)
	x <- coords[, 1]
	y <- coords[, 2]
	
	res <- matrix(
		nrow = length(strains),
		ncol = length(strains),
		dimnames = list(strains, strains)
	)
	
	for (i in 2:nrow(res)) {
		for (j in 1:(i - 1)) {
			x_part <- (x[[j]] - x[[i]])^2
			y_part <- (y[[j]] - y[[i]])^2
			res[[i, j]] <- sqrt(x_part + y_part)
		}
	}
	
	diag(res) <- 0
	res |> as.dist() |> as.matrix()
}

replace_strain_names <- function(
		x,
		virus_info,
		from = "analysis",
		to = "short",
		drop = TRUE
	) {
	from_col <- switch(
		from,
		analysis = "analysis_name",
		full = "genbank_strain_name",
		short = "short_name",
		stop("'from' should be 'analysis', 'full', or 'short'.", call. = FALSE)
	)
	
	to_col <- switch(
		to,
		analysis = "analysis_name",
		full = "genbank_strain_name",
		short = "short_name",
		subtype = "subtype",
		stop("'to' should be 'analysis', 'full', 'short', or 'subtype'.", call. = FALSE)
	)
	
	from_vec <- virus_info[[from_col]]
	if (!all(x %in% from_vec)) {
		missing <- setdiff(x, from_vec)
		stop(
			"'x' includes ",
			from,
			" name(s) absent from virus metadata: ",
			paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
	
	vals <- virus_info[[to_col]][match(x, from_vec)]
	if (isTRUE(drop) && is.factor(vals)) {
		vals <- forcats::fct_drop(vals)
	}
	vals
}

extract_year <- function(analysis_names) {
	stringr::str_extract(analysis_names, "[0-9]{4}$") |>
		as.integer()
}

dist_year <- function(names, virus_info, format = "short") {
	years <- names |>
		replace_strain_names(virus_info = virus_info, from = format, to = "analysis") |>
		extract_year()
	res <- stats::dist(years, method = "manhattan") |> as.matrix()
	short_names <- replace_strain_names(names, virus_info = virus_info, from = format, to = "short")
	colnames(res) <- short_names
	rownames(res) <- short_names
	res
}

dist.year <- function(names, format = "short", virus_info = NULL) {
	if (is.null(virus_info)) {
		stop("`virus_info` is required; global metadata reads are not allowed.", call. = FALSE)
	}
	dist_year(names = names, virus_info = virus_info, format = format)
}

with_seed <- function(seed, code) {
	if (is.null(seed)) {
		return(force(code))
	}
	had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
	on.exit({
		if (had_seed) {
			assign(".Random.seed", old_seed, envir = .GlobalEnv)
		} else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
			rm(".Random.seed", envir = .GlobalEnv)
		}
	}, add = TRUE)
	set.seed(seed)
	force(code)
}

perturb_matrix <- function(mat, mag = 0.01, seed = NULL) {
	validate_distance_matrix(mat)
	with_seed(seed, {
		noise <- matrix(
			stats::runif(nrow(mat) * ncol(mat), -mag, mag),
			ncol = ncol(mat),
			nrow = nrow(mat),
			dimnames = dimnames(mat)
		)
		noise <- (noise + t(noise)) / 2
		diag(noise) <- 0
		out <- mat + noise
		diag(out) <- 0
		out
	})
}

make_analysis_settings <- function(mode = Sys.getenv("FLU_TARGETS_MODE", "test")) {
	mode <- tolower(mode)
	if (!mode %in% c("test", "full")) {
		stop("FLU_TARGETS_MODE must be either 'test' or 'full'.", call. = FALSE)
	}
	
	is_test <- identical(mode, "test")
	list(
		mode = mode,
		seed = 370L,
		alignment_method = "Muscle",
		model_test_gamma_categories = if (is_test) 4L else 20L,
		tree_fit_strategy = if (is_test) "fast_nni" else "pml_bb",
		tree_model_full = "FLU+G(20)+I",
		tree_model_fast = "FLU",
		sh_bootstrap = if (is_test) 100L else 1000000L,
		mantel_permutations = if (is_test) 99L else 9999L,
		mantel_bootstrap_reps = if (is_test) 199L else 4000L,
		subtype_contrast_permutations = if (is_test) 99L else 9999L,
		subtype_contrast_bootstrap_reps = if (is_test) 199L else 4000L,
		correlation_bootstrap_reps = if (is_test) 1000L else 4000L,
		correlation_ci_level = 0.95,
		pepi_perturbation_magnitude = 0.01,
		pepi_perturbation_seed = 370L,
		jitter_seed = 132413L
	)
}

#### END OF FILE ####
