###
# Pairwise strain distance calculations
# Zane Billings
###

as_protein_phydat <- function(alignment_result) {
	phangorn::as.phyDat(alignment_result$protein_msa)
}

extract_aligned_proteins <- function(alignment_result) {
	seqs <- alignment_result$aligned_sequences$pro_aligned
	names(seqs) <- alignment_result$aligned_sequences$short_name
	seqs
}

calculate_hamming_distance <- function(protein_phydat) {
	protein_phydat |>
		phangorn::dist.hamming() |>
		as.matrix()
}

calculate_pepitope_distance <- function(aligned_proteins, subtype) {
	dist.pepi(aligned_proteins, subtype = subtype)
}

read_cartography_map <- function(path) {
	validate_file_exists(path, "cartography .ace file")
	Racmacs::read.acmap(path)
}

calculate_cartography_distance <- function(cartography_map, virus_metadata) {
	dist_cart <- racmacs_map_to_distances(cartography_map)
	colnames(dist_cart) <- replace_strain_names(
		colnames(dist_cart),
		virus_info = virus_metadata,
		from = "analysis",
		to = "short"
	)
	rownames(dist_cart) <- replace_strain_names(
		rownames(dist_cart),
		virus_info = virus_metadata,
		from = "analysis",
		to = "short"
	)
	dist_cart
}

calculate_subtype_distances <- function(
		alignment_result,
		cartography_map,
		virus_metadata
	) {
	subtype <- alignment_result$subtype
	protein_phydat <- as_protein_phydat(alignment_result)
	aligned_proteins <- extract_aligned_proteins(alignment_result)
	
	distances <- list(
		year = dist_year(names(aligned_proteins), virus_info = virus_metadata),
		hamming = calculate_hamming_distance(protein_phydat),
		pepi = calculate_pepitope_distance(aligned_proteins, subtype = subtype),
		cart = calculate_cartography_distance(cartography_map, virus_metadata)
	)
	
	distances <- standardize_distance_order(
		distances,
		reference_names = names(aligned_proteins),
		subtype = subtype
	)
	validate_distance_set(distances, subtype = subtype)
	distances
}

standardize_distance_order <- function(distances, reference_names, subtype = "subtype") {
	purrr::imap(
		distances,
		\(mat, method) {
			validate_distance_matrix(mat, paste(subtype, method))
			if (!setequal(rownames(mat), reference_names)) {
				missing <- setdiff(reference_names, rownames(mat))
				extra <- setdiff(rownames(mat), reference_names)
				stop(
					"Distance matrix ",
					method,
					" for ",
					subtype,
					" has incompatible strain names. Missing: ",
					paste(missing, collapse = ", "),
					"; extra: ",
					paste(extra, collapse = ", "),
					call. = FALSE
				)
			}
			mat[reference_names, reference_names, drop = FALSE]
		}
	)
}

validate_distance_set <- function(distances, subtype = "subtype") {
	if (!is.list(distances) || length(distances) == 0) {
		stop("Distance set for ", subtype, " must be a non-empty list.", call. = FALSE)
	}
	purrr::iwalk(distances, \(mat, name) validate_distance_matrix(mat, paste(subtype, name)))
	
	reference_names <- rownames(distances[[1]])
	name_mismatch <- purrr::keep(distances, \(mat) !identical(rownames(mat), reference_names))
	if (length(name_mismatch) > 0) {
		stop(
			"Distance matrices for ",
			subtype,
			" do not share the same strain names/order: ",
			paste(names(name_mismatch), collapse = ", "),
			call. = FALSE
		)
	}
	invisible(distances)
}

combine_distance_tables <- function(distances_by_subtype, unique_pairs = TRUE) {
	purrr::imap_dfr(
		distances_by_subtype,
		\(distances, subtype) purrr::imap_dfr(
			distances,
			\(mat, method) tidy_dist_mat(
				mat,
				unique_pairs = unique_pairs,
				include_diagonal = !unique_pairs
			) |>
				dplyr::mutate(method = method, .before = 1)
		) |>
			dplyr::mutate(subtype = subtype, .before = 1)
	)
}

normalize_distance_table <- function(distance_table) {
	distance_table |>
		dplyr::group_by(.data$subtype, .data$method) |>
		dplyr::mutate(d = normalize_vector(.data$d)) |>
		dplyr::ungroup()
}

normalize_vector <- function(x) {
	xmax <- max(x, na.rm = TRUE)
	xmin <- min(x, na.rm = TRUE)
	if (isTRUE(all.equal(xmax, xmin))) {
		return(x * 0)
	}
	(x - xmin) / (xmax - xmin)
}

lower_triangle_values <- function(mat) {
	validate_distance_matrix(mat)
	mat[lower.tri(mat)]
}

align_distance_matrices <- function(mat1, mat2) {
	validate_distance_matrix(mat1)
	validate_distance_matrix(mat2)
	shared <- intersect(rownames(mat1), rownames(mat2))
	if (length(shared) < 3) {
		stop("Need at least three shared strains for a matrix comparison.", call. = FALSE)
	}
	mat1 <- mat1[shared, shared, drop = FALSE]
	mat2 <- mat2[shared, shared, drop = FALSE]
	list(mat1 = mat1, mat2 = mat2)
}

mantel_permutation_test <- function(
		mat1,
		mat2,
		permutations,
		seed = NULL,
		method = "pearson"
	) {
	aligned <- align_distance_matrices(mat1, mat2)
	mat1 <- aligned$mat1
	mat2 <- aligned$mat2
	x <- lower_triangle_values(mat1)
	y <- lower_triangle_values(mat2)
	observed <- stats::cor(x, y, method = method, use = "complete.obs")
	
	null <- with_seed(seed, {
		replicate(permutations, {
			perm <- sample(seq_len(nrow(mat2)))
			permuted <- mat2[perm, perm, drop = FALSE]
			stats::cor(
				x,
				lower_triangle_values(permuted),
				method = method,
				use = "complete.obs"
			)
		})
	})
	
	tibble::tibble(
		estimate = observed,
		p_value = (sum(abs(null) >= abs(observed), na.rm = TRUE) + 1) / (permutations + 1),
		permutations = permutations,
		method = paste0("mantel_", method)
	)
}

compare_distance_matrices <- function(distances_by_subtype, settings = make_analysis_settings()) {
	purrr::imap_dfr(
		distances_by_subtype,
		\(distances, subtype) {
			method_pairs <- utils::combn(names(distances), 2, simplify = FALSE)
			results <- purrr::map_dfr(method_pairs, \(pair) {
				m1 <- pair[[1]]
				m2 <- pair[[2]]
				aligned <- align_distance_matrices(distances[[m1]], distances[[m2]])
				descriptive <- tibble::tibble(
					estimate = stats::cor(
						lower_triangle_values(aligned$mat1),
						lower_triangle_values(aligned$mat2),
						method = "pearson",
						use = "complete.obs"
					),
					p_value = NA_real_,
					permutations = NA_integer_,
					method = "descriptive_pearson"
				)
				mantel <- mantel_permutation_test(
					distances[[m1]],
					distances[[m2]],
					permutations = settings$mantel_permutations,
					seed = settings$seed,
					method = "pearson"
				)
				dplyr::bind_rows(descriptive, mantel) |>
					dplyr::mutate(m1 = m1, m2 = m2, .before = 1)
			})
			results |>
				dplyr::mutate(subtype = subtype, .before = 1)
		}
	)
}

#### END OF FILE ####
