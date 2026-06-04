###
# P-epitope calculator
# Zane Billings
# 2024-05-04
# Calculate the p-epitope distance between two aligned influenza sequences
###

# First we need to create a reference data frame of p-epitope site locations
get_pepitope_sites <- function(
		subtype,
		sites = c('full', 'a', 'b', 'c', 'd', 'e', 'all_epitopes')
	) {
	
	# standardize subtype
	subtype <- stringr::str_to_lower(subtype)
	
	# The sites are different for H1 and H3. These are based on Amanda's
	# code where she referenced the actual papers.
	if (startsWith(subtype, 'h1')) {
		full <- 1:326
		
		site_a <- c(
			118, 120:122, 126:129, 132:135, 137, 139:143, 146, 147, 149, 165, 252,
			253
		)
		
		site_b <- c(
			124, 125, 152:157, 160, 162, 183:187, 189:191, 193:196
		)
		
		site_c <- c(
			34:38, 40, 41, 43:45, 269:274, 276:278, 283, 288, 292, 295, 297, 298,
			302, 303, 305:310
		)
		
		site_d <- c(
			89, 94:96, 113, 117, 163, 164, 166:174, 176, 179, 198, 200, 202, 204:216,
			222:227, 235, 237, 241, 243:245
		)
		
		site_e <- c(
			47, 48, 50, 51, 53, 54, 56:58, 66, 68:75, 78:80, 82:86, 102, 257:261,
			263, 267
		)
		
	} else if (startsWith(subtype, 'h3')) {
		full <- 1:328
		
		site_a <- c(
			122, 124, 126, 130:133, 135, 137:138, 140, 142:146, 150, 152,
			168
		)
		
		site_b <- c(
			128, 129, 155:160, 163, 165, 186:190, 192:194, 196:198
		)
		
		site_c <- c(
			44:48, 50, 51, 53, 54, 273, 275:276, 278:280, 294, 297, 299, 300,
			304:305, 307:312
		)
		
		site_d <- c(
			96, 102, 103, 117, 121, 167, 170:177, 179, 182, 201, 203, 207:209,
			212:219, 226:230, 238, 240, 242, 244, 246:248
		)
		
		site_e <- c(
			57, 59, 62, 63, 67, 75, 78, 80:83, 86:88, 91, 92, 94, 109, 260:262,
			265
		)
	} else {
		stop("'subtype' should be 'h1n1' or 'h3n2'.")
	}
	
	# Make a vector containing all epitope residues
	all <- c(site_a, site_b, site_c, site_d, site_e) |>
		unique() |>
		sort()
	
	# Make a list containing all the vectors
	res <- list(full, site_a, site_b, site_c, site_d, site_e, all)
	names(res) <- c('full', 'a', 'b', 'c', 'd', 'e', 'all_epitopes')
	
	# Now return the sites that were requested before
	out <- res[sites]
	return(out)
}

# Utility function for subsetting characters in a single string (in R language,
# getting a substring by position for a character vector of length 1)
extract_string_chars <- function(str, pos) {
	residues <- split_sequence_characters(str)
	if (any(pos > length(residues))) {
		stop("p-epitope site positions exceed the aligned sequence length.", call. = FALSE)
	}
	paste(residues[pos], collapse = "")
}

extract_site_residues <- function(seq, pos) {
	extract_string_chars(seq, pos) |>
		split_sequence_characters() |>
		toupper()
}

pepitope_site_mask <- function(seqs, sites) {
	residue_matrix <- do.call(rbind, purrr::map(seqs, extract_site_residues, pos = sites))
	apply(residue_matrix, 2, \(site) all(is_standard_amino_acid(site)))
}

pepitope_complete_site_masks <- function(seqs, subtype) {
	p_epi_sites <- get_pepitope_sites(subtype, sites = c('a', 'b', 'c', 'd', 'e'))
	purrr::map(p_epi_sites, \(sites) pepitope_site_mask(seqs, sites))
}

pepitope_epitope_distance <- function(seq_1, seq_2, sites, site_mask = NULL) {
	residues_1 <- extract_site_residues(seq_1, sites)
	residues_2 <- extract_site_residues(seq_2, sites)
	comparable <- is_standard_amino_acid(residues_1) & is_standard_amino_acid(residues_2)
	if (!is.null(site_mask)) {
		if (length(site_mask) != length(sites)) {
			stop("p-epitope site mask must match the number of epitope sites.", call. = FALSE)
		}
		comparable <- comparable & site_mask
	}
	if (!any(comparable)) {
		return(NA_real_)
	}
	mean(residues_1[comparable] != residues_2[comparable])
}

# Compute the p-epitope distance between two strings
pepitope <- function(seq_1, seq_2, subtype, site_masks = NULL) {
	# Get the numbers for the residues in each epitope
	p_epi_sites <- get_pepitope_sites(subtype, sites = c('a', 'b', 'c', 'd', 'e'))
	
	epi_dists <-
		purrr::imap_dbl(
			p_epi_sites,
			\(current_sites, site_name) pepitope_epitope_distance(
				seq_1,
				seq_2,
				current_sites,
				site_mask = site_masks[[site_name]]
			)
		)
	
	# Now we get the maximum of those and return it
	if (all(is.na(epi_dists))) {
		return(NA_real_)
	}
	p_epi <- max(epi_dists, na.rm = TRUE)
	return(p_epi)
}

# pepitope(prot_h1[[1]], prot_h1[[2]], 'h1')
# pepitope(prot_h3[[1]], prot_h3[[2]], 'h1')

# Compute a p-epitope distance matrix for all sequences in a character vector
dist.pepi <- function(seqs, subtype, deletion = c("pairwise", "complete")) {
	deletion <- match.arg(deletion)
	if (!is.character(seqs)) {
		stop("seqs must be a character vector of aligned amino-acid sequences.", call. = FALSE)
	}
	if (is.null(names(seqs)) || any(is.na(names(seqs))) || any(names(seqs) == "")) {
		stop("seqs must have non-missing strain names.", call. = FALSE)
	}
	validate_unique_values(names(seqs), "p-epitope sequence names")
	site_masks <- if (identical(deletion, "complete")) {
		pepitope_complete_site_masks(seqs, subtype)
	} else {
		NULL
	}
	# Set up an empty matrix to hold results
	res <- matrix(
		NA_real_,
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	
	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	if (length(seqs) > 1) {
		for (i in 2:nrow(res)) {
			for (j in 1:(i-1)) {
				distance <- pepitope(seqs[[i]], seqs[[j]], subtype, site_masks = site_masks)
				res[[i, j]] <- distance
				res[[j, i]] <- distance
			}
		}
	}
	
	# Set the diagonal to zero
	diag(res) <- 0
	return(res)
}

# test <- dist.pepi(prot_h1, 'h1')

##### END OF FILE ####
