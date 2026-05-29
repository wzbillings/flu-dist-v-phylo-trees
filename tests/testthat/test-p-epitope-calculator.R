test_that("p-epitope site lookup returns expected H1 and H3 structures", {
	skip_if_not_installed("stringr")

	h1_sites <- get_pepitope_sites("H1N1", sites = c("full", "a", "all_epitopes"))
	h3_sites <- get_pepitope_sites("h3n2", sites = c("full", "a", "all_epitopes"))

	expect_equal(names(h1_sites), c("full", "a", "all_epitopes"))
	expect_equal(length(h1_sites$full), 326)
	expect_equal(length(h3_sites$full), 328)
	expect_true(all(h1_sites$a %in% h1_sites$all_epitopes))
	expect_true(all(h3_sites$a %in% h3_sites$all_epitopes))
	expect_error(get_pepitope_sites("b"), "'subtype' should be 'h1n1' or 'h3n2'.", fixed = TRUE)
})

test_that("p-epitope distance is zero for identical aligned sequences", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("purrr")
	skip_if_not_installed("stringr")

	seqs <- c(one = paste(rep("A", 326), collapse = ""), two = paste(rep("A", 326), collapse = ""))

	expect_equal(pepitope(seqs[[1]], seqs[[2]], "h1n1"), 0)

	distances <- dist.pepi(seqs, "h1n1")
	expect_equal(dim(distances), c(2L, 2L))
	expect_equal(rownames(distances), names(seqs))
	expect_equal(colnames(distances), names(seqs))
	expect_equal(unname(diag(distances)), c(0, 0))
	expect_equal(distances, t(distances))
})

test_that("p-epitope distance matrix detects a toy epitope change", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("purrr")
	skip_if_not_installed("stringr")

	seq_one <- paste(rep("A", 326), collapse = "")
	seq_two_chars <- rep("A", 326)
	seq_two_chars[118] <- "B"
	seq_two <- paste(seq_two_chars, collapse = "")
	seqs <- c(one = seq_one, two = seq_two)

	distances <- dist.pepi(seqs, "h1n1")

	expect_gt(distances["one", "two"], 0)
	expect_equal(distances["one", "two"], distances["two", "one"])
})
