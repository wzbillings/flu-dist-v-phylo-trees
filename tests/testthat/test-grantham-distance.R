test_that("Grantham substitution matrix contains canonical values", {
	mat <- grantham_substitution_matrix()

	expect_equal(dim(mat), c(20L, 20L))
	expect_equal(mat["A", "A"], 0)
	expect_equal(mat["A", "V"], 64)
	expect_equal(mat["I", "L"], 5)
	expect_equal(mat["C", "W"], 215)
	expect_equal(mat, t(mat))
})

test_that("pairwise Grantham distance is averaged over comparable sites", {
	expect_equal(grantham_pair_distance("AAAA", "AAAA"), 0)
	expect_equal(grantham_pair_distance("A", "V"), 64)
	expect_equal(grantham_pair_distance("A-DX", "VWEA"), (64 + 45) / 2)
	expect_true(is.na(grantham_pair_distance("----", "XXXX")))
	expect_error(grantham_pair_distance(NA_character_, "AAAA"), "must not be missing")
	expect_error(grantham_pair_distance("AAA", "AAAA"), "same aligned length")
})

test_that("Grantham distance matrices are symmetric and preserve strain names", {
	seqs <- c(one = "AA", two = "AV", three = "VV")

	distances <- dist_grantham(seqs)

	expect_equal(dim(distances), c(3L, 3L))
	expect_equal(rownames(distances), names(seqs))
	expect_equal(colnames(distances), names(seqs))
	expect_equal(unname(diag(distances)), c(0, 0, 0))
	expect_equal(distances, t(distances))
	expect_equal(distances["one", "two"], 32)
	expect_equal(distances["one", "three"], 64)
})

test_that("complete-deletion Grantham distances remove non-comparable sites across all sequences", {
	seqs <- c(
		one = "AAAA",
		two = "AVXA",
		three = "AV-A"
	)

	distances <- dist_grantham(seqs, deletion = "complete")

	expect_equal(distances["one", "two"], 64 / 3)
	expect_equal(distances["one", "three"], 64 / 3)
	expect_equal(distances["two", "three"], 0)
	expect_equal(distances, t(distances))
})
