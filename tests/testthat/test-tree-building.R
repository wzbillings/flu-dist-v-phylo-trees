make_support_toy_trees <- function() {
	ml_tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
	alt_tree <- ape::read.tree(text = "((a:1,c:1):1,(b:1,(d:1,e:1):1):1);")
	bootstrap_trees <- list(ml_tree, ml_tree, alt_tree, ml_tree)
	class(bootstrap_trees) <- "multiPhylo"
	list(ml_tree = ml_tree, bootstrap_trees = bootstrap_trees)
}

test_that("branch-support details and summaries preserve bootstrap context", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	toy <- make_support_toy_trees()
	detail <- calculate_branch_support_detail(
		toy$ml_tree,
		toy$bootstrap_trees,
		subtype = "h1",
		bootstrap_replicates = length(toy$bootstrap_trees)
	)
	summary <- summarise_branch_support_detail(detail)

	expect_equal(nrow(detail), ape::Nnode(toy$ml_tree))
	expect_equal(
		names(detail),
		c("subtype", "node", "bootstrap_replicates", "support_percent", "support_category")
	)
	expect_true(all(detail$support_percent >= 0 & detail$support_percent <= 100))
	expect_equal(summary$subtype, "h1")
	expect_equal(summary$bootstrap_replicates, 4)
	expect_equal(summary$internal_branches, ape::Nnode(toy$ml_tree))
	expect_true(summary$mean_support_percent <= 100)
	expect_true(summary$branches_ge_70_percent <= summary$internal_branches)
})

test_that("branch-support details reject inconsistent bootstrap context", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")

	toy <- make_support_toy_trees()

	expect_error(
		calculate_branch_support_detail(
			toy$ml_tree,
			toy$bootstrap_trees,
			subtype = "h1",
			bootstrap_replicates = length(toy$bootstrap_trees) + 1L
		),
		"must match the number of bootstrap trees"
	)
})

test_that("bootstrap topology stability reports RF-style companion metrics", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	toy <- make_support_toy_trees()

	stability <- calculate_bootstrap_topology_stability(
		toy$ml_tree,
		toy$bootstrap_trees,
		subtype = "h1"
	)
	summary <- summarise_bootstrap_topology_stability(stability)

	expect_equal(nrow(stability), 4)
	expect_equal(
		names(stability),
		c(
			"subtype",
			"bootstrap_replicate",
			"rf_distance",
			"normalized_rf_distance",
			"weighted_rf_distance",
			"branch_score_distance",
			"path_distance",
			"distance_status"
		)
	)
	expect_equal(stability$rf_distance[[1]], 0)
	expect_equal(stability$normalized_rf_distance[[1]], 0)
	expect_equal(summary$subtype, "h1")
	expect_equal(summary$bootstrap_replicates, 4)
	expect_equal(summary$distance_failures, 0)
	expect_true(summary$identical_topology_fraction > 0)
	expect_true(summary$median_normalized_rf_distance >= 0)
})

test_that("bootstrap topology stability records unresolved-tree distance failures", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	toy <- make_support_toy_trees()
	star_tree <- ape::stree(length(toy$ml_tree$tip.label), type = "star")
	star_tree$tip.label <- toy$ml_tree$tip.label
	bootstrap_trees <- list(star_tree)
	class(bootstrap_trees) <- "multiPhylo"

	stability <- suppressMessages(
		calculate_bootstrap_topology_stability(
			toy$ml_tree,
			bootstrap_trees,
			subtype = "h1"
		)
	)
	summary <- summarise_bootstrap_topology_stability(stability)

	expect_true(all(is.na(dplyr::select(stability, -subtype, -bootstrap_replicate, -distance_status))))
	expect_equal(stability$distance_status[[1]], "distance_failed")
	expect_equal(summary$distance_failures, 1)
	expect_equal(summary$usable_topology_replicates, 0)
	expect_true(is.na(summary$identical_topology_fraction))
	expect_false(is.nan(summary$identical_topology_fraction))
})

test_that("bootstrap topology summary ignores failed rows with partial distances", {
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	topology_stability <- tibble::tibble(
		subtype = c("h1", "h1"),
		bootstrap_replicate = c(1L, 2L),
		rf_distance = c(0, 0),
		normalized_rf_distance = c(0, 1),
		weighted_rf_distance = c(0, 10),
		branch_score_distance = c(0, 20),
		path_distance = c(0, 30),
		distance_status = c("ok", "distance_failed")
	)

	summary <- summarise_bootstrap_topology_stability(topology_stability)

	expect_equal(summary$bootstrap_replicates, 2)
	expect_equal(summary$usable_topology_replicates, 1)
	expect_equal(summary$distance_failures, 1)
	expect_equal(summary$identical_topology_fraction, 1)
	expect_equal(summary$median_normalized_rf_distance, 0)
})

test_that("ML support wrapper accepts precomputed bootstrap trees for toy objects", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")
	skip_if_not_installed("tibble")
	skip_if_not_installed("dplyr")

	toy <- make_support_toy_trees()
	tree_analysis <- list(
		subtype = "h1",
		ml_tree = list(tree = toy$ml_tree)
	)

	support <- calculate_ml_tree_support(
		tree_analysis,
		settings = make_analysis_settings("test"),
		bootstrap_trees = toy$bootstrap_trees
	)

	expect_equal(support$subtype, "h1")
	expect_s3_class(support$bootstrap_trees, "multiPhylo")
	expect_s3_class(support$supported_tree, "phylo")
	expect_equal(support$branch_support_summary$bootstrap_replicates, 4)
	expect_equal(support$topology_stability_summary$bootstrap_replicates, 4)
})

test_that("ML support bootstrap runs with evaluated optNni settings", {
	skip_if_not_installed("ape")
	skip_if_not_installed("phangorn")

	base_sequence <- strsplit("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT", "")[[1]]
	mutate_site <- function(site, value) {
		out <- base_sequence
		out[[site]] <- value
		out
	}
	sequences <- list(
		a = base_sequence,
		b = mutate_site(8, "T"),
		c = mutate_site(18, "T"),
		d = mutate_site(28, "T"),
		e = mutate_site(38, "A"),
		f = mutate_site(1, "T")
	)
	dna_phydat <- phangorn::phyDat(sequences, type = "DNA")
	start_tree <- phangorn::NJ(phangorn::dist.ml(dna_phydat))
	fit <- phangorn::pml(start_tree, dna_phydat, model = "JC")
	settings <- make_analysis_settings("test")
	settings$ml_support_bootstrap <- 1L
	settings$ml_support_opt_nni <- TRUE

	bootstrap_trees <- run_ml_support_bootstrap(fit, settings, seed_offset = 1L)

	expect_s3_class(bootstrap_trees, "multiPhylo")
	expect_equal(length(bootstrap_trees), 1)
})
