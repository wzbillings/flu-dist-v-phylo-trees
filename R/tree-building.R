###
# Tree construction and comparison functions
# Zane Billings
###

run_model_test <- function(protein_phydat, settings = make_analysis_settings()) {
	phangorn::modelTest(
		protein_phydat,
		model = "FLU",
		G = TRUE,
		I = TRUE,
		k = settings$model_test_gamma_categories
	) |>
		as.data.frame()
}

fit_ml_tree <- function(protein_phydat, settings = make_analysis_settings()) {
	with_seed(settings$seed, {
		if (identical(settings$tree_fit_strategy, "pml_bb")) {
			return(phangorn::pml_bb(
				protein_phydat,
				model = settings$tree_model_full,
				rearrangement = "stochastic",
				method = "unrooted",
				site.rate = "gamma_quadrature"
			))
		}
		
		start_tree <- protein_phydat |>
			phangorn::dist.ml(model = settings$tree_model_fast) |>
			phangorn::NJ()
		
		start_fit <- phangorn::pml(
			start_tree,
			protein_phydat,
			model = settings$tree_model_fast
		)
		
		phangorn::optim.pml(
			start_fit,
			model = settings$tree_model_fast,
			rearrangement = "NNI",
			optNni = TRUE,
			optBf = FALSE,
			optQ = FALSE,
			optGamma = FALSE,
			optInv = FALSE
		)
	})
}

prepare_distances_for_tree_building <- function(distances, settings = make_analysis_settings()) {
	out <- distances
	if ("pepi" %in% names(out)) {
		out$pepi <- perturb_matrix(
			out$pepi,
			mag = settings$pepi_perturbation_magnitude,
			seed = settings$pepi_perturbation_seed
		)
	}
	out
}

build_neighbor_joining_trees <- function(distances, settings = make_analysis_settings()) {
	distances <- prepare_distances_for_tree_building(distances, settings)
	purrr::map(distances, phangorn::NJ)
}

score_neighbor_joining_trees <- function(
		neighbor_joining_trees,
		ml_tree,
		protein_phydat
	) {
	purrr::imap(
		neighbor_joining_trees,
		\(tree, tree_name) phangorn::pml(
			tree,
			protein_phydat,
			site.rate = "gamma_quadrature",
			bf = ml_tree$bf,
			Q = ml_tree$Q,
			inv = ml_tree$inv,
			k = ml_tree$k,
			shape = ml_tree$shape,
			rate = ml_tree$rate,
			model = ml_tree$model,
			ASC = FALSE
		) |>
			phangorn::optim.pml()
	)
}

run_sh_test <- function(ml_tree, scored_neighbor_trees, settings = make_analysis_settings()) {
	tree_list <- c("ml" = list(ml_tree), scored_neighbor_trees)
	phangorn::SH.test(tree_list, B = settings$sh_bootstrap)
}

extract_ml_tree_distances <- function(ml_tree) {
	ml_tree$tree |> ape::cophenetic.phylo()
}

calculate_delta_scores <- function(protein_phydat) {
	phangorn::delta.score(protein_phydat, arg = "all")
}

as_multiphylo <- function(tree_fits) {
	trees <- purrr::map(tree_fits, "tree")
	class(trees) <- "multiPhylo"
	trees
}

calculate_tree_distance_metrics <- function(ml_tree, scored_neighbor_trees) {
	tree_fits <- c("ml" = list(ml_tree), scored_neighbor_trees)
	multi <- as_multiphylo(tree_fits)
	methods <- list(
		SPR = phangorn::SPR.dist,
		RF = phangorn::RF.dist,
		wRF = phangorn::wRF.dist,
		KF = phangorn::KF.dist,
		path = phangorn::path.dist
	)
	purrr::map(methods, \(f) f(multi) |> as.matrix())
}

calculate_tree_analysis <- function(
		alignment_result,
		distances,
		settings = make_analysis_settings()
	) {
	protein_phydat <- as_protein_phydat(alignment_result)
	model_test <- run_model_test(protein_phydat, settings)
	ml_tree <- fit_ml_tree(protein_phydat, settings)
	neighbor_joining_trees <- build_neighbor_joining_trees(distances, settings)
	scored_neighbor_trees <- score_neighbor_joining_trees(
		neighbor_joining_trees,
		ml_tree,
		protein_phydat
	)
	sh_test <- run_sh_test(ml_tree, scored_neighbor_trees, settings)
	
	list(
		subtype = alignment_result$subtype,
		model_test = model_test,
		ml_tree = ml_tree,
		neighbor_joining_trees = neighbor_joining_trees,
		scored_neighbor_trees = scored_neighbor_trees,
		sh_test = sh_test,
		ml_tree_distances = extract_ml_tree_distances(ml_tree),
		delta_scores = calculate_delta_scores(protein_phydat),
		tree_distance_metrics = calculate_tree_distance_metrics(ml_tree, scored_neighbor_trees)
	)
}

#### END OF FILE ####
