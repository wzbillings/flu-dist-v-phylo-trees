###
# Tree construction and comparison functions
# Zane Billings
###

model_test_gamma_categories <- function(settings = make_analysis_settings()) {
	categories <- as.integer(settings$model_test_gamma_categories)
	categories <- unique(categories)
	if (length(categories) == 0 || any(is.na(categories)) || any(categories < 2L)) {
		stop("Model-test gamma categories must be integer values of at least 2.", call. = FALSE)
	}
	categories
}

run_model_test <- function(protein_phydat, settings = make_analysis_settings()) {
	model_tests <- purrr::map_dfr(
		model_test_gamma_categories(settings),
		\(gamma_categories) {
			phangorn::modelTest(
				protein_phydat,
				model = "FLU",
				G = TRUE,
				I = TRUE,
				k = gamma_categories,
				control = phangorn::pml.control(trace = 0)
			) |>
				as.data.frame() |>
				dplyr::mutate(model_test_gamma_categories = gamma_categories, .before = 1)
		}
	)

	model_tests |>
		dplyr::arrange(.data$AICc, dplyr::desc(.data$logLik)) |>
		dplyr::group_by(.data$Model) |>
		dplyr::slice(1) |>
		dplyr::ungroup() |>
		dplyr::arrange(.data$AICc)
}

calculate_model_test <- function(alignment_result, settings = make_analysis_settings()) {
	alignment_result |>
		as_protein_phydat() |>
		run_model_test(settings)
}

best_model_by_aicc <- function(model_test) {
	check_required_columns(model_test, c("Model", "AICc"), "model-test table")
	model_test$Model[[which.min(model_test$AICc)]]
}

model_log_likelihood_loss_fraction <- function(model_test, candidate_model) {
	check_required_columns(model_test, c("Model", "logLik", "AICc"), "model-test table")
	if (!candidate_model %in% model_test$Model) {
		return(NA_real_)
	}
	best_model <- best_model_by_aicc(model_test)
	best_log_likelihood <- model_test$logLik[match(best_model, model_test$Model)]
	candidate_log_likelihood <- model_test$logLik[match(candidate_model, model_test$Model)]
	loss <- max(0, best_log_likelihood - candidate_log_likelihood)
	loss / max(abs(best_log_likelihood), .Machine$double.eps)
}

summarise_common_model_candidates <- function(model_tests_by_subtype) {
	if (is.null(names(model_tests_by_subtype)) || any(names(model_tests_by_subtype) == "")) {
		stop("Model-test list must be named by subtype.", call. = FALSE)
	}
	common_models <- Reduce(intersect, purrr::map(model_tests_by_subtype, "Model"))
	if (length(common_models) == 0) {
		stop("No common candidate models are present in every subtype model-test table.", call. = FALSE)
	}

	purrr::map_dfr(
		common_models,
		\(candidate_model) {
			losses <- purrr::map_dbl(
				model_tests_by_subtype,
				model_log_likelihood_loss_fraction,
				candidate_model = candidate_model
			)
			combined_aicc <- purrr::map_dbl(
				model_tests_by_subtype,
				\(model_test) model_test$AICc[match(candidate_model, model_test$Model)]
			) |>
				sum()
			tibble::tibble(
				Model = candidate_model,
				combined_AICc = combined_aicc,
				max_log_likelihood_loss_fraction = max(losses, na.rm = TRUE)
			)
		}
	) |>
		dplyr::arrange(.data$combined_AICc)
}

choose_tree_model <- function(
		model_tests_by_subtype,
		performance_tolerance = make_analysis_settings()$model_performance_tolerance
	) {
	candidate_summary <- summarise_common_model_candidates(model_tests_by_subtype)
	acceptable <- candidate_summary |>
		dplyr::filter(.data$max_log_likelihood_loss_fraction <= performance_tolerance) |>
		dplyr::arrange(.data$combined_AICc)

	if (nrow(acceptable) > 0) {
		selected_model <- acceptable$Model[[1]]
		selected_models <- stats::setNames(
			rep(selected_model, length(model_tests_by_subtype)),
			names(model_tests_by_subtype)
		)
		strategy <- "common"
	} else {
		selected_model <- NA_character_
		selected_models <- purrr::map_chr(model_tests_by_subtype, best_model_by_aicc)
		strategy <- "subtype-specific"
	}

	list(
		strategy = strategy,
		selected_model = selected_model,
		selected_models = selected_models,
		performance_tolerance = performance_tolerance,
		candidate_summary = candidate_summary
	)
}

default_tree_model <- function(settings = make_analysis_settings()) {
	if (identical(settings$mode, "full")) {
		settings$tree_model_full
	} else {
		settings$tree_model_fast
	}
}

fit_ml_tree <- function(
		protein_phydat,
		settings = make_analysis_settings(),
		model = NULL
	) {
	if (is.null(model)) {
		model <- default_tree_model(settings)
	}
	with_seed(settings$seed, {
		if (identical(settings$tree_fit_strategy, "pml_bb")) {
			return(phangorn::pml_bb(
				protein_phydat,
				model = model,
				rearrangement = "stochastic",
				method = "unrooted",
				site.rate = "gamma_quadrature",
				control = phangorn::pml.control(trace = 0)
			))
		}
		
		phangorn::pml_bb(
			protein_phydat,
			model = model,
			rearrangement = "NNI",
			method = "unrooted",
			site.rate = "gamma_quadrature",
			control = phangorn::pml.control(trace = 0)
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
			phangorn::optim.pml(control = phangorn::pml.control(trace = 0))
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
		settings = make_analysis_settings(),
		model_test = NULL,
		model_choice = NULL,
		selected_model = NULL
	) {
	protein_phydat <- as_protein_phydat(alignment_result)
	if (is.null(model_test)) {
		model_test <- run_model_test(protein_phydat, settings)
	}
	if (is.null(selected_model)) {
		selected_model <- if (!is.null(model_choice)) {
			model_choice$selected_models[[alignment_result$subtype]]
		} else {
			best_model_by_aicc(model_test)
		}
	}
	ml_tree <- fit_ml_tree(protein_phydat, settings, model = selected_model)
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
		selected_model = selected_model,
		model_choice = model_choice,
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
