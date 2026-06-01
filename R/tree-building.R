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

as_phylo_tree <- function(tree_or_fit) {
	if (inherits(tree_or_fit, "phylo")) {
		return(tree_or_fit)
	}
	if (is.list(tree_or_fit) && inherits(tree_or_fit$tree, "phylo")) {
		return(tree_or_fit$tree)
	}
	stop("Expected a phylo tree or fitted pml-style object with a phylo tree.", call. = FALSE)
}

validate_bootstrap_trees <- function(bootstrap_trees, reference_tree) {
	if (!inherits(bootstrap_trees, "multiPhylo") || length(bootstrap_trees) == 0) {
		stop("Bootstrap trees must be a non-empty multiPhylo object.", call. = FALSE)
	}
	reference_labels <- reference_tree$tip.label
	mismatched <- purrr::keep(
		seq_along(bootstrap_trees),
		\(i) !setequal(bootstrap_trees[[i]]$tip.label, reference_labels)
	)
	if (length(mismatched) > 0) {
		stop(
			"Bootstrap tree tip labels must match the reference ML tree. ",
			"Mismatched replicate(s): ",
			paste(mismatched, collapse = ", "),
			call. = FALSE
		)
	}
	invisible(bootstrap_trees)
}

support_category <- function(support_percent) {
	dplyr::case_when(
		is.na(support_percent) ~ "missing",
		support_percent >= 95 ~ ">=95%",
		support_percent >= 90 ~ "90-94%",
		support_percent >= 70 ~ "70-89%",
		support_percent > 0 ~ "<70%",
		TRUE ~ "0%"
	)
}

calculate_branch_support_detail <- function(
		reference_tree,
		bootstrap_trees,
		subtype,
		bootstrap_replicates = length(bootstrap_trees)
	) {
	reference_tree <- as_phylo_tree(reference_tree)
	validate_bootstrap_trees(bootstrap_trees, reference_tree)
	actual_bootstrap_replicates <- length(bootstrap_trees)
	bootstrap_replicates <- as.integer(bootstrap_replicates)
	if (length(bootstrap_replicates) != 1 || is.na(bootstrap_replicates) || bootstrap_replicates < 1) {
		stop("Bootstrap replicate count must be a positive integer.", call. = FALSE)
	}
	if (bootstrap_replicates != actual_bootstrap_replicates) {
		stop("Bootstrap replicate count must match the number of bootstrap trees.", call. = FALSE)
	}
	
	support_counts <- ape::prop.clades(reference_tree, bootstrap_trees, rooted = FALSE)
	if (length(support_counts) != ape::Nnode(reference_tree)) {
		stop("Branch-support counts did not match the number of internal ML-tree nodes.", call. = FALSE)
	}
	support_percent <- as.numeric(support_counts) / bootstrap_replicates * 100
	
	tibble::tibble(
		subtype = subtype,
		node = length(reference_tree$tip.label) + seq_along(support_percent),
		bootstrap_replicates = bootstrap_replicates,
		support_percent = support_percent,
		support_category = support_category(support_percent)
	)
}

summarise_branch_support_detail <- function(branch_support_detail) {
	check_required_columns(
		branch_support_detail,
		c("subtype", "bootstrap_replicates", "support_percent"),
		"branch-support detail"
	)
	branch_support_detail |>
		dplyr::group_by(.data$subtype, .data$bootstrap_replicates) |>
		dplyr::summarise(
			internal_branches = dplyr::n(),
			mean_support_percent = mean(.data$support_percent, na.rm = TRUE),
			median_support_percent = stats::median(.data$support_percent, na.rm = TRUE),
			min_support_percent = min(.data$support_percent, na.rm = TRUE),
			branches_ge_70_percent = sum(.data$support_percent >= 70, na.rm = TRUE),
			branches_ge_90_percent = sum(.data$support_percent >= 90, na.rm = TRUE),
			branches_ge_95_percent = sum(.data$support_percent >= 95, na.rm = TRUE),
			.groups = "drop"
		)
}

add_branch_support_labels <- function(reference_tree, branch_support_detail) {
	reference_tree <- as_phylo_tree(reference_tree)
	ordered_support <- branch_support_detail |>
		dplyr::arrange(.data$node) |>
		dplyr::pull("support_percent")
	if (length(ordered_support) != ape::Nnode(reference_tree)) {
		stop("Branch-support detail must contain one row per internal ML-tree node.", call. = FALSE)
	}
	reference_tree$node.label <- sprintf("%.0f", ordered_support)
	reference_tree
}

tree_distance_to_reference <- function(reference_tree, tree, distance_function, ...) {
	warned <- FALSE
	value <- tryCatch(
		withCallingHandlers(
			as.numeric(distance_function(reference_tree, tree, ...)),
			warning = function(w) {
				warned <<- TRUE
				invokeRestart("muffleWarning")
			}
		),
		error = function(e) NA_real_
	)
	if (warned) {
		return(NA_real_)
	}
	value
}

max_or_na <- function(x) {
	if (all(is.na(x))) {
		return(NA_real_)
	}
	max(x, na.rm = TRUE)
}

mean_or_na <- function(x) {
	if (length(x) == 0 || all(is.na(x))) {
		return(NA_real_)
	}
	mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
	if (length(x) == 0 || all(is.na(x))) {
		return(NA_real_)
	}
	stats::median(x, na.rm = TRUE)
}

usable_topology_metric <- function(metric, distance_status) {
	dplyr::if_else(distance_status == "ok", metric, NA_real_)
}

calculate_bootstrap_topology_stability <- function(reference_tree, bootstrap_trees, subtype) {
	reference_tree <- as_phylo_tree(reference_tree)
	validate_bootstrap_trees(bootstrap_trees, reference_tree)
	
	purrr::map_dfr(
		seq_along(bootstrap_trees),
		\(replicate_id) {
			tree <- bootstrap_trees[[replicate_id]]
			distances <- c(
				rf_distance = tree_distance_to_reference(
					reference_tree,
					tree,
					phangorn::RF.dist,
					normalize = FALSE,
					rooted = FALSE
				),
				normalized_rf_distance = tree_distance_to_reference(
					reference_tree,
					tree,
					phangorn::RF.dist,
					normalize = TRUE,
					rooted = FALSE
				),
				weighted_rf_distance = tree_distance_to_reference(
					reference_tree,
					tree,
					phangorn::wRF.dist,
					normalize = FALSE,
					rooted = FALSE
				),
				branch_score_distance = tree_distance_to_reference(
					reference_tree,
					tree,
					phangorn::KF.dist,
					rooted = FALSE
				),
				path_distance = tree_distance_to_reference(
					reference_tree,
					tree,
					phangorn::path.dist,
					use.weight = FALSE
				)
			)
			distance_status <- if (any(is.na(distances))) {
				distances[] <- NA_real_
				"distance_failed"
			} else {
				"ok"
			}
			
			tibble::tibble(
				subtype = subtype,
				bootstrap_replicate = as.integer(replicate_id),
				rf_distance = unname(distances[["rf_distance"]]),
				normalized_rf_distance = unname(distances[["normalized_rf_distance"]]),
				weighted_rf_distance = unname(distances[["weighted_rf_distance"]]),
				branch_score_distance = unname(distances[["branch_score_distance"]]),
				path_distance = unname(distances[["path_distance"]]),
				distance_status = distance_status
			)
		}
	)
}

summarise_bootstrap_topology_stability <- function(topology_stability) {
	check_required_columns(
		topology_stability,
		c(
			"subtype",
			"bootstrap_replicate",
			"rf_distance",
			"normalized_rf_distance",
			"weighted_rf_distance",
			"branch_score_distance",
			"path_distance",
			"distance_status"
		),
		"bootstrap topology-stability data"
	)
	topology_stability |>
		dplyr::group_by(.data$subtype) |>
		dplyr::summarise(
			bootstrap_replicates = dplyr::n(),
			usable_topology_replicates = sum(.data$distance_status == "ok", na.rm = TRUE),
			distance_failures = sum(.data$distance_status != "ok", na.rm = TRUE),
			identical_topology_fraction = mean_or_na(dplyr::if_else(
				.data$distance_status == "ok",
				.data$rf_distance == 0,
				NA
			)),
			median_normalized_rf_distance = median_or_na(usable_topology_metric(
				.data$normalized_rf_distance,
				.data$distance_status
			)),
			mean_normalized_rf_distance = mean_or_na(usable_topology_metric(
				.data$normalized_rf_distance,
				.data$distance_status
			)),
			max_normalized_rf_distance = max_or_na(usable_topology_metric(
				.data$normalized_rf_distance,
				.data$distance_status
			)),
			median_weighted_rf_distance = median_or_na(usable_topology_metric(
				.data$weighted_rf_distance,
				.data$distance_status
			)),
			median_branch_score_distance = median_or_na(usable_topology_metric(
				.data$branch_score_distance,
				.data$distance_status
			)),
			median_path_distance = median_or_na(usable_topology_metric(
				.data$path_distance,
				.data$distance_status
			)),
			.groups = "drop"
		)
}

ml_support_replicates <- function(settings = make_analysis_settings()) {
	reps <- as.integer(settings$ml_support_bootstrap)
	if (length(reps) != 1 || is.na(reps) || reps < 1L) {
		stop("ML support bootstrap replicates must be a positive integer.", call. = FALSE)
	}
	reps
}

run_ml_support_bootstrap <- function(ml_tree, settings = make_analysis_settings(), seed_offset = 0L) {
	if (!inherits(ml_tree, "pml")) {
		stop("ML support bootstrap requires a fitted pml object.", call. = FALSE)
	}
	with_seed(settings$seed + 30000L + seed_offset, {
		if (isTRUE(settings$ml_support_opt_nni)) {
			return(phangorn::bootstrap.pml(
				ml_tree,
				bs = ml_support_replicates(settings),
				trees = TRUE,
				optNni = TRUE,
				control = phangorn::pml.control(trace = 0)
			))
		}
		phangorn::bootstrap.pml(
			ml_tree,
			bs = ml_support_replicates(settings),
			trees = TRUE,
			optNni = FALSE,
			control = phangorn::pml.control(trace = 0)
		)
	})
}

calculate_ml_tree_support <- function(
		tree_analysis,
		settings = make_analysis_settings(),
		seed_offset = 0L,
		bootstrap_trees = NULL
	) {
	subtype <- tree_analysis$subtype
	reference_tree <- as_phylo_tree(tree_analysis$ml_tree)
	if (is.null(bootstrap_trees)) {
		bootstrap_trees <- run_ml_support_bootstrap(tree_analysis$ml_tree, settings, seed_offset)
	}
	validate_bootstrap_trees(bootstrap_trees, reference_tree)
	bootstrap_replicates <- length(bootstrap_trees)
	branch_support_detail <- calculate_branch_support_detail(
		reference_tree,
		bootstrap_trees,
		subtype = subtype,
		bootstrap_replicates = bootstrap_replicates
	)
	topology_stability <- calculate_bootstrap_topology_stability(
		reference_tree,
		bootstrap_trees,
		subtype = subtype
	)
	
	list(
		subtype = subtype,
		bootstrap_replicates = bootstrap_replicates,
		bootstrap_trees = bootstrap_trees,
		supported_tree = add_branch_support_labels(reference_tree, branch_support_detail),
		branch_support_detail = branch_support_detail,
		branch_support_summary = summarise_branch_support_detail(branch_support_detail),
		topology_stability = topology_stability,
		topology_stability_summary = summarise_bootstrap_topology_stability(topology_stability)
	)
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
