library(targets)

tar_source("R")

tar_option_set(
	seed = 370,
	format = "rds"
)

list(
	tar_target(analysis_mode, Sys.getenv("FLU_TARGETS_MODE", "test")),
	tar_target(analysis_settings, make_analysis_settings(analysis_mode)),
	
	# Raw/source and first-pass cartography inputs.
	tar_target(sequence_source_file, "data/full-sequences.xlsx", format = "file"),
	tar_target(virus_name_file, "data/UGAFluVac-virus-names.csv", format = "file"),
	tar_target(h1_cartography_file, "data/h1_post_all_2d.ace", format = "file"),
	tar_target(h3_cartography_file, "data/h3_post_all_2d.ace", format = "file"),
	tar_target(
		source_input_audit,
		validate_source_inputs(
			sequence_source_file,
			virus_name_file,
			c(h1 = h1_cartography_file, h3 = h3_cartography_file)
		)
	),
	
	# Virus metadata and sequence cleaning.
	tar_target(virus_metadata, load_virus_metadata(virus_name_file)),
	tar_target(raw_sequences, load_raw_sequences(sequence_source_file)),
	tar_target(clean_sequences, clean_sequence_data(raw_sequences, virus_metadata)),
	tar_target(
		sequence_cleaning_summary,
		sequence_cleaning_audit(raw_sequences, virus_metadata, clean_sequences)
	),
	tar_target(
		clean_sequence_file,
		write_rds_target(clean_sequences, "results/derived/clean-data.rds"),
		format = "file"
	),
	tar_target(
		sequence_cleaning_summary_file,
		write_rds_target(sequence_cleaning_summary, "results/derived/sequence-cleaning-summary.rds"),
		format = "file"
	),
	
	# Multiple sequence alignments.
	tar_target(h1_alignment, align_h1_sequences(clean_sequences, analysis_settings)),
	tar_target(h3_alignment, align_h3_sequences(clean_sequences, analysis_settings)),
	tar_target(
		alignment_summary,
		dplyr::bind_rows(alignment_audit(h1_alignment), alignment_audit(h3_alignment))
	),
	tar_target(
		h1_rna_alignment_file,
		write_rds_target(h1_alignment$nucleotide_msa, "results/derived/h1-rna-alignment.rds"),
		format = "file"
	),
	tar_target(
		h1_protein_alignment_file,
		write_rds_target(h1_alignment$protein_msa, "results/derived/h1-pro-alignment.rds"),
		format = "file"
	),
	tar_target(
		h1_aligned_sequences_file,
		write_rds_target(h1_alignment$aligned_sequences, "results/derived/h1-seqs-aligned.rds"),
		format = "file"
	),
	tar_target(
		h3_rna_alignment_file,
		write_rds_target(h3_alignment$nucleotide_msa, "results/derived/h3-rna-alignment.rds"),
		format = "file"
	),
	tar_target(
		h3_protein_alignment_file,
		write_rds_target(h3_alignment$protein_msa, "results/derived/h3-pro-alignment.rds"),
		format = "file"
	),
	tar_target(
		h3_aligned_sequences_file,
		write_rds_target(h3_alignment$aligned_sequences, "results/derived/h3-seqs-aligned.rds"),
		format = "file"
	),
	tar_target(
		alignment_summary_file,
		write_rds_target(alignment_summary, "results/derived/alignment-summary.rds"),
		format = "file"
	),
	
	# Sequence and cartography distance calculations.
	tar_target(h1_cartography_map, read_cartography_map(h1_cartography_file)),
	tar_target(h3_cartography_map, read_cartography_map(h3_cartography_file)),
	tar_target(
		h1_distances,
		calculate_subtype_distances(h1_alignment, h1_cartography_map, virus_metadata)
	),
	tar_target(
		h3_distances,
		calculate_subtype_distances(h3_alignment, h3_cartography_map, virus_metadata)
	),
	tar_target(distances_by_subtype, list(h1 = h1_distances, h3 = h3_distances)),
	tar_target(distance_table, combine_distance_tables(distances_by_subtype, unique_pairs = TRUE)),
	tar_target(distance_table_normalized, normalize_distance_table(distance_table)),
	tar_target(
		distance_comparisons,
		compare_distance_matrices(distances_by_subtype, analysis_settings)
	),
	tar_target(
		h1_distance_file,
		write_rds_target(h1_distances, "results/derived/h1-dists.rds"),
		format = "file"
	),
	tar_target(
		h3_distance_file,
		write_rds_target(h3_distances, "results/derived/h3-dists.rds"),
		format = "file"
	),
	tar_target(
		distance_table_file,
		write_rds_target(distance_table, "results/derived/dist-df-unnormalized.rds"),
		format = "file"
	),
	tar_target(
		distance_table_normalized_file,
		write_rds_target(distance_table_normalized, "results/derived/dist-df-normalized.rds"),
		format = "file"
	),
	tar_target(
		distance_comparison_file,
		write_rds_target(distance_comparisons, "results/derived/distance-comparisons.rds"),
		format = "file"
	),
	
	# ML trees, distance trees, and tree comparisons.
	tar_target(h1_tree_analysis, calculate_tree_analysis(h1_alignment, h1_distances, analysis_settings)),
	tar_target(h3_tree_analysis, calculate_tree_analysis(h3_alignment, h3_distances, analysis_settings)),
	tar_target(tree_analyses, list(h1 = h1_tree_analysis, h3 = h3_tree_analysis)),
	tar_target(
		distances_with_tree_by_subtype,
		purrr::imap(
			distances_by_subtype,
			\(distances, subtype) c(
				list(cophenetic = tree_analyses[[subtype]]$ml_tree_distances),
				distances
			)
		)
	),
	tar_target(
		distance_comparisons_with_tree,
		compare_distance_matrices(distances_with_tree_by_subtype, analysis_settings)
	),
	tar_target(
		full_distance_table,
		make_full_distance_table(distances_by_subtype, tree_analyses)
	),
	tar_target(
		h1_model_test_file,
		write_rds_target(h1_tree_analysis$model_test, "results/derived/h1-model-test.rds"),
		format = "file"
	),
	tar_target(
		h3_model_test_file,
		write_rds_target(h3_tree_analysis$model_test, "results/derived/h3-model-test.rds"),
		format = "file"
	),
	tar_target(
		h1_ml_tree_file,
		write_rds_target(h1_tree_analysis$ml_tree, "results/derived/h1-ml-tree.rds"),
		format = "file"
	),
	tar_target(
		h3_ml_tree_file,
		write_rds_target(h3_tree_analysis$ml_tree, "results/derived/h3-ml-tree.rds"),
		format = "file"
	),
	tar_target(
		h1_nj_tree_file,
		write_rds_target(h1_tree_analysis$neighbor_joining_trees, "results/derived/nj-trees-h1.rds"),
		format = "file"
	),
	tar_target(
		h3_nj_tree_file,
		write_rds_target(h3_tree_analysis$neighbor_joining_trees, "results/derived/nj-trees-h3.rds"),
		format = "file"
	),
	tar_target(
		h1_sh_test_file,
		write_rds_target(h1_tree_analysis$sh_test, "results/derived/h1-sh-test.rds"),
		format = "file"
	),
	tar_target(
		h3_sh_test_file,
		write_rds_target(h3_tree_analysis$sh_test, "results/derived/h3-sh-test.rds"),
		format = "file"
	),
	tar_target(
		h1_tree_distance_file,
		write_rds_target(h1_tree_analysis$ml_tree_distances, "results/derived/h1-tree-dist.rds"),
		format = "file"
	),
	tar_target(
		h3_tree_distance_file,
		write_rds_target(h3_tree_analysis$ml_tree_distances, "results/derived/h3-tree-dist.rds"),
		format = "file"
	),
	tar_target(
		h1_delta_score_file,
		write_rds_target(h1_tree_analysis$delta_scores, "results/derived/h1-delta-scores.rds"),
		format = "file"
	),
	tar_target(
		h3_delta_score_file,
		write_rds_target(h3_tree_analysis$delta_scores, "results/derived/h3-delta-scores.rds"),
		format = "file"
	),
	tar_target(
		h1_tree_metric_file,
		write_rds_target(h1_tree_analysis$tree_distance_metrics, "results/derived/h1-tree-metrics.rds"),
		format = "file"
	),
	tar_target(
		h3_tree_metric_file,
		write_rds_target(h3_tree_analysis$tree_distance_metrics, "results/derived/h3-tree-metrics.rds"),
		format = "file"
	),
	tar_target(
		distance_comparison_with_tree_file,
		write_rds_target(distance_comparisons_with_tree, "results/derived/distance-comparisons-with-tree.rds"),
		format = "file"
	),
	
	# Manuscript-ready figures and tables.
	tar_target(
		correlation_plot_file,
		write_correlation_plot(
			full_distance_table,
			"results/Figures/corr-plot.png",
			analysis_settings
		),
		format = "file"
	),
	tar_target(
		cophenetic_mantel_summary,
		calculate_cophenetic_mantel_summary(distances_with_tree_by_subtype, analysis_settings)
	),
	tar_target(
		cophenetic_mantel_summary_file,
		write_rds_target(
			cophenetic_mantel_summary,
			"results/derived/cophenetic-mantel-summary.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_mantel_table_file,
		write_cophenetic_mantel_table(
			cophenetic_mantel_summary,
			"results/Tables/distance-mantel-table.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_subtype_contrast_summary,
		calculate_cophenetic_subtype_contrast_summary(
			distances_with_tree_by_subtype,
			analysis_settings
		)
	),
	tar_target(
		cophenetic_subtype_contrast_summary_file,
		write_rds_target(
			cophenetic_subtype_contrast_summary,
			"results/derived/cophenetic-subtype-contrast-summary.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_subtype_contrast_table_file,
		write_subtype_contrast_table(
			cophenetic_subtype_contrast_summary,
			"results/Tables/subtype-mantel-contrast-table.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_subtype_contrast_plot_file,
		write_subtype_contrast_plot(
			cophenetic_subtype_contrast_summary,
			"results/Figures/subtype-mantel-contrast-plot.png"
		),
		format = "file"
	),
	tar_target(
		cophenetic_subtype_contrast_sensitivity,
		calculate_cophenetic_subtype_contrast_sensitivity(
			distances_with_tree_by_subtype,
			analysis_settings
		)
	),
	tar_target(
		cophenetic_subtype_contrast_sensitivity_file,
		write_rds_target(
			cophenetic_subtype_contrast_sensitivity,
			"results/derived/cophenetic-subtype-contrast-sensitivity.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_subtype_contrast_sensitivity_table_file,
		write_subtype_contrast_sensitivity_table(
			cophenetic_subtype_contrast_sensitivity,
			"results/Tables/subtype-mantel-contrast-sensitivity-table.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_correlation_summary,
		calculate_cophenetic_correlation_summary(full_distance_table, analysis_settings)
	),
	tar_target(
		cophenetic_correlation_summary_file,
		write_rds_target(
			cophenetic_correlation_summary,
			"results/derived/descriptive-pearson-cophenetic-correlation-summary.rds"
		),
		format = "file"
	),
	tar_target(
		cophenetic_correlation_table_file,
		write_cophenetic_correlation_table(
			cophenetic_correlation_summary,
			"results/Tables/descriptive-pearson-distance-correlation-table.rds"
		),
		format = "file"
	),
	tar_target(
		h1_ml_tree_plot_file,
		write_ml_tree_plot(h1_tree_analysis, "results/Figures/h1-ml-tree.png", "H1N1"),
		format = "file"
	),
	tar_target(
		h3_ml_tree_plot_file,
		write_ml_tree_plot(h3_tree_analysis, "results/Figures/h3-ml-tree.png", "H3N2"),
		format = "file"
	),
	tar_target(
		combined_ml_tree_plot_file,
		write_combined_ml_tree_plot(
			h1_tree_analysis,
			h3_tree_analysis,
			"results/Figures/ml-trees.png"
		),
		format = "file"
	),
	tar_target(
		stat_table_file,
		write_stat_table(
			h1_tree_analysis,
			h3_tree_analysis,
			"results/Tables/stat-table.rds"
		),
		format = "file"
	),
	
	# Quarto manuscript render.
	tar_target(manuscript_qmd, "products/manuscript.qmd", format = "file"),
	tar_target(supplement_qmd, "products/supplement.qmd", format = "file"),
	tar_target(bibliography_file, "products/project-refs.bib", format = "file"),
	tar_target(csl_file, "products/plos-computational-biology.csl", format = "file"),
	tar_target(
		manuscript_docx,
		render_quarto_manuscript(
			manuscript_qmd,
			correlation_plot_file,
			combined_ml_tree_plot_file,
			cophenetic_mantel_table_file,
			cophenetic_subtype_contrast_table_file,
			cophenetic_subtype_contrast_plot_file,
			stat_table_file,
			bibliography_file,
			csl_file
		),
		format = "file"
	),
	tar_target(
		supplement_docx,
		render_quarto_supplement(
			supplement_qmd,
			cophenetic_mantel_summary_file,
			cophenetic_mantel_table_file,
			cophenetic_subtype_contrast_summary_file,
			cophenetic_subtype_contrast_table_file,
			cophenetic_subtype_contrast_sensitivity_file,
			cophenetic_subtype_contrast_sensitivity_table_file,
			cophenetic_correlation_table_file,
			bibliography_file,
			csl_file
		),
		format = "file"
	)
)
