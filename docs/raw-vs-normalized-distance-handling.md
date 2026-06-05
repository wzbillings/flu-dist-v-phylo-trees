# Raw Versus Normalized Distance Handling Audit

Date: 2026-06-05

Purpose: document which manuscript metrics use raw distance values and which
use min-max normalized values. This is an audit record for the current
targets-based analysis. The table is written so it can be adapted into the
supplement if a reviewer asks for explicit raw-versus-normalized reporting.

## Policy

Use raw distances for analyses where one metric is used at a time, including
neighbor-joining tree construction. For direct cross-metric comparisons, use raw
distances when the statistic is scale-invariant, such as Mantel r, Pearson r, or
rank-based Spearman summaries.
Normalize distances only for explicitly labeled visual or audit outputs where
different units would otherwise make the display or interpretation depend on
arbitrary scale differences.

When normalization is used, min-max normalization is applied by distance metric
across all displayed or audited unique off-diagonal strain pairs. This preserves
H1N1/H3N2 subtype differences within a metric while putting each metric on a
0-1 display scale.

## Audit Table

| Output | Reported statistic | Distance scale | Normalization scope | Reason | Supplement placement |
|---|---|---|---|---|---|
| Source distance matrices | Temporal, Grantham, p-Epitope, cartographic, and ML-tree cophenetic distances | Raw distances | None | These are the metric-specific analysis inputs and preserve each method's original unit or branch-length scale. | Supplement-ready audit |
| Neighbor-joining distance trees | Distance-derived tree topologies and likelihood/topology comparison summaries | Raw distances | None | Each tree is built from one metric at a time, so cross-metric scale differences do not enter tree construction. The p-Epitope matrix may receive the existing tiny seeded tie-breaking perturbation, but it is not min-max normalized. | Main text |
| Main Mantel correlation table | Mantel r, matrix-permutation p-values, and strain-bootstrap intervals | Raw distances | None | Mantel r is a correlation statistic and is invariant to positive linear rescaling of either distance matrix. | Main text |
| Rank-based Mantel sensitivity table | Spearman Mantel r, matrix-permutation p-values, and strain-bootstrap intervals | Raw distances | None | Spearman Mantel r is computed on ranks and is invariant to monotonic transformations; it is a sensitivity statistic rather than the primary matrix-comparison estimand. | Supplement |
| Subtype contrast table and plot | H3N2 minus H1N1 Mantel-r contrasts | Raw distances | None | The contrast is computed on subtype-specific correlation estimates, which are scale-invariant within each matrix comparison. | Main text |
| Descriptive Pearson correlation table | Pearson r and strain-bootstrap intervals | Raw distances | None | Pearson r is scale-invariant under positive linear transformations, so normalization is not needed for the statistic and raw inputs keep provenance clear. | Supplement-ready audit |
| Descriptive Spearman correlation table | Spearman rho and strain-bootstrap intervals | Raw distances | None | Spearman rho is rank-based and is included as a monotonic-association sensitivity summary, not as independent-observation inference. | Supplement |
| Leave-one-strain-out influence | Mantel-r and descriptive Pearson-r point-estimate changes | Raw distances | None | The influence summaries compare correlation estimates after removing one strain from the existing raw matrices. | Supplement-ready audit |
| Alignment and sequence sensitivity | Correlation-estimate changes, distance-matrix correlations, and within-metric absolute distance changes | Raw distances | None | Association changes are correlation-based; scale-sensitive absolute distance changes are reported only within the same distance metric. | Supplement-ready audit |
| Secondary sequence-distance sensitivity | Correlation-estimate changes versus primary sequence-distance comparators | Raw distances | None | Secondary metrics are compared through scale-invariant association estimates and are not used for distance-tree construction in the current analysis. | Supplement-ready audit |
| Correlation plot | Unique-pair scatterplot axes | Min-max normalized distances | By metric across all displayed unique off-diagonal pairs | Normalization is used only to display differently scaled metrics on comparable 0-1 axes; the figure caption and axis labels identify the normalized scale. | Main text |
| Normalized distance-table artifact | `distance_table_normalized` and `dist-df-normalized.rds` | Min-max normalized distances | By metric across all unique off-diagonal pairs | This derived artifact supports labeled visual or audit use only and is not the input for raw-distance analyses. | Supplement-ready audit |

## Pipeline Notes

- `distance_table` and `full_distance_table` are raw unique-pair tables.
- `distance_table_normalized` is a labeled derived artifact for visual or audit
  use only.
- `distance_scale_audit` and `distance_scale_audit_table_file` provide a
  supplement-ready audit table that can be added to `products/supplement.qmd`
  if requested.
- Normalized sensitivity analyses are not currently generated or interpreted.

