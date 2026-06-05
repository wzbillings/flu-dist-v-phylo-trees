# Test Suite Assessment

## Project testing assessment

This appears to be a non-package R analysis project for an influenza phylogenetics and distance-comparison manuscript. It has a Quarto manuscript under `products/`, reusable R helper files under `R/`, protected/source data under `data/`, generated or derived artifacts under `results/`, and an `renv.lock` environment lockfile. There is no top-level `DESCRIPTION`, so the test suite is set up as a lightweight non-package `{testthat}` suite.

The project does use `{targets}` through `_targets.R`. The pipeline sources `R/` and defines data import, sequence cleaning, alignment, distance calculation, tree analysis, manuscript figure/table, and Quarto render targets.

A small `{testthat}` suite is worthwhile because the repository contains stable reusable helper functions that encode assumptions about input validation, strain-name mapping, distance-matrix structure, sequence cleaning, and deterministic transformations. These are exactly the kinds of functions that can silently break downstream analyses if changed. A broad suite that runs the full pipeline would not be appropriate because it would be slow, data-dependent, and likely brittle during scientific revision.

## Testing scope

The implemented tests cover fast deterministic behavior:

- reusable validation and normalization utilities;
- analysis settings for test and full pipeline modes;
- distance-matrix tidying, ordering, validation, normalization helpers, and
  local Grantham amino-acid distance calculations;
- toy-data behavior for protein-sequence cleaning, cartography-overlap
  inclusion, and sequence-cleaning audits;
- strain provenance, accession/source completeness, inclusion/exclusion flow
  summaries, and unique-pair count audit helpers;
- protein-only alignment summaries and alignment-derived full-length checks;
- alignment-method validation and matrix-only alignment sensitivity helpers,
  including alternative-alignment distance estimates and model-selection flags;
- p-epitope site lookup and basic p-epitope distance-matrix structure;
- supplementary secondary sequence-distance metrics, including amino-acid
  Hamming distance, normalized optimal string alignment distance,
  p-all-epitope distance, BLOSUM62-derived dissimilarity, and their
  supplement-ready sensitivity table helper;
- sequence character classification, p-epitope ambiguous-residue exclusion,
  complete-deletion alternatives for sequence-derived distances, and
  complete-sequence-only matrix-sensitivity helpers;
- lightweight manuscript/table helper functions that do not render plots, tables, or documents;
- strain-level bootstrap helper behavior for distance-correlation intervals;
- Mantel matrix-comparison summaries that exclude diagonal and duplicate
  distance-matrix entries, including rank-based Spearman Mantel sensitivity
  summaries;
- raw-versus-normalized distance handling guardrails, including raw
  single-metric tree inputs, raw scale-invariant association summaries, and
  labeled normalized visual/audit outputs;
- leave-one-strain-out influence summaries, flags, and supplement-ready
  display helpers for matrix-association sensitivity analysis;
- antigenic-cartography diagnostics extracted from stored `.ace` map inputs,
  including map dimensions, strain/serum counts, stress metadata, and titer
  missingness summaries;
- tree-comparison table construction from toy tree-analysis objects, including
  delta log likelihood, SH p-values, RF distance, supplemental topology
  distances, and common-versus-subtype model-selection logic.
- ML-tree bootstrap support and topology-stability wrappers using small toy
  phylogenetic objects, plus supplement-ready support table formatting.

The tests intentionally do not cover:

- full `{targets}` pipeline execution;
- multiple sequence alignment runs;
- model fitting, maximum-likelihood tree fitting, SH tests, or bootstrap-heavy
  analyses on real alignments;
- Quarto rendering;
- private/source data contents;
- generated figures, rendered Word documents, or fragile plot snapshots;
- internet access or external data downloads.

Assumptions:

- Tests are run from the project root.
- Project dependencies are restored through `renv` or otherwise installed.
- Scientific results may change during revision; these tests guard helper behavior and structural assumptions, not final manuscript conclusions.

## Implemented test files

- `tests/testthat/helper-project.R`: sources the project helper files needed by the test suite.
- `tests/run-tests.R`: runs each test block and writes a timestamped Markdown log to `tests/logs/`.
- `tests/logs/.gitignore`: keeps generated Markdown test-run logs out of version control by default.
- `tests/testthat/test-utils.R`: tests general validation, normalization, matrix reshaping, strain-name replacement, temporal distances, seeded perturbations, and analysis settings.
- `tests/testthat/test-data-processing.R`: tests protein-only sequence
  cleaning, subtype splitting, duplicate source-record detection,
  cartography-overlap inclusion, and cleaning-audit summaries with inline toy
  data.
- `tests/testthat/test-provenance-audit.R`: tests strain provenance records,
  source/accession completeness labels, inclusion/exclusion flow summaries, and
  unique off-diagonal pair-count audits with inline toy data.
- `tests/testthat/test-alignment.R`: tests protein-only alignment audit
  summaries, full-length checks computed from aligned non-gap protein lengths,
  and validation of supported `msa` alignment methods.
- `tests/testthat/test-alignment-sensitivity.R`: tests matrix-only
  alternative-alignment sensitivity summaries and model-selection sensitivity
  flagging with inline toy distance matrices and model-test tables.
- `tests/testthat/test-distance-calc.R`: tests distance-set validation, distance-matrix ordering, combined distance tables, metric-level distance-table normalization, vector normalization, shared-strain alignment, unique off-diagonal pair extraction, and small deterministic Mantel permutation result structures.
- `tests/testthat/test-grantham-distance.R`: tests canonical Grantham matrix values, per-comparable-site averaging, explicit gap/ambiguous residue exclusion, and named symmetric matrix output.
- `tests/testthat/test-p-epitope-calculator.R`: tests p-epitope site lookup, invalid subtype handling, and basic p-epitope distance-matrix properties.
- `tests/testthat/test-secondary-sequence-distance.R`: tests supplementary amino-acid Hamming, normalized optimal string alignment, p-all-epitope, BLOSUM62-derived distance matrices, and secondary sequence-distance sensitivity table construction.
- `tests/testthat/test-sequence-sensitivity-audit.R`: tests aligned sequence
  character classification, complete-sequence-only matrix sensitivity, and
  pairwise-deletion versus complete-deletion missing-data sensitivity on toy
  aligned sequences.
- `tests/testthat/test-cartography-diagnostics.R`: tests pure helper behavior
  for `.ace`-derived map diagnostics, including titer missingness summaries,
  unavailable optimizer metadata reporting, and supplement-ready table labels.
- `tests/testthat/test-plots-and-tables.R`: tests stable display labels, distance labels, color palette names, normalized plotting data, raw descriptive distance-correlation inputs, supplement-ready raw-versus-normalized audit tables, weighted correlations, Fisher/Wald confidence interval guardrails, strain-level bootstrap weighting helpers, manuscript-ready Mantel summary structure, rank-based Mantel and descriptive Spearman sensitivity summaries, tree-comparison table helpers, supplemental topology-distance table helpers, and model-selection summaries.
- `tests/testthat/test-tree-building.R`: tests ML-tree support and
  topology-stability helper behavior with toy `phylo` and `multiPhylo` objects,
  and verifies that neighbor-joining preparation preserves raw non-p-epitope
  distance scales, avoiding expensive real-alignment bootstrap runs.
- `tests/testthat/test-subtype-contrast.R`: tests toy-matrix subtype contrast estimates, Bayesian-bootstrap contrast interval structure, Fisher-z and permutation sensitivity summaries, and manuscript-ready subtype contrast table/plot constructors.
- `tests/testthat/test-influence-analysis.R`: tests toy-matrix
  leave-one-strain-out influence estimates, influence flag validation, and
  supplement-ready influence table/plot constructors.

## How to run the tests

From the project root:

```r
source("tests/run-tests.R")
```

or from a shell:

```sh
Rscript tests/run-tests.R
```

The runner writes a timestamped Markdown log under `tests/logs/` with the run datetime, pass/fail/skip counts, and a row for each `test_that()` block. It exits with a non-zero status if any test block fails.

For ad hoc development runs that do not need a Markdown log:

```r
testthat::test_dir("tests/testthat")
```

To run one file while developing:

```r
testthat::test_file("tests/testthat/test-utils.R")
```

This is not an R package, so `devtools::test()` is not the primary command unless the project is later converted into a package.

## Limitations and recommendations

The current suite is intentionally small and should run without rebuilding targets or reading source data. Future tests would be useful for alignment summaries, tree metric wrappers, and manuscript table formatting, but those may require either heavier dependencies, more careful toy phylogenetic objects, or source-code refactoring outside `tests/`.

If this analysis grows, consider adding lightweight tests for target factory functions or pipeline metadata assumptions. Avoid tests that require running expensive model fits, modifying `_targets.R`, or depending on local/private data paths.
