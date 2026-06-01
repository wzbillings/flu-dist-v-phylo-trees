# Comparing antigenic distance metrics for influenza

This repository contains the analysis pipeline, manuscript files, and project
records for an in-progress influenza distance-metrics manuscript. The project
compares how several ways of measuring pairwise differences between influenza A
strains relate to maximum-likelihood phylogenetic tree distances for H1N1 and
H3N2 strain panels.

The current goal is to turn an earlier course-project and dissertation-appendix
analysis into an auditable journal manuscript for influenza vaccine scientists,
virologists, immunologists, and collaborators who need to understand when
phylogenetic or sequence distance is a useful proxy for antigenic change, and
when it is not.

## Intended uses

This repository is intended for:

- developing and auditing the influenza distance-metrics manuscript;
- rerunning the project analysis through the `targets` pipeline;
- tracing results from source inputs to derived tables, figures, and rendered
  manuscript outputs;
- reviewing documented analysis decisions and citation-support materials.

It is not a general-purpose influenza analysis package. Functions in `R/` are
project helpers for this manuscript workflow rather than a stable public API.

## Scientific focus

The analysis compares:

- temporal distance between virus isolation years;
- Grantham amino-acid distance;
- p-Epitope distance;
- antigenic cartography distance estimated from HAI response data;
- maximum-likelihood phylogenetic distances and distance-based
  neighbor-joining tree topologies.

The central working claim is cautious: for the current H1N1 and H3N2 strain
panel, sequence-derived and temporal distances can agree with phylogenetic
distance more strongly than cartographic distance does. That mismatch matters
because cartography is based on observed immune response data. Final scientific
claims should be made only after the full analysis mode has been regenerated and
the manuscript text has been reviewed against the current results.

## Contributors

Project contributors include Savannah L. Miller, Murphy H. John,
Amanda L. Skarlupka, Ted M. Ross, Justin Bahl, Andreas Handel, and
W. Zane Billings.

## Repository structure

- `_targets.R`: reproducible analysis pipeline from source inputs to derived
  outputs, manuscript tables, figures, and rendered documents.
- `R/`: reusable analysis functions for data cleaning, alignment, distance
  calculation, tree fitting, plots, tables, and manuscript rendering.
- `data/`: source inputs and legacy derived files. Treat
  `data/full-sequences.xlsx` and `data/UGAFluVac-virus-names.csv` as immutable
  raw/source inputs.
- `results/`: generated analysis outputs, including derived RDS files,
  manuscript-ready tables, and figures.
- `products/`: Quarto manuscript and supplement sources, bibliography files,
  citation style files, and rendered Word outputs.
- `docs/`: analysis decision log, project planning notes, and bibliography
  support materials.
- `tests/`: lightweight non-package `testthat` suite for deterministic helper
  behavior.
- `archive/`: legacy scripts and historical materials retained for provenance.

## Installation

This is an R analysis project managed with `renv`. The lockfile currently
records R 4.6.0 with Bioconductor 3.23 repositories.

From the repository root, restore the recorded package environment:

```r
renv::restore()
```

The project uses a `targets` pipeline, Quarto manuscript files, and several
R/Bioconductor packages recorded in `renv.lock`.

## Usage

Inspect the pipeline without running it:

```r
targets::tar_manifest()
```

Run the development pipeline:

```r
targets::tar_make()
```

By default, the pipeline uses `FLU_TARGETS_MODE=test`, which keeps stochastic
and expensive steps small enough for routine development checks, including
reduced Mantel, SH-test, subtype-contrast, and ML-tree support bootstrap
settings. Publication-scale settings are selected by setting
`FLU_TARGETS_MODE=full` before running the pipeline:

```r
Sys.setenv(FLU_TARGETS_MODE = "full")
targets::tar_make()
```

Rendered manuscript outputs are generated through pipeline targets from the
Quarto sources in `products/`.

## Tests

Run the lightweight helper-function test suite from the repository root:

```sh
Rscript tests/run-tests.R
```

The tests are intentionally smaller than a full pipeline rebuild. They guard
input validation, distance-matrix helpers, deterministic transformations,
manuscript table helpers, and other low-cost assumptions that can otherwise
silently affect downstream analyses.

## Development notes

- Do not manually edit raw inputs, generated tables, figures, rendered
  documents, or derived RDS files. Update the code or pipeline target that
  creates them.
- Keep reusable analysis logic in `R/` and pipeline definitions in `_targets.R`.
- Record consequential research, analysis, and publication decisions in
  `docs/analysis-decision-log.md`.
- Treat test-mode numerical results as pipeline checks, not final manuscript
  evidence.
- Regenerate full-mode outputs before updating final interpretation, submission
  materials, or public-release claims.

## Citation

No formal paper citation is available yet. If citing this work before
publication, cite the repository title, contributors, date accessed, and commit
hash, and confirm details with the project team. Once the manuscript is
published, cite the published article instead.

## License

This repository is licensed under the GNU Affero General Public License v3.0.
See [`LICENSE`](LICENSE) for the full license text.

In brief, AGPL-3.0 permits use, copying, modification, and redistribution under
copyleft terms, and includes source-sharing obligations for modified versions
made available over a network. This summary is provided for orientation only;
the license file controls.

## Disclaimer

Most of the code in this repository was generated with assistance from a large language model. All code and related materials were reviewed by a human, but they remain part of an in-progress influenza manuscript and should not be treated as legal, medical, regulatory, statistical, financial, or other professional advice.

The materials in this repository may be incomplete, incorrect, outdated, or unsuitable for a specific project, institution, journal, funder, dataset, regulatory context, or research purpose.

Users are responsible for independently reviewing, adapting, validating, and maintaining any analysis code, workflow materials, manuscript text, prompts, scripts, documentation, and other repository contents before using them in real research or public-facing work. Appropriate human review should be applied to scientific claims, statistical methods, privacy protections, authorship decisions, citations, compliance obligations, and public releases.

To the fullest extent permitted by applicable law, the maintainer disclaims responsibility for any loss, harm, liability, or other consequence arising from the use, modification, interpretation, or redistribution of this work by others.

Other listed contributors provided inspiration, advice, or feedback for the project, but they are not responsible for the current contents of this repository unless explicitly stated otherwise.

## Maintenance and contact

Maintainer: Zane Billings <wz.billings@gmail.com>

This repository is intended to remain practical and lightweight. Additions
should be clear, reusable, and grounded in actual academic or research workflow
needs.
