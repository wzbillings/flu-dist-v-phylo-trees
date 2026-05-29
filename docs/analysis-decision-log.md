# Analysis Decision Log

This log records consequential research, analysis, reproducibility, and
publication decisions for the influenza distance-metrics manuscript project.
Entries are based on documented human responses unless otherwise stated.

## 2026-05-29 - Manuscript Identity and Rebuild Strategy

**Decision:** The project will be expanded and revised into a publishable
journal manuscript, specifically an original empirical research article. The
existing draft will not be patched into submission form; it will be rebuilt from
scratch using the current draft as an ideas and methods skeleton.

**Rationale:** The current draft captures useful ideas and a rough methods
structure, but the intended output is now a full research article rather than a
term-project or dissertation-appendix artifact.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, section 1.

**Alternatives considered:** Dissertation appendix, reproducible archive,
descriptive report, computational methods note, or abandoning formal
publication.

**Impact:** Future manuscript work should prioritize argument reconstruction,
expanded background, reproducible analyses, and journal-facing reporting rather
than line edits to the existing draft.

## 2026-05-29 - Intended Audience and Article Positioning

**Decision:** The primary intended audience is influenza vaccine scientists. The
paper should be readable to influenza scientists without technical backgrounds
and should avoid overly methods-heavy language.

**Rationale:** The project aims to influence how influenza vaccine and virology
researchers choose antigenic distance metrics, not primarily to present a new
computational method.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, scope and audience section.

**Alternatives considered:** Antigenic cartography researchers, evolutionary
biologists, computational biologists, vaccine modelers, dissertation committee
readers.

**Impact:** Background, figures, methods, and results should emphasize
interpretability and biological meaning. Technical details should be complete
but can be shifted to methods or supplement when appropriate.

## 2026-05-29 - Core Claim Boundaries

**Decision:** The working core claim is that phylogenetic distance is a poor
proxy for antigenic change in influenza. The manuscript should avoid claiming
that results are representative of all flu strains or all populations. It should
also avoid mechanistic claims that cartographic mismatch proves hemagglutinin
sequence is not the only factor in individual immune responses.

**Rationale:** The author wants the paper to make a general argument about H1N1
and H3N2 distance metrics, while still treating panel sampling and
generalizability as important limitations.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, project identity and claim
boundaries sections.

**Alternatives considered:** A narrower claim limited only to the selected
strain panel; broader mechanistic claims about individual immune-response
determinants.

**Impact:** Future drafting must distinguish metric performance from mechanism.
Generalization beyond the study panel requires careful wording and literature
support. Limitation language should be explicit in the abstract and discussion.

## 2026-05-29 - Study Scope and Subtype Framing

**Decision:** The manuscript will focus on the current heterologous influenza
strain panel rather than adding more cohorts or datasets now. H1N1 and H3N2 will
remain in one paper, and differences between subtype-specific results are a
central contribution. Neuraminidase sequence or inhibition data are future work,
not required for this manuscript.

**Rationale:** The current strain panel is treated by the author as a major
human cohort strain panel for this kind of analysis. Larger strain sets may
allow phylogenetic analyses but usually do not support comparable antigenic
cartography.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, scope and target audience
section.

**Alternatives considered:** Expanding to more cohorts, splitting H1N1 and H3N2,
or requiring neuraminidase analyses before publication.

**Impact:** The pipeline should support the current panel first. Future-work
sections can discuss neuraminidase and larger multiplex-assay-enabled panels.

## 2026-05-29 - Data Sharing, Ethics, and Public Release Assumptions

**Decision:** The author reports that all analysis data are permitted to remain
in the repository and be shared; HAI data, Racmacs `.ace` files, cartographic
coordinates, and derived cartographic distances are approved for sharing. The
author reports no human participant data, identifying metadata, site
information, dates, or restricted study details in the repository. The author
reports that the data have been designated exempt nonhuman-subjects data by the
institutional IRB.

**Rationale:** These answers remove several immediate data-safety blockers for
local analysis and documentation work, while still requiring provenance and IRB
details before final manuscript submission or public release.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, data rights/privacy section.

**Alternatives considered:** Controlled data storage, removing files from Git
tracking, replacing data with synthetic examples, or treating the repository as
collaborator-only.

**Impact:** Code and reproducibility work may proceed locally, but final public
release still requires a curated repository and a verified data availability
statement. The methods must include a study ethics subsection once IRB details
are supplied.

## 2026-05-29 - Raw and Derived Data Boundaries

**Decision:** `data/full-sequences.xlsx` and
`data/UGAFluVac-virus-names.csv` are raw source data and must remain immutable.
All other files in `data/`, all files in `results/`, and any `.docx` or `.pdf`
files in `products/` should be treated as derived outputs that need to be
reproducible.

**Rationale:** This protects the stated source data while allowing derived
artifacts to be regenerated by a pipeline.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, reproducibility and pipeline
boundaries section.

**Alternatives considered:** Treating all tracked data as raw or preserving
existing derived outputs as fixed historical artifacts.

**Impact:** Future code work should not edit the two raw source files. Moving,
deleting, or regenerating tracked derived outputs still requires explicit
confirmation because it changes repository contents and provenance.

## 2026-05-29 - Primary Distance-Matrix Analysis Direction

**Decision:** The distance-matrix comparison should be rebuilt using unique
off-diagonal strain pairs only. Diagonal and duplicate pairwise distances should
be excluded from reported correlations and figures. The primary matrix
comparison should use Mantel correlations, while preserving descriptive
correlations where they help a nontechnical audience interpret calibration-like
linear relationships.

**Rationale:** The current full-matrix Pearson correlation approach treats
non-independent duplicate and diagonal distances as if they were independent
observations. Mantel-style matrix comparisons better match the structure of the
data while retaining interpretable correlation summaries.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, analysis decisions section;
prior publishability review,
`reviews/manuscript-publishability-review.md`.

**Alternatives considered:** Ordinary Pearson correlations on all matrix
entries, descriptive-only correlations, Spearman correlations, matrix
regression.

**Impact:** Figure 1 and any associated tables should be regenerated after the
analysis is corrected. P-values and confidence intervals should at least appear
in a supplement if they are straightforward to calculate and defensible.

## 2026-05-29 - Primary Tree-Comparison Direction

**Decision:** Change in log likelihood should be the main tree-comparison metric,
with an LRT considered if constructible. Existing SH-test and RF-distance
outputs remain relevant but should be revisited in the new analysis design.

**Rationale:** The author wants tree comparison framed around likelihood
differences where possible.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, analysis decisions section.

**Alternatives considered:** RF distance as the main comparison metric,
SH-test-centered framing, or topology-distance-only reporting.

**Impact:** The tree-comparison code and results tables should be redesigned
around the selected primary metric after feasibility is checked.

## 2026-05-29 - Reproducible Workflow Direction

**Decision:** The project should be converted to a `targets` pipeline. All
analysis steps should be part of the pipeline from sequence cleaning through
final graphs. Cartography code needs to move into this repository in a future
step. Existing rendered Word documents should be regenerated after workflow
repair.

**Rationale:** A pipeline is needed to make the analysis auditable from source
data through manuscript-ready outputs.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, reproducibility section.

**Alternatives considered:** Makefile, numbered scripts, Quarto-only workflow,
or treating existing RDS files as fixed outputs.

**Impact:** Future code work should prioritize `_targets.R`, reusable functions,
dependency cleanup, and reproducible output generation.

## 2026-05-29 - Publication and Preprint Strategy

**Decision:** Journal submission is desired. The target venue type is influenza
or virology journal, with aspirational targets including PNAS or Proceedings B
if feasible, followed by something around the level of mBio. PLOS Computational
Biology, Journal of Virology, and Influenza and Other Respiratory Viruses are
currently named as backups. A preprint is desired after initial coauthor
approval, and the author will approve and upload to medRxiv manually. The author
strongly prefers open-access journals and journals run by independent societies.

**Rationale:** The author wants the work to influence influenza researchers'
choice of distance metrics while preserving open-access and society-journal
preferences.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`, publication strategy
section.

**Alternatives considered:** No journal submission, methods venue, broad-scope
journal, computational biology venue, or no preprint.

**Impact:** Venue-specific claims, competitiveness, formatting, and policies
must be verified later from official journal sources before submission planning.
No preprint or submission actions should occur without explicit human approval.

## 2026-05-29 - Deferred Human Inputs

**Decision:** The following items are deferred and require future human input:
accession table and citations, IRB details, coauthor table, collaborator approval
timing, final sensitivity analyses, branch-support/model-selection/alignment/
cartography/missing-data check requirements, and final public-release package
contents.

**Rationale:** The author identified these as future-step inputs or areas where
additional explanation is needed.

**Evidence / citation:** Human prerequisites responses,
`docs/2026-05-28-human-prerequisites-responses.md`.

**Alternatives considered:** Proceeding as if these details were resolved.

**Impact:** These items should not block initial exploratory bibliography work
or pipeline scaffolding, but they do block final manuscript claims, public
release, submission, and any definitive data availability/ethics statements.

## 2026-05-29 - Seed Bibliography Parameters

**Decision:** The seed annotated bibliography should target 100 citations if
possible. Output should include Markdown annotations and BibTeX entries.
Verification level should be verified metadata. The preferred time window is
2010 onward, with older sources included only when necessary or seminal, such
as Hobson 1972. Preprints may be included. Retracted papers are excluded.

**Rationale:** The bibliography is intended to expand the manuscript background
from the current introduction and limited bibliography while keeping citation
quality high enough for later manuscript development.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Smaller seed set, unverified candidate list,
post-2010 sources only, peer-reviewed-only list, excluding preprints.

**Impact:** Bibliography work can begin without waiting for the remaining code
or manuscript prerequisites. Preprints should be clearly labeled, and retraction
status should be checked where feasible rather than assumed.

## 2026-05-29 - Code Work Authorization and Off-Limits Areas

**Decision:** The agent is authorized to modify `R/`, create `_targets.R`,
update Quarto/manuscript files, update `renv.lock`, and regenerate outputs.
The following areas are off-limits for code work: `archive/`, raw data files,
`renv/`, `reviews/`, `scripts/`, and housekeeping files. Raw data files are
those previously identified as `data/full-sequences.xlsx` and
`data/UGAFluVac-virus-names.csv`.

**Rationale:** The project needs a reproducible pipeline and code cleanup, but
raw data, archived material, local environment machinery, review artifacts,
script utilities, and repository housekeeping should not be altered during code
work.

**Evidence / citation:** Human response in chat on 2026-05-29; prior raw-data
boundary decision in this log.

**Alternatives considered:** Documentation-only work, broader repository
reorganization, modifying raw inputs, or modifying archived files.

**Impact:** Future code work may proceed in the authorized areas. Any proposed
change outside those boundaries requires fresh human confirmation.

## 2026-05-29 - Derived Output Cleanup Authorization

**Decision:** Tracked derived files may be cleared out now. Derived filenames may
change as long as equivalent files can eventually be regenerated. Historical
versions are considered recoverable from Git history if needed.

**Rationale:** The pipeline should become the source of truth for derived data,
tables, figures, and manuscript outputs rather than preserving obsolete tracked
artifacts as active project state.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Leaving derived files in place while scaffolding
the new pipeline, moving files to archive, or retaining old filenames exactly.

**Impact:** Before deleting or moving derived files, the agent should verify
absolute paths and stay within authorized locations. The first pipeline pass
should make clear which outputs are not yet regenerated.

## 2026-05-29 - Expensive Computation Strategy

**Decision:** Expensive or stochastic computations may be run, but the pipeline
should include a test mode with faster settings that can be turned on or off in
`_targets.R`. Work should proceed sequentially, verifying that the pipeline runs
before launching expensive computations.

**Rationale:** The project needs reproducible expensive steps eventually, but
development should not be blocked by slow ML tree fitting or high-bootstrap
tree-comparison settings.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Scaffolding only, always running full expensive
settings, or retaining precomputed ML objects as fixed outputs.

**Impact:** `_targets.R` should expose a clear test/full mode switch. Early
pipeline checks should use fast deterministic settings where possible.

## 2026-05-29 - Placeholder Inputs for Accession and Ethics Details

**Decision:** Accession/source citation details and IRB details may be treated
as placeholders for now. The author will provide the accession table and IRB
details while other work proceeds.

**Rationale:** These details are required for final manuscript integrity but do
not need to block initial pipeline restructuring, seed bibliography work, or
methods scaffolding.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Blocking all manuscript/code work until tables and
IRB details are supplied.

**Impact:** Any manuscript text involving accession provenance or ethics must be
clearly marked with placeholders until the human supplies verified details.

## 2026-05-29 - Cartography Input Boundary for First Pipeline Pass

**Decision:** The first pipeline pass should treat existing `.ace` cartography
files as inputs. Cartography-generation code will be integrated later after the
current pipeline is remodeled.

**Rationale:** This allows the project to make progress on reproducibility using
available cartography outputs without blocking on cartography-generation code.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Blocking pipeline work until cartography-generation
code is available, or attempting to reconstruct cartography generation from
scratch immediately.

**Impact:** The first `_targets.R` design should explicitly distinguish source
inputs, `.ace` cartography inputs, and derived outputs. Later pipeline revisions
can replace `.ace` inputs with upstream cartography-generation targets.

## 2026-05-29 - Robustness Checks To Propose

**Decision:** The agent should propose a default robustness-check set for human
approval rather than requiring the human to specify all checks in advance.

**Rationale:** Some robustness-check categories require technical translation
before the author can approve or reject them.

**Evidence / citation:** Human response in chat on 2026-05-29.

**Alternatives considered:** Blocking code work until the human independently
specifies all robustness checks.

**Impact:** The proposed default set should be reviewed before being treated as
required analysis.

## 2026-05-29 - Explicit Sequence-Metadata Join

**Decision:** The rebuilt pipeline should join sequence records to virus
metadata by normalized strain name rather than binding rows by positional order.

**Rationale:** The old workflow used an order-dependent `bind_cols()` step,
which made provenance fragile if either source file changed order. An explicit
join makes mismatched, missing, or duplicate strain names auditable.

**Evidence / citation:** Pipeline refactor in `R/data-processing.R` and
`_targets.R`.

**Alternatives considered:** Preserve the old order-based bind or treat the
legacy `data/viruses_used.rds` object as an input.

**Impact:** The pipeline now validates metadata matching before sequence
cleaning. Any future source-data name change should fail loudly instead of
silently recoding strains.

## 2026-05-29 - Pipeline-Derived Output Location

**Decision:** New pipeline-derived intermediate files should be written under
`results/derived/`; manuscript-ready figures and tables remain under
`results/Figures/` and `results/Tables/` for continuity with the existing
manuscript.

**Rationale:** `data/full-sequences.xlsx` and
`data/UGAFluVac-virus-names.csv` are immutable source inputs, so derived
sequence-cleaning and alignment outputs should no longer be regenerated into
`data/`.

**Evidence / citation:** `_targets.R` file-target paths.

**Alternatives considered:** Continue writing derived RDS files into `data/` or
rename all output directories in one large cleanup.

**Impact:** Legacy derived files in `data/` are not used by the first-pass
pipeline and can be removed later after regeneration is verified.

## 2026-05-29 - Test and Full Pipeline Settings

**Decision:** The pipeline defaults to `FLU_TARGETS_MODE=test`, with seed 370,
99 Mantel permutations, 199 Mantel Bayesian-bootstrap draws, 100 SH-test
bootstrap resamples, 4 gamma categories for model-test scaffolding, and a fast
NNI ML-tree fitting strategy. Full mode is selected with
`FLU_TARGETS_MODE=full` and uses 9,999 Mantel permutations, 4,000 Mantel
Bayesian-bootstrap draws, 1,000,000 SH-test resamples, 20 gamma categories, and
`pml_bb` tree fitting.

**Rationale:** The test mode is intended to verify pipeline integrity without
launching publication-scale stochastic or expensive computations.

**Evidence / citation:** `make_analysis_settings()` in `R/utils.R`.

**Alternatives considered:** Always use publication-scale settings or preserve
precomputed ML tree objects as fixed historical outputs.

**Impact:** Test-mode numerical results should not be interpreted as final
manuscript evidence. Full-mode results must be regenerated before final
scientific claims are updated.

## 2026-05-29 - Deterministic p-Epitope Distance Perturbation

**Decision:** The small perturbation used to make p-epitope distance matrices
usable for neighbor joining should be seeded, symmetric, and leave matrix
diagonals at zero.

**Rationale:** The old workflow added random noise to the whole matrix without
local seed control, which could alter diagonal entries and reduce
reproducibility.

**Evidence / citation:** `perturb_matrix()` in `R/utils.R` and
`prepare_distances_for_tree_building()` in `R/tree-building.R`.

**Alternatives considered:** Remove perturbation entirely or preserve the old
unseeded full-matrix perturbation.

**Impact:** Tree-building remains compatible with the existing p-epitope
distance workflow while making the stochastic adjustment auditable.

## 2026-05-29 - Figure 1 Normalization and Correlation Table

**Decision:** Figure 1 should show min-max normalized distance values and use a
colorblind-friendly subtype palette to distinguish H1N1 and H3N2. Correlation
numbers should be removed from the figure and reported in a separate table with
overall and subtype-specific descriptive Pearson correlations for each distance
metric versus ML-tree cophenetic distance.

**Rationale:** Normalization improves visual comparability across distance
metrics without changing correlations under a global positive linear
transformation. Moving correlation estimates to a table keeps the figure
visually cleaner and avoids repeating all pairwise matrix correlations.

**Evidence / citation:** User request on 2026-05-29; implementation in
`R/plots-and-tables.R` and `_targets.R`.

**Alternatives considered:** Keep unnormalized axes and in-panel labels; report
all pairwise correlations among distance metrics.

**Impact:** Figure 1 is a comparative visualization rather than a numeric
results table. The separate table should be used for exact reported
correlations.

## 2026-05-29 - Correlation Confidence Intervals

**Decision:** The distance-versus-cophenetic correlation table should include
95% intervals. The pipeline should use Bayesian bootstrap percentile intervals
through `bayesboot` when available, fall back to BCa intervals through `boot`,
and finally fall back to Fisher-z/Wald intervals if neither bootstrap path is
available.

**Rationale:** The interval calculation should be reproducible and explicit,
while still allowing the pipeline to run if optional bootstrap packages are not
available. The intervals are descriptive for unique off-diagonal strain pairs
and should not be interpreted as fully accounting for all matrix dependence.

**Evidence / citation:** User request on 2026-05-29; implementation in
`R/plots-and-tables.R`, `_targets.R`, and `renv.lock`.

**Alternatives considered:** Report point estimates only, use only Wald-type
intervals, or block the pipeline if `bayesboot` is unavailable.

**Impact:** Table output now records the interval method in the derived summary
object. Manuscript interpretation should describe these as descriptive
correlation intervals unless a stronger matrix-aware inference is later added.

## 2026-05-29 - Strain-Level Bootstrap Correlation Intervals

**Decision:** Treat the H1N1/H3N2 strains as a sample from a broader strain
universe for correlation-interval estimation. Bootstrap uncertainty should be
computed over strain units, not over individual pairwise distance rows. For
overall correlations, the bootstrap is stratified by subtype so the observed
H1N1/H3N2 pair contribution is preserved in each replicate.

**Rationale:** Pairwise distance rows are not independent because each strain
appears in multiple pairs. Drawing Bayesian-bootstrap weights over strains and
propagating them to pairwise distances as endpoint-weight products better
matches the intended sampling interpretation.

**Evidence / citation:** User clarification on 2026-05-29 that the strains
represent a broader universe; implementation in `R/plots-and-tables.R`.

**Alternatives considered:** Keep bootstrapping pairwise rows, use unstratified
strain weights for the overall rows, or report only Mantel-style permutation
tests without confidence intervals.

**Impact:** Correlation intervals now represent uncertainty under a
strain-sampling estimand. Point estimates are unchanged, but interval widths may
differ from row-level pairwise bootstraps because matrix dependence is now
represented in the resampling procedure.

## 2026-05-29 - Replace Primary Hamming Distance with Grantham Distance

**Decision:** Replace the primary whole-sequence Hamming distance metric with a
local Grantham amino-acid distance metric. Keep temporal distance, Grantham
distance, p-epitope distance, cartographic distance, and ML-tree cophenetic
distance as the primary comparison set.

**Rationale:** Prior work in this project indicated that whole-sequence Hamming
distance and p-epitope distance were highly correlated, making Hamming distance a
less informative additional primary metric. Grantham distance preserves a
whole-sequence amino-acid comparison while weighting substitutions by
physicochemical difference. The local implementation avoids adding an
unnecessary dependency and reports a per-comparable-site mean distance so pairs
with different numbers of gaps or ambiguous residues remain comparable.

**Evidence / citation:** Grantham R. 1974. Amino acid difference formula to help
explain protein evolution. Science 185(4154):862-864. DOI:
10.1126/science.185.4154.862. The `ahgroup/agdist` package documentation was
consulted as a reference implementation; both `agdist` and this repository use
AGPL-3.0-compatible licensing, but the needed Grantham functionality was
implemented locally in `R/grantham-distance.R`.

**Alternatives considered:** Keep Hamming distance as a co-primary sequence
metric, add `agdist` as a dependency, or compute summed rather than averaged
Grantham distances. These were rejected because the task is to replace Hamming,
the required functionality is small, and summed distances would be harder to
compare across pairs with different numbers of comparable aligned sites.

**Impact:** Main distance outputs, figures, tables, and manuscript terminology
now use Grantham distance instead of Hamming distance. p-epitope distance remains
in the primary analysis.

## 2026-05-29 - Primary Mantel Table and Supplementary Pearson Reporting

**Decision:** The main manuscript distance-matrix table should report Mantel
correlations for ML-tree cophenetic distance versus temporal, Grantham,
p-epitope, and cartographic distance. The table should include overall,
H1N1-specific, and H3N2-specific rows; unique off-diagonal strain pairs only;
Bayesian bootstrap 95% intervals over strain units; and two-sided
matrix-permutation p-values for the observed strain panel. Descriptive Pearson
correlations should move to the supplement.

**Rationale:** Mantel correlations are the primary matrix-comparison statistic
approved for this analysis, while Pearson correlations remain useful for
interpreting approximate linear agreement between distance scales. Overall rows
are retained because reviewers are likely to ask for them, but subtype-specific
rows remain central to interpretation.

**Evidence / citation:** Human approval in chat on 2026-05-29; implementation
in `R/distance-calc.R`, `R/plots-and-tables.R`, `_targets.R`,
`products/manuscript.qmd`, and `products/supplement.qmd`.

**Alternatives considered:** Main-table Pearson correlations with Mantel
p-values, subtype-only reporting, omitting Mantel intervals, and reporting
Bayesian-bootstrap tail probabilities as p-values. Bootstrap tail probabilities
were not treated as p-values because the permutation p-values are the
pre-specified apparent-sample matrix-correspondence test.

**Impact:** Manuscript claims should refer to Mantel agreement or matrix
association as the primary distance-matrix result. Pearson correlations should
be labeled descriptive and interpreted only as scale-calibration summaries.
