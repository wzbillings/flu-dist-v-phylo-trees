# Analysis Decision Log

This log records consequential research, analysis, reproducibility, and
publication decisions for the influenza distance-metrics manuscript project.
Entries are based on documented human responses unless otherwise stated.

## 2026-06-05 - Strain-Panel Coverage and Generalizability Summary

**Decision:** Add panel-summary outputs that describe the expert-selected
H1N1/H3N2 analysis panel by subtype, isolation year, cohort-season Fluzone
vaccine-component status, unique off-diagonal pair counts, temporal span, and
available cartographic/phylogenetic distance coverage. The panel may be
described as expert-selected to span historical H1N1/H3N2 antigenic space, but
the manuscript must state that it is not a representative survey of all
influenza strains, cohorts, or eras.

Use `data/fluzone_vaccine_strains_by_season.csv` as the canonical source for
Fluzone vaccine components during the cohort-study seasons. One-time permission
was granted to standardize the source short names `SG/16` to `Sing/16` and
`BR/18` to `Bris/18`; after that edit, the file returns to read-only source
status. Seasons with any `excluded_from_high_dose = TRUE` component are treated
as trivalent high-dose Fluzone seasons; seasons with no excluded components are
treated as quadrivalent high-dose Fluzone seasons.

**Rationale:** The manuscript needs an auditable panel description so readers
can see how much of the H1N1/H3N2 historical and vaccine-component space is
covered without implying representative sampling. The vaccine-status definition
is limited to cohort-season Fluzone components because that is the approved
canonical source now available in the repository.

**Evidence / citation:** ZB approvals in chat on 2026-06-05; implementation in
`R/panel-summary.R`, `_targets.R`, `products/manuscript.qmd`,
`products/supplement.qmd`, and related tests.

**Alternatives considered:** Omit vaccine status; infer vaccine status from
short names without a source file; keep source aliases and handle them through
analysis-code aliasing; or use broader wording that implies representativeness.
These alternatives were rejected because ZB supplied a canonical vaccine-source
file, approved short-name standardization in that file, and approved explicit
generalizability limitation language.

**Impact:** Panel coverage, vaccine matching, and space-coverage outputs should
be regenerated through the targets pipeline. Final publication claims should
remain panel-specific unless independent panels or additional cohorts are
analyzed.

## 2026-06-05 - Rank-Based Mantel and Descriptive Spearman Sensitivity

**Decision:** Add rank-based Spearman Mantel correlations and descriptive
Spearman correlations as supplement-reported sensitivity analyses for monotonic
but nonlinear relationships between ML-tree cophenetic distance and the four
primary candidate distance metrics. Keep Pearson-based Mantel correlations as
the primary distance-matrix estimand. Do not reinterpret or strengthen main
manuscript claims from the rank-based outputs unless ZB reviews the regenerated
results and approves a claim change.

For the rank-based Mantel sensitivity, use the same unique off-diagonal
strain-pair representation, subtype-stratified overall matrix-label
permutation scheme, random seed conventions, and test/full Mantel permutation
settings as the primary Mantel table. Estimate intervals with the same
strain-unit bootstrap settings used for the existing Mantel and descriptive
correlation intervals; for Spearman rows, compute the bootstrap statistic as a
weighted correlation of rank-transformed pairwise distances.

**Rationale:** The primary Mantel table evaluates linear matrix association on
the raw distance scales, while rank-based Spearman summaries check whether the
qualitative association pattern is robust to monotonic but nonlinear
relationships. Keeping rank-based outputs supplement-only avoids changing the
primary estimand or over-expanding the main result set.

**Evidence / citation:** ZB approval in chat on 2026-06-05; implementation in
`R/plots-and-tables.R`, `_targets.R`, `R/manuscript.R`,
`products/manuscript.qmd`, `products/supplement.qmd`, and related tests.

**Alternatives considered:** Report Spearman Mantel only without descriptive
Spearman; omit bootstrap intervals for rank-based outputs; promote rank-based
results to the main table or main table footnote; or interpret rank-based
changes as updated manuscript conclusions. ZB requested descriptive Spearman
and bootstrap intervals. Main-table promotion and interpretation changes were
deferred to preserve the current primary analysis framing.

**Impact:** The supplement should include a rank-based Mantel sensitivity table
and a descriptive Spearman correlation table. Full-mode regeneration is needed
before publication-scale interpretation. Rank-based results are sensitivity
analyses and require human review before changing the manuscript's scientific
claims.

## 2026-06-05 - Raw and Normalized Distance Handling

**Decision:** Use raw distances for analyses where one distance metric is used
at a time, including neighbor-joining tree construction. For cross-metric
comparisons, keep raw distances when the statistic is scale-invariant, including
Mantel correlations, descriptive Pearson correlations, subtype Mantel-r
contrasts, leave-one-strain-out correlation changes, and alignment or
secondary-distance correlation-change sensitivity summaries. Use min-max
normalized distances only for explicitly labeled visual or audit outputs where
scale differences would otherwise affect display or interpretation. When
normalization is used, normalize by distance metric across all displayed or
audited unique off-diagonal strain pairs, not separately by subtype.

Normalized sensitivity analyses are not reported in the current manuscript pass.
Instead, maintain an audit record and supplement-ready table so normalized
handling can be documented in the supplement if requested by a reviewer.

**Rationale:** Single-metric analyses do not compare across distance units, so
normalizing before tree construction would change the input scale without
answering the scientific question. Correlation-based cross-metric summaries are
invariant to positive linear rescaling, so raw inputs preserve provenance
without changing the statistic. Normalization is useful for Figure 1-style
visual comparison because temporal, sequence, cartographic, and tree distances
have different units and ranges.

**Evidence / citation:** ZB approval in chat on 2026-06-05; implementation in
`R/distance-calc.R`, `R/plots-and-tables.R`, `_targets.R`,
`docs/raw-vs-normalized-distance-handling.md`, `products/manuscript.qmd`, and
related tests.

**Alternatives considered:** Report normalized sensitivity analyses immediately,
normalize separately by subtype and metric, or use a single global
normalization across all metrics. Normalized sensitivity analyses were deferred
because the approved current need is an audit, not an additional results
family. Subtype-specific normalization was rejected because it would erase
within-metric cross-subtype range differences in displays. Global
all-metric normalization was rejected because larger-unit metrics would dominate
the scale and obscure smaller-range metrics.

**Impact:** Raw and normalized outputs must be clearly labeled. Figure 1 and
the normalized distance-table artifact are visual/audit products, not inputs to
the primary Mantel, Pearson, subtype-contrast, influence, sensitivity, or
neighbor-joining analyses. The audit table can be added to the supplement if a
reviewer requests explicit raw-versus-normalized reporting.

## 2026-06-04 - Supplementary Secondary Sequence-Distance Metrics

**Decision:** Add supplementary-only secondary sequence-distance metrics for
sensitivity analysis: amino-acid Hamming distance, normalized optimal string
alignment distance, $p$-all-epitope distance, and a BLOSUM62-derived
dissimilarity. Keep temporal distance, Grantham distance, $p$-Epitope distance,
and cartographic distance as the primary distance metrics. Do not add
nucleotide sequence information or nucleotide-distance analyses to this paper.
Do not implement LogDet unless a reviewer or coauthor requests it. Do not use
the secondary distances to build neighbor-joining trees unless full-mode
outputs show a major difference or a reviewer/coauthor requests those topology
checks.

For all secondary sequence metrics, use the same pairwise-deletion rule as the
primary sequence-derived distances: compare standard amino-acid residues and
exclude gaps, ambiguous or missing residue codes, and known nonstandard
residue codes from the pair-specific denominator. Compute Hamming distance as
the proportion of comparable aligned sites that differ. Compute optimal string
alignment after pairwise deletion and normalize by the number of comparable
sites. Compute $p$-all-epitope using the same H1N1/H3N2 epitope site sets as
$p$-Epitope, but average the five site-specific Hamming distances rather than
taking the maximum. Compute the BLOSUM62-derived dissimilarity per comparable
site as `S(a,a) + S(b,b) - 2 * S(a,b)` from the BLOSUM62 substitution score
matrix, then average across comparable sites.

Compare each secondary metric's association with the primary ML-tree
cophenetic distance against the closest primary sequence-distance comparator:
Hamming, optimal string alignment, and BLOSUM62 against Grantham; and
$p$-all-epitope against $p$-Epitope. Use the existing 0.10 absolute
correlation-change flag threshold for the supplementary sensitivity table.

**Rationale:** Grantham remains the primary whole-sequence amino-acid metric
because it weights substitutions by physicochemical difference and replaced
Hamming as the primary metric in an earlier decision. The secondary metrics
provide simpler or alternative sequence-distance checks without changing the
primary analysis. Nucleotide analyses are excluded because nucleotide changes
are harder to map directly to antigenicity in this manuscript's framing.
LogDet and secondary distance-tree topology checks add complexity and are
deferred unless the sensitivity outputs or external feedback justify them.

**Evidence / citation:** Human approval in chat on 2026-06-04; implementation
in `R/secondary-sequence-distance.R`, `_targets.R`,
`products/manuscript.qmd`, `products/supplement.qmd`, and related tests.

**Alternatives considered:** Include nucleotide p-distance, add LogDet or
other model-based distances immediately, promote Hamming back to a primary
metric, or build secondary-metric neighbor-joining trees in the first pass.
These alternatives were rejected or deferred for the current manuscript pass.

**Impact:** The supplement should report secondary sequence-distance
sensitivity outputs, but main-text interpretation should remain based on the
primary four metrics until full-mode outputs and coauthor review are complete.

## 2026-06-04 - Matrix-Only Alignment Sensitivity Analysis

**Decision:** Add a matrix-only alignment sensitivity analysis for the current
amino-acid workflow. Keep MUSCLE as the primary alignment method and use
ClustalW and ClustalOmega, through the existing `msa` package, as alternative
protein alignments. Recompute sequence-derived Grantham and p-epitope distance
matrices under each alternative alignment, then compare their matrix
associations with the existing primary MUSCLE ML-tree cophenetic distances.
Use the existing 0.10 absolute correlation-change threshold as the material
alignment-distance flag. Do not include nucleotide alignments or nucleotide
distances. Do not refit alternative ML trees or distance-based trees in this
first pass.

Also compare FLU-family amino-acid model-test results across primary and
alternative alignments. Flag model sensitivity only when the primary selected
tree model's log-likelihood loss under an alternative alignment exceeds the
existing 10% model-performance tolerance, paralleling the common-versus-subtype
model-selection rule.

**Rationale:** The manuscript needs to document whether the primary
sequence-distance conclusions depend on the MUSCLE amino-acid alignment.
ClustalW and ClustalOmega are available through the locked `msa` dependency and
do not require external command-line binaries. A matrix-only first pass is
computationally modest and directly evaluates the alignment dependence of
sequence-derived distance matrices while preserving the current ML tree as the
reference outcome. Alternative tree refitting is scientifically broader and
more expensive, so it is deferred unless a coauthor or reviewer asks for it.

**Evidence / citation:** Human approval in chat on 2026-06-04; implementation
in `R/alignment.R`, `R/alignment-sensitivity.R`, `R/utils.R`, `_targets.R`,
`products/manuscript.qmd`, `products/supplement.qmd`, and related tests. The
`msa` package documents integrated support for MUSCLE, ClustalW, and
ClustalOmega protein alignments.

**Alternatives considered:** Include nucleotide alignments or nucleotide
distance metrics; refit ML trees and distance-based trees under each
alternative alignment immediately; treat any model-name change as a sensitivity
flag even if the primary selected model remains within the approved
performance tolerance; or block the pipeline on unavailable external aligner
binaries. Nucleotide analyses were rejected as out of scope for the current HA
protein-distance question. Alternative tree refits are deferred. Model changes
within tolerance are recorded but not flagged as material. External aligners
are not used in this pass.

**Impact:** The supplement should report alignment-distance and model-selection
sensitivity tables. Alignment sensitivity flags require human interpretation
before final manuscript claims are strengthened or changed. Publication-scale
interpretation still requires full-mode regeneration and review.

## 2026-06-04 - Complete-Sequence and Missing-Data Sensitivity Audit

**Decision:** Add a targeted sequence-completeness and missing-data audit for
the current analysis panel. Treat standard amino-acid residues as comparable
for sequence-derived distances. Exclude alignment gaps, ambiguous or missing
residue codes (`B`, `J`, `X`, `Z`, and `?`), and known nonstandard amino-acid
codes (`O` and `U`) from the relevant pair-specific denominators rather than
counting them as zero-distance substitutions. Surface any other aligned
sequence characters as unexpected audit findings for human source-data review.
Add a complete-sequence-only matrix sensitivity that removes included
non-full-length strains from the existing distance matrices, plus a cheap
pairwise-deletion versus complete-deletion sensitivity for Grantham and
p-epitope distances. Use the existing 0.10 absolute correlation-change
threshold as the material-change flag for the targeted matrix-association
sensitivity summaries. Do not refit alignments, ML trees, distance trees, or
cartography maps in this pass.

**Rationale:** `MS/85` is currently expected to be the only included
non-full-length analysis strain, and the prior leave-one-strain-out influence
analysis is relevant but not sufficient to document sequence completeness and
missing-data assumptions. Pairwise deletion preserves the primary estimand while
avoiding treating unknown residues or gaps as biological matches. Complete
deletion is cheap enough to report as a denominator sensitivity, whereas
complete-sequence-only tree refitting and gap-as-state phylogenetic analyses
would answer broader questions and may require heavier computation.

**Evidence / citation:** Human approval in chat on 2026-06-04; implementation
in `R/grantham-distance.R`, `R/p-epitope-calculator.R`,
`R/sequence-sensitivity-audit.R`, `_targets.R`, `products/manuscript.qmd`, and
`products/supplement.qmd`. Pairwise and complete deletion are documented
treatments for gaps, missing data, and ambiguous states in sequence-distance
software documentation such as `phangorn` and MEGA.

**Alternatives considered:** Leave p-epitope gap/ambiguity handling delegated
to `phangorn::dist.hamming()` without explicit denominator control; treat
ambiguous residues as mismatches; treat gaps as a distinct biological state;
delete non-comparable sites globally for the primary analysis; or refit
complete-sequence-only alignments and trees immediately. The first two options
were rejected for transparency and biological interpretability. The latter
three are deferred unless requested by ZB, coauthors, or reviewers.

**Impact:** Sequence-derived distance calculations now handle gaps and
ambiguous residues consistently across Grantham and p-epitope metrics. The
supplement should report the aligned character audit, complete-sequence matrix
sensitivity, and pairwise-versus-complete deletion sensitivity. Any unexpected
sequence characters or sensitivity flags require human review before
publication-scale interpretation.

## 2026-06-03 - Leave-One-Strain-Out Matrix Influence Analysis

**Decision:** Add a leave-one-strain-out influence analysis for the current
distance-matrix associations. For each subtype, remove one strain at a time
from the existing ML-tree cophenetic and candidate distance matrices, then
recompute the observed association with temporal, Grantham, p-epitope, and
cartographic distance on the remaining unique off-diagonal strain pairs. Use an
absolute change of at least 0.10 correlation units as the flagging threshold.
Retain both Mantel-r and descriptive Pearson point-estimate rows in the derived
summary, focus the supplement display on Mantel r, and include a short
sensitivity-analysis section in the main manuscript. Do not refit alignments,
ML trees, or cartography maps for this pass.

**Rationale:** The approved robustness-check set called for a
leave-one-strain-out influence check for matrix correlations. Subsetting the
existing matrices directly evaluates whether the reported matrix-association
estimates are driven by individual strains while preserving the current
targets pipeline and avoiding a much more expensive, scientifically different
delete-and-refit analysis.

**Evidence / citation:** Human approval in chat on 2026-06-03; implementation
in `R/influence-analysis.R`, `_targets.R`, `products/manuscript.qmd`, and
`products/supplement.qmd`.

**Alternatives considered:** Refit the MSA, ML tree, and cartography-derived
objects after each strain deletion; report only descriptive Pearson estimates;
or omit the analysis from the main manuscript when flags do not materially
affect interpretation. The detailed refit-style leave-one-out analysis is
deferred unless requested by a coauthor or reviewer.

**Impact:** The sensitivity analysis should help diagnose whether individual
strains strongly influence temporal, Grantham, p-epitope, or cartographic
agreement with ML-tree distance. Flagged rows require scientific interpretation
before strengthening or changing manuscript claims. Full-mode outputs should be
regenerated before final publication-scale interpretation.

## 2026-06-01 - Strain Provenance and Inclusion Audit Outputs

**Decision:** Add pipeline-derived strain provenance and inclusion audit
outputs for the current sequence source table. Report all source sequence rows,
analysis inclusion status, rows excluded because they are absent from the stored
cartography maps, source/accession completeness, source and alignment
full-length status, and unique off-diagonal pair counts by subtype and distance
metric. Place detailed provenance and audit tables in the supplement-ready
outputs, with a concise main-manuscript methods note.

**Rationale:** The publishability review identified missing accession,
provenance, exclusion, and pair-count reporting as a reproducibility gap. The
new protein sequence CSV already contains source/accession fields and
cartography-overlap inclusion status can be derived from the pipeline, so these
audits should be generated from targets rather than manually maintained.

**Evidence / citation:** Implementation in `R/provenance-audit.R`,
`_targets.R`, `products/manuscript.qmd`, and `products/supplement.qmd`.
The analysis set remains governed by the 2026-06-01 cartography-overlap
inclusion decision above.

**Alternatives considered:** Put the full provenance table in the main
manuscript, report only aggregate counts, or leave accession/source details as
unrendered internal metadata. Supplement-first reporting keeps the main text
concise while preserving auditability.

**Impact:** The tables make excluded source rows and pair counts explicit
without changing the analysis panel or distance calculations. Final
source-citation wording still requires the deferred human accession/source
citation review before submission.

## 2026-06-01 - Protein Sequence Source Table and Cartography-Overlap Inclusion

**Decision:** Replace the prior sequence workbook plus virus-name metadata CSV
inputs with `data/UGAFluVac-sequences.csv`, a protein-sequence source table
containing strain names, short names, sequence source databases, accession
numbers, source full-length flags, and hemagglutinin protein sequences. Keep the
analysis set restricted to strains present in the stored H1N1 and H3N2
cartography maps. Sequence records not present in the stored maps should be
retained in the source table for provenance but excluded from the current
distance-metric analysis. Remove nucleotide-sequence requirements from the
pipeline because nucleotide sequences are not used in the current distance,
tree, or manuscript outputs. Replace the previous hard-coded full-length rule
with source full-length flags plus post-alignment non-gap protein-length
checks.

**Rationale:** Cartographic distance comparison is a central contribution of
the manuscript, and the available first-pass cartography maps define the strain
panel that can support that comparison. The updated sequence table supplies
accession/source provenance in a simpler auditable format and avoids carrying
unused nucleotide fields. Full-length status should be auditable from the
sequence data and alignment rather than inferred from a hard-coded strain name.

**Evidence / citation:** Human instruction in chat on 2026-06-01; implementation
in `R/data-processing.R`, `R/alignment.R`, `_targets.R`,
`products/manuscript.qmd`, and `products/supplement.qmd`.

**Alternatives considered:** Expand the analysis panel to every sequence row in
the new source table, keep requiring nucleotide sequences, preserve the deleted
virus-name metadata CSV, or keep a strain-name exception for full-length status.

**Impact:** The current analysis panel remains 18 H1N1 strains and 21 H3N2
strains. Seven source sequence rows are excluded from the current analysis
because they are not present in the stored cartography maps. `MS/85` is the
only source and analysis-panel sequence flagged as non-full-length in the new
source table and alignment audit.

## 2026-06-01 - Antigenic Cartography Input Diagnostics

**Decision:** Extract and report antigenic-cartography diagnostics from the
stored H1N1 and H3N2 `.ace` inputs. Mandatory reported diagnostics are map
dimensionality, strain/antigen count, serum count, stored projection count,
stored map stress, measured and missing titer counts, missingness percentage,
lower-bound/censored titer counts, and observed titer range. Optimizer history
beyond stored stress and projection count should be treated as best-effort and
reported as unavailable when absent from the stored `.ace` input rather than
blocking manuscript work. Report the diagnostics in the supplement first, with
a concise main-methods note.

**Rationale:** The cartographic distance methods need to be auditable from the
available first-pass `.ace` inputs while preserving the boundary that full
cartography-generation reproducibility will be addressed later. Supplement-first
placement documents the map and titer metadata without overloading the main
results.

**Evidence / citation:** Human approval in chat on 2026-06-01; implementation
in `R/cartography-diagnostics.R`, `_targets.R`, `products/manuscript.qmd`, and
`products/supplement.qmd`.

**Alternatives considered:** Require complete optimizer history before
publication work can continue, omit unavailable diagnostics entirely, or place
the full diagnostics table in the main text.

**Impact:** The `.ace` files remain read-only first-pass inputs. A separate
paper is planned to describe construction of the human-serum antigenic
cartographies. The current manuscript should acknowledge that using the same HAI
dataset to build the maps and compare cartographic distances is a limitation,
while noting that no comparable alternative is currently available for this
strain panel.

## 2026-05-31 - ML Tree Branch Support and Bootstrap Topology Stability

**Decision:** Add ML-tree robustness diagnostics using nonparametric
site-bootstrap resampling of each fitted ML tree. Use 10 bootstrap replicates in
test mode and 1,000 replicates in full mode. Report branch-support and
bootstrap topology-stability summaries in the supplement only for now; do not
add support labels to the main ML-tree figure unless a later manuscript decision
requires it. Include companion topology-stability metrics such as RF,
normalized RF, weighted RF, branch-score, and path distances from each
bootstrap replicate tree to the optimized ML tree.

**Rationale:** Branch support directly addresses reviewer concerns that the ML
tree topology needs robustness evidence. Supplement-first reporting avoids
overloading the main figure while still making support and stability evidence
auditable. Test-mode settings keep development checks feasible, while full-mode
settings are intended for publication-scale interpretation.

**Evidence / citation:** Human approval in chat on 2026-05-31; implementation
in `R/tree-building.R`, `R/plots-and-tables.R`, `_targets.R`,
`products/manuscript.qmd`, and `products/supplement.qmd`.

**Alternatives considered:** Put bootstrap support labels directly on the main
ML-tree figure, report branch support only without topology-stability companion
metrics, or use topology-stability summaries only if branch-support plotting
proved fragile.

**Impact:** Test-mode outputs should be treated as pipeline checks only.
Publication claims about ML-tree robustness require regenerating the support
targets in `FLU_TARGETS_MODE=full` and reviewing the resulting supplemental
tables.

## 2026-05-31 - Tree Comparison and ML Model Selection

**Decision:** Rebuild tree-comparison reporting around the decrease in log
likelihood from the optimized ML tree. SH-test p-values and Robinson-Foulds
distances remain in the main tree-comparison table as secondary checks. SPR,
weighted RF, Kuhner-Felsenstein/branch-score, and path distances should be
reported in the supplement as topology-distance diagnostics. Do not implement a
standard likelihood-ratio test for distance-tree versus ML-tree topology
comparisons because these fixed topology comparisons are not regular nested
model comparisons with a usual chi-square reference distribution.

**Rationale:** Delta log likelihood directly answers how much sequence-model
support is lost by imposing each distance-derived neighbor-joining topology.
The SH test is an established likelihood-based topology comparison, and RF-type
distances are useful topology checks, but neither should displace the primary
likelihood-loss summary. A standard LRT would overstate inferential precision
because the compared trees are alternative fixed/optimized topologies rather
than nested parameter restrictions.

**Evidence / citation:** Human approval in chat on 2026-05-31; implementation
in `R/tree-building.R`, `R/plots-and-tables.R`, `_targets.R`, and
`products/manuscript.qmd`. Method citations added for the FLU model,
neighbor joining, SH test, RF distance, and AICc model selection.

**Alternatives considered:** Keep RF distance as the main tree-comparison
metric, frame SH test as the primary comparison, implement an LRT, or omit
additional topology-distance diagnostics.

**Impact:** The main manuscript tree table should be generated as
`results/Tables/tree-comparison-table.rds`. The supplemental topology table
should be generated as `results/Tables/tree-topology-distance-table.rds`.

## 2026-05-31 - Common-versus-Subtype Tree Model Rule

**Decision:** Compare FLU-family amino-acid substitution models with and without
invariant sites and discrete gamma rate variation. In full mode, evaluate gamma
category counts 4, 8, 12, 16, 20, 24, 28, and 32. Select the AICc-supported
model within each subtype, but prefer a common model for both subtypes when its
log-likelihood loss is no more than 10% relative to the subtype-specific
AICc-supported model in either subtype. If no common model satisfies this rule,
use subtype-specific models and discuss why subtype-specific model choices may
be appropriate.

**Rationale:** A common model improves comparability between H1N1 and H3N2 tree
comparisons, but forcing a common model would be inappropriate if it materially
degrades fit for one subtype. The 10% loss rule operationalizes the human
preference for a common model without allowing substantial subtype-specific
performance loss.

**Evidence / citation:** Human instruction in chat on 2026-05-31.

**Alternatives considered:** Always use a single preselected FLU + G + I model,
always use subtype-specific AICc minima, or choose based on one subtype only.

**Impact:** Tree-fitting targets now use `tree_model_choice`, and the model
selection table should be regenerated before final tree-comparison claims are
interpreted. Test-mode model-selection results are for pipeline verification
only.

## 2026-05-30 - Subtype Contrast for Matrix-Association Estimates

**Decision:** Quantify subtype differences in distance-metric performance as
the H3N2 Mantel correlation minus the H1N1 Mantel correlation for each distance
metric. The primary uncertainty interval should use a Bayesian bootstrap over
strain units via `bayesboot`. Fisher-z and independent Mantel permutation
contrast calculations should be reported in the supplement as sensitivity
analyses, not as the primary estimand.

**Rationale:** Subtype-specific differences are central to the manuscript, and
the contrast should preserve the already-approved matrix-association framework
without treating individual pairwise distances as independent observations. The
Fisher-z sensitivity is familiar but overconfident for pairwise matrices, while
the permutation sensitivity is matrix-aware but tests a different
no-association null rather than direct equality of subtype-specific
correlations.

**Evidence / citation:** Human approval in chat on 2026-05-30; Mantel and
Bayesian-bootstrap methods cited in the manuscript and supplement.

**Alternatives considered:** Descriptive differences only, Fisher-z
approximation as the primary method, and permutation contrast as the primary
method.

**Impact:** The targets pipeline should produce a main-table and contrast plot
for H3N2-minus-H1N1 Mantel differences, plus supplement-ready sensitivity
summaries. Main-text wording should remain cautious and strain-panel-specific
until full-mode results and manuscript revisions are complete.

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
