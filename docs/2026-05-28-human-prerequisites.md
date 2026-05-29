# Human Prerequisites Before Substantive Project Work

Date: 2026-05-28

Purpose: this file lists the questions and human actions that must be resolved
before an agent should begin substantive manuscript revision, analysis changes,
pipeline restructuring, data release preparation, or journal submission work for
this project.

This is not a revision plan. It is a gatekeeping checklist. The goal is to avoid
making scientific, authorship, privacy, or publication decisions by accident.

## Current Status - Updated 2026-05-29

Human responses have been reviewed in
`docs/2026-05-28-human-prerequisites-responses.md`, and consequential decisions
have been recorded in `docs/analysis-decision-log.md`.

### Resolved Enough to Guide Initial Work

- [x] Project identity: expand and revise into a publishable original empirical
  research article.
- [x] Revision posture: rebuild the manuscript from scratch using the current
  draft as an ideas and methods skeleton.
- [x] Primary audience: influenza vaccine scientists, including readers without
  deep technical/computational backgrounds.
- [x] Scope: focus on the current H1N1/H3N2 heterologous human strain panel.
- [x] Subtype framing: keep H1N1 and H3N2 in one paper, with subtype differences
  treated as a central finding.
- [x] Neuraminidase analyses: future work only; not required for this manuscript.
- [x] Data-sharing posture for current local work: the author reports that the
  analysis data may remain in the repository and be shared, subject to final
  public-release curation.
- [x] Ethics posture for current local work: the author reports that repository
  data contain no human participant data or identifying metadata and that the
  institutional IRB designated the data as exempt nonhuman-subjects data.
- [x] Raw data boundaries: `data/full-sequences.xlsx` and
  `data/UGAFluVac-virus-names.csv` are raw source files and must not be edited.
- [x] Derived output boundaries: other files in `data/`, all files in `results/`,
  and `.docx`/`.pdf` files in `products/` should be treated as derived outputs
  that need to be reproducible.
- [x] Distance-matrix repair direction: use unique off-diagonal strain pairs;
  exclude diagonals and duplicates from reported correlations and figures.
- [x] Main matrix-comparison direction: use Mantel correlations while preserving
  descriptive correlation summaries where useful for nontechnical interpretation.
- [x] Main tree-comparison direction: use change in log likelihood as the main
  tree-comparison metric, with LRT considered if constructible.
- [x] Workflow direction: convert the project to a `targets` pipeline.
- [x] Output direction: regenerate manuscript-ready outputs from the repaired
  workflow rather than treating existing Word documents as final.
- [x] Publication intent: journal submission is desired; current venue family is
  influenza/virology, with specific journal fit to be verified later.
- [x] Preprint intent: desired after initial coauthor approval; author will
  approve and upload to medRxiv manually.
- [x] Seed bibliography settings: target 100 citations if possible; Markdown
  plus BibTeX; verified metadata; 2010 onward except necessary or seminal older
  sources; preprints allowed; retracted papers excluded.
- [x] Code work authorization: agent may modify `R/`, create `_targets.R`,
  update Quarto/manuscript files, update `renv.lock`, and regenerate outputs.
- [x] Code work off-limits areas: `archive/`, raw data files, `renv/`,
  `reviews/`, `scripts/`, and housekeeping files.
- [x] Derived output cleanup: tracked derived files may be cleared out now, and
  derived filenames may change if equivalent files can eventually be
  regenerated.
- [x] Expensive computation strategy: expensive/stochastic steps may run after
  sequential pipeline verification, with a test mode for faster settings in
  `_targets.R`.
- [x] Placeholder accession and IRB inputs: accession/source citation table and
  IRB details can be placeholders for now while work proceeds.
- [x] Cartography boundary for first pipeline pass: treat existing `.ace` files
  as inputs; integrate cartography-generation code later.

### Still Requires Human Confirmation

- [ ] Accession table and source citations for sequences.
- [ ] IRB details for the study ethics subsection.
- [ ] Coauthor/contributor table and approval timing.
- [ ] Final sensitivity-analysis requirements.
- [ ] Branch support, model-selection, alignment, cartography, and missing-data
  checks required before submission.
- [ ] Final public-release plan for a curated repository that does not expose
  internal or unneeded files through commit history.
- [ ] Venue-specific policy and fit verification before submission or preprint.
- [ ] Human approval of the default robustness-check set before treating it as
  required analysis.

### Minimum Answer Needed to Begin Seed Annotated Bibliography

Resolved 2026-05-29. The seed annotated bibliography can begin with these
parameters:

- [x] Target count: 100 citations if possible.
- [x] Output format: Markdown annotations plus BibTeX.
- [x] Verification level: verified metadata.
- [x] Time window: 2010 onward, with older sources only if necessary or seminal.
- [x] Preprints: include, clearly labeled.
- [x] Exclusions: no retracted papers.

### Proposed Default Robustness-Check Set - Pending Human Approval

These checks are proposed for the rebuilt analysis plan. They should not be
treated as required until the human approves, revises, or rejects them.

#### Distance-Matrix Comparison

- [ ] Use only unique off-diagonal strain pairs for descriptive pairwise plots
  and summaries.
- [ ] Report Mantel correlations as the primary matrix-comparison statistic.
- [ ] Include descriptive Pearson correlations for interpretability, clearly
  labeled as descriptive rather than independent-observation inference.
- [ ] Consider Spearman/Mantel rank correlations as a secondary check for
  monotonic but nonlinear relationships.
- [ ] Compare normalized and unnormalized distance matrices where normalization
  changes interpretation.
- [ ] Use a fixed random seed and documented permutation count for Mantel tests,
  with fast test-mode and full-analysis settings.

#### Tree Construction and Tree Comparison

- [ ] Re-run ML tree fitting from pipeline inputs rather than relying on fixed
  precomputed tree objects.
- [ ] Reproduce and document model comparison for each subtype.
- [ ] Decide whether the same substitution model must be used for both subtypes
  or whether subtype-specific best models are acceptable.
- [ ] Report the primary tree-comparison metric as change in log likelihood,
  with LRT considered only if the comparison is statistically valid for the
  fitted models.
- [ ] Retain SH-test and RF-distance summaries as secondary topology-comparison
  checks if still methodologically appropriate.
- [ ] Add branch-support or topology-stability checks if feasible in full mode.

#### Alignment and Sequence Handling

- [ ] Keep the current amino-acid MSA workflow as the primary sequence analysis
  unless revised later.
- [ ] Verify how gaps, ambiguous residues, and incomplete sequences are handled.
- [ ] Run a leave-one-strain-out influence check for matrix correlations, or at
  minimum flag strains with high influence on the reported correlations.
- [ ] Check whether any incomplete strain should be excluded as a sensitivity
  analysis after verifying the exact strain names and sequence completeness.

#### Antigenic Cartography Inputs

- [ ] Treat existing `.ace` files as first-pass cartography inputs.
- [ ] Report available map metadata such as dimensionality, stress, number of
  optimizations, and titer handling if available from the `.ace` files.
- [ ] Confirm strain-name harmonization between cartography, sequence, and tree
  objects.
- [ ] Defer full cartography-generation reproducibility until the upstream code
  is provided or reconstructed.

#### Panel, Generalizability, and Missing Data

- [ ] Document how the strain panel was selected and why it is scientifically
  useful despite not being representative of all influenza strains.
- [ ] Report counts of strains, pairwise comparisons, excluded records, and
  missing values at each pipeline stage.
- [ ] Separate H1N1 and H3N2 primary results, then compare subtype-specific
  patterns explicitly.
- [ ] Avoid combining subtypes unless there is a clearly justified secondary
  analysis.

#### Pipeline and Computational Reproducibility

- [ ] Add a test mode in `_targets.R` with reduced optimization/bootstrap/
  permutation settings.
- [ ] Add a full mode in `_targets.R` for publication-quality computations.
- [ ] Record seeds and computational settings for stochastic steps.
- [ ] Make file targets explicit for manuscript-ready tables and figures.

## Blocking Questions

### 1. Project Identity

- [ ] Should this project be treated as a publishable journal manuscript, a
  dissertation appendix, a reproducible teaching/project archive, or some
  combination?
- [ ] If it should become a journal manuscript, what is the intended article
  type: descriptive empirical paper, computational methods note, virology
  short report, reproducibility report, or another type?
- [ ] What is the one-sentence claim the paper is allowed to make?
- [ ] What claims should the paper explicitly avoid making, even if they sound
  attractive?
- [ ] Is the current goal to salvage the existing draft, rebuild it from
  scratch, or document why it should not be pursued as a formal publication?

### 2. Scope and Target Audience

- [ ] Who is the primary intended audience: influenza virologists, antigenic
  cartography researchers, evolutionary biologists, computational biologists,
  vaccine modelers, or dissertation committee readers?
- [ ] Should the manuscript focus narrowly on the strain panel in this
  repository, or should it be expanded with additional cohorts, strains, or
  datasets?
- [ ] Should H1N1 and H3N2 remain in one paper, or should the analysis be
  reframed around the subtype-specific differences?
- [ ] Should neuraminidase sequence or inhibition data be considered required
  for publication, or only discussed as future work?

### 3. Data Rights, Privacy, and Release Constraints

- [ ] Which sequence records came from GISAID, UniProt, GenBank, collaborators,
  or other sources?
- [ ] Are any GISAID-derived sequences, metadata, alignments, distance matrices,
  or derived files currently tracked in this repository?
- [ ] Are the GISAID-derived materials allowed to remain in the repository if
  the repository is public or pushed to GitHub?
- [ ] Are the HAI data, Racmacs `.ace` files, cartographic coordinates, or
  derived cartographic distances approved for sharing with collaborators,
  reviewers, or the public?
- [ ] Are any human participant data, potentially identifying metadata, site
  information, dates, or restricted study details present in the repository?
- [ ] What exact data availability statement is legally and ethically permitted?
- [ ] Should any data files be removed from Git tracking, moved to controlled
  storage, or replaced with synthetic/toy data before public release?

### 4. Authorship, Collaboration, and Permission

- [ ] Who must approve any manuscript submission?
- [ ] Who must approve any public repository, data release, preprint, or
  supplemental material?
- [ ] Who qualifies as an author, contributor, data provider, or acknowledged
  collaborator?
- [ ] Are there collaborators or data generators whose approval is required
  before analyzing, publishing, or sharing derived results?
- [ ] Are there unpublished, embargoed, dissertation, collaborator, or reviewer
  constraints that limit what can be written or released?

### 5. Scientific Claim Boundaries

- [ ] Is the project allowed to claim that temporal distance is unsuitable in
  general, or only that it performed poorly in this selected H1N1 panel?
- [ ] Is the project allowed to claim that genetic and antigenic evolution do
  not agree, or only that the selected HA phylogenetic distances and HAI-based
  cartographic distances diverged in this dataset?
- [ ] Is the project allowed to make any claim about individual immune response
  determinants, or is that outside the evidence available here?
- [ ] Which analyses were planned before seeing the results, and which were
  exploratory or post hoc?
- [ ] What limitations must be stated plainly in the abstract, discussion, or
  data availability section?

### 6. Analysis Decisions Requiring Human Approval

- [ ] Should the distance-matrix comparison be rebuilt using unique off-diagonal
  strain pairs only?
- [ ] Should the main inference use Mantel/permutation-style tests, descriptive
  correlations only, or another method?
- [ ] Should diagonal and duplicate pairwise distances be excluded from all
  reported correlations and figures?
- [ ] Should Pearson, Spearman, Mantel correlation, matrix regression, or another
  metric be the primary comparison?
- [ ] Should p-values and confidence intervals be reported at all for the
  distance comparisons?
- [ ] What sensitivity analyses are required before the manuscript is worth
  resubmitting or submitting?
- [ ] Should ML tree fitting be rerun from source data, or should existing
  precomputed ML tree objects be treated as fixed historical outputs?
- [ ] What branch support, model-selection, alignment, cartography, and missing
  data checks are required?

### 7. Reproducibility and Pipeline Boundaries

- [ ] Should the project be converted to a `targets` pipeline?
- [ ] If not `targets`, what should be the authoritative workflow: Makefile,
  Quarto project render, numbered scripts, or another system?
- [ ] Which outputs should be reproducibly regenerated, and which should remain
  archived historical artifacts?
- [ ] Which files in `data/` are raw/protected source data and must remain
  immutable?
- [ ] Which files are derived outputs and should be moved or regenerated outside
  protected source-data locations?
- [ ] Should the existing rendered Word documents be kept as historical outputs
  or regenerated from Quarto after the workflow is repaired?

### 8. Publication Strategy

- [ ] Is a journal submission actually desired?
- [ ] If yes, should the target be a broad-scope journal, influenza/virology
  journal, evolutionary biology journal, computational biology journal, or
  methods venue?
- [ ] Is a preprint desired, and if so, who must approve it first?
- [ ] Are APCs, open-access requirements, data-sharing policies, or institutional
  repository requirements constraints on venue choice?
- [ ] Should the work be positioned as a modest descriptive result rather than a
  major computational biology contribution?

## Human Actions Required Before Agent Work Begins

### Must Do First

- [ ] Decide the project identity and write a one-paragraph scope statement.
- [ ] Confirm whether the repository may be public, private, or collaborator-only.
- [ ] Confirm GISAID, UniProt, GenBank, collaborator, and human-study data
  permissions.
- [ ] Confirm which data files are raw/protected and which outputs may be
  regenerated or moved.
- [ ] Confirm who must approve publication, preprint posting, repository release,
  and data/code sharing.
- [ ] Confirm the allowed strength of the main manuscript claim.
- [ ] Confirm whether the agent is authorized to modify analysis scripts,
  manuscript text, pipeline files, bibliography files, and generated outputs.

### Should Do Before Analysis Changes

- [ ] Decide the primary distance-comparison method.
- [ ] Decide whether p-values and confidence intervals should be reported for
  distance-matrix comparisons.
- [ ] Decide the required sensitivity analyses.
- [ ] Decide whether ML trees should be rerun or treated as historical outputs.
- [ ] Decide whether to implement `targets` or another workflow manager.
- [ ] Decide where analysis decisions should be recorded. If no better location
  exists, use `docs/analysis-decision-log.md`.

### Should Do Before Manuscript Rewriting

- [ ] Select the intended audience.
- [ ] Select the intended venue type, even if not a final journal.
- [ ] Decide which claims must be weakened, removed, or retained.
- [ ] Provide or approve the citation areas that must be checked.
- [ ] Confirm whether acknowledgements and authorship information are complete
  and safe to publish.
- [ ] Confirm required ethics, funding, competing-interest, author-contribution,
  code-availability, and data-availability statements.

### Should Do Before Public Release or Submission

- [ ] Review the repository for restricted data, credentials, private URLs, and
  unpublished sensitive material.
- [ ] Confirm that no GISAID-restricted materials are redistributed in violation
  of applicable access terms.
- [ ] Confirm that human-study data sharing complies with consent, IRB, data-use
  agreements, and collaborator expectations.
- [ ] Confirm that all generated figures, tables, and manuscript outputs can be
  regenerated from documented code or are clearly labeled as historical outputs.
- [ ] Confirm that coauthors have reviewed the final manuscript and repository
  release package.
- [ ] Confirm the final venue and verify current journal-specific policies before
  submission.

## Information the Human Should Provide to the Agent

- [ ] A short scope statement for the project.
- [ ] A list of files or directories the agent must not read, modify, commit, or
  discuss.
- [ ] A data provenance note for each source dataset.
- [ ] A decision on whether the repository is intended to be public.
- [ ] A list of required approvers for publication, preprint, data release, and
  code release.
- [ ] A preferred publication path, or explicit instruction to defer publication
  strategy.
- [ ] Approval or rejection of `targets` as the workflow system.
- [ ] Approval or rejection of rerunning expensive phylogenetic computations.
- [ ] Any known prior reviewer comments, decision letters, dissertation committee
  feedback, or collaborator constraints that should guide revisions.

## Work the Agent Should Not Start Until These Are Resolved

- Manuscript claim strengthening or reframing.
- Journal-specific submission preparation.
- Data release preparation.
- Public repository cleanup.
- Raw data modification.
- Removal or relocation of tracked data files.
- Rebuilding the analysis pipeline in a way that changes output provenance.
- Recomputing expensive or stochastic phylogenetic analyses.
- Updating bibliography entries based on unverified assumptions.
- Writing a final data availability statement.

## Notes for Future Work

- Any change to inclusion criteria, distance metrics, model family, alignment
  strategy, sensitivity analyses, data availability, repository release status,
  or publication target should be recorded in a decision log.
- Any unresolved answer should be labeled `Unclear`, `Not found`, `Assumed`, or
  `Requires human confirmation` in future review or revision documents.
- Secrets and credentials should never be pasted into repository files or chat.
  If code requires credentials, document only the environment variable names.
