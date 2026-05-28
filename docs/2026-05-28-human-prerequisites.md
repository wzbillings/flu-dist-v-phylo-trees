# Human Prerequisites Before Substantive Project Work

Date: 2026-05-28

Purpose: this file lists the questions and human actions that must be resolved
before an agent should begin substantive manuscript revision, analysis changes,
pipeline restructuring, data release preparation, or journal submission work for
this project.

This is not a revision plan. It is a gatekeeping checklist. The goal is to avoid
making scientific, authorship, privacy, or publication decisions by accident.

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
