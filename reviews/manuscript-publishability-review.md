# Manuscript Publishability Review

Date: 2026-05-28

Scope: this report preserves the manuscript publishability critique produced
from static inspection of the repository, manuscript files, code, generated
outputs, bibliography, and project documentation. It is a review artifact only.
It does not modify manuscript, analysis, data, pipeline, generated output, or
bibliography files.

## 1. Executive Verdict

Verdict: `Not currently publishable; substantial revision needed`.

Current status: no-go for journal submission. The core idea is salvageable, but
the current manuscript is a short class-project-style report, not a defensible
publishable paper.

Main reason: the central quantitative comparison is under-specified and partly
invalid as reported. The Pearson correlations appear to use full symmetric
distance matrices including duplicate pairs and diagonals, treating
non-independent pairwise distances as independent observations. That makes the
reported confidence intervals and inferential framing unreliable.

Strongest aspect: the project has a real, coherent comparative question and
actual linked sequence/cartography/phylogenetic outputs.

Weakest aspect: methods/reproducibility. There is no `_targets.R`, no
authoritative pipeline, incomplete data provenance, missing method citations,
case-sensitive path risks, and some dependencies used in scripts are absent
from `renv.lock`.

Most realistic publication path: narrow the manuscript into a modest
computational/virology methods note or descriptive comparison, fix the
distance-matrix statistics, rebuild the reproducibility pipeline, and submit to
a broad-scope or evolutionary/virology venue only after data-sharing permissions
are settled.

Checks run during the original review: static repository inspection,
manuscript/Word extraction, citation-key check, generated table/figure
inspection, and read-only RDS summaries via `Rscript --vanilla`. Quarto render
and expensive tree-building commands were not run.

## 2. One-Paragraph Blunt Summary

This is not a publishable manuscript yet. It reads like exactly what the README
says it is: a term-project/dissertation-appendix artifact with interesting
preliminary results. The paper currently overreaches from a small virus panel
and weakly reported distance comparisons to broad claims about temporal
methods, antigenic evolution, and immune response determinants. The biggest
problem is not prose. It is that the statistical comparison of distance
matrices is not publication-safe, the phylogenetic reporting is too thin, the
code is not organized as a reproducible pipeline, and the data-sharing situation
may be nontrivial because the project uses GISAID-derived sequences and human
immunological data. The promising part is that there is a coherent, publishable
kernel: "How well do temporal, sequence, p-epitope, and cartographic distances
agree with HA phylogenetic structure for a defined influenza strain panel?"
Everything else should be cut or rebuilt around that.

## 3. Cognitive Verification Questions and Answers

### A. Contribution

Evidence status: supported by manuscript and repository.

The manuscript tries to compare temporal, Hamming, p-epitope, and
antigenic-cartography distances against ML phylogenetic distances for H1N1 and
H3N2 strains. The contribution is empirical/descriptive with a computational
methods flavor. It is explicit but underdeveloped. One-sentence claim to
existence: "In a selected influenza strain panel, sequence-derived distances
track HA phylogeny more closely than HAI-based antigenic cartography, suggesting
cartography captures immune-response structure not reducible to HA sequence
distance." Intended audience is influenza evolution, antigenic cartography, and
computational vaccinology readers. Why they would care is currently
under-argued.

### B. Fit

Evidence status: supported by README and manuscript; journal fit is expert
judgment requiring final verification.

This is not currently a PLOS Computational Biology-style paper. It lacks a new
method, large dataset, strong biological insight, or mechanistic result. It is
closer to a short computational virology note, reproducibility appendix, or
descriptive comparison. The README explicitly says it was "probably never
published formally in an academic journal," which is accurate for the present
version.

### C. Evidence

Evidence status: supported by manuscript, repository, and generated outputs.

The manuscript has two figures and one table. The table shows strong
differences, e.g. H1N1 temporal distance tree has Delta ll = 2239.4 and RF = 26,
while H1N1 Hamming has Delta ll = 0.6 and RF = 2. That is promising. But
evidence is incomplete: no branch support, no map uncertainty/stress reporting,
no accession/provenance table, no sensitivity analyses, no
unique-pair/permutation-based distance-matrix inference, and no justification
for model-selection compromise.

### D. Reproducibility and Auditability

Evidence status: supported by repository.

Not adequate. There is `renv.lock`, but no `_targets.R`, no Makefile, no Quarto
project file, and scripts are ad hoc. `R/tree-building.R` has
`FIT_ML_TREES <- FALSE`, so the core ML outputs are read from precomputed RDS
files rather than rebuilt. Scripts mix `results` and `Results`, which is fragile
on case-sensitive systems. `R/plots-and-tables.R` uses packages such as
`GGally`, `ggpp`, `ggtree`, `hgp`, `cowplot`, `patchwork`, `magick`, and
`ggpubr`; several are not present in `renv.lock` based on inspection.

### E. Argument

Evidence status: supported by manuscript.

The argument is visible but too thin. The introduction jumps from universal
influenza vaccines to a small distance-comparison exercise. The discussion turns
descriptive strain-panel findings into broad advice: "temporal methods should
be avoided." That is too strong. The paper loses the reader where "antigenic
distance," "phylogenetic distance," "immune response data," and "evolutionary
distance" are treated as interchangeable targets.

### F. Literature and Citation Integrity

Evidence status: supported by repository and external standards.

Citation keys used in the manuscript are present, excluding Quarto
cross-references. But important method citations are missing: FLU substitution
model, neighbor joining, SH test, RF distance, Mantel/permutation approaches for
distance matrices, antigenic cartography foundations, and phylogenetic reporting
standards. The manuscript cites `phangorn` and `msa`, but not enough of the
actual methods being used.

### G. Publishability

Evidence status: inferred from manuscript/repository plus publication standards.

Not currently publishable. Fixable if narrowed. Shortest credible path: reframe
as a descriptive computational comparison, fix the distance-matrix inference,
add reproducible pipeline and data availability/ethics statements, and target a
broad or applied venue. Highest-upside path: expand to more cohorts/strain
panels, include neuraminidase, external validation, antigenic map uncertainty,
and a stronger influenza-evolution framing. Work not worth doing now: polishing
prose, aiming high, or adding speculative universal-vaccine claims before the
core methods are repaired.

## 4. Reviewer Reports

### Editor-in-Chief / Desk-Rejection Reviewer

Likely decision now: desk reject or reject without review. The title and
abstract promise a publishable distance-metric comparison, but the manuscript
body is too short, too informal, and not sufficiently self-contained.
Desk-rejection risks: unclear novelty, weak article identity, inadequate
methods detail, overbroad conclusions, missing data availability, and unclear
journal fit. Prior rejection letter: not found. If this had been rejected
already, the rejection would probably be fair unless it was submitted explicitly
as a class-project appendix or technical note.

### Subject-Matter Expert Reviewer

Strength: comparing cartography, sequence distance, p-epitope, temporal
distance, and ML tree structure is a useful influenza question. Weakness: the
manuscript does not engage deeply enough with antigenic drift, immune imprinting,
antigenic cartography uncertainty, HA vs NA contributions, or why human
HAI-derived cartography might differ from HA phylogeny. The claim that HA
sequence is not the only determinant of immune response is plausible, but this
analysis does not directly test individual immune-response determinants.

### Methods / Statistics / Causal Inference Reviewer

Fatal-to-major issue: pairwise matrix entries are treated as independent in
`R/plots-and-tables.R`, and the code uses all `n^2` entries per subtype/method:
H1N1 has 324 rows per method where only 153 unique off-diagonal pairs exist;
H3N2 has 441 where only 210 unique off-diagonal pairs exist. Use upper/lower
triangle only and use Mantel/permutation or related distance-matrix methods. The
manuscript should not report ordinary Pearson CIs as if each pairwise distance
is independent. Also: no branch support, no sensitivity to alignment/model
choices, no explicit estimand, and no separation between descriptive and
inferential claims. Add citations: Legendre & Fortin on Mantel/distance-matrix
tests, Shimodaira-Hasegawa, Robinson-Foulds, Saitou-Nei, and FLU substitution
model.

### Results and Evidence Reviewer

Convincing: Table 1 has a coherent pattern, especially the H1N1
temporal/cartographic divergence from ML topology and the stronger Hamming
agreement. Weak: Figure 1 is visually useful but statistically overconfident;
duplicate/diagonal points and narrow CIs undermine it. Missing: all eight NJ
trees, accession/strain table, model-test table, branch support, antigenic map
diagnostics, and sensitivity analysis from `archive/sensitivity-analysis.docx`.

### Reproducibility / Open Science Reviewer

Blockers: no pipeline source of truth; key outputs precomputed; missing
dependencies; mixed path case; no render instructions; no data availability
statement; no decision log. Data safety risk: `R/data-processing.R` says raw
sequences came from GISAID and UniProt, and `archive/gisaidr-test.R` uses
environment variables for GISAID credentials. GISAID redistribution terms and
human HAI data permissions need explicit review before public release.

### Writing, Structure, and Argument Reviewer

The manuscript needs a full rebuild, not copyediting. The introduction should
define the concrete scientific problem, not start with universal vaccine
generalities. Methods need enough detail to reproduce the analysis. Results
should be organized around pre-specified comparisons. Discussion should stop
overclaiming and explain what cartography-vs-phylogeny disagreement can and
cannot mean.

### Citation and Literature Reviewer

No missing citation keys found for actual manuscript citations, but the
bibliography is underused and incomplete for methods. `dang2010` appears in the
bibliography and should likely support the FLU model. Add placeholders where not
yet verified: `[ADD CITATION: neighbor-joining method]`,
`[ADD CITATION: Shimodaira-Hasegawa topology test]`,
`[ADD CITATION: Robinson-Foulds tree distance]`,
`[ADD CITATION: distance matrix correlation / Mantel or permutation inference]`,
`[ADD CITATION: antigenic cartography foundational paper]`.

### Publication Strategy Reviewer

Do not submit as-is. Do not aim at PLOS Computational Biology, PLOS Pathogens,
JVI, or Bioinformatics unless substantially expanded. Realistic route is a
broad-scope reproducible descriptive paper after methodological repair. A
preprint is not advisable until the matrix-inference and data-sharing issues are
fixed.

### Project Manager / Revision Planner

Go/no-go: no-go for submission. Go for revision only if the author wants a real
paper rather than a dissertation appendix. First decision: paper identity.
Second: data permissions. Third: statistical repair. Fourth: pipeline.
Polishing comes last.

## 5. Major Problems Ranked by Severity

| Severity | Evidence | Why It Matters | Fix | Effort | Human Confirmation | Citation / Decision Log |
|---|---|---|---|---|---|---|
| Fatal | `R/plots-and-tables.R`; H1 324 and H3 441 matrix rows per method | Treats duplicate/diagonal non-independent distances as independent | Use unique off-diagonal pairs; use Mantel/permutation/distance-matrix methods | Medium | No | Citation yes; log yes |
| Major | `products/manuscript.qmd:181-200` | Conclusions overclaim beyond small descriptive design | Rewrite claims as descriptive and strain-panel-specific | Medium | Yes | Citation yes; log yes |
| Major | No `_targets.R`; ad hoc scripts | Reviewer cannot audit source-to-output chain | Add `targets` or equivalent pipeline | High | No | targets docs; log yes |
| Major | `R/tree-building.R:87`; ML tree fitting disabled | Core results depend on precomputed RDS files | Make reproducible rebuild path; document expensive steps | Medium | No | Log yes |
| Major | Missing branch support/model sensitivity | Phylogenetic conclusions lack robustness | Add bootstrap/support, model-selection table, alignment sensitivity | High | No | Citation yes; log yes |
| Major | GISAID/human HAI data in repo | Public sharing may be restricted | Create data availability and access plan | Medium | Yes | GISAID/PLOS policies; log yes |
| Major | Missing methods detail | Cannot reproduce HAI map/distance calculations | Add sequence IDs, map settings, titer handling, dimensions, optimizations | Medium | Yes | Citation yes |
| Moderate | Mixed `Results`/`results` paths | Render likely fails on case-sensitive systems | Standardize paths | Low | No | No |
| Moderate | Dependencies used but absent from lock | Reproduction may fail after `renv::restore()` | Update dependency record or remove unused deps | Low-medium | No | No |
| Minor | `products/manuscript.qmd:185` missing space | Signals draft immaturity | Copyedit after substantive fixes | Low | No | No |

## 6. Strengths Worth Preserving

1. The comparison question is real: sequence, temporal, p-epitope,
   cartographic, and phylogenetic distances are worth comparing if framed
   narrowly.
2. The generated outputs already show a plausible biological pattern, especially
   H1N1 temporal failure and H3N2 temporal/sequence alignment with tree
   distance.
3. The Quarto manuscript consumes generated figures/tables rather than doing all
   heavy computation inline. Preserve that separation, but put the generation
   under a real pipeline.
4. `renv.lock` exists and records R 4.3.3/Bioconductor 3.18. That is a useful
   starting point, but not sufficient.
5. The project already has a sensitivity-analysis outline in
   `archive/sensitivity-analysis.docx`; convert that from a list into actual
   planned sensitivity checks.

## 7. Missing Work Required for Publication

Conceptual: define whether this is a descriptive influenza paper, a
computational methods note, or a reproducibility appendix.

Literature: add influenza antigenic evolution, antigenic cartography
foundations, distance-matrix inference, phylogenetic reporting, SH/RF/NJ/model
citations.

Methods/analysis: redo correlations; add permutation/Mantel-style inference;
report unique pair counts; justify model choice; add branch support; add
alignment/model/cartography sensitivity analyses.

Results: add all NJ trees or supplement; add model-test table; add map
diagnostics; add strain/accession/provenance table.

Reproducibility/pipeline: create `_targets.R` or Makefile; fix path case;
document render command; ensure `renv.lock` covers every package; avoid derived
outputs under protected `data/`.

Data safety/privacy: verify GISAID redistribution, UniProt/GenBank accessions,
HAI data permissions, IRB/consent constraints, and coauthor approval.

Writing/structure: rewrite title, abstract, intro, methods, results,
discussion, limitations.

Submission materials: data availability statement, code availability statement,
ethics statement, author contribution statement, competing interests, funding,
reporting checklist.

## 8. Journal and Venue Recommendations

| Venue | Fit | Reframing | Risks | Competitiveness | Target |
|---|---|---|---|---|---|
| PLOS ONE | Broad rigorous descriptive study | "Validated descriptive comparison of distance metrics" | Data availability and methods rigor | Medium | Realistic after fixes |
| BMC Ecology and Evolution | Phylogenetics/evolution scope | Emphasize evolutionary-method comparison | Needs stronger phylogenetic reporting | Medium | Realistic/stretch |
| Infection, Genetics and Evolution | Pathogen evolution + bioinformatics | Infectious disease evolutionary genetics angle | Too narrow unless biological insight strengthened | Medium-high | Stretch |
| Virus Evolution | Virus evolution audience | Stronger viral evolution contribution | Current novelty too weak | High | Stretch |
| Archives of Virology | Broad virology incl. phylogeny/immunity | Applied virology comparison | May view as too computational/descriptive | Medium | Realistic/fallback |
| Influenza and Other Respiratory Viruses | Influenza-specific audience | Influenza antigenic characterization | Official scope needs verification | Medium | Realistic/fallback |
| BMC Bioinformatics | Computational biological data methods | Only if turned into reusable method/workflow | Current paper has no new algorithm/software | High | Stretch |
| Bioinformatics Advances | Computational biology venue | Only if packaged as reusable analysis/software note | Current work too application-specific | High | Stretch |

Avoid for now: PLOS Computational Biology, PLOS Pathogens, Journal of Virology,
Nature/Cell/Science-family venues, and methods journals requiring genuinely new
algorithms. JVI explicitly warns against phylogenetic analyses without clear
biological/mechanistic insight; that describes the current risk.

Preprint: not yet. Post only after statistical repair and data-permission
review. Registered report: no. Software paper: only if a reusable
package/workflow is built. Data paper: no unless permissions allow curated data
release.

## 9. Revision Roadmap

### Phase 1: Decide the Paper's Core Identity

Action: choose "descriptive computational influenza comparison" or
"dissertation reproducibility appendix." Owner: author/coauthors. Priority:
highest. Effort: low. Output: one-paragraph scope statement. Decision log: yes.

### Phase 2: Fix Fatal and Major Problems

Action: redo distance-matrix correlations using unique off-diagonal pairs and
permutation/Mantel-style inference; remove invalid CIs. Owner: analyst.
Priority: highest. Effort: medium. Dependency: identity decision. Output:
revised Figure 1/table. Decision log: yes.

### Phase 3: Rebuild the Argument

Action: rewrite title/abstract/intro/discussion around the narrowed claim.
Owner: author/editor/domain expert. Priority: high. Effort: medium. Output: new
manuscript skeleton. Decision log: yes for claim changes.

### Phase 4: Strengthen Evidence and Robustness

Action: add branch support, model-selection justification, all NJ trees,
alignment/model/cartography sensitivity checks. Owner: analyst/domain expert.
Priority: high. Effort: high. Output: supplemental methods/results. Decision
log: yes.

### Phase 5: Improve Reproducibility and Auditability

Action: add `_targets.R` or Makefile, fix paths, complete `renv`, document
render/rebuild commands, separate raw/derived data. Owner: analyst. Priority:
high. Effort: high. Output: clean reproducibility workflow. Decision log: yes.

### Phase 6: Prepare for Submission

Action: pick venue, verify scope/current policies, prepare data/code
availability, ethics, author contributions, cover letter, reporting checklist.
Owner: author/coauthors. Priority: medium after fixes. Effort: medium. Output:
submission package. Decision log: yes.

## 10. Recommended Next Actions

1. Decide whether this should be a real manuscript or only a reproducible
   dissertation appendix.
2. Confirm whether GISAID-derived sequences and HAI/cartography data can be
   shared or must be restricted.
3. Recompute Figure 1 using unique off-diagonal distances only.
4. Replace ordinary Pearson CI/p-value framing with permutation/Mantel-style
   distance-matrix inference.
5. Add a model-selection and phylogenetic-support section.
6. Convert the sensitivity-analysis outline into actual checks.
7. Build `_targets.R` or another authoritative pipeline.
8. Fix `Results`/`results` path inconsistencies.
9. Update `renv.lock` so all used packages are captured.
10. Rewrite the abstract and conclusion after the analysis is repaired, not
    before.

## Sources and Standards Cited in the Original Review

- [Legendre & Fortin on distance-matrix tests](https://pubmed.ncbi.nlm.nih.gov/21565094/)
- [MIAPA phylogenetic reporting](https://pmc.ncbi.nlm.nih.gov/articles/PMC3167193/)
- [PLOS data availability policy](https://journals.plos.org/plosmedicine/s/data-availability)
- [GISAID Terms of Use](https://gisaid.org/terms-of-use/)
- [targets manual](https://books.ropensci.org/targets/walkthrough.html)
- [renv documentation](https://rstudio.github.io/renv/)
- [Quarto freeze docs](https://quarto.org/docs/projects/code-execution.html)
- [PLOS ONE criteria](https://journals.plos.org/plosone/s/criteria-for-publication)
- [PLOS Computational Biology scope](https://journals.plos.org/ploscompbiol/s/journal-information)
- [Virus Evolution scope](https://academic.oup.com/VE/pages/About)
- [BMC Bioinformatics scope](https://bmcbioinformatics.biomedcentral.com/)
- [BMC Ecology and Evolution scope](https://link.springer.com/journal/12862)
- [Infection, Genetics and Evolution scope](https://www.sciencedirect.com/journal/infection-genetics-and-evolution)
- [JVI scope](https://journals.asm.org/journal/jvi/scope)
