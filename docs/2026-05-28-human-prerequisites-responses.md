# Human Prerequisites Responses

Date started: 2026-05-28
Responder: ZB

## 1. Project Identity

### Q: Should this project be treated as a publishable journal manuscript, dissertation appendix, reproducible archive, or combination?
Response: we will expand and revise it to a publishable journal manuscript.
Status: Decided
Decision-log entry needed: yes
Notes:

### Q: If it should become a journal manuscript, what is the intended article type: descriptive empirical paper, computational methods note, virology short report, reproducibility report, or another type?
Response: original empirical research article
Status: Decided
Log: no
Notes: this is an original research article that provides expansion on the potential mechanisms of the disagreement in antigenic distance metrics that we found in our previous work.

### Q: What is the one-sentence claim the paper is allowed to make?

Response: Phylogenetic distance is a poor proxy for antigenic change in influenza.
Status: Decided
Log: no
Notes:

### Q: What claims should the paper explicitly avoid making, even if they sound
  attractive?

Response: That our claims are representative of all flu strains, or populations. cartographic mismatch means hemagglutinin sequence is not the only factor in immune responses

### Q: Is the current goal to salvage the existing draft, rebuild it from
  scratch, or document why it should not be pursued as a formal publication?

Response: we will rebuild the manuscript from scratch using this draft as a template.
The manuscript will need a lot of expansion and revision, but we just want to
use the ideas and methods skeleton presented here.

### Q: Who is the primary intended audience: influenza virologists, antigenic
  cartography researchers, evolutionary biologists, computational biologists,
  vaccine modelers, or dissertation committee readers?

Response: the primary intended audience will be influenza vaccine scientists.
It should be readable to influenza scientists without technical backgrounds and
avoid overly methods-heavy language.

### Q: Should the manuscript focus narrowly on the strain panel in this
  repository, or should it be expanded with additional cohorts, strains, or
  datasets?

Response: this is actually the largest heterologous influenza strain panel in a
human cohort study so far. So we will focus on this specific strain panel, but
note that it was chosen by experts and analyzing it in this way still provides
valuable unknown information about antigenic cartography. Studies that use
many different influenza strains can typically only perform phylogenetic analyses,
not antigenic cartography.

### Q: Should H1N1 and H3N2 remain in one paper, or should the analysis be
  reframed around the subtype-specific differences?

Response: H1N1 and H3N2 should be in the same paper. The discussion about the similarities
and differences in results between the subtype is a crucial new finding.

### Q: Should neuraminidase sequence or inhibition data be considered required
  for publication, or only discussed as future work?

Response: This data does not currently exist, so it is for future work only.

### Q: Which sequence records came from GISAID, UniProt, GenBank, collaborators,
  or other sources?

Response: We have this information compiled for another paper and I will provide you
with the accession table and citations in a future step.

### Q: Are any GISAID-derived sequences, metadata, alignments, distance matrices,
  or derived files currently tracked in this repository?

Response: sequences are stored in the dataset in this repo, which is permissable
by their sharing standards. No other information and no files are from GISAID.

### Q: Are the GISAID-derived materials allowed to remain in the repository if
  the repository is public or pushed to GitHub?

Response: yes.

### Q: Are the HAI data, Racmacs `.ace` files, cartographic coordinates, or
  derived cartographic distances approved for sharing with collaborators,
  reviewers, or the public?

Response: yes.

### Q: Are any human participant data, potentially identifying metadata, site
  information, dates, or restricted study details present in the repository?

Response: no.
Notes: our data has been designated exempt nonhuman subjects data by our institutional IRB.
The methods section should include a "Study ethics" subsection that explains this
and states that the original study was IRB-approved, I will provide the IRB details
in a future step.

### Q: What exact data availability statement is legally and ethically permitted?

Response: All data used for the analysis are provided in our repository. Our institution's
IRB has designated the data as non-human subjects data, and all sequences are obtained
from accessible sources.

### Q: Should any data files be removed from Git tracking, moved to controlled
  storage, or replaced with synthetic/toy data before public release?

Response: no. When we are ready to release we will curate a new public-facing
repo that does not contain internal or unneeded files, and does not allow
extraction of these things from commit history.

### Q: Who must approve any manuscript submission?

Response: right now, only me. I will pass the work to coauthors at the appropriate stages.

### Q: Who must approve any public repository, data release, preprint, or
  supplemental material?

Response: see above.

### Q: Who qualifies as an author, contributor, data provider, or acknowledged
  collaborator?

Response: I will provide a table of coauthors in a future step.

### Q: Are there collaborators or data generators whose approval is required
  before analyzing, publishing, or sharing derived results?

Response: I will handle this manually at the correct time.

### Q: Are there unpublished, embargoed, dissertation, collaborator, or reviewer
  constraints that limit what can be written or released?

Response: No.

### Q: Is the project allowed to claim that temporal distance is unsuitable in
  general, or only that it performed poorly in this selected H1N1 panel?

Response: Yes.
Notes: The argument for H3N2 is based on the fact that while the match is better
it still isn't public, but it's really bad for H1N1 and we should be comparing
subtypes on the same metric, so there's no good reason to use temporal for
any subtype.

### Q: Is the project allowed to claim that genetic and antigenic evolution do
  not agree, or only that the selected HA phylogenetic distances and HAI-based
  cartographic distances diverged in this dataset?

Response: Of course the latter is true, but this is a limitation of our study
and any possible study of this type. We believe that our results generalize
to all H1N1 and H3N2 strains in general and this is the conclusion we want to make,
but the generalizability and panel sampling is an important limitation to mention.

### Q: Is the project allowed to make any claim about individual immune response
  determinants, or is that outside the evidence available here?

Response: we can only make general claims. We know that cartography is an imperfect
but relatively good measure of how differently strains appear to the immune system.
So if we can see differents in cartography and evolution, we know that evolutionary
patterns and sequences are not optimized for immune evolution. We just can't
say anything about the mechanisms of the trends we see.

### Q: Which analyses were planned before seeing the results, and which were
  exploratory or post hoc?

Response: All of the analyses in the current draft were planned from the start, although
we did not formally preregister this manuscript.

### Q: What limitations must be stated plainly in the abstract, discussion, or
  data availability section?


Response: the limitations that are currently in the document and that I've noted
in previous answers. Additionally, we want to note that the panels are not representative
of all influenza strains and more will be better and easier to do when multiplex
assays are available in the near future. We also know that while the panel was chosen
by experts based on standard lab strains they do not evenly cover the antigenic distance
space for all metrics.

### Q: Should the distance-matrix comparison be rebuilt using unique off-diagonal
  strain pairs only?

Response: yes.

### Q: Should the main inference use Mantel/permutation-style tests, descriptive
  correlations only, or another method?

Response: we want to include descriptive correlations because they are easy
to interpret for a nontechnical audience and linear relationships have a clear
meaning here (measurement calibration). However, we should use the Mantel version of the statistic.

### Q: Should diagonal and duplicate pairwise distances be excluded from all
  reported correlations and figures?

Response: yes.

### Q: Should Pearson, Spearman, Mantel correlation, matrix regression, or another
  metric be the primary comparison?

The primary metric should be Mantel correlations for the matrix comparison and change in loglikelihood
(with LRT if it is constructible) should be the main tree comparison metric.

### Q: Should p-values and confidence intervals be reported at all for the
  distance comparisons?

Response: they at least need to be in a table in the supplement if they are straightforward
to calculate.

### Q: What sensitivity analyses are required before the manuscript is worth
  resubmitting or submitting?

Response: return to this after initial result draft is done.

### Q: Should ML tree fitting be rerun from source data, or should existing
  precomputed ML tree objects be treated as fixed historical outputs?

Response: We should do all analysis steps as part of our pipeline, from sequence cleaning
through to final graphs. This will mean cartography code needs to move here too in a future step.

### Q: What branch support, model-selection, alignment, cartography, and missing
  data checks are required?

Response: not sure what this means enough to provide an answer.

### Q: Should the project be converted to a `targets` pipeline?

Response: yes

### Q: If not `targets`, what should be the authoritative workflow: Makefile,
  Quarto project render, numbered scripts, or another system?

Response: N/A

### Q: Which outputs should be reproducibly regenerated, and which should remain
  archived historical artifacts?

Response: all artifacts are preserved in the git history if there are issues, so
treat all as needing to be reproducible.

### Q: Which files in `data/` are raw/protected source data and must remain
  immutable?

Response: `full-sequences.xlsx` and `UGAFluVac-virus-names.csv` are the raw data
files, you should never change them. I will only chnage them if I identify
issues in the source data.

### Q: Which files are derived outputs and should be moved or regenerated outside
  protected source-data locations?

Response: all other files in `data/`, all files in `results/`, and any `.docx` or
`pdf` files in `products/`.

### Q: Should the existing rendered Word documents be kept as historical outputs
  or regenerated from Quarto after the workflow is repaired?

Response: regenerated.

### Q: Is a journal submission actually desired?

Response: yes

### Q: If yes, should the target be a broad-scope journal, influenza/virology
  journal, evolutionary biology journal, computational biology journal, or
  methods venue?

Response: influenza or virology journal. Our top target is probably PNAS or Proc Soc B
if you think either of those is achievable, followed by something on the level of
mBio. PLOS CB and Journal of Virology have liked our papers in the past but they,
along with Influenza and Other Respiratory Viruses, are our backups.

### Q: Is a preprint desired, and if so, who must approve it first?

Response: Yes, once we get an initial round of coauthor approval. I will
approve and upload to MedRxiv myself.

### Q: Are APCs, open-access requirements, data-sharing policies, or institutional
  repository requirements constraints on venue choice?

Response: Strongly prefer journals with open-access requirements and journals
run by idependent societies.

### Q: Should the work be positioned as a modest descriptive result rather than a
  major computational biology contribution?

Response: this is not groundbreaking, but it is very suggestive of what phenomena
underlie our previous findings. Many papers use bad distance metrics and we want
to influence them to use the best distance metrics in their studies.

