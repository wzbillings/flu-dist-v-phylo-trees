# AGENTS.md

## Purpose and scope

This repository contains an in-progress research, data science, manuscript, technical report, analysis, or drafting project. Treat it as a research record, not just a codebase.

Apply these instructions unless a more specific `AGENTS.md` in a subdirectory overrides them. Follow the repository's existing conventions over generic examples in this file.

Priorities, in order:

1. Protect data, privacy, and secrets.
2. Preserve reproducibility and provenance.
3. Keep claims, analyses, outputs, and citations traceable.
4. Make small, reviewable changes.
5. Prefer correctness and auditability over speed or cleverness.

## Use available research skills

When the task involves literature review, manuscript planning, manuscript review, journal targeting, reproducible research planning, citation integrity, or academic project triage, use `$academic-research-suite` as appropriate.

Do not invoke the skill for ordinary coding, formatting, file management, or mechanical cleanup unless the task has a substantive research component.

## Start-of-task checklist

Before editing, inspect the repository structure and existing conventions. Look for:

- `README.md`, `AGENTS.md`, `CONTRIBUTING.md`, project notes, protocols, preregistrations, or decision logs.
- `.gitignore` and other ignore files.
- Workflow files such as `_targets.R`, `_quarto.yml`, `Makefile`, `Snakefile`, GitHub Actions, or notebooks.
- Dependency files such as `renv.lock`, `pyproject.toml`, `requirements.txt`, `environment.yml`, `DESCRIPTION`, or `package.json`.
- Common project directories such as `R/`, `src/`, `analysis/`, `manuscript/`, `reports/`, `docs/`, `data/`, `outputs/`, or `figures/`.

For complex, ambiguous, or high-risk tasks, make a short plan before editing. For straightforward low-risk tasks, proceed directly.

Do not reorganize the repository unless explicitly asked or unless the current structure blocks reproducibility.

## Data, privacy, and secrets

Treat raw, source, restricted, and confidential data as protected research artifacts.

Do not edit, overwrite, reformat, relocate, or delete protected data unless explicitly instructed. This applies especially to:

- `data/`, `data/raw/`, `data/input/`
- `private-data/`, `restricted-data/`
- `inst/extdata/`
- any directory labeled raw, source, confidential, restricted, or private

If data transformation is needed, implement it as code and write derived outputs to an appropriate derived-data, cache, or analysis-output location. Raw data should remain immutable.

Before adding or modifying files, check for:

- personally identifiable information or protected health information
- confidential participant, collaborator, business, or reviewer information
- embargoed or unpublished results
- credentials, tokens, private URLs, keys, or secrets
- identifying metadata in exported files

Never move sensitive information into public documentation, examples, tests, rendered outputs, logs, or generated artifacts. Use synthetic, simulated, de-identified, or toy data for examples.

Never hard-code secrets in code, documentation, YAML, shell scripts, GitHub Actions, notebooks, or Quarto documents. Use environment variables, local secret files, credential managers, or CI/CD secret stores. If code needs a secret, document only the variable name.

Do not read from or write to `.Renviron`, `.Rprofile`, `.secrets`, SSH keys, service-account keys, or similar local configuration files unless explicitly instructed.

## Git and generated files

Respect `.gitignore`. Do not add generated outputs, caches, large intermediates, private data, temporary files, local IDE files, or machine-specific artifacts unless the repository clearly tracks them.

Files that usually should not be tracked include:

- `.Rhistory`, `.RData`, `.Renviron`, `.Rproj.user/`
- `.quarto/`, `_targets/`, `_targets-r/`
- cache directories and temporary files
- logs containing data or private paths
- private data and credentials

Do not commit files between 50 MB and 95 MB without explicit user confirmation. If a requested change would add or modify a file in this size range, pause before committing it, identify the file and its size, explain why it may be inappropriate for Git, and ask whether to proceed or use an alternative such as Git LFS, external storage, data regeneration, or `.gitignore`.

Do not commit files larger than 95 MB. If a requested change would require committing such a file, treat it as a blocking issue: stop, explain which file is too large, and recommend an alternative such as Git LFS, external storage, data regeneration, or adding the file to `.gitignore`.

Do not create commits, push to remotes, rewrite history, delete branches, or change repository visibility unless explicitly asked.

## Reproducible workflow rules

Prefer the project's established workflow. If the project uses a pipeline tool, treat that pipeline as the source of truth for derived data, models, tables, figures, diagnostics, and manuscript-ready outputs.

For R projects using `targets`:

- Put reusable functions in `R/`.
- Put target definitions in `_targets.R` or files sourced by it.
- Avoid ad hoc scripts that duplicate pipeline logic.
- Avoid hidden global state and interactive-session dependencies.
- Use descriptive target names.
- Use file targets for generated files.
- Control random number generation with `tar_option_set(seed = ...)` or explicit documented seeds.

Quarto, R Markdown, notebooks, and manuscript files should render or run from a clean session. They should usually consume pipeline outputs rather than recreate major analysis steps inline.

Do not manually edit generated tables, figures, model outputs, or derived datasets. Change the code that generates them.

## Dependencies and environment

Do not casually alter dependency management.

If the project uses `renv`, `pak`, `requirements.txt`, `uv`, `conda`, `poetry`, `DESCRIPTION`, or another environment system, use the project's normal workflow.

When adding a dependency:

- Prefer established, maintained packages.
- Avoid adding packages for small tasks already handled by existing dependencies or base language features.
- Document why the dependency is needed when the reason is not obvious.
- Update lock files only when appropriate.

## Research integrity and manuscript standards

Substantive scientific, statistical, methodological, or domain-specific claims must be supported by one of:

- project results
- cited academic literature
- cited official documentation
- explicitly labeled expert judgment
- clearly stated assumptions

Do not invent citations, DOIs, PMIDs, author lists, journal details, URLs, or bibliography entries. Use clear placeholders when citation information is missing, for example:

```text
[ADD CITATION: methodological reference for clustered standard errors]
```

Do not strengthen claims beyond the evidence. Do not convert exploratory findings into confirmatory claims. Do not remove uncertainty to make a manuscript sound stronger.

Clearly distinguish:

- prespecified analyses
- primary outcomes
- secondary outcomes
- exploratory analyses
- post hoc analyses
- sensitivity analyses
- robustness checks

Use stronger causal language only when the design and assumptions justify it.

When editing manuscripts, improve precision, transparency, argument structure, and evidentiary support. Avoid filler, exaggerated novelty claims, vague causal language, and unsupported claims of practical importance.

Where relevant, consider reporting guidelines such as CONSORT, STROBE, PRISMA, TRIPOD, RECORD, CHEERS, ARRIVE, CARE, SQUIRE, SAMPL, or discipline-specific standards. Do not claim compliance unless the checklist has actually been completed.

## Decision logging

Important research decisions should be documented in an existing decision log or in `docs/analysis-decision-log.md` if no better location exists.

Log decisions involving:

- inclusion or exclusion criteria
- outcome definitions
- covariate selection
- transformations
- model family or prior distributions
- missing data or outlier handling
- sensitivity analyses
- deviations from protocols or preregistrations
- changes after reviewer comments
- journal targeting decisions

Suggested format:

```markdown
## YYYY-MM-DD — Decision Title

**Decision:** What was decided.

**Rationale:** Why this decision was made.

**Evidence / citation:** Relevant manuscript section, code location, protocol, reviewer comment, or literature citation.

**Alternatives considered:** Other plausible options and why they were not chosen.

**Impact:** How this affects interpretation, reproducibility, or publication strategy.
```

## Coding and analysis standards

Use clear, boring, maintainable code.

Prefer explicit function arguments, descriptive names, small functions, readable pipelines, clear error messages, comments explaining why, and project-root-aware paths.

Avoid hidden global state, unnecessary metaprogramming, unexplained magic constants, copy-pasted analysis blocks, monolithic scripts, `setwd()`, absolute local paths, and model selection based on favorable results.

Data cleaning code should make clear what changed, why it changed, and how many records were affected when relevant. Do not silently drop rows, recode values, or impute missing data.

Analysis code should make the inferential target clear, including the estimand, target population, unit of analysis, model, assumptions, diagnostics, and sensitivity analyses when relevant.

Tables and figures should be generated reproducibly and should include clear captions, units, populations, transformations, exclusions, and uncertainty intervals where appropriate.

## Change management and final response

Prefer small, coherent edits. Avoid mixing unrelated changes. Do not perform large rewrites, major restructuring, or sweeping formatting changes unless explicitly requested.

When completing a task, summarize:

- what changed
- why it changed
- files affected
- checks run or not run
- whether outputs need regeneration
- whether citations or decision logs need follow-up

## Checks before finishing

Before finalizing changes, consider whether the task affected:

- data privacy or secret exposure
- reproducibility or pipeline integrity
- manuscript claims or citations
- generated outputs or decision logs
- dependency files or README instructions
- Quarto rendering, tests, or public-release readiness

For code changes, run the smallest relevant checks supported by the project. Examples:

```r
targets::tar_manifest()
targets::tar_make()
quarto::quarto_render("manuscript/manuscript.qmd")
```

Use the project's actual commands when they differ. Do not run expensive, destructive, publication, export, or submission commands without confirming that they are appropriate.

## Stop and ask before high-risk changes

Ask for human confirmation before:

- modifying raw or restricted data
- deleting files
- changing primary outcomes, estimands, inclusion criteria, exclusion criteria, or main models
- changing interpretation of the main result
- adding or removing major claims
- adding substantial dependencies
- changing licenses
- changing public/private release assumptions
- uploading, exporting, submitting, or publishing outputs
- making irreversible repository changes

For low-risk cleanup, documentation improvements, small reproducibility fixes, and obvious typo corrections, proceed without asking unless project instructions say otherwise.
