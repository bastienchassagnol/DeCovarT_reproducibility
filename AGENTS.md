# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code): copies of the ubuntu-genesis bundle rules, refreshed by `agent-sync install-project`, plus the project-only `bibliography` rule (section below).

Division of content: the DeCovarT vignettes hold the theory (models, assumptions, literature, reasons for each choice); this book holds the practice (Suppinger data, algorithms as implemented, settings, code, outputs) and links to vignette anchors instead of restating theory.

## Book chapters (`docs/`)

- **Scope.** The scope section of a chapter is one or two sentences: what the
  chapter builds, from which data, and the vignette it applies.
- **Code goes to the appendix.** The body does not walk through functions,
  scripts or output files. Each chapter ends with an appendix section
  (`## Code {#sec-<chapter>-code .appendix}`) holding one table: function or
  script, file, purpose, and the section it supports (`@sec-…`).
- **Folded code.** Code chunks are folded by default (`code-fold: true`,
  `code-summary: "Show the code"`) and every page has the Code menu
  (`code-tools: true`), both set once in `docs/_quarto.yml`; a chunk that
  must stay open sets `#| code-fold: false`.
- **Verbatim names.** Function, method, argument, package, file, column and
  layer names are always in backticks: `presence_gate()`, `min_cells`,
  `glmGamPoi`, `scored_slots.rds`.
- **Options not taken.** Report the filters and thresholds that were
  applied. Filters considered but not applied go to a table footnote
  (`tinytable::tt(notes = ...)`) or to the vignette, not to the body.
- **Figures.** Every saved figure has, in the same folder and with the same
  base name, an `.rds` holding the minimal data to redraw it: a tibble, or a
  named list of tibbles (data, annotations, settings). `save_figure(data = )`
  in `R/02_00_utils.R` writes it. Bars, lollipops and box plots are ordered by
  decreasing value. A figure that needs the full page width uses
  [the full-width trick](https://www.productive-r-workflow.com/quarto-tricks#full-width-graph):
  `#| column: screen` and `#| out-width: 100%` on the chunk, or
  `::: {.column-screen}` around a static image (`.column-page` keeps the
  margins).
- **Dictionaries.** Every variable read from the Suppinger data or from
  `data/dictionaries/` is documented (meaning, units, allowed values, missing
  coding, source) in `data/dictionaries/suppinger_metadata_dictionary.csv`.
  A new column or dictionary file adds its rows there in the same change.

## Rules and skills: where they live

Source of truth: `Config/.cursor/` in
[ubuntu-genesis](https://github.com/bastienchassagnol/ubuntu-genesis).

```yaml
# read by agent-sync install-project
packs:
  - quarto-r-analysis
project_rules: copy
```

| Level | Rules | Skills | Written by |
| --- | --- | --- | --- |
| 1. Project | `.cursor/rules/*.mdc`, `.claude/rules/*.md` | none | `agent-sync install-project` (bundle copies) + `bibliography.mdc` (this repository) |
| 2. User | `~/.cursor/rules`, `~/.claude/rules` | `~/.cursor/skills`, `~/.claude/skills` | `agent-sync install-user`; do not hand-edit |
| 3. Machine-local | none | `firecrawl*`, `~/.agents/skills` | not vendored |
| 4. Host built-in | none | `~/.cursor/skills-cursor`, plugins, `~/.claude/skills/synced/` | the tool itself |

- Project and user copies of a bundle rule are the same file. Claude Code loads
  both levels together (no winner), Cursor prefers the project copy; keeping
  them identical avoids any conflict. Refresh after a bundle change:
  `"$AGENT_BUNDLE_ROOT/agent-sync/agent-sync.sh" install-project .`
  (copies the pack into both folders, prunes rules removed from the bundle,
  keeps `bibliography.mdc` and mirrors it to `.claude/rules/bibliography.md`).
- `bibliography.mdc` is project-only: the curation scripts live in the sibling
  `../DeCovarT/scripts/curate_bibliography/`, and `check_bibliography.py --clean`
  never runs on the book bibliography.
- Claude Code reads `paths:` in a rule's frontmatter; Cursor reads `globs:` and
  `alwaysApply:`. Bundle rules carry both.
- Symbolic verification lives in `verification/` (`run_all.py`).
- Do not commit a symlink that leaves the repository; the only tracked symlink
  is `slides/_extensions/pandoc-ext/diagram`, relative and inside the repository.

## Code graph (CodeGraph first, Graphify second)

CodeGraph is the default for code. Graphify is secondary and applies to mixed material. Routing follows the `code-intelligence` rule; this section adds the repository paths.

- **Code questions** (`R/`, `scripts/`, `verification/`, Python helpers: callers, callees, impact, how a call path reaches a function): where `.codegraph/` exists at the repository root, call `codegraph_explore` (MCP) or `codegraph explore "<question>"` before Grep, Glob or Read. Treat the returned source as already read. The index auto-syncs while the MCP server runs. `.codegraph/` is not committed; build it once per checkout with `codegraph init`. `codegraph.json` excludes `renv/`, `docs/slides-deck/`, `docs/_extensions/`, `docs/_freeze/`, `docs/_book/` and `slides/_extensions/`.
- **Cross-document questions** (links across book chapters in `docs/`, the DeCovarT vignettes, papers and design notes): use `graphify query "<q>"` only where `graphify-out/graph.json` exists. The folder is not committed and is rebuilt at milestones (annotated tag or version bump, via `Config/.cursor/code-knowledge/scripts/graphify-milestone.sh` in the ubuntu-genesis bundle), not on every commit. Do not run `graphify hook install` for code sync: CodeGraph covers it. `.graphifyignore` already excludes generated and vendored paths (`renv/`, `docs/slides-deck/`, `docs/_extensions/`).
- **Neither index exists**: for a code question, tell the user to run `codegraph init` and wait; for a cross-document question, tell them to run the `/graphify` skill from the repository root. Do not substitute a long manual exploration.
