# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code). Files shared with the DeCovarT package are copied from it; the others are specific to this repository.

Division of content: the DeCovarT vignettes hold the theory (models, assumptions, literature, reasons for each choice); this book holds the practice (Suppinger data, algorithms as implemented, settings, code, outputs) and links to vignette anchors instead of restating theory.

## Book chapters (`docs/`)

- **Scope.** The scope section of a chapter is one or two sentences: what the
  chapter builds, from which data, and the vignette it applies.
- **Code goes to the appendix.** The body does not walk through functions,
  scripts or output files. Each chapter ends with an appendix section
  (`## Code {#sec-<chapter>-code .appendix}`) holding one table: function or
  script, file, purpose, and the section it supports (`@sec-…`).
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

## Skills and rules: where they live, which one wins

Source of truth: `Config/.cursor/` in
[ubuntu-genesis](https://github.com/bastienchassagnol/ubuntu-genesis), cloned to
`~/repos/ubuntu-genesis`. `agent-sync.sh` copies it to the user level of both tools.

| Level | Rules | Skills | Written by |
| --- | --- | --- | --- |
| 1. Project | `.cursor/rules/*.mdc`, `.claude/rules/*.md` | `.cursor/skills`, `.claude/skills` (none yet) | this repository, committed |
| 2. User | `~/.cursor/rules`, `~/.claude/rules` | `~/.cursor/skills`, `~/.claude/skills` | `agent-sync.sh install-user`; do not hand-edit |
| 3. Machine-local | none | `firecrawl*`, `~/.agents/skills` (`npx skills`) | not vendored, not in git |
| 4. Host built-in | none | `~/.cursor/skills-cursor`, plugins, `~/.claude/skills/synced/` | the tool itself |

The tools disagree on a name clash, so a name is defined at one level only:

- Cursor merges rules with the project over the user level (Team, Project, User).
- Claude Code loads user and project rules together and does not define a winner.
  For skills of the same name, a personal skill (`~/.claude/skills`) beats a project skill.
- Claude Code reads only `paths:` in a rule's frontmatter; Cursor reads `globs:` and
  `alwaysApply:`. Bundle rules carry both. A rule without either loads in every session.
- Claude Code does not read `.cursor/` or `.agents/`.

Project copies in `.claude/rules/` and `.cursor/rules/` have the same names as bundle
rules and may differ from them. Read the project copy as the authority for this
repository. Do not add a project skill whose name the bundle already uses.

### Keeping the levels in step

- Edit in the bundle, commit there, `git pull`; the `post-merge` hook runs `install-user`.
- Edit under `~/.cursor` or `~/.claude`: the `afterFileEdit` (Cursor) and `PostToolUse`
  (Claude Code) hook, `~/.cursor/hooks/agent-sync-after-edit.sh`, runs `publish`
  and redeploys. Review `git diff` in the bundle, then commit.
- Log: `~/.cache/agent-sync/hook.log`. Check: `$AGENT_BUNDLE_ROOT/agent-sync/agent-sync.sh verify`,
  or open a shell with `CLAUDE_AGENT_SYNC=1` (set in `~/.zshrc`; output in `verify.log`).
- `AGENT_BUNDLE_ROOT` is exported in `~/.zshrc` and in the `env` block of
  `~/.claude/settings.json`. Hook shells do not read `~/.zshrc`.

### Symlinks

- Do not commit a symlink that leaves the repository (for example `.cursor/rules`
  pointing at a bundle outside it): it breaks on clone, in CI and on Windows, and a
  symlinked `skills/` puts every bundle skill at project level. Use copies.
  A relative link inside the repository is acceptable.
- `~/.bashrc` is a separate file that hands interactive shells to zsh. Do not link it
  to `~/.zshrc`: that file uses zsh-only syntax (`setopt`, `typeset -g`).
- Audit of `dtoo_project/` on 2026-11-18: the only tracked symlink is
  `slides/_extensions/pandoc-ext/diagram`, relative and inside the repository.

## Code graph (graphify)

The `graphify` rule applies only where `graphify-out/` exists, and the folder is not committed. If `graphify-out/graph.json` is missing, tell the user to run the `/graphify` skill from the repository root before exploring the code, and wait for the build. Do not substitute a long manual exploration. `.graphifyignore` already excludes generated and vendored paths (`renv/`, `docs/slides-deck/`, `docs/_extensions/`).
