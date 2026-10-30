# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code). Files shared with the DeCovarT package are copied from it; the others are specific to this repository.

Division of content: the DeCovarT vignettes hold the theory (models, assumptions, literature, reasons for each choice); this book holds the practice (Suppinger data, algorithms as implemented, settings, code, outputs) and links to vignette anchors instead of restating theory.

Rules already state repository scope, R style, scientific writing, and statistical modelling. This file only records how they combine here:

- `decovart-analysis.mdc` overrides any generic R **package** workflow. This tree is an analysis project.
- `quarto-book.mdc`, `bibliography.mdc`, and `nature-methods-chapters.mdc` attach when editing `docs/` or `bibliography/`.
- `literature-intake.mdc` applies when converting office files, slides, or PDFs.
