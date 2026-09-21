# DeCovarT_reproducibility

DeCovarT's companion analysis repository. All datasets are freely
available via GEO. If you use the method in your own work, please cite
the package and this repository.

## How to render the Quarto book

From the repository root, after `renv::restore()`:

```bash
quarto render docs
# first time only, once gh-pages is allowed:
# quarto publish gh-pages docs
```

The book format is **HTML only** for now (`pdf` is commented in
`docs/_quarto.yml`). CI publishes HTML to the `gh-pages` branch.

R chunks in `docs/` stay `eval: false` until artefacts exist. Mermaid
and static Markdown still execute; [`freeze: auto`](https://quarto.org/docs/projects/code-execution.html#freeze)
caches those results in `docs/_freeze`.
