# sctransform_panel

This folder holds the first gene panel of the feature-selection chapter: `SCTransform` v2 refitted per slot, genes ranked by a type-equal between/within ratio of Pearson residuals, and a 500-gene export. It was moved here from `scripts/legacy/` on 2026-11-18. The panel was dropped because the live selection fits a negative-binomial GLM with `glmGamPoi::glm_gp()` and ranks genes by the cell-weighted deviance R² (`scripts/02_02_cross_cell_consistency.R`); none of the current scripts sources these files.

- `02_90_stats_prepare_panel.R`: loads GSE229513, refits `SCTransform` per slot, writes the ranking tables to `output/naive_marker_selection/tables/` and the slim Seurat object to `data/intermediate/suppinger_sct_allowlist_seurat.rds`. Sources `02_90_stats_rank_sctransform.R`.
- `02_90_stats_rank_sctransform.R`: the between/within split of Pearson residuals and the rank helpers.
- `02_90_vis_legacy_panel.R`: figures of the legacy panel, written under `output/naive_marker_selection/`.

Outputs still read: the tables in `output/naive_marker_selection/tables/` are read by `docs/01-data-suppinger.qmd`, `docs/02-feature-selection.qmd` (legacy SCTransform sections, `tbl-fs-sct-settings`) and `docs/06-sensitivity-misspecification.qmd` (rank A, `ranks_global_equal_type.csv`). Do not delete them.

Run from the repository root, as in the header of each script; paths are relative to `getwd()`.
