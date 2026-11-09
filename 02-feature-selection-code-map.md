# Feature-selection code map

Read this before editing the Suppinger gene and cell-type selection.
The live pipeline has been written and has **not been run**. Do not
launch the scripts in `scripts/02_*.R` unless the task is to produce
outputs. The current 500-gene panel on disk still comes from the
archived `SCTransform` rank.

Repository: `DeCovarT_reproducibility` (this directory). Theory and
proofs stay in the DeCovarT vignette
`vignettes/03-1-feature-selection.qmd`. This book chapter records what
is implemented:
`docs/02-feature-selection.qmd`.

## Naming

`02` is the feature-selection chapter. The second number is the tier.

| Second number | Tier |
|---|---|
| `00` | Shared helpers |
| `01` | Quality control and prior knowledge |
| `02` | Cross-cell consistency |
| `03` | Cross-platform consistency |
| `05` | Classifier perspective (not a tier) |

A third token is used only when that tier has both a computation file
and a plot file. Computation is `stats_`. Plots are `vis_`. Tiers (ii)
and (iii) are computation only, so they have no role token.

Covariance files stay unnumbered (`R/cv_covariance_selection.R`,
`R/ebic_regularisation_path.R`, `R/plot_grn_regularisation.R`). They
belong to chapter 03. Tier (iv) does not select genes and has no file.

## Two models

Live fits call `glmGamPoi::glm_gp()` with an explicit design:

- M0: `log m_gc = beta_0g + log s_c`
- M1: `log m_gc = beta_0g + beta_g,j(c) + log s_c`
- `offset = log(s_c)`, `size_factors = 1`
- `s_c` is the full RNA library size divided by its median in the slot
- overdispersion is estimated gene by gene
- residuals are not clipped

`Seurat::SCTransform()` is not on this path. It is archived under
`scripts/legacy/` because it produced the panel already on disk. New
drivers must not source that folder.

`DeCovarT/R/01_02_generate_synthetic_networks.R` is package API for
synthetic networks. It is not a vignette toy and it stays in the package.

## Source order

1. `R/02_00_utils.R` (`%||%`, `assay_matrix()`, `seurat_metadata_chr()`, `set_overlap()`)
2. `R/02_01_stats_quality_control.R` (needs step 1; defines `presence_gate()` and `map_marker_symbols()`)
3. `R/02_01_vis_projections.R` (plots only)
4. `R/02_02_cross_cell.R` (needs steps 1 and 2)
5. Tier (iii), after step 1: `R/02_03_influence.R`, `R/02_03_bulk_consistency.R`, `R/02_03_plugin_composition.R`, then `R/02_03_tail_bounds.R` (the tail file calls the plug-in)
6. `R/02_05_classifiers.R` (needs step 1)

`glmGamPoi` is not in `renv.lock`. `org.Mm.eg.db` is used offline for
marker symbols (no `biomaRt`). `mixOmics` is used for PLS-DA and is
checked with `requireNamespace()`. `quantreg` is the quantile band.
None of these is installed by the scripts.

## Tier (i) — `R/02_01_stats_quality_control.R`

| Function | Role |
|---|---|
| `keep_bulk_matched_cells()` | Time in `48h`/`72h`/`96h`, non-missing label, library size > 0 |
| `scored_types()` | Drop artefact labels and types with fewer than `min_cells` cells |
| `read_bulk_htseq()` | `data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz`; row names are Ensembl ids; column `gene_name` |
| `standard_wnt_samples()` | Three columns whose names contain the time token and `_CTL_` (standard Wnt). `_Early_` is the early-Wnt arm and is not the matched pair |
| `shared_bulk_single_cell_genes()` | Intersect RNA symbols with bulk `gene_name`. `keep` is TRUE only when the symbol has exactly one Ensembl id after stripping a version suffix (`.16`). `row_id` is the first bulk row of that id |
| `bulk_symbol_counts()` | Bulk count matrix on those symbols |
| `is_mitochondrial_gene()` | `^mt-`, case insensitive |
| `drop_mitochondrial_genes()` | Removes those symbols. Ribosomal and sex-chromosome symbols stay |
| `cell_mitochondrial_fraction()` | Share of counts on `mt-` symbols. Stored, not thresholded |
| `map_marker_symbols()` | Offline symbol map. Needs `org.Mm.eg.db` |
| `presence_gate()` | Gene detected and varying in every scored type |
| `prepare_scored_counts()` | Runs the cell filter, records the mitochondrial fraction, then keeps shared non-mitochondrial genes and scored types |

`presence_gate(min_frac)` defaults to `0.1`. Every call in the scripts
and in the chapter passes `min_frac = 0.05`. Report the value that was
used; do not silently switch to the default.

Driver: `scripts/02_01_stats_quality_control.R`.

Fixed settings: `min_cells = 20`, `min_frac = 0.05`, times
`48h`/`72h`/`96h`, metadata column `celltypeannotation`, allowed types
are `include_in_ranking == TRUE` in
`data/dictionaries/suppinger_celltype_crosswalk.csv`.

Outputs (not produced): `output/feature_selection/tier_01/`, including
`scored_slots.rds`, `shared_genes.csv`, `mito_fraction.csv`,
`dropped_mt_genes.csv`, `cell_counts.csv`, `presence_<time>.csv`,
`marker_symbol_map.csv`.

The script stops if either raw file is missing:

- `data/raw/GSE229513_gastruloidsobject.rds` (DVC)
- `data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz`

## Tier (i) plots — `R/02_01_vis_projections.R`

`plot_type_pca_plsda()` draws two scatters of cells, coloured by scored
cell type:

- PCA: `log1p(counts / size_factor)`, genes centred and scaled, `prcomp`
- PLS-DA: `mixOmics::plsda()` on the same log matrix, `ncomp = 2`, response = cell type

`size_factor` in the driver is the full RNA library divided by its
median in the slot (the same scale as the GLM offset). The function
does not save files. No UMAP, no Venn, no barcode distribution.

Driver: `scripts/02_01_vis_projections.R`. Figures would be written to
`output/feature_selection/tier_01/figures/` (`pca_cell_type_<time>.pdf`,
`plsda_cell_type_<time>.pdf`). Not run. Colours come from
`data/dictionaries/suppinger_celltype_colours.csv`.

## Tier (ii) — `R/02_02_cross_cell.R`

| Function | Role |
|---|---|
| `compare_type_effect_tests()` | `glm_gp()` for M0 and M1; four tests: `lr`, `anova_vst`, `ql_f`, `welch_perm` |
| `type_test_kendall()` | Kendall tau-b of the log p-value ranks |
| `type_test_calibration()` | Share of null p-values below 0.01 and 0.05 |
| `rank_genes_omnibus()` | BH gate on one test (default `ql_f`), then rank by deviance R^2 |

`glmGamPoi` reports overdispersion `1 / psi`. The likelihood-ratio
refits M0 at the M1 overdispersion. The quasi-likelihood test is
`glmGamPoi::test_de(fit, reduced_design = ~1)`. The permutation is
Freedman-Lane on the Pearson residuals of M0, Welch's G, one shared
permutation of labels per draw, `B = 2000`.

Driver: `scripts/02_02_cross_cell_consistency.R`.

| Argument | Value in the driver |
|---|---|
| `n_perm` | 2000 |
| `fdr` | 0.05 |
| `min_frac` | 0.05 |
| `n_select` | 500 |
| `gate_test` | `ql_f` |
| `seed` (fit) | 20261101 |
| seed for null labels | 20261102 |
| seed for the null fit | 20261103 |

Reads `output/feature_selection/tier_01/scored_slots.rds`. Writes
`output/feature_selection/tier_02/` (`tests_<time>.csv` has all four
tests; `ranked_<time>.csv` is the gate-and-rank table;
`cross_cell.rds` holds both and the selected-gene vector `panel`).
Not run.

## Tier (iii)

| File | Functions |
|---|---|
| `R/02_03_influence.R` | `cook_cells()` (M1 via `glm_gp()`), `cook_genes()` |
| `R/02_03_bulk_consistency.R` | `bulk_consistency_filter()` (range and MAD; heuristic) |
| `R/02_03_plugin_composition.R` | `nnls_lawson_hanson()`, `plugin_composition_demixsc()` |
| `R/02_03_tail_bounds.R` | one function per tail statement, plus the dispatcher |

Tail statements. `tail_threshold()` checks the unit, then calls the
matching function. `quantile_band_flags()` is separate (model-free).

| Function | Method | Statement | Unit |
|---|---|---|---|
| `tail_t1_clt()` | `clt` | T1, Lindeberg-Feller | gene and independent-gene module |
| `tail_t1_log()` | `clt_log` | T1-log | gene only |
| `tail_t2_clt_dep()` | `clt_dep` | T2 | module only; needs `adjacency` and `cell_cor` |
| `tail_t3_poisson()` | `poisson` | T3 | gene and module |
| `tail_t4_nb()` | `nb` | T4 | gene and module |
| `tail_t5_nb_dep()` | `nb_dep` | T5 | module only; needs `adjacency` |
| `quantile_band_flags()` | quantile regression | QR | genes in one replicate |

`level` defaults to `0.01` per side. A stored threshold is the positive
distance used by the flag: upper flag `statistic >= threshold`, lower
flag `statistic <= -threshold`. On a gene with zero variance of `D`,
`tail_threshold()` returns `NA` thresholds and warns once.

Plug-in: deMixSC-style weighted NNLS, `c_reg = 2`, empty G2 (no
same-aliquot benchmark), `a_i = sum(beta)`, `p = beta / a_i`. The
Lawson-Hanson solver is inlined; the `nnls` package is not required.

Cook flags: leverage above `h_mult` times its mean (default 3) and
`|standardised residual| > r_max` (default 3). Both conditions. Nothing
is clipped.

Driver: `scripts/02_03_cross_platform_consistency.R`.

| Argument | Value |
|---|---|
| `z_range`, `mad_k` | 3, 3 |
| `level` | 0.01 |
| `c_reg` | 2 |
| `max_degree`, `mean_degree` | 10, 3 |
| `h_mult`, `r_max` | 3, 3 |

The panel is the tier (ii) selected genes united with mapped marker
symbols. Bulk columns at each time are the three `_CTL_` samples from
`standard_wnt_samples()`. Modules are the lineage groups in
`data/dictionaries/mouse_gastruloid_signaling_markers.csv`.

Writes `output/feature_selection/tier_03/`. Not run.

## Perspective — `R/02_05_classifiers.R`

`classifier_markers()`: stability selection with `glmnet` (grouped
multinomial), `mixOmics::splsda()`, and `ranger`. Defaults in the
chapter call: `n_subsample = 100`, `frac = 0.5`, `n_top = 50`,
`threshold = 0.6`, `seed = 20261104`. No driver script in this pass.
The call is the `eval: false` chunk `fs-classifiers-call` in the chapter.

## Legacy archive — `scripts/legacy/`

Do not import these from a live driver.

| File | Role |
|---|---|
| `scripts/legacy/02_90_stats_rank_sctransform.R` | `rank_genes_sctransform()`, `anova_ss()`, `anova_ms_equal_type()`, `rank_ratio()`. Calls `Seurat::SCTransform()` |
| `scripts/legacy/02_90_stats_prepare_panel.R` | The run that built the current panel (outputs 2026-09-28 to 2026-10-01). `n_genes_global = 500`, `n_genes_celltype = 50`, `min_cells = 20`, `sct_seed = 1` |
| `scripts/legacy/02_90_vis_legacy_panel.R` | Redraws `output/naive_marker_selection/` tables. No barcode distribution |

The chapter appendix `Legacy panel: SCTransform ratio` (`#sec-fs-legacy`)
is `.appendix .unnumbered .unlisted`, as are `#sec-fs-sct` and
`#sec-fs-ms-ratio`. It documents that historical fit. It is not the
procedure to rerun for a new panel.

## What not to change by accident

- Do not move `DeCovarT/R/01_02_generate_synthetic_networks.R`.
- Do not put `SCTransform` files back in `R/` or in the root of `scripts/`.
- Do not set `eval: true` on the chapter chunks.
- Do not threshold the mitochondrial fraction unless a cutoff is supplied on purpose. None is fixed.
- Do not treat the legacy 500-gene RDS (`data/intermediate/equal_type_top500_by_time.rds`) as the output of tiers (ii) or (iii). Those tiers have not been run.
