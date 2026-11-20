# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_05_vis_enrichment.R \
#   > "logs/02_05_vis_enrichment_$(date +%F).log" 2>&1 &
#
# Enrichment figures (appendix of chapter 2): summary lollipop over the
# curated themes (data/dictionaries/enrichment_themes.csv, a draft
# interpretation written by a language model and checked against the
# tables), Manhattan plots of the ORA, GSEA volcano, semantic map of the
# significant GO BP terms.

root <- here::here()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_05_vis_enrichment.R"))

enr_dir <- file.path(root, "output", "feature_selection", "enrichment")
fig_dir <- file.path(enr_dir, "figures")
png_dir <- file.path(root, "docs", "images", "feature_selection")
enr <- utils::read.csv(file.path(enr_dir, "enrichment.csv"), stringsAsFactors = FALSE)
sig <- utils::read.csv(file.path(enr_dir, "enrichment_significant.csv"), stringsAsFactors = FALSE)
themes_path <- file.path(root, "data", "dictionaries", "enrichment_themes.csv")

if (file.exists(themes_path)) {
  themes <- utils::read.csv(themes_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  save_figure(
    plot_enrichment_lollipop(
      sig,
      themes,
      title = "Enriched themes: metabolism, signalling and lineage programmes per type",
      subtitle = "Curated, non-redundant terms (draft interpretation, model-written, checked against the tables); BH 0.05"
    ),
    "fs_enrichment_lollipop",
    fig_dir,
    png_dir = png_dir,
    data = list(themes = themes, terms = merge(themes[, c("theme", "time_point", "cell_type", "method", "database", "id")], sig)),
    width = 11,
    height = 0.3 * nrow(themes) + 2.4
  )
}
n_sig <- stats::aggregate(cbind(n = padj <= 0.05) ~ method + database, enr, sum)
save_figure(
  plot_enrichment_manhattan(
    enr,
    title = "Over-representation of the genes confirmed above the average of the other types",
    subtitle = "One point per term; dashed: largest p passing BH 0.05 in the panel; universe: genes tested in the slot"
  ),
  "fs_enrichment_manhattan",
  fig_dir,
  png_dir = png_dir,
  data = enr[enr$method == "ORA", c("time_point", "cell_type", "database", "id", "description", "pval", "padj")],
  width = 12,
  height = 8.5
)
save_figure(
  plot_gsea_volcano(
    enr,
    title = "Gene-set enrichment along the one-versus-average ranking",
    subtitle = "Positive NES: set concentrated among genes above the average of the other types; dashed: BH 0.05"
  ),
  "fs_enrichment_gsea_volcano",
  fig_dir,
  png_dir = png_dir,
  data = enr[enr$method == "GSEA", c("time_point", "cell_type", "database", "id", "description", "nes", "padj")],
  width = 12,
  height = 8.5
)
clusters_path <- file.path(enr_dir, "go_clusters.csv")
if (file.exists(clusters_path)) {
  clusters <- utils::read.csv(clusters_path, stringsAsFactors = FALSE)
  save_figure(
    plot_semantic_map(
      clusters,
      sig,
      title = sprintf("%s significant GO BP terms fall into %s semantic clusters", fmt_n(nrow(clusters)), fmt_n(length(unique(clusters$cluster)))),
      subtitle = "MDS of 1 - semantic similarity (binary cut at 0.7); size: type x time sets with the term"
    ),
    "fs_enrichment_semantic_map",
    fig_dir,
    png_dir = png_dir,
    data = clusters,
    width = 9,
    height = 8
  )
}
message("Wrote ", fig_dir)
