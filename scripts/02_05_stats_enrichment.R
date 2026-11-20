# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_05_stats_enrichment.R \
#   > "logs/02_05_stats_enrichment_$(date +%F).log" 2>&1 &
#
# ORA and GSEA per type and slot (appendix of chapter 2), from the post hoc
# contrasts of scripts/02_04_stats_deg_contrasts.R. Primary family: the
# permuted Welch t with the graphical within-gene procedure. Writes
# output/feature_selection/enrichment/enrichment.csv (all terms),
# enrichment_significant.csv (BH <= 0.05) and go_clusters.csv
# (simplifyEnrichment clusters of the significant GO BP terms).

root <- here::here()
source(file.path(root, "R", "02_05_stats_enrichment.R"))

fdr <- 0.05
primary_test <- "welch_perm"
seed <- 20261121L
deg_dir <- file.path(root, "output", "feature_selection", "deg")
out_dir <- file.path(root, "output", "feature_selection", "enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg <- readRDS(file.path(deg_dir, "deg.rds"))
dec <- deg$decisions[deg$decisions$test == primary_test, , drop = FALSE]
t2g <- hallmark_term2gene()

rows <- list()
for (tm in names(deg$results)) {
  res <- deg$results[[tm]]
  genes <- rownames(res$log2fc)
  entrez <- symbol_to_entrez(genes)
  universe <- unique(stats::na.omit(entrez))
  message(tm, ": ", length(universe), " of ", length(genes), " genes mapped to Entrez")
  for (tp in res$types) {
    d <- dec[dec$time_point == tm & dec$type == tp, , drop = FALSE]
    confirmed <- unique(stats::na.omit(entrez[d$gene[d$a_graph %in% TRUE]]))
    ranks <- stats::setNames(res$p$welch_perm$t_a[, tp], entrez[genes])
    e <- enrich_one_type(confirmed, universe, ranks, t2g, seed = seed)
    if (!is.null(e) && nrow(e) > 0L) {
      e$time_point <- tm
      e$cell_type <- tp
      e$n_confirmed <- length(confirmed)
      rows[[paste(tm, tp)]] <- e
    }
    message("  ", tp, ": ", length(confirmed), " confirmed genes, ",
            sum(e$padj <= fdr, na.rm = TRUE), " terms at BH ", fdr)
  }
}
enr <- do.call(rbind, rows)
rownames(enr) <- NULL
utils::write.csv(enr, file.path(out_dir, "enrichment.csv"), row.names = FALSE, fileEncoding = "UTF-8")
sig <- enr[enr$padj <= fdr & !is.na(enr$padj), , drop = FALSE]
utils::write.csv(sig, file.path(out_dir, "enrichment_significant.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Semantic clusters of the significant GO BP terms (all types, both methods).
go_ids <- unique(sig$id[sig$database == "GO BP"])
clusters <- NULL
if (length(go_ids) >= 5L) {
  sim <- simplifyEnrichment::GO_similarity(go_ids, ont = "BP", db = "org.Mm.eg.db")
  # cutoff 0.7 instead of the default 0.85: the default splits the 1,902
  # terms of this run into 2 clusters only.
  cl <- simplifyEnrichment::cluster_terms(sim, method = "binary_cut", control = list(cutoff = 0.7))
  mds <- stats::cmdscale(stats::as.dist(1 - sim), k = 2)
  clusters <- data.frame(id = rownames(sim), cluster = cl, mds1 = mds[, 1], mds2 = mds[, 2])
  utils::write.csv(clusters, file.path(out_dir, "go_clusters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}
saveRDS(list(enrichment = enr, significant = sig, go_clusters = clusters, fdr = fdr, primary_test = primary_test),
        file.path(out_dir, "enrichment.rds"))
utils::capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
message("Wrote ", out_dir)
