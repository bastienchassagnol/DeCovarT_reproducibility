#' Enrichment of the genes that drive the type effect (appendix of chapter 2).
#'
#' Over-representation (ORA) of the genes confirmed as enriched in one type
#' (hypothesis A of `R/02_04_stats_deg_contrasts.R`, primary family:
#' permuted Welch t, graphical procedure), against the genes tested in the
#' slot as universe; and gene-set enrichment (GSEA) of all tested genes
#' ranked by the one-versus-average Welch t. Databases: GO biological
#' process (`clusterProfiler`), Reactome (`ReactomePA`) and the MSigDB
#' Hallmark collection (`msigdbr`).

#' Symbols to Entrez identifiers (mouse), first match kept.
#'
#' @param symbols Gene symbols.
#' @return Named character vector of Entrez ids (names = symbols), `NA`
#'   when unmapped.
symbol_to_entrez <- function(symbols) {
  map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(
    unique(symbols),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db"
  )))
  map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
  stats::setNames(map$ENTREZID[match(symbols, map$SYMBOL)], symbols)
}

#' Hallmark gene sets for mouse as a TERM2GENE table.
#'
#' @return Data frame `term`, `gene` (Entrez ids).
hallmark_term2gene <- function() {
  h <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
  unique(data.frame(term = h$gs_name, gene = as.character(h$ncbi_gene)))
}

#' Tidy one clusterProfiler result into the common long format.
#'
#' @param res An `enrichResult` or `gseaResult`, or `NULL`.
#' @param method `"ORA"` or `"GSEA"`.
#' @param database Database label.
#' @return Data frame (possibly empty) with `method`, `database`, `id`,
#'   `description`, `set_size`, `count`, `nes`, `pval`, `padj`, `genes`.
tidy_enrichment <- function(res, method, database) {
  empty <- data.frame(
    method = character(), database = character(), id = character(),
    description = character(), set_size = numeric(), count = numeric(),
    nes = numeric(), pval = numeric(), padj = numeric(), genes = character()
  )
  if (is.null(res)) {
    return(empty)
  }
  d <- as.data.frame(res)
  if (nrow(d) == 0L) {
    return(empty)
  }
  if (method == "ORA") {
    data.frame(
      method = method,
      database = database,
      id = d$ID,
      description = d$Description,
      set_size = as.numeric(sub("/.*", "", d$BgRatio)),
      count = d$Count,
      nes = NA_real_,
      pval = d$pvalue,
      padj = d$p.adjust,
      genes = d$geneID
    )
  } else {
    data.frame(
      method = method,
      database = database,
      id = d$ID,
      description = d$Description,
      set_size = d$setSize,
      count = lengths(strsplit(d$core_enrichment, "/")),
      nes = d$NES,
      pval = d$pvalue,
      padj = d$p.adjust,
      genes = d$core_enrichment
    )
  }
}

#' ORA and GSEA of one type in one slot, three databases.
#'
#' @param confirmed Entrez ids of the confirmed genes of the type.
#' @param universe Entrez ids of the tested genes of the slot.
#' @param ranks Named numeric vector (Entrez ids), the ranking statistic.
#' @param t2g_hallmark TERM2GENE table from [hallmark_term2gene()].
#' @param min_genes Smallest number of confirmed genes for an ORA.
#'   Default `10`.
#' @param seed Seed for the GSEA permutations.
#' @return Data frame in the format of [tidy_enrichment()]; all terms are
#'   kept (`pvalueCutoff = 1`), the BH filter is applied downstream.
enrich_one_type <- function(confirmed, universe, ranks, t2g_hallmark, min_genes = 10L, seed = 1L) {
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  ranks <- sort(ranks[!is.na(names(ranks)) & is.finite(ranks)], decreasing = TRUE)
  ranks <- ranks[!duplicated(names(ranks))]
  out <- list()
  if (length(confirmed) >= min_genes) {
    out$ora_go <- tidy_enrichment(safe(clusterProfiler::enrichGO(
      confirmed, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID", ont = "BP",
      universe = universe, pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 500, readable = TRUE
    )), "ORA", "GO BP")
    out$ora_reactome <- tidy_enrichment(safe(ReactomePA::enrichPathway(
      confirmed, organism = "mouse", universe = universe, pvalueCutoff = 1,
      qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = TRUE
    )), "ORA", "Reactome")
    out$ora_hallmark <- tidy_enrichment(safe(clusterProfiler::enricher(
      confirmed, universe = universe, TERM2GENE = t2g_hallmark,
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500
    )), "ORA", "Hallmark")
  }
  withr::with_seed(seed, {
    out$gsea_go <- tidy_enrichment(safe(clusterProfiler::gseGO(
      ranks, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID", ont = "BP",
      minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, eps = 0,
      verbose = FALSE, seed = TRUE
    )), "GSEA", "GO BP")
    # gsePathway() takes no `eps` or `seed`; the outer seed fixes the draw.
    out$gsea_reactome <- tidy_enrichment(safe(ReactomePA::gsePathway(
      ranks, organism = "mouse", minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 1, verbose = FALSE
    )), "GSEA", "Reactome")
    out$gsea_hallmark <- tidy_enrichment(safe(clusterProfiler::GSEA(
      ranks, TERM2GENE = t2g_hallmark, minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 1, eps = 0, verbose = FALSE, seed = TRUE
    )), "GSEA", "Hallmark")
  })
  do.call(rbind, out)
}
