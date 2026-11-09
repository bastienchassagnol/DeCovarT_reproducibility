#' Tier (i) of the gene selection: cells, genes, markers, presence.
#'
#' Live fits do not call `Seurat::SCTransform()`. Cell filters, the
#' bulk-single-cell symbol intersection, the mitochondrial-symbol drop
#' and the presence gate live here. `glmGamPoi::glm_gp()` is the mean
#' model of tiers (ii) and (iii), in `R/02_02_cross_cell.R` and
#' `R/02_03_influence.R`.
#'
#' Requires `R/02_00_utils.R` (`%||%`) to be sourced first.
#' `org.Mm.eg.db` is used only by [map_marker_symbols()]. Written, not run.

#' Map dictionary gene symbols onto mouse symbols in a count matrix.
#'
#' Resolution is offline: case-fold against `universe`, then
#' `org.Mm.eg.db` `SYMBOL` and `ALIAS`. There is no `biomaRt` call.
#'
#' @param markers Character vector of symbols from the signalling
#'   dictionary (`gene_marker`).
#' @param universe Character vector of symbols present in the RNA
#'   count matrix.
#' @return A data frame with `gene_marker`, `mouse_symbol`, and
#'   `map_method`. `mouse_symbol` is `NA` when the marker is
#'   unmapped or ambiguous.
#' @details
#'   `map_method` is one of `casefold_universe`,
#'   `ambiguous_universe`, `org.Mm.eg.db_symbol`,
#'   `org.Mm.eg.db_symbol_absent`, `org.Mm.eg.db_alias`,
#'   `org.Mm.eg.db_alias_absent`, `ambiguous_alias`, or `unmapped`.
map_marker_symbols <- function(markers, universe) {
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop(
      "Package org.Mm.eg.db is required for offline mouse symbol ",
      "mapping. Install it with BiocManager::install(\"org.Mm.eg.db\")."
    )
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Package AnnotationDbi is required for offline mouse symbol mapping."
    )
  }
  markers <- unique(as.character(markers))
  universe <- unique(as.character(universe))
  upper_hits <- split(universe, toupper(universe))
  symbols <- AnnotationDbi::keys(
    org.Mm.eg.db::org.Mm.eg.db,
    keytype = "SYMBOL"
  )
  alias_tbl <- suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = symbols,
      columns = "ALIAS",
      keytype = "SYMBOL"
    )
  )
  alias_tbl <- alias_tbl[
    !is.na(alias_tbl$ALIAS) & nzchar(alias_tbl$ALIAS),
    ,
    drop = FALSE
  ]
  alias_upper <- split(alias_tbl$SYMBOL, toupper(alias_tbl$ALIAS))
  symbol_upper <- split(symbols, toupper(symbols))

  resolve_one <- function(marker) {
    key <- toupper(marker)
    in_matrix <- upper_hits[[key]] %||% character()
    if (length(unique(in_matrix)) == 1L) {
      return(list(mouse = unique(in_matrix), method = "casefold_universe"))
    }
    if (length(unique(in_matrix)) > 1L) {
      return(list(mouse = NA_character_, method = "ambiguous_universe"))
    }
    db_hit <- unique(symbol_upper[[key]] %||% character())
    in_universe <- intersect(db_hit, universe)
    if (length(in_universe) == 1L) {
      return(list(mouse = in_universe, method = "org.Mm.eg.db_symbol"))
    }
    if (length(db_hit) == 1L && length(in_universe) == 0L) {
      return(list(mouse = db_hit, method = "org.Mm.eg.db_symbol_absent"))
    }
    alias_hit <- unique(alias_upper[[key]] %||% character())
    alias_in <- intersect(alias_hit, universe)
    if (length(alias_in) == 1L) {
      return(list(mouse = alias_in, method = "org.Mm.eg.db_alias"))
    }
    if (length(alias_hit) == 1L && length(alias_in) == 0L) {
      return(list(
        mouse = alias_hit,
        method = "org.Mm.eg.db_alias_absent"
      ))
    }
    if (length(alias_in) > 1L || length(db_hit) > 1L) {
      return(list(mouse = NA_character_, method = "ambiguous_alias"))
    }
    list(mouse = NA_character_, method = "unmapped")
  }

  resolved <- lapply(markers, resolve_one)
  data.frame(
    gene_marker = markers,
    mouse_symbol = vapply(resolved, function(x) x$mouse, character(1)),
    map_method = vapply(resolved, function(x) x$method, character(1)),
    stringsAsFactors = FALSE
  )
}

#' Cells in a bulk-matched time slot with a label and a library.
#'
#' Keeps a cell when its time is one of `time_levels`, its cell-type
#' label is non-missing, and its library size is positive. Artefact
#' labels and the minimum cell count are applied later, by
#' [scored_types()].
#'
#' @param time_point Time label per cell, for example `"48h"`.
#' @param cell_type Cell-type label per cell.
#' @param library_size Full RNA library size per cell (column sum of
#'   the unfiltered count matrix).
#' @param time_levels Times that match the standard-Wnt bulk.
#'   Default `c("48h", "72h", "96h")`.
#' @return A logical vector, one value per cell.
keep_bulk_matched_cells <- function(
  time_point,
  cell_type,
  library_size,
  time_levels = c("48h", "72h", "96h")
) {
  n <- length(time_point)
  if (length(cell_type) != n || length(library_size) != n) {
    stop(
      "`time_point`, `cell_type` and `library_size` ",
      "must have one value per cell."
    )
  }
  time_point <- as.character(time_point)
  cell_type <- as.character(cell_type)
  ok_time <- time_point %in% time_levels
  ok_type <- !is.na(cell_type) & nzchar(cell_type)
  ok_lib <- is.finite(library_size) & library_size > 0
  ok_time & ok_type & ok_lib
}

#' Cell types that enter the scores of one slot.
#'
#' A type is scored when it is in `allowed` and the slot has at least
#' `min_cells` cells of that type. Artefact labels stay out of `allowed`
#' (`include_in_ranking` is FALSE in the crosswalk).
#'
#' @param cell_type Cell-type labels of one time slot.
#' @param allowed Types that may be scored.
#' @param min_cells Minimum cells of a type in the slot. Default `20`.
#' @return Character vector of scored type labels.
scored_types <- function(cell_type, allowed, min_cells = 20L) {
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 1) {
    stop("`min_cells` must be one number at least 1.")
  }
  cell_type <- as.character(cell_type)
  allowed <- unique(as.character(allowed))
  tab <- table(cell_type, useNA = "no")
  keep <- names(tab)[as.integer(tab) >= min_cells]
  intersect(keep, allowed)
}

#' Read the Suppinger bulk HTSeq table.
#'
#' Row names are Ensembl ids. `gene_name` holds the mouse symbol.
#' Count columns are every other column. The file is
#' `data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz`.
#'
#' @param path Path to the gzipped table.
#' @return A list with `genes` (columns `row_id`, `ensembl_id` with the
#'   version suffix removed, `gene_name`) and `counts` (numeric matrix,
#'   row names equal to `row_id`).
read_bulk_htseq <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing ", path, ". Place ",
      "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz in data/raw/."
    )
  }
  tab <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!"gene_name" %in% colnames(tab)) {
    stop("Column gene_name is absent from the bulk HTSeq table.")
  }
  row_id <- rownames(tab)
  genes <- data.frame(
    row_id = row_id,
    ensembl_id = sub("\\..*$", "", row_id),
    gene_name = as.character(tab$gene_name),
    stringsAsFactors = FALSE
  )
  count_cols <- setdiff(colnames(tab), "gene_name")
  counts <- as.matrix(tab[, count_cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- row_id
  list(genes = genes, counts = counts)
}

#' Standard-Wnt bulk columns at one time.
#'
#' Column names follow the Suppinger HTSeq grammar used in
#' GastroDeconv2FateMap: a token `48h` / `72h` / `96h`, the arm marker
#' `_CTL_` (standard Wnt) or `_Early_` (early Wnt), and `rep` plus a
#' digit. The matched pair uses the three `_CTL_` columns at `time_point`.
#'
#' @param sample_names Column names of the HTSeq count matrix.
#' @param time_point One of `"48h"`, `"72h"`, `"96h"`.
#' @return The three matching column names.
standard_wnt_samples <- function(sample_names, time_point) {
  sample_names <- as.character(sample_names)
  has_time <- grepl("\\d+h", sample_names)
  time_of <- ifelse(
    has_time,
    sub(".*?(\\d+h).*", "\\1", sample_names),
    NA_character_
  )
  keep <- !is.na(time_of) &
    time_of == time_point &
    grepl("_CTL_", sample_names, fixed = TRUE)
  if (sum(keep) != 3L) {
    stop(
      "Expected 3 standard-Wnt (_CTL_) samples at ",
      time_point,
      "; matched ",
      sum(keep),
      ". Columns: ",
      paste(sample_names, collapse = ", "),
      "."
    )
  }
  sample_names[keep]
}

#' Symbols shared by the single-cell RNA assay and the bulk HTSeq table.
#'
#' The intersection is on mouse symbols. A symbol is dropped when the
#' bulk rows that carry it have more than one Ensembl id after the
#' version suffix is removed. That is the alignment used for Suppinger
#' in GastroDeconv2FateMap `scripts/02_01_deconvolution_preprocessing.R`.
#' Symbols with exactly one Ensembl id are kept. The first bulk row of
#' that id is recorded in `row_id`.
#'
#' @param sc_genes Symbols in the single-cell RNA matrix.
#' @param bulk_genes `gene_name` of [read_bulk_htseq()], one per bulk row.
#' @param bulk_ensembl Stripped Ensembl id per bulk row.
#' @param bulk_row_id Row name per bulk row. Default `NULL` uses the
#'   position as a character id.
#' @return A data frame with `gene`, `n_ensembl`, `keep` and `row_id`
#'   (`NA` when `keep` is FALSE), in the order of `sc_genes`.
shared_bulk_single_cell_genes <- function(
  sc_genes,
  bulk_genes,
  bulk_ensembl,
  bulk_row_id = NULL
) {
  sc_genes <- unique(as.character(sc_genes))
  bulk_genes <- as.character(bulk_genes)
  bulk_ensembl <- sub("\\..*$", "", as.character(bulk_ensembl))
  if (length(bulk_genes) != length(bulk_ensembl)) {
    stop("`bulk_genes` and `bulk_ensembl` must have one entry per bulk row.")
  }
  if (is.null(bulk_row_id)) {
    bulk_row_id <- as.character(seq_along(bulk_genes))
  }
  if (length(bulk_row_id) != length(bulk_genes)) {
    stop("`bulk_row_id` must have one entry per bulk row.")
  }
  shared <- intersect(sc_genes, unique(bulk_genes))
  empty <- data.frame(
    gene = character(),
    n_ensembl = integer(),
    keep = logical(),
    row_id = character(),
    stringsAsFactors = FALSE
  )
  if (length(shared) == 0L) {
    return(empty)
  }
  in_shared <- bulk_genes %in% shared
  gene_v <- bulk_genes[in_shared]
  ens_v <- bulk_ensembl[in_shared]
  row_v <- as.character(bulk_row_id[in_shared])
  n_map <- tapply(ens_v, gene_v, function(ids) length(unique(ids)))
  first_row <- tapply(seq_along(gene_v), gene_v, function(idx) {
    ids <- unique(ens_v[idx])
    if (length(ids) != 1L) {
      return(NA_character_)
    }
    hit <- row_v[idx][ens_v[idx] == ids]
    hit[[1L]]
  })
  data.frame(
    gene = shared,
    n_ensembl = as.integer(n_map[shared]),
    keep = as.integer(n_map[shared]) == 1L,
    row_id = as.character(first_row[shared]),
    stringsAsFactors = FALSE
  )
}

#' Bulk counts on unambiguous symbols.
#'
#' @param htseq Output of [read_bulk_htseq()].
#' @param shared Output of [shared_bulk_single_cell_genes()].
#' @param genes Symbols to extract. Default: every `keep` symbol.
#' @return A numeric matrix, symbols as row names, bulk samples as
#'   columns. The row is the `row_id` stored for that symbol.
bulk_symbol_counts <- function(htseq, shared, genes = NULL) {
  kept <- shared[shared$keep, , drop = FALSE]
  if (is.null(genes)) {
    genes <- kept$gene
  }
  genes <- as.character(genes)
  missing <- setdiff(genes, kept$gene)
  if (length(missing) > 0L) {
    stop(
      length(missing),
      " genes are absent from the unambiguous shared symbols."
    )
  }
  row_id <- kept$row_id[match(genes, kept$gene)]
  if (anyNA(row_id) || any(!row_id %in% rownames(htseq$counts))) {
    stop("A kept symbol has no row in the bulk count matrix.")
  }
  out <- htseq$counts[row_id, , drop = FALSE]
  rownames(out) <- genes
  out
}

#' Mouse mitochondrial symbols (`mt-`, either case).
#'
#' @param genes Character vector of gene symbols.
#' @return A logical vector, TRUE for symbols matching `^mt-`.
is_mitochondrial_gene <- function(genes) {
  grepl("^mt-", as.character(genes), ignore.case = TRUE)
}

#' Drop mouse mitochondrial symbols.
#'
#' The Suppinger bulk alignment does not drop these symbols. This
#' function is the explicit drop. Ribosomal and sex-chromosome symbols
#' are left in place.
#'
#' @param genes Character vector of gene symbols.
#' @return `genes` with mitochondrial symbols removed.
drop_mitochondrial_genes <- function(genes) {
  genes <- as.character(genes)
  genes[!is_mitochondrial_gene(genes)]
}

#' Share of counts on mitochondrial symbols, per cell.
#'
#' Call this before [drop_mitochondrial_genes()] subsets the matrix.
#' The fraction is stored. No cutoff is applied: the Suppinger
#' preprocessing does not fix one.
#'
#' @param counts Genes by cells. Row names are gene symbols.
#' @return A numeric vector of length `ncol(counts)`, named by the
#'   cell barcodes. `NA` when the cell has library size zero.
cell_mitochondrial_fraction <- function(counts) {
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene symbols as row names.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package Matrix is required to sum count columns.")
  }
  mt <- is_mitochondrial_gene(rownames(counts))
  total <- Matrix::colSums(counts)
  mt_count <- if (any(mt)) {
    Matrix::colSums(counts[mt, , drop = FALSE])
  } else {
    numeric(ncol(counts))
  }
  frac <- ifelse(total > 0, as.numeric(mt_count) / as.numeric(total), NA_real_)
  names(frac) <- colnames(counts)
  frac
}

#' Detection in every cell type.
#'
#' DeCovarT needs a non-singular within-type covariance for every type,
#' so a gene must be detected in a minimum fraction of the cells of
#' **every** type, not only in one, and must vary in every type. The gate
#' also removes perfect markers (zero in some type), which break
#' identifiability once gene-wise platform ratios are allowed.
#'
#' @param counts Genes by cells (dense or sparse `Matrix`).
#' @param cell_type Cell-type label per column of `counts`.
#' @param min_frac Minimum fraction of cells with a non-zero count, in
#'   each type. Default `0.1`; the tier (ii) calls pass `0.05`.
#' @return A data frame with `gene`, `min_detection` (smallest fraction
#'   over types), `weakest_type`, `min_variance` (smallest within-type
#'   variance of the counts), and `pass` (`min_detection >= min_frac` and
#'   `min_variance > 0`), plus the attribute `detection` (genes by types
#'   matrix of fractions).
presence_gate <- function(counts, cell_type, min_frac = 0.1) {
  cell_type <- droplevels(factor(cell_type))
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (anyNA(cell_type)) {
    stop("`cell_type` has missing labels; drop those cells first.")
  }
  types <- levels(cell_type)
  detection <- matrix(
    NA_real_,
    nrow = nrow(counts),
    ncol = length(types),
    dimnames = list(rownames(counts), types)
  )
  variance <- detection
  for (tp in types) {
    idx <- which(cell_type == tp)
    x_j <- counts[, idx, drop = FALSE]
    detection[, tp] <- Matrix::rowSums(x_j > 0) / length(idx)
    m1 <- Matrix::rowMeans(x_j)
    m2 <- Matrix::rowMeans(x_j * x_j)
    n_j <- length(idx)
    variance[, tp] <- if (n_j > 1L) {
      pmax(m2 - m1^2, 0) * n_j / (n_j - 1)
    } else {
      0
    }
  }
  min_detection <- apply(detection, 1L, min)
  min_variance <- apply(variance, 1L, min)
  out <- data.frame(
    gene = rownames(counts) %||% as.character(seq_len(nrow(counts))),
    min_detection = min_detection,
    weakest_type = types[max.col(-detection, ties.method = "first")],
    min_variance = min_variance,
    pass = min_detection >= min_frac & min_variance > 0,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "detection") <- detection
  out
}

#' Scored cells and shared genes of each bulk-matched time.
#'
#' Applies [keep_bulk_matched_cells()], records
#' [cell_mitochondrial_fraction()] on that cell set, restricts genes to
#' the unambiguous shared symbols that are not mitochondrial, then
#' splits the cells by time and [scored_types()].
#'
#' `library_size` is the full RNA library, computed before any gene
#' subset. The returned slot counts are already subset to the scored
#' types and to the genes that remain.
#'
#' @param counts Genes by cells, symbols as row names.
#' @param time_point,cell_type,library_size One value per column.
#' @param allowed_types Types with `include_in_ranking` TRUE.
#' @param bulk_genes,bulk_ensembl,bulk_row_id Columns of the HTSeq
#'   gene table from [read_bulk_htseq()].
#' @param time_levels Bulk-matched times.
#' @param min_cells Passed to [scored_types()]. Default `20`.
#' @return A list with `slots` (one list per time: `counts`,
#'   `cell_type`, `library_size`, `time_point`, `scored_types`),
#'   `shared`, `mito` (barcode, time, type, fraction) and `dropped_mt`.
prepare_scored_counts <- function(
  counts,
  time_point,
  cell_type,
  library_size,
  allowed_types,
  bulk_genes,
  bulk_ensembl,
  bulk_row_id = NULL,
  time_levels = c("48h", "72h", "96h"),
  min_cells = 20L
) {
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene symbols as row names.")
  }
  keep <- keep_bulk_matched_cells(
    time_point,
    cell_type,
    library_size,
    time_levels
  )
  counts <- counts[, keep, drop = FALSE]
  time_point <- as.character(time_point[keep])
  cell_type <- as.character(cell_type[keep])
  library_size <- library_size[keep]
  fraction <- cell_mitochondrial_fraction(counts)
  mito <- data.frame(
    barcode = colnames(counts),
    time_point = time_point,
    cell_type = cell_type,
    mito_fraction = fraction,
    stringsAsFactors = FALSE
  )
  shared <- shared_bulk_single_cell_genes(
    rownames(counts),
    bulk_genes,
    bulk_ensembl,
    bulk_row_id
  )
  dropped_mt <- rownames(counts)[is_mitochondrial_gene(rownames(counts))]
  kept_symbols <- drop_mitochondrial_genes(shared$gene[shared$keep])
  counts <- counts[intersect(kept_symbols, rownames(counts)), , drop = FALSE]
  slots <- lapply(time_levels, function(tm) {
    in_time <- time_point == tm
    types <- scored_types(cell_type[in_time], allowed_types, min_cells)
    in_slot <- in_time & cell_type %in% types
    list(
      counts = counts[, in_slot, drop = FALSE],
      cell_type = cell_type[in_slot],
      library_size = library_size[in_slot],
      time_point = tm,
      scored_types = types
    )
  })
  names(slots) <- time_levels
  list(
    slots = slots,
    shared = shared,
    mito = mito,
    dropped_mt = dropped_mt
  )
}
