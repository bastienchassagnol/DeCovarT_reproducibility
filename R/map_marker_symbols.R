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
    stop("Package AnnotationDbi is required for offline mouse symbol mapping.")
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
