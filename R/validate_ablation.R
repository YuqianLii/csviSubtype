#' Summarize ablation stability
#'
#' @param feature_sets A named list of character vectors or data.frames.
#' @param footprint_sets A named list of character vectors or data.frames.
#' @param groups Character vector of group names.
#'
#' @return A list of summary tables.
#' @export
csvi_summarize_ablation <- function(
    feature_sets,
    footprint_sets,
    groups = c("A1", "B1")
) {
  normalize_tokens <- function(x) {
    if (length(x) == 0L) return(character(0))
    x <- unlist(strsplit(as.character(x), "[;,]", perl = TRUE))
    x <- trimws(x)
    x <- toupper(x)
    x <- x[!is.na(x) & x != ""]
    unique(x)
  }

  required_keys <- as.vector(unlist(lapply(
    groups,
    function(g) paste0(c("full", "noCSS", "noCIS", "noLRS"), "_", g)
  )))

  miss_feature <- setdiff(required_keys, names(feature_sets))
  miss_fp <- setdiff(required_keys, names(footprint_sets))

  if (length(miss_feature) > 0L) {
    stop(
      "feature_sets is missing keys: ",
      paste(miss_feature, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(miss_fp) > 0L) {
    stop(
      "footprint_sets is missing keys: ",
      paste(miss_fp, collapse = ", "),
      call. = FALSE
    )
  }

  get_genes <- function(x) {
    if (is.null(x)) return(character(0))

    if (is.character(x)) {
      return(normalize_tokens(x))
    }

    x <- data.table::as.data.table(x)
    if (nrow(x) == 0L) return(character(0))

    candidate_cols <- c("gene", "mouse_gene", "gene_key", "human_gene_any")
    gcol <- candidate_cols[candidate_cols %in% names(x)][1]

    if (is.na(gcol)) {
      stop("No gene-like column found. Expected one of: gene, mouse_gene, gene_key, human_gene_any", call. = FALSE)
    }

    normalize_tokens(x[[gcol]])
  }

  feature_summary <- data.table::rbindlist(lapply(groups, function(g) {
    full_genes <- get_genes(feature_sets[[paste0("full_", g)]])

    data.table::rbindlist(lapply(c("noCSS", "noCIS", "noLRS"), function(m) {
      gm <- get_genes(feature_sets[[paste0(m, "_", g)]])

      data.table::data.table(
        groupAB = g,
        mode = m,
        size_feat_full = length(full_genes),
        size_feat_mode = length(gm),
        jacc_feat_vs_full = jaccard_similarity(full_genes, gm)
      )
    }))
  }))

  footprint_summary <- data.table::rbindlist(lapply(groups, function(g) {
    full_genes <- get_genes(footprint_sets[[paste0("full_", g)]])

    data.table::rbindlist(lapply(c("noCSS", "noCIS", "noLRS"), function(m) {
      gm <- get_genes(footprint_sets[[paste0(m, "_", g)]])

      data.table::data.table(
        groupAB = g,
        mode = m,
        size_fp_full = length(full_genes),
        size_fp_mode = length(gm),
        jacc_fp_vs_full = jaccard_similarity(full_genes, gm),
        dropped = paste(setdiff(full_genes, gm), collapse = ";"),
        added = paste(setdiff(gm, full_genes), collapse = ";")
      )
    }))
  }))

  merged <- merge(
    feature_summary,
    footprint_summary,
    by = c("groupAB", "mode"),
    all = TRUE
  )

  list(
    feature_summary = feature_summary,
    footprint_summary = footprint_summary,
    merged = merged
  )
}
