#' Build scRNA-based subtype gene features
#'
#' This function computes subtype-oriented scRNA feature scores from
#' CSS, CIS, LRS, and DEG tables.
#'
#' @param css A data.frame/data.table for CSS summary.
#' @param cis A data.frame/data.table for CIS summary.
#' @param lrs A data.frame/data.table for ligand-receptor pairs.
#' @param deg A data.frame/data.table for DEG summary.
#' @param tau_vascular Numeric. CSS gate threshold.
#' @param top_k_sender Integer. Number of top senders retained per group.
#' @param gamma Numeric. Exponent used to transform LRS weights.
#' @param eps Numeric. Small constant for numerical stability.
#' @param mode Character. Currently only `"receptor"` is implemented in v0.1.
#' @param scrna_q_high Numeric. High threshold for q_sc tier.
#' @param scrna_q_mid Numeric. Mid threshold for q_sc tier.
#' @param split_complex_receptor Logical. Whether receptor complexes are split.
#'
#' @return A named list with group-level features and intermediate tables.
#' @export
csvi_build_scrna_features <- function(
    css,
    cis,
    lrs,
    deg,
    tau_vascular = 0.10,
    top_k_sender = 3,
    gamma = 1,
    eps = 1e-12,
    mode = c("receptor", "pathway"),
    scrna_q_high = 0.80,
    scrna_q_mid = 0.50,
    split_complex_receptor = TRUE
) {
  mode <- match.arg(mode)

  if (mode == "pathway") {
    stop(
      "pathway mode is not yet implemented in csviSubtype v0.1. ",
      "Use mode = 'receptor' for the package core, and keep the manuscript-level ",
      "pathway workflow in inst/scripts.",
      call. = FALSE
    )
  }

  css <- data.table::as.data.table(css)
  cis <- data.table::as.data.table(cis)
  lrs <- data.table::as.data.table(lrs)
  deg <- data.table::as.data.table(deg)

  if (!("gene" %in% names(deg))) {
    stop("deg must contain 'gene'", call. = FALSE)
  }

  deg[, gene := toupper(as.character(gene))]

  if (!("logFC" %in% names(deg))) {
    if ("avg_log2FC" %in% names(deg)) deg[, logFC := avg_log2FC]
    if ("avg_logFC" %in% names(deg)) deg[, logFC := avg_logFC]
  }

  if (!("logFC" %in% names(deg))) {
    stop("deg must contain 'logFC' or equivalent", call. = FALSE)
  }

  if (!("direction" %in% names(deg))) {
    deg[, direction := data.table::fifelse(logFC > 0, "up_in_B1", "up_in_A1")]
  }

  deg[, M := abs(logFC)]

  css_grp <- standardize_css(css)
  css_grp[, Cp := css_gate_fun(R_v, tau_vascular)]

  cis_dt <- prep_cis(cis, css_grp, eps = eps)

  cis_ordered <- cis_dt[order(groupAB, -S_int)]
  cis_ordered[, sender_rank := seq_len(.N), by = "groupAB"]
  cis_top <- cis_ordered[sender_rank <= top_k_sender]
  cis_top[, sender_rank := NULL]

  lrs_dt <- prep_lrs(lrs, gamma = gamma)

  if (nrow(lrs_dt) > 0L) {
    lrs_dt <- merge(
      lrs_dt,
      cis_top[, c("groupAB", "sender", "S_int"), with = FALSE],
      by.x = c("groupAB", "source"),
      by.y = c("groupAB", "sender"),
      all = FALSE
    )
  }

  if (nrow(lrs_dt) == 0L) {
    lrs_use <- data.table::data.table(
      groupAB = character(),
      source = character(),
      ligand = character(),
      receptor = character(),
      A0 = numeric(),
      S_int = numeric()
    )
  } else if (isTRUE(split_complex_receptor)) {
    lrs_use <- data.table::rbindlist(
      lapply(seq_len(nrow(lrs_dt)), function(i) {
        row <- lrs_dt[i]
        rs <- split_complex(row$receptor)
        if (length(rs) == 0L) return(NULL)

        w <- row$A0 / length(rs)

        data.table::data.table(
          groupAB = row$groupAB,
          source = row$source,
          ligand = row$ligand,
          receptor = rs,
          A0 = w,
          S_int = row$S_int
        )
      }),
      fill = TRUE
    )
  } else {
    lrs_use <- lrs_dt[, c("groupAB", "source", "ligand", "receptor", "A0", "S_int"), with = FALSE]
  }

  if (nrow(lrs_use) > 0L) {
    lrs_use[, A0_hat := safe_div(A0, sum(A0, na.rm = TRUE), eps), by = c("groupAB", "source")]
    lrs_use[!is.finite(A0_hat), A0_hat := 0]
    lrs_use[, contrib := S_int * A0_hat]

    lrs_tmp <- data.table::copy(lrs_use)
    lrs_tmp[, gene := receptor]

    gsum <- lrs_tmp[
      ,
      list(sum_contrib = sum(contrib, na.rm = TRUE)),
      by = c("groupAB", "gene")
    ]
  } else {
    gsum <- data.table::data.table(
      groupAB = character(),
      gene = character(),
      sum_contrib = numeric()
    )
  }

  out <- lapply(c("A1", "B1"), function(g) {
    tmp <- merge(
      deg[, c("gene", "logFC", "direction", "M"), with = FALSE],
      gsum[groupAB == g, c("gene", "sum_contrib"), with = FALSE],
      by = "gene",
      all.x = TRUE
    )

    tmp[is.na(sum_contrib), sum_contrib := 0]
    tmp[, groupAB := g]
    tmp[, F_gene := M * sum_contrib]

    tmp
  })

  out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)

  out[, q_M := percentile_rank_desc(M), by = "groupAB"]
  out[, q_F := 0]
  out[F_gene > 0, q_F := percentile_rank_desc(F_gene), by = "groupAB"]
  out[, q_sc := pmax(q_M, q_F)]
  out[, scrna_tier := assign_scrna_tier(q_sc, high = scrna_q_high, mid = scrna_q_mid)]
  out[, scrna_percentile := q_sc]
  out[, chainSupport := (sum_contrib > 0)]

  list(
    css_group = css_grp,
    cis_all = cis_dt,
    cis_top = cis_top,
    lrs_use = lrs_use,
    features = out,
    features_by_group = list(
      A1 = out[groupAB == "A1"],
      B1 = out[groupAB == "B1"]
    )
  )
}
