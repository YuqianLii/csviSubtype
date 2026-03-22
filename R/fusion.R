#' Fuse GWAS and scRNA evidence into subtype barcodes
#'
#' @param gwas A data.frame/data.table of GWAS gene-level evidence.
#' @param scrna_a1 A data.frame/data.table of A1 scRNA feature table.
#' @param scrna_b1 A data.frame/data.table of B1 scRNA feature table.
#' @param gwas_high Numeric threshold for high GWAS tier.
#' @param gwas_mid Numeric threshold for mid GWAS tier.
#' @param scrna_high_q Numeric threshold for high scRNA tier.
#' @param scrna_mid_q Numeric threshold for mid scRNA tier.
#'
#' @return A named list containing GWAS summary and A1/B1 full/binary/strict ledgers.
#' @export
csvi_fuse_barcode <- function(
    gwas,
    scrna_a1,
    scrna_b1,
    gwas_high = 0.80,
    gwas_mid = 0.50,
    scrna_high_q = 0.80,
    scrna_mid_q = 0.50
) {
  gwas <- data.table::as.data.table(gwas)
  scrna_a1 <- data.table::as.data.table(scrna_a1)
  scrna_b1 <- data.table::as.data.table(scrna_b1)

  tier_gwas <- function(pp4, high = 0.8, mid = 0.5) {
    if (is.na(pp4)) return("none")
    if (pp4 >= high) return("high")
    if (pp4 >= mid) return("mid")
    "low"
  }

  zone_rule_strict <- function(gwas_tier, scrna_tier) {
    if (gwas_tier == "high" && scrna_tier == "high") return("red")
    if ((gwas_tier == "high" && scrna_tier == "mid") ||
        (gwas_tier == "mid" && scrna_tier == "high")) return("orange")
    if (gwas_tier == "mid" && scrna_tier == "mid") return("yellow")
    "grey"
  }

  prep_scrna <- function(dt) {
    dt <- data.table::copy(dt)

    if (!("gene" %in% names(dt))) {
      stop("Each scRNA table must contain a 'gene' column.", call. = FALSE)
    }

    dt[, gene := toupper(as.character(gene))]
    dt[, gene_key := gene]

    if (!("q_sc" %in% names(dt))) {
      if (!("M" %in% names(dt))) stop("scRNA table must contain 'M' when q_sc is absent.", call. = FALSE)
      if (!("F_gene" %in% names(dt))) stop("scRNA table must contain 'F_gene' when q_sc is absent.", call. = FALSE)

      if (!("q_M" %in% names(dt))) {
        dt[, q_M := percentile_rank_desc(M)]
      }

      if (!("q_F" %in% names(dt))) {
        dt[, q_F := 0]
        dt[F_gene > 0, q_F := percentile_rank_desc(F_gene)]
      }

      dt[, q_sc := pmax(q_M, q_F)]
    }

    dt[, scrna_tier := assign_scrna_tier(q_sc, high = scrna_high_q, mid = scrna_mid_q)]
    data.table::setorder(dt, -q_sc, -F_gene, -M)
    dt[!duplicated(gene_key)]
  }

  req_gwas <- c("trait", "gene", "mouse_gene_1to1", "gwas_pp4")
  if (!all(req_gwas %in% names(gwas))) {
    stop(
      "gwas must contain columns: ", paste(req_gwas, collapse = ", "),
      call. = FALSE
    )
  }

  gwas[, gene := as.character(gene)]
  gwas[, mouse_gene_1to1 := toupper(as.character(mouse_gene_1to1))]
  gwas <- gwas[!is.na(mouse_gene_1to1) & mouse_gene_1to1 != ""]
  gwas[, gene_key := mouse_gene_1to1]

  gwas_gene <- gwas[
    ,
    list(
      mouse_gene = mouse_gene_1to1[1],
      human_gene_any = paste(unique(gene), collapse = ";"),
      gwas_pp4_max = max(gwas_pp4, na.rm = TRUE),
      gwas_traits = paste(unique(trait), collapse = ";")
    ),
    by = "gene_key"
  ]

  gwas_gene[is.infinite(gwas_pp4_max), gwas_pp4_max := NA_real_]
  gwas_gene[, gwas_tier := vapply(gwas_pp4_max, tier_gwas, character(1), high = gwas_high, mid = gwas_mid)]

  A1 <- prep_scrna(scrna_a1)
  B1 <- prep_scrna(scrna_b1)

  merge_one <- function(scrna_dt) {
    full <- merge(gwas_gene, scrna_dt, by = "gene_key", all = TRUE)

    full[, has_gwas := !is.na(gwas_pp4_max)]
    full[, has_scrna := !is.na(groupAB) & groupAB != ""]

    full[, gwas_tier := data.table::fifelse(is.na(gwas_tier), "none", gwas_tier)]
    full[, scrna_tier := data.table::fifelse(is.na(scrna_tier), "none", scrna_tier)]

    full[, zone_binary := data.table::fifelse(has_gwas & has_scrna, "red", "grey")]
    full[, include_binary := (zone_binary == "red")]

    full[, zone_strict := mapply(zone_rule_strict, gwas_tier, scrna_tier)]
    full[, include_strict := zone_strict %in% c("red", "orange")]

    if (!("mouse_gene" %in% names(full))) {
      full[, mouse_gene := NA_character_]
    }
    full[is.na(mouse_gene) & !is.na(gene), mouse_gene := gene]

    list(
      full = full,
      binary = full[include_binary == TRUE],
      strict = full[include_strict == TRUE]
    )
  }

  list(
    gwas_gene_summary = gwas_gene,
    A1 = merge_one(A1),
    B1 = merge_one(B1)
  )
}
