# =========================
# baseline_scoreboard_v1.R
# Build baselines + scoreboard from:
#   - footprint_merge_A1_full.tsv / footprint_merge_B1_full.tsv
#   - feature_A1_qsc_top50.tsv / feature_B1_qsc_top50.tsv
#   - (optional) feature_gwas_gene_vascular_strict.tsv (for reference)
# =========================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Config (edit paths if needed) ----
cfg <- list(
  fp_A1 = "footprint_merge_A1_full.tsv",
  fp_B1 = "footprint_merge_B1_full.tsv",
  scr_A1_top50 = "feature_A1_qsc_top50.tsv",
  scr_B1_top50 = "feature_B1_qsc_top50.tsv",
  gwas_strict_ref = "feature_gwas_gene_vascular_strict.tsv",  # optional
  outdir = "output/baselines_v1"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

jaccard <- function(a, b) {
  a <- unique(na.omit(a)); b <- unique(na.omit(b))
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

score_vs_barcode <- function(gene_list, barcode_set, gwas_universe = NULL) {
  gene_list <- unique(na.omit(gene_list))
  barcode_set <- unique(na.omit(barcode_set))
  inter <- intersect(gene_list, barcode_set)
  
  precision <- if (length(gene_list) == 0) NA_real_ else length(inter) / length(gene_list)
  coverage  <- if (length(barcode_set) == 0) NA_real_ else length(inter) / length(barcode_set)
  
  out <- list(
    n_genes = length(gene_list),
    n_overlap_barcode = length(inter),
    precision_vs_barcode = precision,
    coverage_vs_barcode  = coverage
  )
  
  if (!is.null(gwas_universe)) {
    gwas_universe <- unique(na.omit(gwas_universe))
    in_gwas <- intersect(gene_list, gwas_universe)
    out$n_in_gwas <- length(in_gwas)
    out$frac_in_gwas <- if (length(gene_list) == 0) NA_real_ else length(in_gwas) / length(gene_list)
  }
  
  out
}

# ---- Load ----
fpA <- fread(cfg$fp_A1)
fpB <- fread(cfg$fp_B1)

# Columns expected in footprint_merge_*:
# human_gene_any, gwas_pp4_max, gwas_tier, q_sc, chainSupport, scrna_tier, include_binary
# (If your column names differ, edit here.)
stopifnot(all(c("human_gene_any","gwas_pp4_max","include_binary") %in% names(fpA)))
stopifnot(all(c("human_gene_any","gwas_pp4_max","include_binary") %in% names(fpB)))

# scRNA top50 features
scrA <- fread(cfg$scr_A1_top50)
scrB <- fread(cfg$scr_B1_top50)
stopifnot("gene" %in% names(scrA))
stopifnot("gene" %in% names(scrB))

# optional GWAS strict reference
HAS_GWAS_REF <- file.exists(cfg$gwas_strict_ref)
if (HAS_GWAS_REF) {
  gref <- fread(cfg$gwas_strict_ref)
  if (!"gene" %in% names(gref)) HAS_GWAS_REF <- FALSE
}

# ---- Define: GWAS-backed universe (from merged footprints) ----
gwas_universe_A <- fpA[!is.na(gwas_pp4_max) & !is.na(human_gene_any),
                       .(gene = clean_gene(human_gene_any),
                         gwas_pp4_max = max(gwas_pp4_max, na.rm = TRUE),
                         gwas_tier = gwas_tier[which.max(gwas_pp4_max)])]
gwas_universe_A <- unique(gwas_universe_A, by = "gene")[order(-gwas_pp4_max)]
gwas_universe_A_genes <- gwas_universe_A$gene

gwas_universe_B <- fpB[!is.na(gwas_pp4_max) & !is.na(human_gene_any),
                       .(gene = clean_gene(human_gene_any),
                         gwas_pp4_max = max(gwas_pp4_max, na.rm = TRUE),
                         gwas_tier = gwas_tier[which.max(gwas_pp4_max)])]
gwas_universe_B <- unique(gwas_universe_B, by = "gene")[order(-gwas_pp4_max)]
gwas_universe_B_genes <- gwas_universe_B$gene

# sanity: they should be identical universes across A1/B1 in your current design
gwas_universe_genes <- sort(unique(c(gwas_universe_A_genes, gwas_universe_B_genes)))

# ---- Define: final binary barcode (from merged footprints) ----
barcodeA <- fpA[include_binary == TRUE & !is.na(human_gene_any), unique(clean_gene(human_gene_any))]
barcodeB <- fpB[include_binary == TRUE & !is.na(human_gene_any), unique(clean_gene(human_gene_any))]

# also export annotated barcode tables
barcode_tbl_A <- fpA[include_binary == TRUE & !is.na(human_gene_any),
                     .(gene = clean_gene(human_gene_any),
                       mouse_gene, gwas_traits, gwas_pp4_max, gwas_tier,
                       q_sc, scrna_tier, chainSupport, logFC, direction)]
barcode_tbl_A <- unique(barcode_tbl_A, by = "gene")[order(-gwas_pp4_max)]
barcode_tbl_B <- fpB[include_binary == TRUE & !is.na(human_gene_any),
                     .(gene = clean_gene(human_gene_any),
                       mouse_gene, gwas_traits, gwas_pp4_max, gwas_tier,
                       q_sc, scrna_tier, chainSupport, logFC, direction)]
barcode_tbl_B <- unique(barcode_tbl_B, by = "gene")[order(-gwas_pp4_max)]

# ---- Define: scRNA-only Top50 baselines ----
scrA_genes <- unique(clean_gene(scrA$gene))
scrB_genes <- unique(clean_gene(scrB$gene))

# ---- Export baseline lists ----
fwrite(gwas_universe_A, file.path(cfg$outdir, "baseline_A1_GWASonly_universe.tsv"), sep = "\t")
fwrite(gwas_universe_B, file.path(cfg$outdir, "baseline_B1_GWASonly_universe.tsv"), sep = "\t")
fwrite(data.table(gene = scrA_genes), file.path(cfg$outdir, "baseline_A1_scRNAonly_top50.tsv"), sep = "\t")
fwrite(data.table(gene = scrB_genes), file.path(cfg$outdir, "baseline_B1_scRNAonly_top50.tsv"), sep = "\t")
fwrite(barcode_tbl_A, file.path(cfg$outdir, "final_A1_binary_barcode.tsv"), sep = "\t")
fwrite(barcode_tbl_B, file.path(cfg$outdir, "final_B1_binary_barcode.tsv"), sep = "\t")

if (HAS_GWAS_REF) {
  fwrite(data.table(gene = unique(clean_gene(gref$gene))),
         file.path(cfg$outdir, "reference_GWAS_strict_40genes.tsv"), sep = "\t")
}

# ---- Scoreboard ----
rows <- list()

# A1
rows[["A1_GWASonly_universe"]] <- c(
  groupAB = "A1", method = "GWAS-only (universe from merge)",
  score_vs_barcode(gwas_universe_A_genes, barcodeA, gwas_universe_genes)
)

rows[["A1_scRNAonly_top50"]] <- c(
  groupAB = "A1", method = "scRNA-only (top50 q_sc)",
  score_vs_barcode(scrA_genes, barcodeA, gwas_universe_genes)
)

rows[["A1_final_barcode"]] <- c(
  groupAB = "A1", method = "Final binary barcode",
  score_vs_barcode(barcodeA, barcodeA, gwas_universe_genes)
)

# B1
rows[["B1_GWASonly_universe"]] <- c(
  groupAB = "B1", method = "GWAS-only (universe from merge)",
  score_vs_barcode(gwas_universe_B_genes, barcodeB, gwas_universe_genes)
)

rows[["B1_scRNAonly_top50"]] <- c(
  groupAB = "B1", method = "scRNA-only (top50 q_sc)",
  score_vs_barcode(scrB_genes, barcodeB, gwas_universe_genes)
)

rows[["B1_final_barcode"]] <- c(
  groupAB = "B1", method = "Final binary barcode",
  score_vs_barcode(barcodeB, barcodeB, gwas_universe_genes)
)

scoreboard <- rbindlist(lapply(rows, as.data.table), fill = TRUE)

# subtype-specificity within each method (A1 vs B1 Jaccard)
subtype_j <- data.table(
  method = c("GWAS-only (universe from merge)", "scRNA-only (top50 q_sc)", "Final binary barcode"),
  jaccard_A1_vs_B1 = c(
    jaccard(gwas_universe_A_genes, gwas_universe_B_genes),
    jaccard(scrA_genes, scrB_genes),
    jaccard(barcodeA, barcodeB)
  )
)

fwrite(scoreboard, file.path(cfg$outdir, "scoreboard.tsv"), sep = "\t")
fwrite(subtype_j, file.path(cfg$outdir, "subtype_specificity_jaccard.tsv"), sep = "\t")

cat("\n[OK] Wrote outputs to: ", cfg$outdir, "\n", sep = "")
cat("[Info] A1 barcode n=", length(barcodeA), " | B1 barcode n=", length(barcodeB), "\n", sep = "")
cat("[Info] GWAS universe n=", length(gwas_universe_genes), "\n", sep = "")





# =========================
# compare_discovery_vs_validation.R
# Compare barcode stability + PP4 replication
# Inputs:
#   footprint_merge_A1_full_discovery.tsv           (discovery)
#   footprint_merge_B1_full_discovery.tsv           (discovery)
#   footprint_merge_A1_full_validation.tsv          (validation)
#   footprint_merge_B1_full_validation.tsv          (validation)
# Outputs (output/validation_check/):
#   FigV1_barcode_stability.(png/pdf)
#   FigV2_PP4_replication_scatter.(png/pdf)
#   TableV0_barcode_diff.tsv
#   TableV1_barcode_PP4_disc_vs_val.tsv
#   TableV2_barcode_stability.tsv
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cfg <- list(
  disc_A1 = "footprint_merge_A1_full_discovery.tsv",
  disc_B1 = "footprint_merge_B1_full_discovery.tsv",
  val_A1  = "footprint_merge_A1_full_validation.tsv",
  val_B1  = "footprint_merge_B1_full_validation.tsv",
  outdir  = "output/validation_check",
  dpi     = 300
)
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(x))
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}
jaccard <- function(a, b) {
  a <- unique(na.omit(a)); b <- unique(na.omit(b))
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

read_merge <- function(path, group) {
  dt <- fread(path)
  need <- c("human_gene_any","gwas_pp4_max","include_binary")
  stopifnot(all(need %in% names(dt)))
  dt[, gene := clean_gene(human_gene_any)]
  dt[, gwas_pp4_max := as_num(gwas_pp4_max)]
  dt[, include_binary := as.logical(include_binary)]
  dt[, groupAB := group]
  dt[!is.na(gene)]
}

# ---- read ----
dA <- read_merge(cfg$disc_A1, "A1")
dB <- read_merge(cfg$disc_B1, "B1")
vA <- read_merge(cfg$val_A1,  "A1")
vB <- read_merge(cfg$val_B1,  "B1")

disc <- rbindlist(list(dA, dB), fill = TRUE)
val  <- rbindlist(list(vA, vB), fill = TRUE)

# ---- barcode sets (one row per groupAB; genes stored as list) ----
barcode_disc <- disc[include_binary == TRUE & !is.na(gene),
                     .(gene_disc = list(sort(unique(gene)))), by = groupAB]

barcode_val  <- val [include_binary == TRUE & !is.na(gene),
                     .(gene_val  = list(sort(unique(gene)))), by = groupAB]

stab <- merge(barcode_disc, barcode_val, by = "groupAB", all = TRUE)

stab[, `:=`(
  n_disc  = lengths(gene_disc),
  n_val   = lengths(gene_val),
  overlap = mapply(function(a,b) length(intersect(a,b)), gene_disc, gene_val),
  jaccard = mapply(jaccard, gene_disc, gene_val)
)]

stab_plot <- stab[, .(groupAB, n_disc, n_val, overlap, jaccard)][order(groupAB)]
fwrite(stab_plot, file.path(cfg$outdir, "TableV2_barcode_stability.tsv"), sep = "\t")

# ---- FigV1: barcode stability ----
p1 <- ggplot(stab_plot, aes(x = groupAB, y = jaccard)) +
  geom_col() +
  geom_text(aes(label = sprintf("J=%.2f\novl=%d", jaccard, overlap)),
            vjust = -0.3, size = 3) +
  ylim(0, 1.05) +
  labs(title = "Barcode stability: discovery-defined vs validation-derived",
       x = NULL, y = "Jaccard(barcode_disc, barcode_val)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face="bold", size=14))

ggsave(file.path(cfg$outdir, "FigV1_barcode_stability.pdf"), p1, width = 7, height = 3.8)
ggsave(file.path(cfg$outdir, "FigV1_barcode_stability.png"), p1, width = 7, height = 3.8, dpi = cfg$dpi)

# ---- per-gene PP4 replication (discovery barcode genes) ----
disc_bar_genes <- unique(disc[include_binary == TRUE & !is.na(gene), .(groupAB, gene)])

pp4_disc <- disc[, .(pp4_disc = max(gwas_pp4_max, na.rm = TRUE)), by = .(groupAB, gene)]
pp4_val  <- val [, .(pp4_val  = max(gwas_pp4_max, na.rm = TRUE)), by = .(groupAB, gene)]

# handle all-NA -> -Inf
pp4_disc[is.infinite(pp4_disc), pp4_disc := NA_real_]
pp4_val [is.infinite(pp4_val),  pp4_val  := NA_real_]

pp4_tbl <- merge(disc_bar_genes, pp4_disc, by = c("groupAB","gene"), all.x = TRUE)
pp4_tbl <- merge(pp4_tbl, pp4_val,  by = c("groupAB","gene"), all.x = TRUE)

fwrite(pp4_tbl[order(groupAB, -pp4_disc)],
       file.path(cfg$outdir, "TableV1_barcode_PP4_disc_vs_val.tsv"), sep = "\t")




# ---- FigV2: PP4 replication scatter (barcode union; missing -> 0; show zeros) ----
barcode_union <- merge(barcode_disc, barcode_val, by = "groupAB", all = TRUE)

barcode_union[, gene_union := mapply(function(a, b) {
  a <- if (is.null(a) || length(a) == 0) character() else a
  b <- if (is.null(b) || length(b) == 0) character() else b
  sort(unique(c(a, b)))
}, gene_disc, gene_val)]

barcode_union_long <- barcode_union[, .(gene = unlist(gene_union)), by = groupAB]
barcode_union_long <- unique(barcode_union_long[!is.na(gene)])

pp4_fig <- merge(barcode_union_long, pp4_disc, by = c("groupAB","gene"), all.x = TRUE)
pp4_fig <- merge(pp4_fig, pp4_val,  by = c("groupAB","gene"), all.x = TRUE)

# status label
pp4_fig[, status := fifelse(!is.na(pp4_disc) & !is.na(pp4_val), "Both",
                            fifelse(!is.na(pp4_disc) &  is.na(pp4_val), "Discovery-only",
                                    fifelse( is.na(pp4_disc) & !is.na(pp4_val), "Validation-only",
                                             "Missing")))]

# plot values: missing -> 0
pp4_fig[, pp4_disc_plot := ifelse(is.na(pp4_disc), 0, pp4_disc)]
pp4_fig[, pp4_val_plot  := ifelse(is.na(pp4_val),  0, pp4_val)]

pp4_fig[, groupAB := factor(groupAB, levels = c("A1","B1"))]
pp4_fig[, status := factor(status, levels = c("Both","Discovery-only","Validation-only","Missing"))]

# mark zero-border points
pp4_fig[, is_zero := (pp4_disc_plot == 0 | pp4_val_plot == 0)]

p2 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  # non-zero points
  geom_point(
    data = pp4_fig[is_zero == FALSE],
    aes(x = pp4_disc_plot, y = pp4_val_plot, shape = status),
    size = 2.8
  ) +
  # zero-border points: jitter so they don't sit on the panel border
  geom_point(
    data = pp4_fig[is_zero == TRUE],
    aes(x = pp4_disc_plot, y = pp4_val_plot, shape = status),
    position = position_jitter(width = 0.012, height = 0.012),
    size = 3.2
  ) +
  facet_wrap(~ groupAB, nrow = 1) +
  # give a small margin so x=0/y=0 points are fully visible
  scale_x_continuous(limits = c(-0.03, 1.03), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(-0.03, 1.03), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Replication of GWAS colocalization support for barcode genes",
    subtitle = "Barcode union (discovery U validation); missing PP4 plotted as 0",
    x = "PP4 (discovery)", y = "PP4 (validation)", shape = "Gene status"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# labels (prefer ggrepel)
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p2 <- p2 + ggrepel::geom_text_repel(
    data = pp4_fig,
    aes(x = pp4_disc_plot, y = pp4_val_plot, label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15
  )
} else {
  p2 <- p2 + geom_text(
    data = pp4_fig,
    aes(x = pp4_disc_plot, y = pp4_val_plot, label = gene),
    vjust = -0.6, size = 3
  )
}

ggsave(file.path(cfg$outdir, "FigV2_PP4_replication_scatter.pdf"), p2, width = 10, height = 4.6)
ggsave(file.path(cfg$outdir, "FigV2_PP4_replication_scatter.png"), p2, width = 10, height = 4.6, dpi = cfg$dpi)


# ---- identify which genes changed (A1/B1) ----
diff_list <- rbindlist(lapply(stab$groupAB, function(g) {
  a <- stab[groupAB == g]$gene_disc
  b <- stab[groupAB == g]$gene_val
  
  a <- if (length(a) == 0 || is.null(a[[1]])) character() else unique(na.omit(a[[1]]))
  b <- if (length(b) == 0 || is.null(b[[1]])) character() else unique(na.omit(b[[1]]))
  
  data.table(
    groupAB = g,
    disc_only = paste(sort(setdiff(a, b)), collapse = ","),
    val_only  = paste(sort(setdiff(b, a)), collapse = ","),
    shared    = paste(sort(intersect(a, b)), collapse = ",")
  )
}), fill = TRUE)

fwrite(diff_list, file.path(cfg$outdir, "TableV0_barcode_diff.tsv"), sep = "\t")

# ---- small console summary ----
cat("[OK] Outputs -> ", cfg$outdir, "\n", sep = "")
print(stab_plot)
print(diff_list)
if (nrow(miss_dt) > 0) {
  cat("[Warn] Missing PP4 for:\n")
  print(miss_dt)
}














# =========================
# plot_baselines_v1.R
# Visualize baseline vs final barcode
# Inputs:
#   baseline_A1_GWASonly_universe.tsv
#   baseline_A1_scRNAonly_top50.tsv
#   baseline_B1_GWASonly_universe.tsv
#   baseline_B1_scRNAonly_top50.tsv
#   final_A1_binary_barcode.tsv
#   final_B1_binary_barcode.tsv
#   scoreboard.tsv
#   subtype_specificity_jaccard.tsv
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------- Config ----------
cfg <- list(
  in_dir  = "output/baselines_v1",                       # folder where TSVs are
  out_dir = "output/baselines_figs_v1",  # figures output
  dpi     = 300
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

fpath <- function(x) file.path(cfg$in_dir, x)

# ---------- Load ----------
scoreboard <- fread(fpath("scoreboard.tsv"))
subtype_j  <- fread(fpath("subtype_specificity_jaccard.tsv"))

bA_gwas <- fread(fpath("baseline_A1_GWASonly_universe.tsv"))
bB_gwas <- fread(fpath("baseline_B1_GWASonly_universe.tsv"))
bA_scr  <- fread(fpath("baseline_A1_scRNAonly_top50.tsv"))
bB_scr  <- fread(fpath("baseline_B1_scRNAonly_top50.tsv"))

barA <- fread(fpath("final_A1_binary_barcode.tsv"))
barB <- fread(fpath("final_B1_binary_barcode.tsv"))

# ---------- Helpers ----------
as_num <- function(x) suppressWarnings(as.numeric(x))
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

# Make sure core columns are numeric
num_cols <- intersect(
  c("n_genes","n_overlap_barcode","precision_vs_barcode","coverage_vs_barcode",
    "n_in_gwas","frac_in_gwas"),
  names(scoreboard)
)
for (cc in num_cols) scoreboard[[cc]] <- as_num(scoreboard[[cc]])

# Order methods nicely
method_levels <- c(
  "GWAS-only (universe from merge)",
  "scRNA-only (top50 q_sc)",
  "Final binary barcode"
)
scoreboard[, method := factor(method, levels = method_levels)]
scoreboard[, groupAB := factor(groupAB, levels = c("A1","B1"))]

# ---------- FIG 1: Precision & Coverage vs final barcode ----------
dt1 <- melt(
  scoreboard,
  id.vars = c("groupAB","method","n_genes","n_overlap_barcode"),
  measure.vars = c("precision_vs_barcode","coverage_vs_barcode"),
  variable.name = "metric",
  value.name = "value"
)
dt1[, metric := factor(metric,
                       levels = c("precision_vs_barcode","coverage_vs_barcode"),
                       labels = c("Precision vs barcode","Coverage of barcode"))]

p1 <- ggplot(dt1, aes(x = method, y = value)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 3) +
  facet_grid(metric ~ groupAB, scales = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Baseline performance against final binary barcode",
       subtitle = "Precision = overlap / n_genes; Coverage = overlap / n_barcode") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$out_dir, "Fig1_scoreboard_precision_coverage.pdf"),
       p1, width = 9, height = 5)
ggsave(file.path(cfg$out_dir, "Fig1_scoreboard_precision_coverage.png"),
       p1, width = 9, height = 5, dpi = cfg$dpi)

# ---------- FIG 2: GWAS-anchored fraction ----------
# (scRNA-only is expected ~0; final barcode expected high)
dt2 <- scoreboard[!is.na(frac_in_gwas), .(groupAB, method, frac_in_gwas, n_in_gwas, n_genes)]
p2 <- ggplot(dt2, aes(x = method, y = frac_in_gwas)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d/%d", n_in_gwas, n_genes)), vjust = -0.3, size = 3) +
  facet_wrap(~ groupAB, nrow = 1) +
  labs(x = NULL, y = "Fraction in GWAS-backed universe",
       title = "GWAS anchoring of each method output",
       subtitle = "Label shows n_in_gwas / n_genes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$out_dir, "Fig2_gwas_anchored_fraction.pdf"),
       p2, width = 9, height = 3.8)
ggsave(file.path(cfg$out_dir, "Fig2_gwas_anchored_fraction.png"),
       p2, width = 9, height = 3.8, dpi = cfg$dpi)

# ---------- FIG 3 (REPLACE): Subtype specificity as lollipop ----------
subtype_j[, jaccard_A1_vs_B1 := as_num(jaccard_A1_vs_B1)]
subtype_j[, method := factor(method, levels = method_levels)]

p3 <- ggplot(subtype_j, aes(x = method, y = jaccard_A1_vs_B1)) +
  geom_segment(aes(xend = method, y = 0, yend = jaccard_A1_vs_B1), linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", jaccard_A1_vs_B1)), vjust = -0.6, size = 3) +
  ylim(-0.05, 1.05) +
  labs(x = NULL, y = "Jaccard(A1, B1)",
       title = "Subtype specificity (A1 vs B1) by method",
       subtitle = "Lower Jaccard = stronger subtype specificity") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(5.5, 5.5, 20, 5.5))

ggsave(file.path(cfg$out_dir, "Fig3_subtype_specificity_jaccard.pdf"),
       p3, width = 8.8, height = 3.8)
ggsave(file.path(cfg$out_dir, "Fig3_subtype_specificity_jaccard.png"),
       p3, width = 8.8, height = 3.8, dpi = cfg$dpi)

# ---------- FIG S1: Pairwise Jaccard heatmap among sets (per subtype) ----------
jaccard <- function(a, b) {
  a <- unique(na.omit(a)); b <- unique(na.omit(b))
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

make_setlist <- function(group) {
  if (group == "A1") {
    list(
      "GWAS-universe" = clean_gene(bA_gwas$gene),
      "scRNA-top50"   = clean_gene(bA_scr$gene),
      "barcode"       = clean_gene(barA$gene)
    )
  } else {
    list(
      "GWAS-universe" = clean_gene(bB_gwas$gene),
      "scRNA-top50"   = clean_gene(bB_scr$gene),
      "barcode"       = clean_gene(barB$gene)
    )
  }
}

calc_pairwise <- function(sets, group) {
  nms <- names(sets)
  out <- rbindlist(lapply(nms, function(i) {
    rbindlist(lapply(nms, function(j) {
      data.table(
        groupAB = group,
        set_i = i, set_j = j,
        jaccard = jaccard(sets[[i]], sets[[j]]),
        overlap = length(intersect(unique(na.omit(sets[[i]])), unique(na.omit(sets[[j]]))))
      )
    }))
  }))
  out
}

pw <- rbindlist(list(
  calc_pairwise(make_setlist("A1"), "A1"),
  calc_pairwise(make_setlist("B1"), "B1")
))

pw[, set_i := factor(set_i, levels = c("GWAS-universe","scRNA-top50","barcode"))]
pw[, set_j := factor(set_j, levels = c("GWAS-universe","scRNA-top50","barcode"))]
pw[, groupAB := factor(groupAB, levels = c("A1","B1"))]

pS1 <- ggplot(pw, aes(x = set_i, y = set_j, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f\n(n=%d)", jaccard, overlap)), size = 3) +
  facet_wrap(~ groupAB, nrow = 1) +
  labs(x = NULL, y = NULL,
       title = "Pairwise overlap among GWAS-universe, scRNA-top50, and final barcode",
       subtitle = "Cell label: Jaccard and overlap count") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$out_dir, "FigS1_pairwise_jaccard_heatmap.pdf"),
       pS1, width = 10, height = 4)
ggsave(file.path(cfg$out_dir, "FigS1_pairwise_jaccard_heatmap.png"),
       pS1, width = 10, height = 4, dpi = cfg$dpi)

# ---------- FIG S2 (REPLACE): Barcode PP4 lollipop (free y per facet) ----------
prep_barcode <- function(dt, group) {
  dt <- copy(dt)
  dt[, gene := clean_gene(gene)]
  dt[, gwas_pp4_max := as_num(gwas_pp4_max)]
  dt[, groupAB := group]
  dt <- dt[!is.na(gene) & !is.na(gwas_pp4_max)]
  dt[order(-gwas_pp4_max)]
}

barA2 <- prep_barcode(barA, "A1")
barB2 <- prep_barcode(barB, "B1")
barAll <- rbindlist(list(barA2, barB2), fill = TRUE)

# reorder genes within each facet
barAll[, gene_f := {
  gg <- gene
  gg[order(-gwas_pp4_max)]
}, by = groupAB]
barAll[, gene_f := factor(gene_f, levels = unique(gene_f)), by = groupAB]

pS2 <- ggplot(barAll, aes(x = gene_f, y = gwas_pp4_max)) +
  geom_segment(aes(xend = gene_f, y = 0, yend = gwas_pp4_max)) +
  geom_point(size = 2.2) +
  facet_wrap(~ groupAB, scales = "free_y", nrow = 1) +
  coord_flip() +
  labs(x = NULL, y = "GWAS PP4 (max)",
       title = "Final binary barcode genes ranked by GWAS colocalization support") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$out_dir, "FigS2_barcode_lollipop_PP4.pdf"),
       pS2, width = 9, height = 4.2)
ggsave(file.path(cfg$out_dir, "FigS2_barcode_lollipop_PP4.png"),
       pS2, width = 9, height = 4.2, dpi = cfg$dpi)




# =========================
# plot_scatter_pp4_qsc.R
# Scatter: x=GWAS PP4 (max), y=q_sc, shape=chainSupport
# Highlight: include_binary == TRUE
# Inputs:
#   footprint_merge_A1_full.tsv
#   footprint_merge_B1_full.tsv
# Outputs:
#   output/baseline_figs_v1/Fig4_scatter_PP4_vs_qsc_chainSupport.(png/pdf)
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cfg <- list(
  fp_A1 = "footprint_merge_A1_full.tsv",
  fp_B1 = "footprint_merge_B1_full.tsv",
  outdir = "output/baselines_figs_v1",
  dpi = 300
)
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(x))
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

# ---- load ----
A1 <- fread(cfg$fp_A1)
B1 <- fread(cfg$fp_B1)

# required columns (edit here if your colnames differ)
need <- c("human_gene_any","gwas_pp4_max","q_sc","chainSupport","include_binary")
stopifnot(all(need %in% names(A1)))
stopifnot(all(need %in% names(B1)))

A1[, groupAB := "A1"]
B1[, groupAB := "B1"]
dt <- rbindlist(list(A1, B1), fill = TRUE)

# ---- clean ----
dt[, gene := clean_gene(human_gene_any)]
dt[, gwas_pp4_max := as_num(gwas_pp4_max)]
dt[, q_sc := as_num(q_sc)]
dt[, chainSupport := as_num(chainSupport)]
dt[, include_binary := as.logical(include_binary)]

# keep GWAS-anchored candidates with scRNA score
dt <- dt[!is.na(gene) & !is.na(gwas_pp4_max) & !is.na(q_sc)]

# ---- shape binning for chainSupport (robust) ----
# If chainSupport is in [0,1], cut into Low/Mid/High; otherwise use quantiles
if (all(dt$chainSupport >= 0 & dt$chainSupport <= 1, na.rm = TRUE)) {
  dt[, chain_bin := cut(chainSupport,
                        breaks = c(-Inf, 1/3, 2/3, Inf),
                        labels = c("Low", "Mid", "High"),
                        right = TRUE)]
} else {
  qs <- quantile(dt$chainSupport, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  qs <- unique(as.numeric(qs))
  if (length(qs) < 4) {
    dt[, chain_bin := factor("All")]
  } else {
    dt[, chain_bin := cut(chainSupport,
                          breaks = qs,
                          include.lowest = TRUE,
                          labels = c("Low", "Mid", "High"))]
  }
}
dt[, chain_bin := as.factor(chain_bin)]
dt[, groupAB := factor(groupAB, levels = c("A1","B1"))]

# ---- plot ----
p <- ggplot(dt, aes(x = gwas_pp4_max, y = q_sc, shape = chain_bin)) +
  geom_point(size = 2, alpha = 0.85) +
  # highlight barcode genes
  geom_point(data = dt[include_binary == TRUE],
             size = 3.2, alpha = 1) +
  facet_wrap(~ groupAB, nrow = 1) +
  labs(
    x = "GWAS PP4 (max)",
    y = "scRNA support (q_sc)",
    shape = "chainSupport",
    title = "Genetic anchoring vs scRNA support",
    subtitle = "Points = GWAS-anchored candidates; larger points = final binary barcode"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# Optional: label barcode genes if ggrepel installed
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(
    data = dt[include_binary == TRUE],
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  )
} else {
  # fallback label without ggrepel (may overlap)
  p <- p + geom_text(
    data = dt[include_binary == TRUE],
    aes(label = gene),
    size = 3,
    vjust = -0.6
  )
}

ggsave(file.path(cfg$outdir, "Fig4_scatter_PP4_vs_qsc_chainSupport.pdf"),
       p, width = 10, height = 4.2)
ggsave(file.path(cfg$outdir, "Fig4_scatter_PP4_vs_qsc_chainSupport.png"),
       p, width = 10, height = 4.2, dpi = cfg$dpi)

cat("[OK] Wrote: ", file.path(cfg$outdir, "Fig4_scatter_PP4_vs_qsc_chainSupport.(pdf/png)"), "\n", sep = "")





