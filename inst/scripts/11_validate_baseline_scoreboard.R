# =========================
# baseline_build_v1.R
# Build manuscript baselines from existing outputs:
#   GWAS-only:  output/summary_v1/feature_gwas_gene_vascular_strict_v2.tsv
#   scRNA-only: output/scRNA_feature_v2/feature_*_qsc_top50.tsv
#   coloc-only: output/coloc_v1/*__strict/coloc_gene_results.tsv
#   MR-only:    output/coloc_v1/*__strict/mr_results.tsv (may be empty)
#   cSVI barcode: output/footprint_v1/{discovery,validation}_footprint_*_binary.tsv
# + naive intersection/union and stability (discovery vs validation)
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cfg <- list(
  root_output = "output",
  outdir = "output/baseline_bench_v1",
  
  # scRNA q_sc Top50 files (A1/B1)
  scrna_top50 = list(
    A1 = "output/scRNA_feature_v2/feature_A1_qsc_top50.tsv",
    B1 = "output/scRNA_feature_v2/feature_B1_qsc_top50.tsv"
  ),
  
  # cSVI barcode (binary footprint) files
  csvi_barcode_files = list(
    A1 = list(
      discovery  = "output/footprint_v1/discovery_footprint_A1_binary.tsv",
      validation = "output/footprint_v1/validation_footprint_A1_binary.tsv"
    ),
    B1 = list(
      discovery  = "output/footprint_v1/discovery_footprint_B1_binary.tsv",
      validation = "output/footprint_v1/validation_footprint_B1_binary.tsv"
    )
  ),
  
  # GWAS gene feature summary (v2 includes discovery + validation)
  gwas_gene_file = "output/summary_v1/feature_gwas_gene_vascular_strict_v2.tsv",
  
  # coloc folder root
  coloc_root = "output/coloc_v1",
  
  # tags we care about
  tags = c("WMH_discovery","WMH_validation","LAC_discovery","LAC_validation"),
  
  # coloc threshold
  pp4_thr = 0.8,
  
  # MR p threshold (keep strict; empty is OK)
  mr_p_thr = 0.05,
  
  # plot helper
  show_count_label = TRUE
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
guess_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

as_genes <- function(x) unique(na.omit(trimws(as.character(x))))

# ---------- read gene list ----------
read_gene_list <- function(file) {
  if (!file.exists(file)) stop("Missing file: ", file)
  dt <- fread(file)
  gene_col <- guess_col(dt, c(
    "human_gene_any","human_gene","gene","Gene","SYMBOL","symbol",
    "gene_symbol","hgnc_symbol","best_gene","best_gene_symbol"
  ))
  if (is.na(gene_col)) stop("No gene-like column in: ", file, "\nCols=", paste(names(dt), collapse=","))
  as_genes(dt[[gene_col]])
}

# scRNA-only gene sets
scrna_sets <- lapply(cfg$scrna_top50, read_gene_list)

# cSVI barcode gene sets (discovery/validation)
barcode_sets <- list(
  A1 = list(
    discovery  = read_gene_list(cfg$csvi_barcode_files$A1$discovery),
    validation = read_gene_list(cfg$csvi_barcode_files$A1$validation)
  ),
  B1 = list(
    discovery  = read_gene_list(cfg$csvi_barcode_files$B1$discovery),
    validation = read_gene_list(cfg$csvi_barcode_files$B1$validation)
  )
)

# ---------- GWAS-only gene sets from feature_gwas_gene_vascular_strict_v2.tsv ----------
read_gwas_sets <- function(file, tags) {
  if (!file.exists(file)) stop("Missing file: ", file)
  dt <- fread(file)
  
  # guess tag column
  tag_col <- guess_col(dt, c("tag","gwas_tag","GWAS","Trait","trait","gwas"))
  if (is.na(tag_col)) {
    stop("Cannot find tag column in GWAS gene feature file.\nCols=", paste(names(dt), collapse=","))
  }
  
  # guess gene column
  gene_col <- guess_col(dt, c("human_gene_any","human_gene","gene","Gene","SYMBOL","symbol","best_gene","best_gene_symbol"))
  if (is.na(gene_col)) {
    stop("Cannot find gene column in GWAS gene feature file.\nCols=", paste(names(dt), collapse=","))
  }
  
  out <- list()
  for (tg in tags) {
    out[[tg]] <- as_genes(dt[get(tag_col) == tg, get(gene_col)])
  }
  out
}

gwas_sets <- read_gwas_sets(cfg$gwas_gene_file, cfg$tags)

# ---------- coloc-only / MR-only: scan strict folders ----------
parse_folder <- function(folder_name) {
  # <TAG>__vs__<QTL>__strict
  m <- regexec("^(.+?)__vs__(.+?)__strict$", folder_name)
  r <- regmatches(folder_name, m)[[1]]
  if (length(r) == 0) return(NULL)
  list(tag = r[2], qtl = r[3])
}

scan_strict_runs <- function(root, tags) {
  dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  runs <- rbindlist(lapply(dirs, function(d) {
    nm <- basename(d)
    info <- parse_folder(nm)
    if (is.null(info)) return(NULL)
    if (!(info$tag %in% tags)) return(NULL)
    data.table(
      dir = d,
      folder = nm,
      tag = info$tag,
      qtl = info$qtl,
      coloc_gene = file.path(d, "coloc_gene_results.tsv"),
      mr = file.path(d, "mr_results.tsv")
    )
  }), fill = TRUE)
  runs
}

runs <- scan_strict_runs(cfg$coloc_root, cfg$tags)
fwrite(runs, file.path(cfg$outdir, "index_strict_runs.tsv"), sep="\t")

read_coloc_genes <- function(file, pp4_thr) {
  if (!file.exists(file)) return(character())
  dt <- fread(file)
  gene_col <- guess_col(dt, c("best_gene","best_gene_symbol","human_gene_any","human_gene","gene","Gene","SYMBOL","symbol"))
  pp4_col  <- guess_col(dt, c("PP4","pp4","pp4_abf","PP4.abf","pp4_max","gwas_pp4_max"))
  if (is.na(gene_col) || is.na(pp4_col)) {
    stop("coloc_gene_results missing gene/PP4 cols: ", file, "\nCols=", paste(names(dt), collapse=","))
  }
  as_genes(dt[!is.na(get(pp4_col)) & get(pp4_col) >= pp4_thr, get(gene_col)])
}

read_mr_genes <- function(file, p_thr) {
  if (!file.exists(file)) return(character())
  dt <- tryCatch(fread(file), error=function(e) data.table())
  if (nrow(dt) == 0) return(character())
  gene_col <- guess_col(dt, c("exposure","gene","Gene","SYMBOL","symbol","best_gene","best_gene_symbol","human_gene"))
  p_col    <- guess_col(dt, c("pval","p","p_value","pval_ivw","pval.wald","pval_wald_ratio","pval_IVW"))
  if (is.na(gene_col) || is.na(p_col)) {
    stop("mr_results missing gene/p cols: ", file, "\nCols=", paste(names(dt), collapse=","))
  }
  as_genes(dt[!is.na(get(p_col)) & get(p_col) < p_thr, get(gene_col)])
}

# build coloc baselines:
# 1) blood-only 2) best-tissue (upper-bound proxy by union after threshold)
build_coloc_sets <- function(runs, pp4_thr) {
  out <- list()
  for (tg in unique(runs$tag)) {
    sub <- runs[tag == tg]
    
    sub_blood <- sub[grepl("eQTLGen_blood", qtl, fixed = TRUE)]
    genes_blood <- unique(unlist(lapply(sub_blood$coloc_gene, read_coloc_genes, pp4_thr=pp4_thr)))
    out[[paste0(tg, "__coloc_eqtlgen_blood")]] <- genes_blood
    
    genes_best <- unique(unlist(lapply(sub$coloc_gene, read_coloc_genes, pp4_thr=pp4_thr)))
    out[[paste0(tg, "__coloc_best_tissue")]] <- genes_best
  }
  out
}

build_mr_sets <- function(runs, p_thr) {
  out <- list()
  for (tg in unique(runs$tag)) {
    sub <- runs[tag == tg]
    
    sub_blood <- sub[grepl("eQTLGen_blood", qtl, fixed = TRUE)]
    genes_blood <- unique(unlist(lapply(sub_blood$mr, read_mr_genes, p_thr=p_thr)))
    out[[paste0(tg, "__mr_eqtlgen_blood")]] <- genes_blood
    
    genes_best <- unique(unlist(lapply(sub$mr, read_mr_genes, p_thr=p_thr)))
    out[[paste0(tg, "__mr_best_tissue")]] <- genes_best
  }
  out
}

coloc_sets <- build_coloc_sets(runs, cfg$pp4_thr)
mr_sets    <- build_mr_sets(runs, cfg$mr_p_thr)

# ---------- Jaccard ----------
jaccard2 <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a)==0 && length(b)==0) return(NA_real_)  # empty-empty => NA
  if (length(a)==0 || length(b)==0) return(0)
  length(intersect(a,b)) / length(union(a,b))
}

# derive WMH/LAC pairs
pairs <- data.table(
  trait = c("WMH","LAC"),
  disc = c("WMH_discovery","LAC_discovery"),
  val  = c("WMH_validation","LAC_validation")
)

methods <- c(
  "csvi_barcode",
  "gwas_only","scrna_only",
  "coloc_eqtlgen_blood","coloc_best_tissue",
  "mr_eqtlgen_blood","mr_best_tissue",
  "naive_intersection_gwas_scrna","naive_union_gwas_scrna"
)

collect_set <- function(groupAB, tag, method) {
  scrna <- scrna_sets[[groupAB]]
  gwas  <- gwas_sets[[tag]]
  
  if (method == "csvi_barcode") {
    if (grepl("_discovery$", tag)) return(barcode_sets[[groupAB]]$discovery)
    if (grepl("_validation$", tag)) return(barcode_sets[[groupAB]]$validation)
    stop("Unknown tag split for csvi_barcode: ", tag)
  }
  
  if (method == "gwas_only") return(gwas)
  if (method == "scrna_only") return(scrna)
  
  if (method == "coloc_eqtlgen_blood") return(coloc_sets[[paste0(tag,"__coloc_eqtlgen_blood")]])
  if (method == "coloc_best_tissue")   return(coloc_sets[[paste0(tag,"__coloc_best_tissue")]])
  
  if (method == "mr_eqtlgen_blood") return(mr_sets[[paste0(tag,"__mr_eqtlgen_blood")]])
  if (method == "mr_best_tissue")   return(mr_sets[[paste0(tag,"__mr_best_tissue")]])
  
  if (method == "naive_intersection_gwas_scrna") return(intersect(gwas, scrna))
  if (method == "naive_union_gwas_scrna")        return(union(gwas, scrna))
  
  stop("Unknown method: ", method)
}

# summary: size + stability (disc vs val) for each trait and groupAB
sum_dt <- rbindlist(lapply(c("A1","B1"), function(g) {
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    tr <- pairs$trait[i]; disc <- pairs$disc[i]; val <- pairs$val[i]
    
    rbindlist(lapply(methods, function(m) {
      S1 <- collect_set(g, disc, m)
      S2 <- collect_set(g, val,  m)
      data.table(
        groupAB = g,
        trait = tr,
        method = m,
        n_discovery = length(unique(S1)),
        n_validation = length(unique(S2)),
        jaccard = jaccard2(S1, S2)
      )
    }))
  }))
}), fill=TRUE)

# status flag
sum_dt[, status := fifelse(n_discovery==0 & n_validation==0, "empty_both",
                           fifelse(n_discovery>0 & n_validation==0, "empty_validation",
                                   fifelse(n_discovery==0 & n_validation>0, "empty_discovery", "ok")))]
sum_dt[, n_label := paste0(n_discovery, "/", n_validation)]

fwrite(sum_dt, file.path(cfg$outdir, "baseline_summary.tsv"), sep="\t")

# plot stability (drop NA bars)
sum_plot <- sum_dt[!is.na(jaccard)]
sum_plot[, method := factor(method, levels=methods)]

p <- ggplot(sum_plot, aes(x=method, y=jaccard)) +
  geom_col() +
  facet_grid(groupAB ~ trait) +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle=30, hjust=1)) +
  labs(title="Baseline stability: discovery vs validation (Jaccard)", x=NULL, y="Jaccard")

if (isTRUE(cfg$show_count_label)) {
  p <- p + geom_text(aes(label=n_label), vjust=-0.2, size=3)
}

ggsave(file.path(cfg$outdir, "Fig_baseline_stability.png"), p, width=13, height=6, dpi=300)

message("[OK] Done. Outputs -> ", cfg$outdir)

# quick sanity prints
cat("\n[Sanity] cSVI barcode sizes:\n")
cat("A1 disc/val:", length(barcode_sets$A1$discovery), "/", length(barcode_sets$A1$validation), "\n")
cat("B1 disc/val:", length(barcode_sets$B1$discovery), "/", length(barcode_sets$B1$validation), "\n")
cat("[Sanity] barcode Jaccard:\n")
cat("A1:", jaccard2(barcode_sets$A1$discovery, barcode_sets$A1$validation), "\n")
cat("B1:", jaccard2(barcode_sets$B1$discovery, barcode_sets$B1$validation), "\n")







# =========================
# build_feature_gwas_gene_vascular_strict_v2.R
# Rebuild GWAS vascular features from strict coloc runs (discovery + validation)
# Input:  output/coloc_v1/*__vs__*__strict/coloc_gene_results.tsv
# Output: output/summary_v1/feature_gwas_gene_vascular_strict_v2.tsv
#         output/summary_v1/gwas_only_genes_<TAG>.tsv
# =========================

suppressPackageStartupMessages({
  library(data.table)
})

cfg <- list(
  coloc_root = "output/coloc_v1",
  outdir     = "output/summary_v1",
  tags       = c("WMH_discovery","WMH_validation","LAC_discovery","LAC_validation"),
  # vascular QTL backgrounds to include
  keep_qtl_regex = "^(eQTLGen_blood|GTExv7_Artery_)",
  pp4_keep  = 0.5
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

parse_folder <- function(folder_name) {
  # <TAG>__vs__<QTL>__strict
  m <- regexec("^(.+?)__vs__(.+?)__strict$", folder_name)
  r <- regmatches(folder_name, m)[[1]]
  if (length(r) == 0) return(NULL)
  list(tag = r[2], qtl = r[3])
}

guess_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

read_coloc_gene_tbl <- function(file, tag, qtl) {
  if (!file.exists(file)) return(NULL)
  dt <- fread(file)
  
  gene_col <- guess_col(dt, c(
    "best_gene","best_gene_symbol","human_gene_any","human_gene",
    "gene","Gene","SYMBOL","symbol","hgnc_symbol"
  ))
  pp4_col <- guess_col(dt, c(
    "PP4","pp4","pp4_abf","PP4.abf","pp4_max","gwas_pp4_max"
  ))
  
  if (is.na(gene_col) || is.na(pp4_col)) {
    stop("Missing gene/PP4 columns in: ", file, "\nCols=", paste(names(dt), collapse=","))
  }
  
  # Optional: if locus id exists, keep best per locus (more conservative)
  locus_col <- guess_col(dt, c(
    "locus_id","locus","region","window","lead_snp","top_snp","rsid","snp"
  ))
  
  dt <- dt[!is.na(get(pp4_col))]
  if (!is.na(locus_col)) {
    setorder(dt, -get(pp4_col))
    dt <- dt[, .SD[1], by = locus_col]
  }
  
  dt <- dt[get(pp4_col) >= cfg$pp4_keep]
  
  out <- data.table(
    tag  = tag,
    qtl  = qtl,
    gene = trimws(as.character(dt[[gene_col]])),
    pp4  = as.numeric(dt[[pp4_col]])
  )
  out <- out[!is.na(gene) & gene != ""]
  out
}

# scan strict dirs
dirs <- list.dirs(cfg$coloc_root, full.names = TRUE, recursive = FALSE)

rows <- rbindlist(lapply(dirs, function(d) {
  nm <- basename(d)
  info <- parse_folder(nm)
  if (is.null(info)) return(NULL)
  if (!(info$tag %in% cfg$tags)) return(NULL)
  if (!grepl(cfg$keep_qtl_regex, info$qtl)) return(NULL)
  
  f <- file.path(d, "coloc_gene_results.tsv")
  read_coloc_gene_tbl(f, info$tag, info$qtl)
}), fill = TRUE)

if (nrow(rows) == 0) stop("No rows collected. Check folder names and keep_qtl_regex.")

# For each (tag, gene) keep max PP4 across vascular QTL backgrounds (upper bound within vascular)
feat <- rows[, .(pp4_max = max(pp4, na.rm=TRUE),
                 qtl_best = qtl[which.max(pp4)]),
             by = .(tag, gene)]
setorder(feat, tag, -pp4_max, gene)

out_file <- file.path(cfg$outdir, "feature_gwas_gene_vascular_strict_v2.tsv")
fwrite(feat, out_file, sep = "\t")

# export per-tag gene sets (GWAS-only baseline in your current definition)
for (tg in cfg$tags) {
  genes <- unique(feat[tag == tg]$gene)
  fwrite(data.table(gene = genes),
         file.path(cfg$outdir, paste0("gwas_only_genes_", tg, ".tsv")),
         sep = "\t")
}

cat("[OK] Written:\n  ", out_file, "\n")
cat("[OK] Tag counts:\n")
print(feat[, .N, by=tag][order(tag)])






suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

infile <- "output/baseline_bench_v1/baseline_summary.tsv"
dt <- fread(infile)

# 只画有定义的 Jaccard
dt <- dt[!is.na(jaccard)]

# 用 discovery/validation 的平均集合大小做 x 轴（你也可以改成 pmax）
dt[, n_mean := (n_discovery + n_validation)/2]
dt[, n_max  := pmax(n_discovery, n_validation)]
dt[, n_min  := pmin(n_discovery, n_validation)]

p <- ggplot(dt, aes(x = n_mean, y = jaccard, label = method)) +
  geom_point() +
  geom_text(vjust = -0.6, size = 3, check_overlap = TRUE) +
  facet_grid(groupAB ~ trait) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(
    title = "Candidate set size vs stability (Jaccard)",
    x = "Candidate set size (mean of discovery/validation, log10)",
    y = "Jaccard (discovery vs validation)"
  )

ggsave("output/baseline_bench_v1/Fig_size_vs_jaccard.png", p, width = 13, height = 6, dpi = 300)