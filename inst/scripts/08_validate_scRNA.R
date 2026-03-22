suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

# ---- 0) 必须确保这些列存在：sample_id / group4 / study ----
md <- pbmc@meta.data

# sample_id：优先用 library；没有就用 orig.ident
if (!("sample_id" %in% colnames(md))) {
  if ("library" %in% colnames(md)) md$sample_id <- as.character(md$library)
  else if ("orig.ident" %in% colnames(md)) md$sample_id <- as.character(md$orig.ident)
  else stop("No sample_id/library/orig.ident in meta.data")
}

# group4：A1/A2/B1/B2（从 sample_id 前缀推断）
if (!("group4" %in% colnames(md))) {
  md$group4 <- NA_character_
  md$group4[grepl("^A1_", md$sample_id)] <- "A1"
  md$group4[grepl("^A2_", md$sample_id)] <- "A2"
  md$group4[grepl("^B1_", md$sample_id)] <- "B1"
  md$group4[grepl("^B2_", md$sample_id)] <- "B2"
}

# study：按 group4 推断
if (!("study" %in% colnames(md))) {
  md$study <- NA_character_
  md$study[md$group4 %in% c("A1","A2")] <- "GSE204803"
  md$study[md$group4 %in% c("B1","B2")] <- "GSE263191"
}

pbmc@meta.data <- md

# ---- 1) 只取 vascular 细胞 ----
stopifnot("vascular_gate" %in% colnames(pbmc@meta.data))
pbmc_v <- subset(pbmc, subset = vascular_gate == TRUE)

# ---- 2) barcode gene list（用你冻结的那版；B1建议用 pooled union=10）----
A1_genes <- c("Sh3pxd2a","Gnb2","Sfr1","Mrpl38")
B1_genes <- c("Zcchc14","Nmt1","Nbeal1","Celf1","Gigyf1","Banp","Slk","Sptbn1","Tab2","Plekhg1")

A1_genes <- intersect(A1_genes, rownames(pbmc_v))
B1_genes <- intersect(B1_genes, rownames(pbmc_v))

# 0) 确保用的是 RNA（或你想用的 assay）
assay_use <- if ("RNA" %in% Assays(pbmc_v)) "RNA" else DefaultAssay(pbmc_v)
DefaultAssay(pbmc_v) <- assay_use

# 1) 把 v5 多 layers 合并（关键一步）
pbmc_v <- SeuratObject::JoinLayers(pbmc_v, assay = assay_use)

# 2) 再跑 module score
pbmc_v <- AddModuleScore(pbmc_v, features = list(A1_genes), name = "A1_score_")
pbmc_v <- AddModuleScore(pbmc_v, features = list(B1_genes), name = "B1_score_")

dt <- as.data.table(pbmc_v@meta.data)

# ---- 4) sample-level 汇总（避免把细胞当独立样本）----
score_by_sample <- dt[, .(
  n_cells_vascular = .N,
  A1_score = mean(A1_score_1, na.rm=TRUE),
  B1_score = mean(B1_score_1, na.rm=TRUE)
), by=.(study, group4, sample_id)]

# ---- 5) 两个 study 内部对照 + 交叉负对照 ----
do_test <- function(d, case, ctrl, col){
  x <- d[group4==case, get(col)]
  y <- d[group4==ctrl, get(col)]
  if (length(x) < 2 || length(y) < 2) return(list(p=NA_real_, d=NA_real_))
  p <- tryCatch(wilcox.test(x, y)$p.value, error=function(e) NA_real_)
  d_eff <- (mean(x)-mean(y)) / sqrt(((sd(x)^2)+(sd(y)^2))/2)
  list(p=p, d=d_eff)
}

A <- score_by_sample[study=="GSE204803"]
B <- score_by_sample[study=="GSE263191"]

res <- rbindlist(list(
  data.table(test="GSE204803: A1_score A1vsA2", do_test(A,"A1","A2","A1_score")),
  data.table(test="GSE204803: B1_score A1vsA2 (neg)", do_test(A,"A1","A2","B1_score")),
  data.table(test="GSE263191: B1_score B1vsB2", do_test(B,"B1","B2","B1_score")),
  data.table(test="GSE263191: A1_score B1vsB2 (neg)", do_test(B,"B1","B2","A1_score"))
), fill=TRUE)

print(res)
print(score_by_sample)

# ---- 6) 简单作图（可直接进 Supplement）----
plot_one <- function(d, y, title){
  ggplot(d, aes(x=group4, y=.data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.12, height=0, size=2) +
    theme_bw() + labs(title=title, x=NULL, y=y)
}
print(plot_one(A, "A1_score", "GSE204803 vascular: A1_score (A1 vs A2)"))
print(plot_one(A, "B1_score", "GSE204803 vascular: B1_score (neg control)"))
print(plot_one(B, "B1_score", "GSE263191 vascular: B1_score (B1 vs B2)"))
print(plot_one(B, "A1_score", "GSE263191 vascular: A1_score (neg control)"))


library(data.table)

out_dir <- "output/scRNA_validation_v1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 如果你有 score_by_sample（每个 sample 一行）
if (exists("score_by_sample")) {
  fwrite(score_by_sample, file.path(out_dir, "barcode_score_by_sample.tsv"), sep = "\t")
}

# 如果你有 res（那张 test/p/cohens_d 的表）
if (exists("res")) {
  fwrite(res, file.path(out_dir, "barcode_stats.tsv"), sep = "\t")
}

cat("[OK] Saved to:", normalizePath(out_dir), "\n")
