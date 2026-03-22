
getwd()
.libPaths(c("E:/R-4.4.0/librarygwas", setdiff(.libPaths(), "E:/R-4.4.0/librarygwas")))

suppressPackageStartupMessages({
  library(data.table)
})

root <- "output/coloc_v1"
outd <- "output/gwas_validation_result"
dir.create(outd, recursive = TRUE, showWarnings = FALSE)

dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("__vs__", basename(dirs))]

parse_nm <- function(nm){
  # <TAG>__vs__<QTL>__(strict|relaxed)
  m <- regexec("^(.+?)__vs__(.+?)__(strict|relaxed)$", nm)
  r <- regmatches(nm, m)[[1]]
  if (length(r) == 0) return(NULL)
  data.table(tag = r[2], qtl = r[3], suffix = r[4])
}

read_one <- function(d){
  nm <- basename(d)
  meta <- parse_nm(nm)
  if (is.null(meta)) return(NULL)
  f <- file.path(d, "coloc_results.tsv")
  if (!file.exists(f)) return(NULL)
  
  dt <- fread(f)
  if (!all(c("chr","start","end","PP4","best_gene") %in% names(dt))) return(NULL)
  
  dt[, locus_key := paste(chr, start, end, sep=":")]
  
  tag <- meta$tag
  dt[, `:=`(
    tag    = tag,
    qtl    = meta$qtl,
    suffix = meta$suffix,
    role   = ifelse(grepl("_validation$", tag), "validation", "discovery"),
    anchor = ifelse(grepl("_validation$", tag), sub("_validation$", "_discovery", tag), tag)
  )]
  dt
}

all <- rbindlist(lapply(dirs, read_one), fill = TRUE)
all <- all[suffix == "strict"]

# 只保留 discovery+validation 都存在的组合
pair_ok <- all[, .N, by=.(anchor, qtl, suffix, role)][, .N, by=.(anchor, qtl, suffix)][N==2]
all <- all[pair_ok, on=.(anchor, qtl, suffix)]

# 宽表：PP4
pp4_w <- dcast(all, anchor + qtl + suffix + locus_key ~ role, value.var = "PP4")
setnames(pp4_w, c("discovery","validation"), c("PP4_discovery","PP4_validation"))

# 宽表：best_gene
bg_w <- dcast(all, anchor + qtl + suffix + locus_key ~ role, value.var = "best_gene")
setnames(bg_w, c("discovery","validation"), c("best_gene_discovery","best_gene_validation"))

cmp <- merge(pp4_w, bg_w, by=c("anchor","qtl","suffix","locus_key"), all=TRUE)

cmp[, gene_match := best_gene_discovery == best_gene_validation]
cmp[, `:=`(
  disc_red    = PP4_discovery >= 0.8,
  disc_orange = PP4_discovery >= 0.5,
  val_red     = PP4_validation >= 0.8,
  val_orange  = PP4_validation >= 0.5
)]

# 你项目里的“红/橙/黄/无复现”口径（可自行调整阈值）
cmp[, tier := fifelse(disc_red & val_red & gene_match, "RED_rep",
                      fifelse(disc_red & val_orange & gene_match, "ORANGE_rep",
                              fifelse(disc_orange & val_orange & gene_match, "YELLOW_rep", "NO_rep")))]

fwrite(cmp[order(anchor, qtl, suffix, -PP4_discovery)], file.path(outd, "PP4_compare_all.tsv"), sep="\t")

summ <- cmp[, .(
  n_locus       = .N,
  n_disc_red    = sum(disc_red, na.rm=TRUE),
  n_val_red     = sum(val_red, na.rm=TRUE),
  n_rep_orange  = sum(tier %in% c("RED_rep","ORANGE_rep"), na.rm=TRUE),
  n_rep_red     = sum(tier == "RED_rep", na.rm=TRUE),
  gene_match_rate = mean(gene_match & disc_red & val_orange, na.rm=TRUE)
), by=.(anchor, qtl, suffix)]

fwrite(summ[order(anchor, qtl, suffix)], file.path(outd, "PP4_summary.tsv"), sep="\t")

stable <- cmp[tier %in% c("RED_rep","ORANGE_rep"),
              .(anchor, qtl, suffix, locus_key,
                best_gene = best_gene_discovery,
                PP4_discovery, PP4_validation)]
fwrite(stable[order(anchor, qtl, suffix, -PP4_validation)], file.path(outd, "stable_loci.tsv"), sep="\t")

stable_genes <- unique(stable[, .(anchor, best_gene)])
fwrite(stable_genes[order(anchor, best_gene)], file.path(outd, "stable_genes.tsv"), sep="\t")

cat("[OK] Wrote files to:", outd, "\n")