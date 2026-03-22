# =========================
# ortholog_v2.R (MINIMAL OUTPUT)
# Only outputs 3 files:
#  1) master_coloc_strict_annotated.tsv
#  2) gene_map_human_mouse_one2one.tsv
#  3) feature_gwas_gene_vascular_strict.tsv
# =========================
getwd()
.libPaths(c("E:/R-4.4.0/librarygwas", setdiff(.libPaths(), "E:/R-4.4.0/librarygwas")))
suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

in_root <- "output/coloc_v1"
out_dir <- "output/summary_v1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[Info] WD:     ", getwd())
message("[Info] InRoot: ", normalizePath(in_root, winslash = "/", mustWork = FALSE))
message("[Info] OutDir: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

if (!dir.exists(in_root)) {
  stop("[Error] Not found: ", in_root,
       "\nMake sure you run this script from the project root where output/coloc_v1 exists.")
}

# -------- helpers --------
parse_from_folder <- function(folder) {
  mode <- NA_character_
  if (grepl("__strict$", folder))  mode <- "strict"
  if (grepl("__relaxed$", folder)) mode <- "relaxed"
  
  base <- sub("__(strict|relaxed)$", "", folder)
  
  trait <- NA_character_
  qtl   <- NA_character_
  if (grepl("__vs__", base, fixed = TRUE)) {
    trait <- sub("__vs__.*$", "", base)
    qtl   <- sub("^.*__vs__", "", base)
  }
  
  qtl_source <- NA_character_
  tissue <- NA_character_
  if (!is.na(qtl) && qtl != "") {
    if (startsWith(qtl, "GTExv7_")) {
      qtl_source <- "GTExv7"
      tissue <- sub("^GTExv7_", "", qtl)
    } else if (startsWith(qtl, "eQTLGen_")) {
      qtl_source <- "eQTLGen"
      tissue <- sub("^eQTLGen_", "", qtl)
    } else {
      qtl_source <- sub("_.*$", "", qtl)
      tissue <- sub("^[^_]+_", "", qtl)
      if (!is.na(tissue) && tissue == qtl) tissue <- NA_character_
    }
  }
  
  list(trait=trait, qtl_source=qtl_source, tissue=tissue, mode=mode)
}

safe_fread <- function(path) {
  tryCatch(
    fread(path, sep="\t", data.table=TRUE, showProgress=FALSE),
    error = function(e) {
      message("[Warn] fread failed: ", path, " | ", e$message)
      NULL
    }
  )
}

# =========================
# 1) Read coloc_results.tsv -> master_coloc
# =========================
coloc_files <- list.files(in_root, pattern="coloc_results\\.tsv$", recursive=TRUE, full.names=TRUE)
message("[Info] Found coloc files: ", length(coloc_files))
if (length(coloc_files) == 0) stop("[Error] No coloc_results.tsv under ", in_root)

read_coloc_one <- function(f) {
  folder <- basename(dirname(f))
  meta <- parse_from_folder(folder)
  if (is.na(meta$mode)) return(NULL)
  
  dt <- safe_fread(f)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  if (!("best_gene" %in% names(dt))) {
    message("[Warn] Skip(no best_gene col): ", f)
    return(NULL)
  }
  
  need_cols <- c("locus","PP4","top_snp_h4","top_snp_h4_pp","nsnps")
  for (cc in need_cols) if (!(cc %in% names(dt))) dt[, (cc) := NA]
  
  dt[, .(
    trait      = meta$trait,
    qtl_source = meta$qtl_source,
    tissue     = meta$tissue,
    mode       = meta$mode,
    locus      = as.character(locus),
    best_gene  = as.character(best_gene),
    best_gene_ens = sub("\\..*$", "", as.character(best_gene)),
    PP4        = suppressWarnings(as.numeric(PP4)),
    top_snp_h4 = as.character(top_snp_h4),
    top_snp_h4_pp = suppressWarnings(as.numeric(top_snp_h4_pp)),
    nsnps      = suppressWarnings(as.integer(nsnps)),
    folder     = folder,
    file       = f
  )]
}

master_coloc <- rbindlist(lapply(coloc_files, read_coloc_one), fill=TRUE)

if (nrow(master_coloc) == 0) {
  dbg <- data.table(
    file = coloc_files,
    folder = basename(dirname(coloc_files)),
    mode = vapply(basename(dirname(coloc_files)), function(x) parse_from_folder(x)$mode, character(1))
  )
  fwrite(head(dbg, 300), file.path(out_dir, "debug_coloc_files_first300.tsv"), sep="\t")
  stop("[Error] master_coloc is empty. See debug_coloc_files_first300.tsv in output/summary_v1.")
}

message("[Info] Master coloc rows: ", nrow(master_coloc))

# =========================
# 2) BioMart annotation: HGNC symbol + mouse 1:1 ortholog
# =========================
genes <- unique(master_coloc$best_gene_ens)
genes <- genes[!is.na(genes) & genes != ""]
message("[Info] Unique best_gene_ens for annotation: ", length(genes))

connect_mart <- function() {
  for (m in c("useast","uswest","asia")) {
    ans <- tryCatch(
      useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", mirror=m),
      error = function(e) NULL
    )
    if (!is.null(ans)) {
      message("[Info] BioMart connected (mirror=", m, ").")
      return(ans)
    }
  }
  ans <- tryCatch(
    useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org"),
    error = function(e) NULL
  )
  if (is.null(ans)) stop("[Error] Cannot connect BioMart (network/Ensembl server).")
  message("[Info] BioMart connected (fallback host=useast.ensembl.org).")
  ans
}
mart <- connect_mart()

# split attribute pages
attrs_gene <- c("ensembl_gene_id","hgnc_symbol","external_gene_name","gene_biotype")
attrs_hom  <- c("ensembl_gene_id",
                "mmusculus_homolog_ensembl_gene",
                "mmusculus_homolog_associated_gene_name",
                "mmusculus_homolog_orthology_type")

fetch_chunk <- function(vals, attrs, max_retry=3) {
  for (k in seq_len(max_retry)) {
    ans <- tryCatch(
      getBM(attributes=attrs, filters="ensembl_gene_id", values=vals, mart=mart),
      error = function(e) e
    )
    if (!inherits(ans, "error")) return(as.data.table(ans))
    message("[Warn] BioMart failed (try ", k, "/", max_retry, "): ", ans$message)
    Sys.sleep(1.2 * k)
  }
  stop("[Error] BioMart failed after retries. Try rerun later.")
}

chunk_size <- 400
chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))

gene_map <- rbindlist(lapply(chunks, fetch_chunk, attrs=attrs_gene), fill=TRUE)
hom_map  <- rbindlist(lapply(chunks, fetch_chunk, attrs=attrs_hom),  fill=TRUE)

gene_map <- unique(gene_map)
hom_map  <- unique(hom_map)
map <- merge(gene_map, hom_map, by="ensembl_gene_id", all=TRUE)
setDT(map)

# human symbol preference
map[, human_symbol := hgnc_symbol]
map[is.na(human_symbol) | human_symbol=="", human_symbol := external_gene_name]
map[is.na(human_symbol) | human_symbol=="", human_symbol := ensembl_gene_id]

# mouse 1:1 only
map[, mouse_symbol_one2one := fifelse(
  mmusculus_homolog_orthology_type == "ortholog_one2one" &
    !is.na(mmusculus_homolog_associated_gene_name) &
    mmusculus_homolog_associated_gene_name != "",
  mmusculus_homolog_associated_gene_name,
  NA_character_
)]

# output #2
map_one2one <- unique(map[!is.na(mouse_symbol_one2one), .(
  ensembl_gene_id,
  human_symbol,
  gene_biotype,
  mouse_ensembl = mmusculus_homolog_ensembl_gene,
  mouse_symbol  = mouse_symbol_one2one
)])
fwrite(map_one2one, file.path(out_dir, "gene_map_human_mouse_one2one.tsv"), sep="\t")

# dedupe mapping before merge (FIXED)
anno_cols <- map[, .(ensembl_gene_id, human_symbol, gene_biotype, mouse_symbol_one2one)]
anno_cols[, has_one2one := as.integer(!is.na(mouse_symbol_one2one) & mouse_symbol_one2one != "")]
setorder(anno_cols, ensembl_gene_id, -has_one2one)
anno_cols <- unique(anno_cols, by = "ensembl_gene_id")
anno_cols[, has_one2one := NULL]

master_anno <- merge(
  master_coloc,
  anno_cols,
  by.x = "best_gene_ens",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

# =========================
# 3) Output #3: minimal GWAS feature for scRNA fusion (vascular, strict)
# =========================
master_anno[, gene := human_symbol]
master_anno[is.na(gene) | gene == "", gene := best_gene_ens]

master_anno[, tissue_group := fifelse(grepl("^Artery_", tissue), "artery",
                                      fifelse(tissue == "blood", "blood", "other"))]

dt_vasc <- master_anno[
  mode == "strict" &
    tissue_group %in% c("artery","blood") &
    !is.na(PP4)
][order(-PP4)]

feat_vasc <- dt_vasc[, .SD[1], by = .(trait, gene)]

feat_vasc_out <- feat_vasc[, .(
  trait,
  gene,
  mouse_gene_1to1 = mouse_symbol_one2one,
  gwas_pp4 = PP4
)]
fwrite(feat_vasc_out, file.path(out_dir, "feature_gwas_gene_vascular_strict.tsv"), sep="\t")

# =========================
# 4) Output #1: strict annotated master for trace-back
# =========================
fwrite(master_anno[mode == "strict"],
       file.path(out_dir, "master_coloc_strict_annotated.tsv"),
       sep="\t")

message("[Done] Outputs:")
message("  - ", file.path(out_dir, "feature_gwas_gene_vascular_strict.tsv"))
message("  - ", file.path(out_dir, "gene_map_human_mouse_one2one.tsv"))
message("  - ", file.path(out_dir, "master_coloc_strict_annotated.tsv"))