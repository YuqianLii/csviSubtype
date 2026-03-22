#GWAS-pre data    ####
# =========================
# gwas_pipeline_v1.R
# One-entry pipeline: raw -> clean -> loci (+ logs + summaries)
# =========================
# NOTE (2026-02-28): EAF/MAF decoupled. Never map MAF -> EAF. Added MAF column; palindromic filtering uses EAF else MAF.
getwd()
.libPaths(c("E:/R-4.4.0/librarygwas", setdiff(.libPaths(), "E:/R-4.4.0/librarygwas")))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -------------------------
# CONFIG
# -------------------------
# 建议 setwd 到项目根目录（里面有 gwas_manifest.tsv 和 gwas_raw 文件夹）
# setwd("D:/GWAS_project")  # 你自己改

MANIFEST <- "gwas_manifest.tsv"

OUT_DIR <- "output/gwas_v1"
DIR_CLEAN <- file.path(OUT_DIR, "clean")
DIR_LOCI  <- file.path(OUT_DIR, "loci")
DIR_LOG   <- file.path(OUT_DIR, "logs")
DIR_ERR   <- file.path(OUT_DIR, "errors")
DIR_SUM   <- file.path(OUT_DIR, "summary")

# loci 参数（先用距离窗口，后续你要换 LD clumping/fine-mapping 再改）
P_THR <- 5e-8
WINDOW_KB <- 250

# 输出格式：推荐 tsv（更快、更稳）
WRITE_AS <- "tsv"  # "tsv" or "csv"

# 跑之前先不处理全量：可用于 debug（比如只读前 2e6 行）
# NOTE: fread 的 nrows 只影响读取量，真正跑全量就设成 Inf / NULL
READ_NROWS <- -1  # NULL 表示全量；例如 2000000 表示只读前200万行

# 如果输出已存在是否跳过（利于断点续跑）
SKIP_IF_EXISTS <- TRUE

# 回文位点过滤策略（建议保留默认）
DROP_PALINDROMIC_IF_EAF_MID <- TRUE
PAL_MID_LOW <- 0.42
PAL_MID_HIGH <- 0.58

# -------------------------
# LOGGING
# -------------------------
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOG, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(..., level = "INFO") {
  msg <- paste0("[", level, "] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# -------------------------
# 00_utils_gwas (merged)
# -------------------------
.upper_allele <- function(x) toupper(trimws(as.character(x)))

.standardize_cols <- function(dt) {
  nms <- names(dt)
  nms_low <- tolower(nms)
  
  pick <- function(cands) {
    idx <- which(nms_low %in% tolower(cands))
    if (length(idx) == 0) return(NA_character_)
    nms[idx[1]]
  }
  
  col_chr  <- pick(c("chr", "chrom", "chromosome"))
  col_bp   <- pick(c("bp", "pos", "position", "base_pair_location"))
  col_snp  <- pick(c("snp", "rsid", "rs_id", "markername", "variant_id", "id", "rs"))
  col_a1   <- pick(c("a1", "alt", "effect_allele", "ea", "allele1", "tested_allele", "t_allele"))
  col_a2   <- pick(c("a2", "ref", "other_allele", "oa", "allele2", "non_tested_allele", "nt_allele"))
  col_beta <- pick(c("beta", "b", "effect", "estimate"))
  col_or   <- pick(c("or", "odds_ratio"))
  col_se   <- pick(c("se", "stderr", "standard_error"))
  col_p    <- pick(c("p", "pval", "p_value", "p-value", "pvalue"))
  col_n    <- pick(c("n", "samplesize", "n_total", "neff", "n_eff", "obs_ct"))
  col_eaf  <- pick(c("eaf", "af", "effect_allele_freq", "freq1", "a1freq", "a1_freq"))
  col_maf  <- pick(c("maf", "minor_allele_freq", "minor_af"))
  if (is.na(col_chr) || is.na(col_bp) || is.na(col_a1) || is.na(col_a2) || is.na(col_p)) {
    stopf("Missing required columns. Need at least CHR/BP + alleles + P. Found columns: %s",
          paste(nms, collapse = ", "))
  }
  
  out <- data.table(
    CHR = suppressWarnings(as.integer(dt[[col_chr]])),
    BP  = suppressWarnings(as.integer(dt[[col_bp]])),
    SNP = if (!is.na(col_snp)) as.character(dt[[col_snp]]) else NA_character_,
    EA  = .upper_allele(dt[[col_a1]]),
    OA  = .upper_allele(dt[[col_a2]]),
    P   = suppressWarnings(as.numeric(dt[[col_p]]))
  )
  
  if (!is.na(col_beta)) {
    out[, BETA := suppressWarnings(as.numeric(dt[[col_beta]]))]
  } else if (!is.na(col_or)) {
    out[, BETA := log(suppressWarnings(as.numeric(dt[[col_or]])))]
  } else {
    out[, BETA := NA_real_]
  }
  
  if (!is.na(col_se)) {
    out[, SE := suppressWarnings(as.numeric(dt[[col_se]]))]
  } else {
    out[, SE := NA_real_]
  }
  
  out[, N   := if (!is.na(col_n)) suppressWarnings(as.numeric(dt[[col_n]])) else NA_real_]
  out[, EAF := if (!is.na(col_eaf)) suppressWarnings(as.numeric(dt[[col_eaf]])) else NA_real_]
  out[, MAF := if (!is.na(col_maf)) suppressWarnings(as.numeric(dt[[col_maf]])) else NA_real_]
  
  out[is.na(SNP) | SNP == "", SNP := paste0(CHR, ":", BP, ":", EA, ":", OA)]
  out[]
}

.qc_sumstats <- function(dt,
                         drop_palindromic_if_eaf_mid = TRUE,
                         pal_mid_low = 0.42,
                         pal_mid_high = 0.58) {
  dt <- copy(dt)
  
  dt <- dt[!is.na(CHR) & !is.na(BP) & CHR > 0 & BP > 0]
  dt <- dt[!is.na(P) & P > 0 & P <= 1]
  dt <- dt[EA %in% c("A","C","G","T") & OA %in% c("A","C","G","T")]
  dt <- dt[EA != OA]
  
  dt[!is.na(SE) & SE <= 0, SE := NA_real_]
  
  # sanitize EAF/MAF
  if (!("MAF" %in% names(dt))) dt[, MAF := NA_real_]
  dt[, EAF := suppressWarnings(as.numeric(EAF))]
  dt[, MAF := suppressWarnings(as.numeric(MAF))]
  dt[!(is.finite(EAF) & EAF > 0 & EAF < 1), EAF := NA_real_]
  dt[!(is.finite(MAF) & MAF > 0 & MAF < 1), MAF := NA_real_]
  dt[is.na(MAF) & is.finite(EAF), MAF := pmin(EAF, 1 - EAF)]

  
  # de-dup: keep smallest P
  setorder(dt, SNP, P)
  dt <- dt[!duplicated(SNP)]
  
  dt[, is_pal := (EA == "A" & OA == "T") | (EA == "T" & OA == "A") |
       (EA == "C" & OA == "G") | (EA == "G" & OA == "C")]
  
  if (drop_palindromic_if_eaf_mid) {
    # Prefer EAF for palindromic disambiguation; if EAF missing, fallback to MAF.
    dt[, PAL_FREQ := fifelse(is.finite(EAF), EAF, fifelse(is.finite(MAF), MAF, NA_real_))]
    # If frequency is missing for a palindromic SNP, drop it (conservative).
    dt <- dt[!(is_pal == TRUE & (is.na(PAL_FREQ) | (PAL_FREQ >= pal_mid_low & PAL_FREQ <= pal_mid_high)))]
    dt[, PAL_FREQ := NULL]
  }
  
  dt[, Z := ifelse(!is.na(BETA) & !is.na(SE), BETA / SE, NA_real_)]
  setorder(dt, CHR, BP)
  dt[]
}

read_and_standardize_sumstats <- function(path, nrows = NULL) {
  if (!file.exists(path)) stopf("File not found: %s", path)
  if (is.numeric(nrows) && length(nrows) == 1L && !is.na(nrows) && nrows > 0) {
    dt_raw <- fread(path, nrows = as.integer(nrows), showProgress = TRUE)
  } else {
    dt_raw <- fread(path, showProgress = TRUE)  # 全量读取（nrows 为 NULL / -1 / 0 都走这里）
  }
  dt_std <- .standardize_cols(dt_raw)
  dt_qc  <- .qc_sumstats(dt_std,
                         drop_palindromic_if_eaf_mid = DROP_PALINDROMIC_IF_EAF_MID,
                         pal_mid_low = PAL_MID_LOW,
                         pal_mid_high = PAL_MID_HIGH)
  dt_qc
}

write_table <- function(dt, path, fmt = "tsv") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (fmt == "csv") {
    fwrite(dt, path, sep = ",")
  } else {
    fwrite(dt, path, sep = "\t")
  }
  invisible(path)
}

# -------------------------
# 01 loci (merged)
# -------------------------
extract_loci_by_distance <- function(dt, p_thr = 5e-8, window_kb = 250, max_loci_per_chr = Inf) {
  sig <- dt[P <= p_thr]
  if (nrow(sig) == 0) {
    return(list(lead = data.table(), loci = data.table()))
  }
  sig <- sig[order(CHR, P, BP)]
  
  lead_list <- list()
  loci_list <- list()
  
  for (chr in sort(unique(sig$CHR))) {
    x <- sig[CHR == chr]
    loci_n <- 0L
    
    while (nrow(x) > 0 && loci_n < max_loci_per_chr) {
      loci_n <- loci_n + 1L
      lead <- x[1]
      
      start_bp <- max(1L, lead$BP - window_kb * 1000L)
      end_bp   <- lead$BP + window_kb * 1000L
      locus_id <- sprintf("chr%s_locus%03d", chr, loci_n)
      
      loci_list[[length(loci_list) + 1L]] <- data.table(
        locus_id = locus_id,
        CHR = chr,
        start_bp = start_bp,
        end_bp = end_bp,
        lead_SNP = lead$SNP,
        lead_BP = lead$BP,
        lead_P = lead$P,
        lead_BETA = lead$BETA,
        lead_SE = lead$SE
      )
      
      lead_list[[length(lead_list) + 1L]] <- data.table(
        locus_id = locus_id,
        CHR = chr,
        BP = lead$BP,
        SNP = lead$SNP,
        EA = lead$EA,
        OA = lead$OA,
        P = lead$P,
        BETA = lead$BETA,
        SE = lead$SE,
        Z = lead$Z
      )
      
      x <- x[!(BP >= start_bp & BP <= end_bp)]
    }
  }
  
  list(
    lead = rbindlist(lead_list, fill = TRUE),
    loci = rbindlist(loci_list, fill = TRUE)
  )
}

# -------------------------
# PROCESS ONE DATASET
# -------------------------
process_one <- function(trait_id, role, ancestry, file_path) {
  log_msg("Start | ", trait_id, " | role=", role, " | ancestry=", ancestry, " | file=", file_path)
  
  clean_path <- file.path(DIR_CLEAN, paste0(trait_id, ".clean.", ifelse(WRITE_AS == "csv", "csv", "tsv")))
  lead_path  <- file.path(DIR_LOCI,  paste0(trait_id, ".lead.",  ifelse(WRITE_AS == "csv", "csv", "tsv")))
  loci_path  <- file.path(DIR_LOCI,  paste0(trait_id, ".loci.",  ifelse(WRITE_AS == "csv", "csv", "tsv")))
  sum_path   <- file.path(DIR_SUM,   paste0(trait_id, ".summary.tsv"))
  
  if (SKIP_IF_EXISTS && file.exists(clean_path) && file.exists(lead_path) && file.exists(loci_path)) {
    log_msg("Skip (outputs exist) | ", trait_id, level = "WARN")
    return(invisible(list(trait_id = trait_id, skipped = TRUE)))
  }
  
  dt <- read_and_standardize_sumstats(file_path, nrows = READ_NROWS)
  
  # quick summary
  sum_dt <- data.table(
    trait_id = trait_id,
    role = role,
    ancestry = ancestry,
    n_variants_clean = nrow(dt),
    n_missing_beta = sum(is.na(dt$BETA)),
    n_missing_se = sum(is.na(dt$SE)),
    n_missing_eaf = sum(is.na(dt$EAF)),
    n_missing_maf = sum(is.na(dt$MAF)),
    n_palindromic = sum(dt$is_pal, na.rm = TRUE),
    min_p = suppressWarnings(min(dt$P, na.rm = TRUE))
  )
  
  # write clean
  write_table(dt, clean_path, fmt = WRITE_AS)
  
  # loci
  loci_res <- extract_loci_by_distance(dt, p_thr = P_THR, window_kb = WINDOW_KB)
  write_table(loci_res$lead, lead_path, fmt = WRITE_AS)
  write_table(loci_res$loci, loci_path, fmt = WRITE_AS)
  
  sum_dt[, n_gws_sig := sum(dt$P <= P_THR, na.rm = TRUE)]
  sum_dt[, n_lead := nrow(loci_res$lead)]
  sum_dt[, n_loci := nrow(loci_res$loci)]
  
  dir.create(dirname(sum_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(sum_dt, sum_path, sep = "\t")
  
  log_msg("Done | ", trait_id,
          " | clean=", nrow(dt),
          " | GWS=", sum_dt$n_gws_sig,
          " | lead=", sum_dt$n_lead,
          " | loci=", sum_dt$n_loci)
  
  invisible(list(trait_id = trait_id, skipped = FALSE))
}

# -------------------------
# MAIN
# -------------------------
dir.create(DIR_CLEAN, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOCI, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ERR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SUM,  recursive = TRUE, showWarnings = FALSE)

if (!file.exists(MANIFEST)) stopf("Manifest not found: %s", MANIFEST)

manifest <- fread(MANIFEST)
need_cols <- c("trait_id", "role", "ancestry", "file")
miss <- setdiff(need_cols, names(manifest))
if (length(miss) > 0) stopf("Manifest missing columns: %s", paste(miss, collapse = ", "))

# normalize
manifest[, role := tolower(role)]
if (any(!manifest$role %in% c("discovery", "validation"))) {
  stopf("role must be discovery/validation. Found: %s", paste(unique(manifest$role), collapse = ", "))
}

# run
log_msg("Pipeline start. n=", nrow(manifest))
results <- list()

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  trait_id <- row$trait_id
  role <- row$role
  ancestry <- row$ancestry
  file_path <- row$file
  
  tryCatch({
    results[[length(results) + 1L]] <- process_one(trait_id, role, ancestry, file_path)
  }, error = function(e) {
    err_file <- file.path(DIR_ERR, paste0(trait_id, ".error.txt"))
    writeLines(as.character(e), err_file)
    log_msg("ERROR | ", trait_id, " | ", as.character(e), level = "ERROR")
  })
}

# merge summaries
sum_files <- list.files(DIR_SUM, pattern = "\\.summary\\.tsv$", full.names = TRUE)
if (length(sum_files) > 0) {
  sum_all <- rbindlist(lapply(sum_files, fread), fill = TRUE)
  fwrite(sum_all, file.path(DIR_SUM, "ALL.summary.tsv"), sep = "\t")
  log_msg("Wrote summary: ", file.path(DIR_SUM, "ALL.summary.tsv"))
}

log_msg("Pipeline finished.")


