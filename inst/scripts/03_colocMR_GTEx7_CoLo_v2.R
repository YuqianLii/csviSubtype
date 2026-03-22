# =========================
# coloc_pipeline_gtexv7_v1.R
# ONE-SHOT: GWAS (with/without EAF) + GTEx v7 allpairs QTL (NO rsid) + MR enabled; FORCE MR optional
# LENIENT VERSION: expand candidate pool + always run strict & relaxed + export QTL to output/coloc_v1
# =========================
getwd()
.libPaths(c("E:/R-4.4.0/librarygwas", setdiff(.libPaths(), "E:/R-4.4.0/librarygwas")))


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(coloc)
  library(TwoSampleMR)
})

# ========= 0) Config =========
cfg <- list(
  gwas_clean_dir = "output/gwas_v1/clean",
  loci_dir       = "output/gwas_v1/loci",
  
  # 用 eQTLGen AF ref 来补 GWAS 的 EAF/MAF（你现成就有，继续复用）
  eqtlgen_af_file  = "input/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
  
  # ---- GTEx v7 allpairs (your inputs) ----
  gtex_allpairs = list(
    Artery_Aorta             = "input/Artery_Aorta.allpairs.txt",
    Artery_Tibial            = "input/Artery_Tibial.allpairs.txt",
    Artery_Coronary          = "input/Artery_Coronary.allpairs.txt",
    Brain_Cortex             = "input/Brain_Cortex.allpairs.txt",
    Brain_Frontal_Cortex_BA9 = "input/Brain_Frontal_Cortex_BA9.allpairs.txt"
  ),
  
  # GTEx clean 输出目录（每个 tissue 一份 clean subset）
  gtex_clean_dir = "input/qtl/GTEx_v7",
  
  gwas_tags = c("LAC_discovery", "WMH_discovery"),
  
  gwas_run_tags = c("LAC_discovery", "WMH_discovery", "LAC_validation", "WMH_validation"),
  
  loci_source_map = list(
    LAC_validation = "LAC_discovery",
    WMH_validation = "WMH_discovery"
  ),  
  # QTL 预处理：用最大窗口一次覆盖（用于从 allpairs 里筛 loci 区域）
  prep_window_kb = 1500,
  
  # ---- 更宽松的参数（扩展候选池）----
  strict_window_kb = 1000,
  strict_min_overlap_gene = 3,
  strict_min_overlap_locus = 5,
  
  relaxed_window_kb = 1500,
  relaxed_min_overlap_gene = 2,
  relaxed_min_overlap_locus = 3,
  
  run_relaxed_always = TRUE,
  topk_gene_per_locus = 5,
  
  # N：GWAS trait meta
  trait_meta = list(
    LAC_discovery  = list(type="cc",    N=225419, s=6030/225419),
    LAC_validation = list(type="cc",    N=262136, s=7338/262136),
    WMH_discovery  = list(type="quant", N=42310,  s=NA_real_),
    WMH_validation = list(type="quant", N=11226,  s=NA_real_),
    QTL            = list(type="quant", N=NA_real_, s=NA_real_)
  ),
  
  coloc_priors = list(p1=1e-4, p2=1e-4, p12=1e-5),
  
  # MR
  run_mr = TRUE,
  pp4_for_mr = 0.8,
  mr_pval_iv = 5e-8,
  mr_run_without_pp4_gate = FALSE,

  # ---- EAF QC & auto-fix (GWAS) ----
  # Some GWAS files provide an 'EAF' column that is actually MAF / not aligned to the effect allele.
  # This QC will detect that pattern and drop bad GWAS EAF, then re-impute EAF from eqtlgen AF ref by EA.
  eaf_qc_enable        = TRUE,
  eaf_qc_n_sample      = 200000L,
  eaf_qc_tol           = 0.02,
  eaf_qc_mismatch_thr  = 0.30,
  eaf_qc_eq_thr        = 0.95,
  export_eaf_qc        = TRUE,

  
  # GTEx QTL：allpairs 没有总样本量列，就用这个常数（你可以之后改成 tissue-specific）
  gtex_N_default = 200,
  
  # ---- big GTEx allpairs robustness ----
  # > gtex_stream_gb: file larger than this -> use streaming (avoid mmap / RAM explosion)
  # chunk_lines: each chunk lines read & processed (streaming mode)
  gtex_stream_gb   = 6,
  gtex_chunk_lines = 200000L,
  
  out_base = "output"
)

# ========= Run only selected tissues (comma-separated) =========
RUN_TISSUES <- Sys.getenv("RUN_TISSUES", "")
if (nzchar(RUN_TISSUES)) {
  keep <- trimws(strsplit(RUN_TISSUES, ",")[[1]])
  keep <- keep[keep %in% names(cfg$gtex_allpairs)]
  if (length(keep) == 0) stop("RUN_TISSUES not found in cfg$gtex_allpairs: ", RUN_TISSUES)
  cfg$gtex_allpairs <- cfg$gtex_allpairs[keep]  # keep your specified order
  message("[RunTissue] Only run: ", paste(names(cfg$gtex_allpairs), collapse=", "))
}

dir.create(cfg$gtex_clean_dir, recursive = TRUE, showWarnings = FALSE)

# ========= helpers =========
find_col <- function(nms, candidates) {
  hit <- intersect(tolower(nms), tolower(candidates))
  if (length(hit) == 0) return(NA_character_)
  nms[match(hit[1], tolower(nms))]
}

drop_ambiguous <- function(dt) {
  amb <- (dt$EA %in% c("A","T") & dt$NEA %in% c("T","A")) |
    (dt$EA %in% c("C","G") & dt$NEA %in% c("G","C"))
  dt[!amb]
}

# ---- GTEx allpairs has NO rsid: use variant_id as SNP key ----
norm_vid <- function(x) gsub("^chr", "", x, ignore.case = TRUE)

comp_base <- function(x) {
  x <- toupper(as.character(x))
  chartr("ACGT", "TGCA", x)
}

# ---- robust parse variant_id (handles X/Y/MT; suppress NA warnings) ----
parse_variant_id <- function(variant_id_vec) {
  parts <- tstrsplit(variant_id_vec, "_", fixed=TRUE)
  if (length(parts) < 4) {
    stop("variant_id parsing failed: expect at least 4 '_' separated fields like chr1_pos_ref_alt_...")
  }
  chr_raw <- parts[[1]]
  pos_raw <- parts[[2]]
  ref_raw <- parts[[3]]
  alt_raw <- parts[[4]]
  
  chr <- gsub("^chr", "", chr_raw, ignore.case=TRUE)
  chr <- toupper(chr)
  chr <- fifelse(chr=="X","23",
                 fifelse(chr=="Y","24",
                         fifelse(chr %in% c("MT","M","MITO"),"25", chr)))
  
  list(
    CHR = suppressWarnings(as.integer(chr)),
    POS = suppressWarnings(as.integer(pos_raw)),
    REF = toupper(as.character(ref_raw)),
    ALT = toupper(as.character(alt_raw))
  )
}

# 在一个 locus 内，把 GWAS rsid 行映射成 GTEx variant_id（通过 CHR+POS+alleles）
map_gwas_region_to_vid <- function(gwas_region_dt, var_lookup_region) {
  # var_lookup_region must have: CHR POS REF ALT variant_id MAF
  if (nrow(gwas_region_dt) == 0 || nrow(var_lookup_region) == 0) return(gwas_region_dt[0])
  
  dt <- copy(gwas_region_dt)
  dt[, row_id := .I]
  dt[, rsid_orig := SNP]
  
  # join by chr/pos
  setkey(var_lookup_region, CHR, POS)
  tmp <- merge(
    dt[, .(row_id, CHR, POS, EA, NEA, MAF)],
    var_lookup_region[, .(CHR, POS, REF, ALT, variant_id, MAF_lookup=MAF)],
    by = c("CHR","POS"),
    allow.cartesian = TRUE
  )
  
  if (nrow(tmp) == 0) return(dt[0])
  
  # direct allele match
  direct <- (tmp$EA == tmp$ALT & tmp$NEA == tmp$REF) | (tmp$EA == tmp$REF & tmp$NEA == tmp$ALT)
  direct[is.na(direct)] <- FALSE
  
  # complement allele match (strand flip)
  EA_c  <- comp_base(tmp$EA)
  NEA_c <- comp_base(tmp$NEA)
  compm <- (EA_c == tmp$ALT & NEA_c == tmp$REF) | (EA_c == tmp$REF & NEA_c == tmp$ALT)
  compm[is.na(compm)] <- FALSE
  
  # ✅ FIX: create columns BEFORE filtering so length always matches
  tmp[, used_comp := compm]
  tmp[, score := fifelse(direct, 2L, 1L)]
  
  # filter to matches
  tmp <- tmp[direct | compm]
  if (nrow(tmp) == 0) return(dt[0])
  
  # prefer direct > complement, then keep first per GWAS row
  setorder(tmp, row_id, -score)
  tmp <- tmp[!duplicated(tmp$row_id)]
  tmp[, score := NULL]
  
  # apply mapping
  dt[, SNP := NA_character_]
  dt[tmp$row_id, SNP := tmp$variant_id]
  
  # if complemented, flip alleles to complemented (beta stays for effect allele, no sign flip needed)
  comp_rows <- tmp$row_id[tmp$used_comp %in% TRUE]
  if (length(comp_rows) > 0) {
    dt[comp_rows, `:=`(EA = comp_base(EA), NEA = comp_base(NEA))]
  }
  
  # fill MAF if missing
  dt[tmp$row_id, MAF := fifelse(is.na(MAF), tmp$MAF_lookup, MAF)]
  
  dt[, row_id := NULL]
  dt <- dt[!is.na(SNP)]
  dt
}

# --- 读取 GWAS/QTL：自动列名映射；允许 EAF/N 缺失；统一生成 MAF（coloc 用） ---
read_sumstats <- function(path) {
  dt <- fread(path)
  
  map <- list(
    SNP  = c("SNP","rsid","RSID","MarkerName","variant_id"),
    CHR  = c("CHR","chr","chrom","chromosome","hg19_chr"),
    POS  = c("POS","pos","BP","bp","position","base_pair_location","hg19_pos","SNPChrPos","snpchrpos"),
    EA   = c("EA","effect_allele","A1","allele1","ALT","alt","Allele1","AlleleAssessed","AlleleA"),
    NEA  = c("NEA","other_allele","A2","allele2","REF","ref","Allele2","OA","AlleleB"),
    BETA = c("BETA","beta","Beta","effect","Estimate","slope","logOR"),
    SE   = c("SE","se","StdErr","stderr","SEbeta","slope_se"),
    P    = c("P","p","PVAL","pval","pvalue","Pvalue","p_value","pval_nominal"),
    EAF  = c("EAF","eaf","AF","af","freq","effect_allele_freq","ALT_AF","alt_af"),
    MAF  = c("MAF","maf"),
    N    = c("N","n","samplesize","sample_size","N_total","Neff","n_samples","num_samples","ma_samples"),
    GENE = c("GENE","gene","gene_id","GENE_ID","gene_name")
  )
  
  rename_if_found <- function(std, cands) {
    col <- find_col(names(dt), cands)
    if (!is.na(col) && col != std) setnames(dt, col, std)
  }
  
  rename_if_found("SNP",  map$SNP)
  rename_if_found("CHR",  map$CHR)
  rename_if_found("POS",  map$POS)
  rename_if_found("EA",   map$EA)
  rename_if_found("NEA",  map$NEA)
  rename_if_found("BETA", map$BETA)
  rename_if_found("SE",   map$SE)
  rename_if_found("P",    map$P)
  rename_if_found("EAF",  map$EAF)
  rename_if_found("MAF",  map$MAF)
  rename_if_found("N",    map$N)
  rename_if_found("GENE", map$GENE)
  
  if (!("EAF" %in% names(dt))) dt[, EAF := NA_real_]
  if (!("MAF" %in% names(dt))) dt[, MAF := NA_real_]
  if (!("N"   %in% names(dt))) dt[, N   := NA_real_]
  
  need <- c("SNP","CHR","POS","EA","NEA","BETA","SE","P")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("Missing columns in ", path, ": ", paste(miss, collapse=", "),
         "\nAvailable columns are: ", paste(names(dt), collapse=", "))
  }
  
  dt[, `:=`(
    CHR = as.integer(CHR),
    POS = as.integer(POS),
    EA  = toupper(as.character(EA)),
    NEA = toupper(as.character(NEA)),
    EAF = as.numeric(EAF),
    MAF = as.numeric(MAF),
    N   = as.numeric(N),
    SE  = as.numeric(SE),
    BETA= as.numeric(BETA),
    P   = as.numeric(P)
  )]
  
  dt[!(is.finite(MAF) & MAF > 0 & MAF < 1), MAF := NA_real_]
  dt[is.na(MAF) & is.finite(EAF), MAF := pmin(EAF, 1 - EAF)]
  dt[!(is.finite(EAF) & EAF > 0 & EAF < 1), EAF := NA_real_]
  
  dt[, varbeta := SE^2]
  as.data.table(dt)
}

choose_loci_path <- function(loci_dir, tag) {
  loci_lead <- file.path(loci_dir, paste0(tag, ".lead.tsv"))
  loci_tsv  <- file.path(loci_dir, paste0(tag, ".loci.tsv"))
  loci_raw  <- file.path(loci_dir, paste0(tag, ".loci"))
  
  if (file.exists(loci_lead)) return(loci_lead)
  if (file.exists(loci_tsv))  return(loci_tsv)
  if (file.exists(loci_raw))  return(loci_raw)
  
  stop("Missing loci file for tag=", tag,
       "\nTried: ", loci_lead, "\n", loci_tsv, "\n", loci_raw)
}

infer_regions <- function(loci_path, window_kb) {
  loc <- fread(loci_path)
  nms0 <- names(loc)
  nms <- tolower(gsub("^#", "", nms0))
  names(loc) <- nms
  
  if (all(c("chr","start","end") %in% nms)) {
    return(as.data.table(loc[, .(chr=as.integer(chr), start=as.integer(start), end=as.integer(end))]))
  }
  if (all(c("chr","bp1","bp2") %in% nms)) {
    return(as.data.table(loc[, .(chr=as.integer(chr), start=as.integer(bp1), end=as.integer(bp2))]))
  }
  
  chr_col <- intersect(nms, c("chr","chrom","chromosome","chrom"))
  pos_col <- intersect(nms, c("pos","bp","position","snpchrpos"))
  if (length(chr_col) >= 1 && length(pos_col) >= 1) {
    chr <- as.integer(loc[[chr_col[1]]])
    pos <- as.integer(loc[[pos_col[1]]])
    return(data.table(chr=chr,
                      start=pos - window_kb*1000L,
                      end=pos + window_kb*1000L))
  }
  
  stop("Cannot infer regions from loci file: ", loci_path,
       "\nColumns are: ", paste(nms0, collapse=", "))
}

extract_region <- function(dt, chr, start, end) {
  dt[CHR == chr & POS >= start & POS <= end]
}

harmonise_two <- function(d1, d2, min_n=20) {
  m <- inner_join(as.data.frame(d1), as.data.frame(d2), by="SNP", suffix=c(".t1",".t2"))
  if (nrow(m) < min_n) return(NULL)
  
  m <- m[!is.na(m$EA.t1) & !is.na(m$NEA.t1) & !is.na(m$EA.t2) & !is.na(m$NEA.t2), ]
  if (nrow(m) < min_n) return(NULL)
  
  same_set <- (m$EA.t1 == m$EA.t2 & m$NEA.t1 == m$NEA.t2) |
    (m$EA.t1 == m$NEA.t2 & m$NEA.t1 == m$EA.t2)
  same_set[is.na(same_set)] <- FALSE
  m <- m[same_set, ]
  if (nrow(m) < min_n) return(NULL)
  
  need_flip <- (m$EA.t1 == m$NEA.t2 & m$NEA.t1 == m$EA.t2)
  m$BETA.t2[need_flip] <- -m$BETA.t2[need_flip]
  if ("EAF.t2" %in% names(m)) m$EAF.t2[need_flip]  <- 1 - m$EAF.t2[need_flip]
  m$EA.t2[need_flip]   <- m$EA.t1[need_flip]
  m$NEA.t2[need_flip]  <- m$NEA.t1[need_flip]
  
  maf1 <- if ("MAF.t1" %in% names(m)) m$MAF.t1 else rep(NA_real_, nrow(m))
  maf2 <- if ("MAF.t2" %in% names(m)) m$MAF.t2 else rep(NA_real_, nrow(m))
  
  maf1 <- ifelse(is.finite(maf1), maf1, ifelse(is.finite(m$EAF.t1), pmin(m$EAF.t1, 1 - m$EAF.t1), NA_real_))
  maf2 <- ifelse(is.finite(maf2), maf2, ifelse(is.finite(m$EAF.t2), pmin(m$EAF.t2, 1 - m$EAF.t2), NA_real_))
  
  t1 <- as.data.table(m %>% transmute(
    SNP, CHR=CHR.t1, POS=POS.t1, EA=EA.t1, NEA=NEA.t1,
    BETA=BETA.t1, SE=SE.t1, P=P.t1, EAF=EAF.t1, N=N.t1,
    varbeta=(SE.t1^2), MAF=maf1
  ))
  
  t2 <- as.data.table(m %>% transmute(
    SNP, CHR=CHR.t2, POS=POS.t2, EA=EA.t2, NEA=NEA.t2,
    BETA=BETA.t2, SE=SE.t2, P=P.t2, EAF=EAF.t2, N=N.t2,
    varbeta=(SE.t2^2), MAF=maf2
  ))
  
  list(t1=drop_ambiguous(t1), t2=drop_ambiguous(t2))
}

# ===== PATCH (A): run_coloc_abf now supports min_snps and N fallback =====
run_coloc_abf <- function(t1, t2, meta1, meta2, priors, min_snps=3) {
  t1 <- t1[order(P)][!duplicated(SNP)]
  t2 <- t2[order(P)][!duplicated(SNP)]
  t1 <- t1[is.finite(BETA) & is.finite(varbeta) & is.finite(MAF) & MAF > 0 & MAF < 1]
  t2 <- t2[is.finite(BETA) & is.finite(varbeta) & is.finite(MAF) & MAF > 0 & MAF < 1]
  if (nrow(t1) < min_snps || nrow(t2) < min_snps) return(NULL)
  
  N1 <- median(t1$N, na.rm=TRUE)
  if (!is.finite(N1) && !is.null(meta1$N) && is.finite(meta1$N)) N1 <- meta1$N
  
  N2 <- median(t2$N, na.rm=TRUE)
  if (!is.finite(N2) && !is.null(meta2$N) && is.finite(meta2$N)) N2 <- meta2$N
  
  if (!is.finite(N1) || !is.finite(N2)) return(NULL)
  
  d1 <- list(
    snp     = t1$SNP,
    beta    = t1$BETA,
    varbeta = t1$varbeta,
    MAF     = t1$MAF,
    N       = N1,
    type    = meta1$type
  )
  if (meta1$type == "cc") d1$s <- meta1$s
  
  d2 <- list(
    snp     = t2$SNP,
    beta    = t2$BETA,
    varbeta = t2$varbeta,
    MAF     = t2$MAF,
    N       = N2,
    type    = meta2$type
  )
  if (meta2$type == "cc") d2$s <- meta2$s
  
  cres <- NULL
  capture.output({
    cres <- coloc.abf(d1, d2, p1=priors$p1, p2=priors$p2, p12=priors$p12)
  })
  cres
}

# ========= AF ref (for GWAS EAF/MAF imputation; rsid-based) =========
load_eqtlgen_af_ref <- function(af_file) {
  af <- fread(af_file)
  need <- c("SNP","AlleleA","AlleleB","AlleleB_all")
  miss <- setdiff(need, names(af))
  if (length(miss) > 0) stop("AF file missing columns: ", paste(miss, collapse=", "))
  
  af <- af[, .(
    SNP = as.character(SNP),
    A = toupper(as.character(AlleleA)),
    B = toupper(as.character(AlleleB)),
    FREQ_B = as.numeric(AlleleB_all)
  )]
  af <- af[is.finite(FREQ_B) & FREQ_B > 0 & FREQ_B < 1]
  setkey(af, SNP)
  af
}
af_ref <- load_eqtlgen_af_ref(cfg$eqtlgen_af_file)

fill_missing_eaf <- function(gwas_dt, af_ref) {
  if (!("EAF" %in% names(gwas_dt))) gwas_dt[, EAF := NA_real_]
  if (!("MAF" %in% names(gwas_dt))) gwas_dt[, MAF := NA_real_]
  if (!("EAF_SRC" %in% names(gwas_dt))) gwas_dt[, EAF_SRC := "missing"]
  
  gwas_dt[, `:=`(EAF=as.numeric(EAF), MAF=as.numeric(MAF))]
  
  gwas_dt[, EAF := fifelse(is.finite(EAF) & EAF > 0 & EAF < 1, EAF, NA_real_)]
  gwas_dt[, MAF := fifelse(is.finite(MAF) & MAF > 0 & MAF < 1, MAF, NA_real_)]
  gwas_dt[, EAF_SRC := fifelse(!is.na(EAF), "gwas", "missing")]
  
  setkey(gwas_dt, SNP)
  setkey(af_ref, SNP)
  
  # 先补 MAF（不依赖 EA）
  gwas_dt[af_ref, on="SNP",
          MAF := fifelse(is.na(MAF), pmin(i.FREQ_B, 1 - i.FREQ_B), MAF)]
  
  # 再补 EAF（依赖 EA 与 A/B 匹配）
  miss_idx <- which(is.na(gwas_dt$EAF))
  if (length(miss_idx) > 0) {
    tmp <- gwas_dt[miss_idx, .(SNP, EA)]
    tmp[, EA := toupper(as.character(EA))]
    tmp <- merge(tmp, af_ref, by="SNP", all.x=TRUE)
    
    tmp[, EAF_imp := fifelse(EA == B, FREQ_B,
                             fifelse(EA == A, 1 - FREQ_B, NA_real_))]
    setkey(tmp, SNP)
    
    gwas_dt[tmp, on="SNP", EAF := fifelse(is.na(EAF), i.EAF_imp, EAF)]
    gwas_dt[tmp, on="SNP",
            EAF_SRC := fifelse(EAF_SRC=="missing" & !is.na(i.EAF_imp), "imputed_eqtlgenAF", EAF_SRC)]
  }
  
  gwas_dt[is.na(MAF) & is.finite(EAF), MAF := pmin(EAF, 1 - EAF)]
  gwas_dt
}

fill_missing_n <- function(gwas_dt, N_const) {
  if (!("N" %in% names(gwas_dt))) gwas_dt[, N := NA_real_]
  gwas_dt[, N := as.numeric(N)]
  if (!is.null(N_const) && is.finite(N_const)) {
    gwas_dt[!is.finite(N), N := as.numeric(N_const)]
  }
  gwas_dt
}

# ========= EAF QC + GWAS preparation cache (avoid re-reading 8M SNPs for strict/relaxed) =========
.GWAS_CACHE <- new.env(parent = emptyenv())
.GWAS_QC    <- new.env(parent = emptyenv())

.eaf_stats <- function(dt) {
  eaf_ok <- is.finite(dt$EAF)
  maf_ok <- is.finite(dt$MAF)
  list(
    n = nrow(dt),
    eaf_missing_rate = mean(!eaf_ok),
    eaf_min = if (any(eaf_ok)) min(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    eaf_median = if (any(eaf_ok)) median(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    eaf_max = if (any(eaf_ok)) max(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    prop_eaf_gt_0_5 = if (any(eaf_ok)) mean(dt$EAF[eaf_ok] > 0.5, na.rm=TRUE) else NA_real_,
    prop_eaf_eq_maf = {
      ok <- eaf_ok & maf_ok
      if (any(ok)) mean(abs(dt$EAF[ok] - dt$MAF[ok]) < 1e-6) else NA_real_
    },
    prop_eaf_imputed = if ("EAF_SRC" %in% names(dt)) mean(dt$EAF_SRC == "imputed_eqtlgenAF", na.rm=TRUE) else NA_real_
  )
}

.check_gwas_eaf_vs_ref <- function(gwas_dt, af_ref, n_sample=200000L, tol=0.02) {
  # Compare GWAS-provided EAF to expected EAF computed from AF ref by EA (only for SNPs that align to A/B).
  if (!all(c("SNP","EA","EAF") %in% names(gwas_dt))) return(list(mismatch=NA_real_, mismatch_flip=NA_real_, matched_n=0L))
  idx <- which(is.finite(gwas_dt$EAF))
  if (length(idx) < 10000L) return(list(mismatch=NA_real_, mismatch_flip=NA_real_, matched_n=0L))

  set.seed(1)
  if (length(idx) > n_sample) idx <- sample(idx, n_sample)

  tmp <- gwas_dt[idx, .(SNP=as.character(SNP), EA=toupper(as.character(EA)), EAF=as.numeric(EAF))]
  tmp <- merge(tmp, af_ref, by="SNP", all.x=TRUE)
  tmp[, EAF_expected := fifelse(EA == B, FREQ_B,
                                fifelse(EA == A, 1 - FREQ_B, NA_real_))]
  tmp <- tmp[is.finite(EAF_expected)]
  if (nrow(tmp) < 10000L) return(list(mismatch=NA_real_, mismatch_flip=NA_real_, matched_n=as.integer(nrow(tmp))))

  mismatch <- mean(abs(tmp$EAF - tmp$EAF_expected) > tol, na.rm=TRUE)
  mismatch_flip <- mean(abs(tmp$EAF - (1 - tmp$EAF_expected)) > tol, na.rm=TRUE)
  list(mismatch=mismatch, mismatch_flip=mismatch_flip, matched_n=as.integer(nrow(tmp)))
}

.qc_and_fix_gwas_eaf <- function(gwas_dt, af_ref, tag, cfg) {
  # If GWAS EAF behaves like MAF and mismatches AF ref by EA, drop GWAS EAF and re-impute later from AF ref.
  if (!isTRUE(cfg$eaf_qc_enable)) return(list(gwas=gwas_dt, qc=NULL))

  if (!("EAF" %in% names(gwas_dt))) gwas_dt[, EAF := NA_real_]
  if (!("MAF" %in% names(gwas_dt))) gwas_dt[, MAF := NA_real_]
  gwas_dt[, `:=`(EAF=as.numeric(EAF), MAF=as.numeric(MAF))]
  gwas_dt[!(is.finite(EAF) & EAF > 0 & EAF < 1), EAF := NA_real_]
  gwas_dt[!(is.finite(MAF) & MAF > 0 & MAF < 1), MAF := NA_real_]
  gwas_dt[is.na(MAF) & is.finite(EAF), MAF := pmin(EAF, 1 - EAF)]

  st <- .eaf_stats(gwas_dt)
  qc <- data.table(
    tag = tag,
    eaf_missing_rate_raw = st$eaf_missing_rate,
    eaf_min_raw = st$eaf_min,
    eaf_median_raw = st$eaf_median,
    eaf_max_raw = st$eaf_max,
    prop_eaf_gt_0_5_raw = st$prop_eaf_gt_0_5,
    prop_eaf_eq_maf_raw = st$prop_eaf_eq_maf,
    mismatch = NA_real_,
    mismatch_if_flipped = NA_real_,
    matched_n = 0L,
    fixed_drop_eaf = FALSE
  )

  max_eaf <- st$eaf_max
  eq_rate <- st$prop_eaf_eq_maf

  looks_like_maf <- is.finite(max_eaf) && (max_eaf <= 0.5 + 1e-12) && is.finite(eq_rate) && (eq_rate >= cfg$eaf_qc_eq_thr)
  if (isTRUE(looks_like_maf)) {
    chk <- .check_gwas_eaf_vs_ref(gwas_dt, af_ref, n_sample=cfg$eaf_qc_n_sample, tol=cfg$eaf_qc_tol)
    qc[, `:=`(mismatch = chk$mismatch, mismatch_if_flipped = chk$mismatch_flip, matched_n = as.integer(chk$matched_n))]

    flip_better <- is.finite(chk$mismatch_flip) && (chk$mismatch_flip + 0.05 < chk$mismatch)

    if (is.finite(chk$mismatch) && chk$mismatch > cfg$eaf_qc_mismatch_thr && isTRUE(flip_better)) {
      message("[FixEAF][", tag, "] Drop GWAS EAF (looks like MAF; mismatch=", signif(chk$mismatch,4),
              ", flip=", signif(chk$mismatch_flip,4), ", matched=", chk$matched_n, "). Re-impute EAF from AF ref by EA.")
      gwas_dt[, EAF := NA_real_]
      gwas_dt[, EAF_SRC := "missing"]
      qc[, fixed_drop_eaf := TRUE]
    } else if (is.finite(chk$mismatch) && chk$mismatch > cfg$eaf_qc_mismatch_thr) {
      message("[FixEAF][", tag, "] EAF looks like MAF but mismatch is high without clear flip improvement -> keep GWAS EAF (avoid false positive).")
    }
  }

  list(gwas=gwas_dt, qc=qc)
}

.get_prepared_gwas <- function(tag) {
  if (exists(tag, envir=.GWAS_CACHE, inherits=FALSE)) return(get(tag, envir=.GWAS_CACHE, inherits=FALSE))

  gwas_file <- file.path(cfg$gwas_clean_dir, paste0(tag, ".clean.tsv"))
  if (!file.exists(gwas_file)) stop("Missing GWAS file: ", gwas_file)

  gwas0 <- read_sumstats(gwas_file)
  qcres <- .qc_and_fix_gwas_eaf(gwas0, af_ref, tag, cfg)
  gwas1 <- qcres$gwas
  assign(tag, qcres$qc, envir=.GWAS_QC)

  gwas1 <- fill_missing_eaf(gwas1, af_ref)
  gwas1 <- fill_missing_n(gwas1, cfg$trait_meta[[tag]]$N)

  # Keep candidate pool; heavy filters happen later after mapping/harmonisation
  gwas1 <- gwas1[order(P)][!duplicated(SNP)]

  # Attach after-imputation stats to QC
  qc <- get(tag, envir=.GWAS_QC, inherits=FALSE)
  if (!is.null(qc) && nrow(qc)==1) {
    st2 <- .eaf_stats(gwas1)
    qc[, `:=`(
      eaf_missing_rate_after = st2$eaf_missing_rate,
      eaf_min_after = st2$eaf_min,
      eaf_median_after = st2$eaf_median,
      eaf_max_after = st2$eaf_max,
      prop_eaf_gt_0_5_after = st2$prop_eaf_gt_0_5,
      prop_eaf_eq_maf_after = st2$prop_eaf_eq_maf,
      prop_eaf_imputed_after = st2$prop_eaf_imputed
    )]
    assign(tag, qc, envir=.GWAS_QC)
  }

  assign(tag, gwas1, envir=.GWAS_CACHE)
  gwas1
}

.get_gwas_qc <- function(tag) {
  if (exists(tag, envir=.GWAS_QC, inherits=FALSE)) return(get(tag, envir=.GWAS_QC, inherits=FALSE))
  NULL
}



# ========= MR helpers =========
prep_for_tsmr <- function(dt, name, is_exposure=TRUE) {
  out <- dt %>% dplyr::transmute(
    SNP = SNP,
    beta = as.numeric(BETA),
    se = as.numeric(SE),
    effect_allele = EA,
    other_allele = NEA,
    eaf = as.numeric(EAF),
    pval = as.numeric(P),
    samplesize = as.numeric(N)
  ) %>% as.data.frame()
  
  if (is_exposure) {
    out$exposure <- name; out$id.exposure <- name
  } else {
    out$outcome <- name; out$id.outcome <- name
  }
  out
}

format_tsmr <- function(df, type=c("exposure","outcome")) {
  type <- match.arg(type)
  TwoSampleMR::format_data(
    df, type = type,
    snp_col="SNP",
    beta_col="beta",
    se_col="se",
    effect_allele_col="effect_allele",
    other_allele_col="other_allele",
    eaf_col="eaf",
    pval_col="pval",
    samplesize_col="samplesize"
  )
}

run_mr_one <- function(t1, t2, tag, chr, start, end, best_gene, pp4, eaf_imp_frac, qtl_name) {
  tryCatch({
    exp_dt <- t2[t2$P < cfg$mr_pval_iv, ]
    if (nrow(exp_dt) < 1) return(NULL)
    
    out_dt <- t1[t1$SNP %in% exp_dt$SNP, ]
    if (nrow(out_dt) < 1) return(NULL)
    
    exp <- prep_for_tsmr(exp_dt, paste0(qtl_name, ":", best_gene), TRUE)
    out <- prep_for_tsmr(out_dt, tag, FALSE)
    
    exp <- format_tsmr(exp, "exposure")
    out <- format_tsmr(out, "outcome")
    
    dat <- suppressWarnings(TwoSampleMR::harmonise_data(exp, out, action=2))
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) < 1) return(NULL)
    
    method_list <- if (nrow(dat) == 1) c("mr_wald_ratio") else c("mr_ivw")
    mr_res <- TwoSampleMR::mr(dat, method_list = method_list)
    
    mr_res$tag <- tag
    mr_res$chr <- chr; mr_res$start <- start; mr_res$end <- end
    mr_res$best_gene <- best_gene
    mr_res$PP4 <- pp4
    mr_res$gwas_eaf_imputed_frac <- eaf_imp_frac
    mr_res$n_iv <- nrow(dat)
    mr_res
  }, error=function(e){
    message("[MR] failed: ", conditionMessage(e))
    NULL
  })
}

# ========= GTEx v7 allpairs -> clean subset (by loci regions) =========
build_merged_regions <- function(loci_tags, window_kb) {
  regions <- rbindlist(lapply(loci_tags, function(tag){
    loci_path <- choose_loci_path(cfg$loci_dir, tag)
    infer_regions(loci_path, window_kb)
  }), fill=TRUE)
  
  setorder(regions, chr, start, end)
  merged <- regions[, {
    s <- start; e <- end
    o <- order(s); s <- s[o]; e <- e[o]
    cur_s <- s[1]; cur_e <- e[1]
    out <- list()
    if (length(s) > 1) {
      for (i in 2:length(s)) {
        if (s[i] <= cur_e) cur_e <- max(cur_e, e[i])
        else { out[[length(out)+1]] <- c(cur_s, cur_e); cur_s <- s[i]; cur_e <- e[i] }
      }
    }
    out[[length(out)+1]] <- c(cur_s, cur_e)
    data.table(start=sapply(out, `[`, 1), end=sapply(out, `[`, 2))
  }, by=.(chr)]
  setnames(merged, c("chr","start","end"), c("CHR","start","end"))
  merged[, CHR := as.integer(CHR)]
  setkey(merged, CHR, start, end)
  merged
}

# ---- GTEx allpairs without rsid, use variant_id as SNP ----
prepare_gtex_clean_one <- function(allpairs_file, tissue, out_file, merged_regions, N_default=200,
                                   stream_gb=6, chunk_lines=200000L) {
  message("[Prep][GTEx] Tissue=", tissue, " | file=", allpairs_file)
  
  # header: avoid fread(nrows=0) (can trigger mmap issues on big files / Windows)
  hdr_line <- readLines(allpairs_file, n=1L, warn=FALSE)
  if (length(hdr_line) == 0) stop("[GTEx] empty file: ", allpairs_file)
  hdr <- strsplit(hdr_line, "\t", fixed=TRUE)[[1]]
  
  c_variant <- find_col(hdr, c("variant_id","variant","variantid"))
  c_geneid  <- find_col(hdr, c("gene_id","geneid"))
  c_genenm  <- find_col(hdr, c("gene_name","genesymbol","gene","GENE"))
  c_gene    <- if (!is.na(c_genenm)) c_genenm else c_geneid
  
  c_p       <- find_col(hdr, c("pval_nominal","pvalue","pval","P","p"))
  c_beta    <- find_col(hdr, c("slope","beta","BETA","effect"))
  c_se      <- find_col(hdr, c("slope_se","se","SE","stderr","StdErr"))
  c_maf     <- find_col(hdr, c("maf","MAF"))
  c_n       <- find_col(hdr, c("n_samples","num_samples","samplesize","N","n","ma_samples"))
  
  need <- c(c_variant, c_gene, c_p, c_beta, c_se, c_maf)
  if (any(is.na(need))) {
    stop("[GTEx] Cannot map required columns in ", allpairs_file,
         "\nHeader:\n", paste(hdr, collapse=", "),
         "\nNeed at least: variant_id, gene_id/gene_name, pval_nominal, slope, slope_se, maf")
  }
  
  idx_variant <- match(c_variant, hdr)
  idx_gene    <- match(c_gene, hdr)
  idx_p       <- match(c_p, hdr)
  idx_beta    <- match(c_beta, hdr)
  idx_se      <- match(c_se, hdr)
  idx_maf     <- match(c_maf, hdr)
  idx_n       <- if (!is.na(c_n)) match(c_n, hdr) else NA_integer_
  
  idx <- c(idx_variant, idx_gene, idx_p, idx_beta, idx_se, idx_maf)
  nm  <- c("variant_id","GENE","P","BETA","SE","MAF")
  if (!is.na(idx_n)) { idx <- c(idx, idx_n); nm <- c(nm, "N_raw") }
  
  fgb <- file.info(allpairs_file)$size / (1024^3)
  do_stream <- is.finite(fgb) && (fgb >= stream_gb)
  
  if (file.exists(out_file)) file.remove(out_file)
  wrote_header <- FALSE
  
  write_chunk <- function(dt) {
    if (nrow(dt) == 0) return(invisible(NULL))
    
    dt[, variant_id := norm_vid(as.character(variant_id))]
    dt[, SNP := variant_id]
    
    dt[, `:=`(
      GENE = as.character(GENE),
      P = as.numeric(P),
      BETA = as.numeric(BETA),
      SE = as.numeric(SE),
      MAF = as.numeric(MAF)
    )]
    
    pv <- parse_variant_id(dt$variant_id)
    dt[, `:=`(CHR = pv$CHR, POS = pv$POS, NEA = pv$REF, EA = pv$ALT)]
    dt[, EAF := NA_real_]
    
    if ("N_raw" %in% names(dt)) {
      dt[, N := as.numeric(N_raw)]
      dt[, N_raw := NULL]
    } else {
      dt[, N := as.numeric(N_default)]
    }
    dt[!is.finite(N), N := as.numeric(N_default)]
    
    dt <- dt[is.finite(CHR) & is.finite(POS)]
    dt <- dt[is.finite(P) & is.finite(BETA) & is.finite(SE)]
    dt <- dt[is.finite(MAF) & MAF > 0 & MAF < 1]
    if (nrow(dt) == 0) return(invisible(NULL))
    
    v <- dt[, .(CHR=as.integer(CHR), start=as.integer(POS), end=as.integer(POS), .row=.I)]
    setkey(v, CHR, start, end)
    ov <- foverlaps(v, merged_regions, nomatch=0L)
    if (nrow(ov) == 0) return(invisible(NULL))
    dt <- dt[ov$.row]
    
    out <- dt[, .(SNP, CHR, POS, EA, NEA, BETA, SE, P, EAF, MAF, N, GENE)]
    fwrite(out, out_file, sep="\t",
           append = wrote_header,
           col.names = !wrote_header)
    wrote_header <<- TRUE
    invisible(NULL)
  }
  
  # fast path (small files)
  if (!do_stream) {
    ok_fast <- TRUE
    tryCatch({
      dt_all <- fread(allpairs_file, select=idx, showProgress=TRUE)
      setnames(dt_all, nm)
      write_chunk(dt_all)
      rm(dt_all); gc(verbose=FALSE)
    }, error=function(e){
      ok_fast <<- FALSE
      message("[Prep][GTEx] fast fread failed -> fallback to streaming. Reason: ", conditionMessage(e))
    })
    if (ok_fast && file.exists(out_file)) {
      message("[Prep][GTEx] Saved clean subset -> ", out_file, " | (fast path)")
      return(invisible(NULL))
    }
  }
  
  # streaming path (big files)
  message("[Prep][GTEx] Streaming mode ON (chunk_lines=", chunk_lines, ")")
  con <- file(allpairs_file, open="r")
  on.exit(try(close(con), silent=TRUE), add=TRUE)
  
  # skip header
  invisible(readLines(con, n=1L, warn=FALSE))
  
  repeat {
    lines <- readLines(con, n=chunk_lines, warn=FALSE)
    if (length(lines) == 0) break
    
    txt <- paste(lines, collapse="\n")
    dt_chunk <- fread(text=txt, sep="\t", header=FALSE, select=idx, showProgress=FALSE)
    setnames(dt_chunk, nm)
    write_chunk(dt_chunk)
    
    rm(dt_chunk, lines, txt)
    gc(verbose=FALSE)
  }
  
  if (!file.exists(out_file)) {
    stop("[GTEx] Streaming finished but produced empty subset. Check build (b37/hg19 vs hg38) or loci windows.")
  }
  
  message("[Prep][GTEx] Saved clean subset -> ", out_file, " | (streaming path)")
}

# ========= 1) Build merged regions once =========
merged_regions <- build_merged_regions(cfg$gwas_tags, cfg$prep_window_kb)

# ========= 2) Prepare GTEx clean subset for each tissue (auto rebuild if missing) =========
need_qtl_cols <- c("SNP","CHR","POS","EA","NEA","BETA","SE","P","EAF","MAF","N","GENE")

gtex_clean_files <- list()
for (tissue in names(cfg$gtex_allpairs)) {
  in_file <- cfg$gtex_allpairs[[tissue]]
  out_file <- file.path(cfg$gtex_clean_dir, paste0(tissue, ".clean.subset.tsv"))
  gtex_clean_files[[tissue]] <- out_file
  
  rebuild <- FALSE
  if (!file.exists(out_file)) {
    rebuild <- TRUE
  } else {
    hdr_line2 <- readLines(out_file, n=1L, warn=FALSE)
    hdr2 <- strsplit(hdr_line2, "\t", fixed=TRUE)[[1]]
    if (!all(need_qtl_cols %in% hdr2)) rebuild <- TRUE
  }
  
  if (rebuild) {
    prepare_gtex_clean_one(in_file, tissue, out_file, merged_regions,
                           N_default = cfg$gtex_N_default,
                           stream_gb = cfg$gtex_stream_gb,
                           chunk_lines = cfg$gtex_chunk_lines)
  } else {
    message("[Prep][GTEx] Found existing clean file (ok), skip: ", out_file)
  }
  gc(verbose=FALSE)
}

# ========= 2.5) Export all QTL used into output/coloc_v1/_inputs =========
coloc_root <- file.path(cfg$out_base, "coloc_v1")
dir.create(coloc_root, recursive=TRUE, showWarnings=FALSE)
inputs_dir <- file.path(coloc_root, "_inputs")
dir.create(inputs_dir, recursive=TRUE, showWarnings=FALSE)

for (tissue in names(gtex_clean_files)) {
  src <- gtex_clean_files[[tissue]]
  dst <- file.path(inputs_dir, paste0("GTExv7_", tissue, ".clean.subset.used.tsv"))
  ok <- file.copy(src, dst, overwrite=TRUE)
  if (!isTRUE(ok)) stop("[ExportQTL] file.copy failed: ", src, " -> ", dst)
  message("[ExportQTL] OK -> ", dst)
}

# ========= 3) Main runner (parameterized by QTL) =========
run_tag_one_qtl <- function(tag, qtl_dt, qtl_name, window_kb, min_overlap_gene, min_overlap_locus, suffix) {
  
  gwas_file <- file.path(cfg$gwas_clean_dir, paste0(tag, ".clean.tsv"))
  if (!file.exists(gwas_file)) stop("Missing GWAS file: ", gwas_file)
  
  loci_tag  <- if (tag %in% names(cfg$loci_source_map)) cfg$loci_source_map[[tag]] else tag

  loci_tag  <- if (!is.null(cfg$loci_source_map) && tag %in% names(cfg$loci_source_map)) cfg$loci_source_map[[tag]] else tag
  loci_path <- choose_loci_path(cfg$loci_dir, loci_tag)
  
  message("== Running: ", tag, " vs ", qtl_name, " | loci: ", basename(loci_path),
          " | window_kb=", window_kb, " | min_gene=", min_overlap_gene)
  
  gwas <- .get_prepared_gwas(tag)

  # ---- EAF QC summary (optional) ----
  qc0 <- .get_gwas_qc(tag)
  if (!is.null(qc0) && nrow(qc0)==1) {
    message(sprintf("[EAFQC][%s] raw_max=%.4f raw_eqMAF=%.4f mismatch=%.4f fixed_drop_eaf=%s | after_missing=%.4f after_prop_gt0.5=%.4f",
                    tag,
                    qc0$eaf_max_raw, qc0$prop_eaf_eq_maf_raw, qc0$mismatch, as.character(qc0$fixed_drop_eaf),
                    qc0$eaf_missing_rate_after, qc0$prop_eaf_gt_0_5_after))
  }
  
  if (!is.finite(median(gwas$N, na.rm=TRUE)) && !(is.finite(cfg$trait_meta[[tag]]$N))) {
    message("[Warn] ", tag, ": N is missing (both GWAS column N and cfg$trait_meta[[tag]]$N). This tag may yield no coloc results.")
  }
  
  regions <- infer_regions(loci_path, window_kb)
  
  out_dir <- file.path(cfg$out_base, "coloc_v1",
                       paste0(tag, "__vs__", qtl_name, "__", suffix))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  # 每个 out_dir 保存一份 QTL used（HARD CHECK）
  qtl_copy2 <- file.path(out_dir, "qtl_used.tsv")
  used_in_inputs <- file.path(inputs_dir, paste0(qtl_name, ".clean.subset.used.tsv"))
  if (file.exists(used_in_inputs)) {
    ok2 <- file.copy(used_in_inputs, qtl_copy2, overwrite=TRUE)
    if (!isTRUE(ok2)) stop("[ExportQTL] file.copy failed: ", used_in_inputs, " -> ", qtl_copy2)
  } else {
    fwrite(qtl_dt, qtl_copy2, sep="\t")
  }
  message("[ExportQTL] OK -> ", qtl_copy2)


  # ---- export EAF QC report (same for strict/relaxed; saved per out_dir for traceability) ----
  if (isTRUE(cfg$export_eaf_qc)) {
    qc_row <- .get_gwas_qc(tag)
    if (!is.null(qc_row) && nrow(qc_row) >= 1) {
      fwrite(qc_row, file.path(out_dir, "eaf_qc.tsv"), sep="\t")
    }
  }
  
  coloc_rows <- list()
  mr_rows <- list()
  gene_rows <- list()
  
  meta1 <- cfg$trait_meta[[tag]]
  meta2 <- cfg$trait_meta$QTL
  
  for (i in seq_len(nrow(regions))) {
    chr <- regions$chr[i]; start <- regions$start[i]; end <- regions$end[i]
    
    # QTL region
    r2 <- extract_region(qtl_dt, chr, start, end)
    
    # GWAS region (rsid key in file)
    r1_raw <- extract_region(gwas, chr, start, end)
    
    if (nrow(r2) >= min_overlap_locus && nrow(r1_raw) >= min_overlap_locus) {
      var_sub <- unique(r2[, .(
        CHR = as.integer(CHR),
        POS = as.integer(POS),
        REF = toupper(as.character(NEA)),
        ALT = toupper(as.character(EA)),
        variant_id = norm_vid(as.character(SNP)),
        MAF = as.numeric(MAF)
      )], by="variant_id")
      
      r1 <- map_gwas_region_to_vid(r1_raw, var_sub)
    } else {
      next
    }
    
    if (nrow(r1) < min_overlap_locus || nrow(r2) < min_overlap_locus) next
    
    eaf_imp_frac <- if ("EAF_SRC" %in% names(r1)) mean(r1$EAF_SRC == "imputed_eqtlgenAF") else NA_real_
    r1 <- r1[order(P)][!duplicated(SNP)]
    
    if (!("GENE" %in% names(r2))) stop("QTL has no GENE column; cannot run gene-level coloc.")
    genes <- unique(r2$GENE)
    if (length(genes) == 0) next
    
    best_summ <- NULL
    best_pp4 <- -Inf
    best_t1 <- NULL
    best_t2 <- NULL
    
    locus_gene_summ <- list()
    
    for (g in genes) {
      r2g <- r2[GENE == g]
      if (nrow(r2g) < min_overlap_gene) next
      r2g <- r2g[order(P)][!duplicated(SNP)]
      if (length(intersect(r1$SNP, r2g$SNP)) < min_overlap_gene) next
      
      harm <- harmonise_two(r1, r2g, min_n = min_overlap_gene)
      if (is.null(harm)) next
      
      t1 <- harm$t1[order(P)][!duplicated(SNP)]
      t2 <- harm$t2[order(P)][!duplicated(SNP)]
      if (nrow(t1) < min_overlap_gene || nrow(t2) < min_overlap_gene) next
      
      cres <- run_coloc_abf(t1, t2, meta1, meta2, cfg$coloc_priors,
                            min_snps = max(2L, as.integer(min_overlap_gene)))
      if (is.null(cres)) next
      
      summ <- as.data.frame(t(cres$summary))
      pp4 <- as.numeric(summ$PP.H4.abf)
      if (is.na(pp4)) next
      
      summ_all <- summ
      summ_all$tag <- tag
      summ_all$chr <- chr; summ_all$start <- start; summ_all$end <- end
      summ_all$nsnps <- nrow(t1)
      summ_all$PP4 <- pp4
      summ_all$gene <- g
      summ_all$window_kb <- window_kb
      summ_all$min_overlap_gene <- min_overlap_gene
      summ_all$gwas_eaf_imputed_frac <- eaf_imp_frac
      
      top <- cres$results[which.max(cres$results$SNP.PP.H4), ]
      summ_all$top_snp_h4 <- top$snp
      summ_all$top_snp_h4_pp <- top$SNP.PP.H4
      
      locus_gene_summ[[length(locus_gene_summ)+1]] <- summ_all
      
      if (pp4 > best_pp4) {
        best_pp4 <- pp4
        summ_best <- summ_all
        summ_best$best_gene <- g
        best_summ <- summ_best
        best_t1 <- t1
        best_t2 <- t2
      }
    }
    
    if (length(locus_gene_summ) > 0) gene_rows <- c(gene_rows, locus_gene_summ)
    if (is.null(best_summ)) next
    
    coloc_rows[[length(coloc_rows) + 1]] <- best_summ
    
    gate_ok <- isTRUE(cfg$mr_run_without_pp4_gate) || (!is.na(best_summ$PP4) && best_summ$PP4 >= cfg$pp4_for_mr)
    if (isTRUE(cfg$run_mr) && isTRUE(gate_ok) && !is.null(best_t1) && !is.null(best_t2)) {
      mr_res <- run_mr_one(best_t1, best_t2,
                           tag=tag, chr=chr, start=start, end=end,
                           best_gene=best_summ$best_gene,
                           pp4=best_summ$PP4,
                           eaf_imp_frac=eaf_imp_frac,
                           qtl_name=qtl_name)
      if (!is.null(mr_res)) mr_rows[[length(mr_rows) + 1]] <- mr_res
    }
  }
  
  # ===== gene-level 输出 =====
  gene_res_path <- file.path(out_dir, "coloc_gene_results.tsv")
  gene_topk_path <- file.path(out_dir, "coloc_gene_topk.tsv")
  
  if (length(gene_rows) > 0) {
    gene_df <- bind_rows(gene_rows)
    fwrite(gene_df, gene_res_path, sep="\t")
    
    setDT(gene_df)
    gene_df[, locus_id := paste0(tag, ":", chr, ":", start, "-", end)]
    setorder(gene_df, -PP4)
    
    gene_topk <- gene_df[, head(.SD, cfg$topk_gene_per_locus), by=.(locus_id)]
    fwrite(gene_topk, gene_topk_path, sep="\t")
  } else {
    fwrite(data.table(), gene_res_path, sep="\t")
    fwrite(data.table(), gene_topk_path, sep="\t")
  }
  
  if (length(coloc_rows) == 0) {
    message("[Warn] No loci passed for ", tag, " in ", suffix, " vs ", qtl_name)
    mr_empty <- data.table(method=character(), nsnp=integer(), b=numeric(), se=numeric(), pval=numeric(),
                           tag=character(), chr=integer(), start=integer(), end=integer(),
                           best_gene=character(), PP4=numeric(), gwas_eaf_imputed_frac=numeric(), n_iv=integer())
    fwrite(mr_empty, file.path(out_dir, "mr_results.tsv"), sep="\t")
    fwrite(data.table(), file.path(out_dir, "coloc_results.tsv"), sep="\t")
    return(list(n=0, mr_n=0, out_dir=out_dir))
  }
  
  coloc_df <- bind_rows(coloc_rows)
  fwrite(coloc_df, file.path(out_dir, "coloc_results.tsv"), sep="\t")
  
  mr_out <- file.path(out_dir, "mr_results.tsv")
  if (length(mr_rows) > 0) {
    mr_df <- bind_rows(mr_rows)
  } else {
    mr_df <- data.table(method=character(), nsnp=integer(), b=numeric(), se=numeric(), pval=numeric(),
                        tag=character(), chr=integer(), start=integer(), end=integer(),
                        best_gene=character(), PP4=numeric(), gwas_eaf_imputed_frac=numeric(), n_iv=integer())
  }
  fwrite(mr_df, mr_out, sep="\t")
  
  message("Saved -> ", out_dir, " (coloc_n=", nrow(coloc_df), ", mr_n=", nrow(mr_df), ")")
  list(n=nrow(coloc_df), mr_n=nrow(mr_df), out_dir=out_dir)
}

# ========= run all =========
message("[RunMode] run_relaxed_always = ", cfg$run_relaxed_always)

for (tissue in names(gtex_clean_files)) {
  qtl_file <- gtex_clean_files[[tissue]]
  qtl_name <- paste0("GTExv7_", tissue)
  
  message("#############################")
  message("### QTL Tissue: ", qtl_name)
  message("#############################")
  
  qtl <- read_sumstats(qtl_file)
  
  for (tag in cfg$gwas_run_tags) {
    message(">>> RUN STRICT: ", tag, " vs ", qtl_name)
    run_tag_one_qtl(tag, qtl, qtl_name,
                    window_kb = cfg$strict_window_kb,
                    min_overlap_gene = cfg$strict_min_overlap_gene,
                    min_overlap_locus = cfg$strict_min_overlap_locus,
                    suffix = "strict")
    
    if (isTRUE(cfg$run_relaxed_always)) {
      message(">>> RUN RELAXED: ", tag, " vs ", qtl_name)
      run_tag_one_qtl(tag, qtl, qtl_name,
                      window_kb = cfg$relaxed_window_kb,
                      min_overlap_gene = cfg$relaxed_min_overlap_gene,
                      min_overlap_locus = cfg$relaxed_min_overlap_locus,
                      suffix = "relaxed")
    }
  }
  
  gc(verbose=FALSE)
}