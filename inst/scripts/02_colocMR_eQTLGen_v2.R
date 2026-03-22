# =========================
# coloc_pipline_v1.R  (ONE-SHOT: supports GWAS with/without EAF + MR enabled; FORCE MR optional)
# LENIENT VERSION: expand candidate pool + always run strict & relaxed + export QTL to output/coloc_v1
# =========================
getwd()
.libPaths(c("E:/R-4.4.0/librarygwas", setdiff(.libPaths(), "E:/R-4.4.0/librarygwas")))
.libPaths()

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
  
  eqtlgen_af_file  = "input/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
  eqtlgen_cis_file = "input/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
  
  qtl_file = "input/qtl/eQTLGen.clean.subset.tsv",
  qtl_name = "eQTLGen_blood",
  
  gwas_tags = c("LAC_discovery", "WMH_discovery"),
  
  # A 模式：真正要跑的（含 validation）
  gwas_run_tags = c("LAC_discovery", "WMH_discovery", "LAC_validation", "WMH_validation"),
  
  # A 模式：validation 用 discovery 的 loci 窗口
  loci_source_map = list(
    LAC_validation = "LAC_discovery",
    WMH_validation = "WMH_discovery"
  ),
  
  # QTL 预处理：用最大窗口一次覆盖
  prep_window_kb = 1500,
  
  # ---- 更宽松的参数（扩展候选池）----
  strict_window_kb = 1000,
  strict_min_overlap_gene = 3,
  strict_min_overlap_locus = 5,
  
  relaxed_window_kb = 1500,
  relaxed_min_overlap_gene = 2,
  relaxed_min_overlap_locus = 3,
  
  # 不再 fallback：每个 tag 都跑 strict + relaxed
  run_relaxed_always = TRUE,
  
  # 候选池输出：每个 locus 输出 topK 基因（同时也保留全量 gene-level 表）
  topk_gene_per_locus = 5,
  
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
  pp4_for_mr = 0.8,       # 你要更宽松也可改成 0.75
  mr_pval_iv = 5e-8,
  mr_run_without_pp4_gate = FALSE,
  
  # ---- EAF QC & auto-fix (GWAS) ----
  # If GWAS-provided EAF is actually MAF / not aligned to EA, MR harmonisation can be corrupted.
  # This QC will (optionally) drop bad GWAS EAF/MAF and re-impute from eqtlgen AF ref by EA.
  eaf_qc_enable = TRUE,
  eaf_qc_n_sample = 200000L,
  eaf_qc_tol = 0.02,
  eaf_qc_mismatch_thr = 0.30,
  eaf_qc_eq_thr = 0.95,
  export_eaf_qc = TRUE,

  eqtlgen_N_default = 31684,
  out_base = "output"
)

dir.create(dirname(cfg$qtl_file), recursive = TRUE, showWarnings = FALSE)

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
    P    = c("P","p","PVAL","pval","pvalue","Pvalue","p_value"),
    EAF  = c("EAF","eaf","AF","af","freq","effect_allele_freq","ALT_AF","alt_af"),
    MAF  = c("MAF","maf"),
    N    = c("N","n","samplesize","sample_size","N_total","Neff"),
    GENE = c("GENE","gene","gene_id","GENE_ID")
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
  
  # NEA 缺失时：用 A1/A2 推断
  if (!("NEA" %in% names(dt))) {
    a1 <- find_col(names(dt), c("A1","allele1","AlleleA"))
    a2 <- find_col(names(dt), c("A2","allele2","AlleleB"))
    if (!is.na(a1) && !is.na(a2) && ("EA" %in% names(dt))) {
      dt[, NEA := fifelse(toupper(get(a1)) == toupper(EA), toupper(get(a2)),
                          fifelse(toupper(get(a2)) == toupper(EA), toupper(get(a1)), NA_character_))]
    }
  }
  
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
  
  # MAF 优先：已有 MAF > EAF 推导
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
  
  # 防止 NA 导致奇怪 subset
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

run_coloc_abf <- function(t1, t2, meta1, meta2, priors) {
  t1 <- t1[order(P)][!duplicated(SNP)]
  t2 <- t2[order(P)][!duplicated(SNP)]
  t1 <- t1[is.finite(BETA) & is.finite(varbeta) & is.finite(MAF) & MAF > 0 & MAF < 1]
  t2 <- t2[is.finite(BETA) & is.finite(varbeta) & is.finite(MAF) & MAF > 0 & MAF < 1]
  if (nrow(t1) < 3 || nrow(t2) < 3) return(NULL)
  
  N1 <- median(t1$N, na.rm=TRUE)
  N2 <- median(t2$N, na.rm=TRUE)
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

# ========= AF ref =========
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

# ========= fill EAF/MAF =========
fill_missing_eaf <- function(gwas_dt, af_ref) {
  if (!("EAF" %in% names(gwas_dt))) gwas_dt[, EAF := NA_real_]
  if (!("MAF" %in% names(gwas_dt))) gwas_dt[, MAF := NA_real_]
  if (!("EAF_SRC" %in% names(gwas_dt))) gwas_dt[, EAF_SRC := "missing"]
  
  gwas_dt[, `:=`(EAF=as.numeric(EAF), MAF=as.numeric(MAF))]
  
  # 清洗后重新标注来源（避免非法EAF变NA但EAF_SRC没变导致无法impute）
  gwas_dt[, EAF := fifelse(is.finite(EAF) & EAF > 0 & EAF < 1, EAF, NA_real_)]
  gwas_dt[, MAF := fifelse(is.finite(MAF) & MAF > 0 & MAF < 1, MAF, NA_real_)]
  gwas_dt[, EAF_SRC := fifelse(!is.na(EAF), "gwas", "missing")]
  
  setkey(gwas_dt, SNP)
  setkey(af_ref, SNP)
  
  # 1) 先补 MAF（不依赖 EA）
  gwas_dt[af_ref, on="SNP",
          MAF := fifelse(is.na(MAF), pmin(i.FREQ_B, 1 - i.FREQ_B), MAF)]
  
  # 2) 再补 EAF（依赖 EA 与 A/B 匹配）
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
  
  # 若仍没有 MAF，用 EAF 推一下
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


# ========= EAF QC & GWAS preparation cache =========
# One-shot preparation per GWAS tag to avoid double-reading for strict/relaxed.
.GWAS_CACHE <- new.env(parent = emptyenv())
.GWAS_QC    <- new.env(parent = emptyenv())

.eaf_stats <- function(dt) {
  eaf_ok <- is.finite(dt$EAF)
  maf_ok <- is.finite(dt$MAF)
  out <- list(
    n = nrow(dt),
    eaf_missing_rate = mean(!eaf_ok),
    eaf_min = if (any(eaf_ok)) min(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    eaf_median = if (any(eaf_ok)) median(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    eaf_max = if (any(eaf_ok)) max(dt$EAF[eaf_ok], na.rm=TRUE) else NA_real_,
    prop_eaf_gt_0_5 = if (any(eaf_ok)) mean(dt$EAF[eaf_ok] > 0.5, na.rm=TRUE) else NA_real_,
    prop_eaf_eq_maf = {
      ok <- eaf_ok & maf_ok
      if (any(ok)) mean(abs(dt$EAF[ok] - dt$MAF[ok]) < 1e-6) else NA_real_
    }
  )
  if ("EAF_SRC" %in% names(dt)) {
    out$prop_eaf_imputed <- mean(dt$EAF_SRC == "imputed_eqtlgenAF", na.rm=TRUE)
  } else {
    out$prop_eaf_imputed <- NA_real_
  }
  out
}

.check_gwas_eaf_vs_ref <- function(gwas_dt, af_ref, n_sample=200000L, tol=0.02) {
  # Returns mismatch stats for GWAS-provided EAF vs expected EAF from AF ref by EA.
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
  if (nrow(tmp) < 10000L) return(list(mismatch=NA_real_, mismatch_flip=NA_real_, matched_n=nrow(tmp)))

  mismatch <- mean(abs(tmp$EAF - tmp$EAF_expected) > tol, na.rm=TRUE)
  mismatch_flip <- mean(abs(tmp$EAF - (1 - tmp$EAF_expected)) > tol, na.rm=TRUE)
  list(mismatch=mismatch, mismatch_flip=mismatch_flip, matched_n=nrow(tmp))
}

.qc_and_fix_gwas_eaf <- function(gwas_dt, af_ref, tag, cfg) {
  # If GWAS EAF looks like MAF and is not aligned to EA (vs AF ref), drop it and re-impute later.
  if (!isTRUE(cfg$eaf_qc_enable)) return(list(gwas=gwas_dt, qc=NULL))

  if (!("EAF" %in% names(gwas_dt))) gwas_dt[, EAF := NA_real_]
  if (!("MAF" %in% names(gwas_dt))) gwas_dt[, MAF := NA_real_]
  gwas_dt[, `:=`(EAF=as.numeric(EAF), MAF=as.numeric(MAF))]
  gwas_dt[!(is.finite(EAF) & EAF > 0 & EAF < 1), EAF := NA_real_]
  gwas_dt[!(is.finite(MAF) & MAF > 0 & MAF < 1), MAF := NA_real_]
  gwas_dt[is.na(MAF) & is.finite(EAF), MAF := pmin(EAF, 1 - EAF)]

  st <- .eaf_stats(gwas_dt)
  max_eaf <- st$eaf_max
  eq_rate <- st$prop_eaf_eq_maf

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
    fixed_drop_eaf = FALSE,
    fixed_drop_maf = FALSE
  )

  # Only do deep mismatch check when EAF looks like MAF (<=0.5 and equals MAF nearly always).
  looks_like_maf <- is.finite(max_eaf) && (max_eaf <= 0.5 + 1e-12) && is.finite(eq_rate) && (eq_rate >= cfg$eaf_qc_eq_thr)
  if (isTRUE(looks_like_maf)) {
    chk <- .check_gwas_eaf_vs_ref(gwas_dt, af_ref, n_sample=cfg$eaf_qc_n_sample, tol=cfg$eaf_qc_tol)
    qc[, `:=`(mismatch = chk$mismatch, mismatch_if_flipped = chk$mismatch_flip, matched_n = as.integer(chk$matched_n))]

    if (is.finite(chk$mismatch) && chk$mismatch > cfg$eaf_qc_mismatch_thr) {
      message("[FixEAF][", tag, "] GWAS EAF behaves like MAF and mismatches AF ref (mismatch=", signif(chk$mismatch, 4),
              ", flip=", signif(chk$mismatch_flip, 4), ", matched=", chk$matched_n, "). ",
              "=> drop GWAS EAF and re-impute from eqtlgen AF by EA.")
      # If MAF equals EAF (most likely derived), drop MAF too to avoid propagating the same bad frequency.
      gwas_dt[, EAF := NA_real_]
      gwas_dt[, EAF_SRC := "missing"]
      qc[, fixed_drop_eaf := TRUE]

      if (is.finite(eq_rate) && eq_rate >= cfg$eaf_qc_eq_thr) {
        gwas_dt[, MAF := NA_real_]
        qc[, fixed_drop_maf := TRUE]
      }
    } else {
      message("[FixEAF][", tag, "] EAF looks like MAF but mismatch not high (mismatch=", signif(chk$mismatch, 4),
              "). Keep as-is.")
    }
  }

  list(gwas=gwas_dt, qc=qc)
}

.get_prepared_gwas <- function(tag) {
  if (exists(tag, envir=.GWAS_CACHE, inherits=FALSE)) {
    return(get(tag, envir=.GWAS_CACHE, inherits=FALSE))
  }

  gwas_file <- file.path(cfg$gwas_clean_dir, paste0(tag, ".clean.tsv"))
  if (!file.exists(gwas_file)) stop("Missing GWAS file: ", gwas_file)

  gwas0 <- read_sumstats(gwas_file)
  qcres <- .qc_and_fix_gwas_eaf(gwas0, af_ref, tag, cfg)
  gwas1 <- qcres$gwas
  assign(tag, qcres$qc, envir=.GWAS_QC)

  # Impute EAF/MAF and N
  gwas1 <- fill_missing_eaf(gwas1, af_ref)
  gwas1 <- fill_missing_n(gwas1, cfg$trait_meta[[tag]]$N)

  # Filter for coloc (MAF and N are required); EAF can be NA.
  gwas1 <- gwas1[is.finite(MAF) & MAF > 0 & MAF < 1 & is.finite(N)]
  gwas1 <- gwas1[order(P)][!duplicated(SNP)]

  # After-imputation stats
  st2 <- .eaf_stats(gwas1)
  qc <- get(tag, envir=.GWAS_QC, inherits=FALSE)
  if (!is.null(qc) && nrow(qc) == 1) {
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

run_mr_one <- function(t1, t2, tag, chr, start, end, best_gene, pp4, eaf_imp_frac) {
  tryCatch({
    exp_dt <- t2[t2$P < cfg$mr_pval_iv, ]
    if (nrow(exp_dt) < 1) return(NULL)
    
    out_dt <- t1[t1$SNP %in% exp_dt$SNP, ]
    if (nrow(out_dt) < 1) return(NULL)
    
    exp <- prep_for_tsmr(exp_dt, paste0(cfg$qtl_name, ":", best_gene), TRUE)
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

# ========= 2) eQTLGen raw -> clean subset =========
pick_col <- function(nms, candidates) {
  hit <- intersect(tolower(nms), tolower(candidates))
  if (length(hit) == 0) return(NA_character_)
  nms[match(hit[1], tolower(nms))]
}

prepare_eqtlgen_clean <- function(cis_file, af_file, loci_tags, out_file, window_kb, N_default) {
  message("[Prep] Reading loci to build region filter...")
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
  setkey(merged, CHR, start, end)
  
  message("[Prep] Reading eQTLGen cis file header...")
  cis_header <- names(fread(cis_file, nrows=0))
  c_snp  <- pick_col(cis_header, c("SNPName","SNP","rsid","RSID"))
  c_gene <- pick_col(cis_header, c("ProbeName","GeneSymbol","HUGO","gene","gene_id"))
  c_p    <- pick_col(cis_header, c("Pvalue","PValue","pvalue","pval"))
  c_z    <- pick_col(cis_header, c("OverallZScore","Zscore","ZScore","z"))
  c_ea   <- pick_col(cis_header, c("AlleleAssessed","AssessedAllele","effect_allele"))
  c_n    <- pick_col(cis_header, c("N","Samplesize","SampleSize","n"))
  
  need <- c(c_snp,c_gene,c_p,c_z,c_ea)
  if (any(is.na(need))) stop("Cannot map columns in cis file.\nHeader:\n", paste(cis_header, collapse=", "))
  
  message("[Prep] Reading eQTLGen cis file (selected columns only)...")
  cis <- fread(cis_file, select = unique(na.omit(c(c_snp,c_gene,c_p,c_z,c_ea,c_n))))
  setnames(cis, c_snp, "SNP")
  setnames(cis, c_gene, "GENE")
  setnames(cis, c_p, "P")
  setnames(cis, c_z, "Z")
  setnames(cis, c_ea, "EA_raw")
  if (!is.na(c_n)) setnames(cis, c_n, "N_raw")
  
  cis[, EA := toupper(as.character(EA_raw))]
  cis[, P := as.numeric(P)]
  cis[, Z := as.numeric(Z)]
  cis[, GENE := as.character(GENE)]
  cis[, SNP := as.character(SNP)]
  cis[, N := if ("N_raw" %in% names(cis)) as.numeric(N_raw) else as.numeric(N_default)]
  cis[, "EA_raw" := NULL]
  if ("N_raw" %in% names(cis)) cis[, "N_raw" := NULL]
  
  message("[Prep] Reading eQTLGen AF file header...")
  af_header <- names(fread(af_file, nrows=0))
  
  a_snp   <- pick_col(af_header, c("SNPName","SNP","rsid","RSID"))
  a_chr   <- pick_col(af_header, c("SNPChr","CHR","chr","hg19_chr"))
  a_pos   <- pick_col(af_header, c("SNPChrPos","POS","pos","hg19_pos"))
  a_a1    <- pick_col(af_header, c("AlleleA","allelea"))
  a_a2    <- pick_col(af_header, c("AlleleB","alleleb"))
  a_freqB <- pick_col(af_header, c("AlleleB_all","alleleb_all"))
  a_allAB <- pick_col(af_header, c("allAB_total","allab_total"))
  a_allB  <- pick_col(af_header, c("allB_total","allb_total"))
  
  need2 <- c(a_snp,a_chr,a_pos,a_a1,a_a2)
  if (any(is.na(need2))) stop("Cannot map columns in AF file.\nHeader:\n", paste(af_header, collapse=", "))
  
  message("[Prep] Reading eQTLGen AF file (selected columns only)...")
  sel <- unique(na.omit(c(a_snp,a_chr,a_pos,a_a1,a_a2,a_freqB,a_allAB,a_allB)))
  af <- fread(af_file, select = sel)
  
  setnames(af, a_snp, "SNP")
  setnames(af, a_chr, "CHR")
  setnames(af, a_pos, "POS")
  setnames(af, a_a1, "A1")
  setnames(af, a_a2, "A2")
  if (!is.na(a_freqB)) setnames(af, a_freqB, "FREQ_B")
  if (!is.na(a_allAB)) setnames(af, a_allAB, "ALL_AB")
  if (!is.na(a_allB))  setnames(af, a_allB,  "ALL_B")
  
  af[, `:=`(
    SNP = as.character(SNP),
    CHR = as.integer(CHR),
    POS = as.integer(POS),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2))
  )]
  
  if ("FREQ_B" %in% names(af)) {
    af[, FREQ_B := as.numeric(FREQ_B)]
  } else if (all(c("ALL_B","ALL_AB") %in% names(af))) {
    af[, `:=`(ALL_B=as.numeric(ALL_B), ALL_AB=as.numeric(ALL_AB))]
    af[, FREQ_B := ALL_B / ALL_AB]
  } else {
    stop("AF file has no AlleleB_all and cannot derive allele frequency from counts.")
  }
  
  af <- af[!is.na(FREQ_B) & FREQ_B > 0 & FREQ_B < 1]
  af[, MAF := pmin(FREQ_B, 1 - FREQ_B)]
  af <- af[!is.na(MAF) & MAF > 0 & MAF < 1]
  
  message("[Prep] Joining cis with AF and filtering to loci regions...")
  setkey(cis, SNP); setkey(af, SNP)
  dt <- cis[af, nomatch=0L]
  
  dt[, NEA := fifelse(A1 == EA, A2, fifelse(A2 == EA, A1, NA_character_))]
  dt <- dt[!is.na(NEA)]
  
  v <- dt[, .(CHR=CHR, start=POS, end=POS, .row=.I)]
  setkey(v, CHR, start, end)
  ov <- foverlaps(v, merged, nomatch=0L)
  if (nrow(ov) == 0) stop("No eQTLGen SNPs overlap your loci regions. Check build consistency (hg19 vs hg38).")
  dt <- dt[ov$.row]
  
  dt[, SE := 1 / sqrt(2 * MAF * (1 - MAF) * N)]
  dt[, BETA := Z * SE]
  
  # EAF 必须是 effect allele (EA) 的频率
  dt[, EAF := fifelse(EA == A2, FREQ_B,
                      fifelse(EA == A1, 1 - FREQ_B, NA_real_))]
  
  dt <- dt[is.finite(EAF) & EAF > 0 & EAF < 1]
  dt <- dt[is.finite(MAF) & MAF > 0 & MAF < 1]
  
  out <- dt[, .(SNP, CHR, POS, EA, NEA, BETA, SE, P, EAF, MAF, N, GENE)]
  fwrite(out, out_file, sep="\t")
  message("[Prep] Saved clean QTL subset -> ", out_file)
}

# ===== QTL 预处理（自动重建）=====
need_qtl_cols <- c("SNP","CHR","POS","EA","NEA","BETA","SE","P","EAF","MAF","N","GENE")
rebuild_qtl <- FALSE
if (!file.exists(cfg$qtl_file)) {
  rebuild_qtl <- TRUE
} else {
  hdr <- names(fread(cfg$qtl_file, nrows=0))
  if (!all(need_qtl_cols %in% hdr)) {
    message("[Prep] Existing QTL file missing cols -> rebuild: ", cfg$qtl_file)
    rebuild_qtl <- TRUE
  } else {
    samp <- fread(cfg$qtl_file, nrows=20000, showProgress=FALSE)
    if (all(c("EAF","MAF") %in% names(samp))) {
      eq_rate <- mean(abs(as.numeric(samp$EAF) - as.numeric(samp$MAF)) < 1e-12, na.rm=TRUE)
      if (is.finite(eq_rate) && eq_rate > 0.995) {
        message(sprintf("[Prep] Detected possible old EAF bug (EAF≈MAF rate=%.3f) -> rebuild QTL", eq_rate))
        rebuild_qtl <- TRUE
      }
    }
  }
}

if (rebuild_qtl) {
  prepare_eqtlgen_clean(cfg$eqtlgen_cis_file, cfg$eqtlgen_af_file,
                        cfg$gwas_tags, cfg$qtl_file, cfg$prep_window_kb, cfg$eqtlgen_N_default)
} else {
  message("[Prep] Found existing QTL clean file (ok), skip: ", cfg$qtl_file)
}

# ========= 2.5) Export QTL used into output/coloc_v1/_inputs (HARD CHECK) =========
coloc_root <- file.path(cfg$out_base, "coloc_v1")
dir.create(coloc_root, recursive=TRUE, showWarnings=FALSE)

inputs_dir <- file.path(coloc_root, "_inputs")
dir.create(inputs_dir, recursive=TRUE, showWarnings=FALSE)

qtl_used_path <- file.path(inputs_dir, paste0(cfg$qtl_name, ".clean.subset.used.tsv"))
ok_copy <- file.copy(cfg$qtl_file, qtl_used_path, overwrite=TRUE)
if (!isTRUE(ok_copy)) {
  stop("[ExportQTL] file.copy failed: ", cfg$qtl_file, " -> ", qtl_used_path)
}
message("[ExportQTL] OK -> ", qtl_used_path)

# ========= 3) Main runner =========
qtl <- read_sumstats(cfg$qtl_file)

run_tag <- function(tag, window_kb, min_overlap_gene, min_overlap_locus, suffix) {
  
  gwas_file <- file.path(cfg$gwas_clean_dir, paste0(tag, ".clean.tsv"))
  if (!file.exists(gwas_file)) stop("Missing GWAS file: ", gwas_file)
  
  loci_tag  <- if (!is.null(cfg$loci_source_map) && tag %in% names(cfg$loci_source_map)) cfg$loci_source_map[[tag]] else tag
  loci_path <- choose_loci_path(cfg$loci_dir, loci_tag)
  
  message("== Running: ", tag, " vs ", cfg$qtl_name, " | loci: ", basename(loci_path),
          " | window_kb=", window_kb, " | min_gene=", min_overlap_gene)
  
    # ---- Prepare GWAS once (EAF QC + (re)impute EAF/MAF from AF ref + fill N) ----
  gwas <- .get_prepared_gwas(tag)
  qc0 <- .get_gwas_qc(tag)
  if (!is.null(qc0) && nrow(qc0)==1) {
    message(sprintf("[EAFQC][%s] raw_missing=%.4f raw_max=%.4f raw_eqMAF=%.4f mismatch=%.4f fixed_drop_eaf=%s | after_missing=%.4f after_prop_gt0.5=%.4f",
                    tag,
                    qc0$eaf_missing_rate_raw, qc0$eaf_max_raw, qc0$prop_eaf_eq_maf_raw,
                    qc0$mismatch, as.character(qc0$fixed_drop_eaf),
                    qc0$eaf_missing_rate_after, qc0$prop_eaf_gt_0_5_after))
  }

  
  regions <- infer_regions(loci_path, window_kb)
  
  out_dir <- file.path(cfg$out_base, "coloc_v1",
                       paste0(tag, "__vs__", cfg$qtl_name, "__", suffix))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  # 每个 out_dir 也保存一份 QTL used（HARD CHECK）
  qtl_copy2 <- file.path(out_dir, "qtl_used.tsv")
  ok2 <- file.copy(cfg$qtl_file, qtl_copy2, overwrite=TRUE)
  if (!isTRUE(ok2)) stop("[ExportQTL] file.copy failed: ", cfg$qtl_file, " -> ", qtl_copy2)
  message("[ExportQTL] OK -> ", qtl_copy2)
  

  # ---- Save EAF QC report (one row) ----
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
    
    r1 <- extract_region(gwas, chr, start, end)
    r2 <- extract_region(qtl,  chr, start, end)
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
      
      cres <- run_coloc_abf(t1, t2, meta1, meta2, cfg$coloc_priors)
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
                           eaf_imp_frac=eaf_imp_frac)
      if (!is.null(mr_res)) mr_rows[[length(mr_rows) + 1]] <- mr_res
    }
  }
  
  # ===== gene-level 输出（修复：先 setDT 再 :=）=====
  gene_res_path <- file.path(out_dir, "coloc_gene_results.tsv")
  gene_topk_path <- file.path(out_dir, "coloc_gene_topk.tsv")
  
  if (length(gene_rows) > 0) {
    gene_df <- bind_rows(gene_rows)
    fwrite(gene_df, gene_res_path, sep="\t")
    
    setDT(gene_df)  # ✅ 关键修复
    gene_df[, locus_id := paste0(tag, ":", chr, ":", start, "-", end)]
    setorder(gene_df, -PP4)
    
    gene_topk <- gene_df[, head(.SD, cfg$topk_gene_per_locus), by=.(locus_id)]
    fwrite(gene_topk, gene_topk_path, sep="\t")
  } else {
    fwrite(data.table(), gene_res_path, sep="\t")
    fwrite(data.table(), gene_topk_path, sep="\t")
  }
  
  if (length(coloc_rows) == 0) {
    message("[Warn] No loci passed for ", tag, " in ", suffix)
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

for (tag in cfg$gwas_run_tags) {
  message(">>> RUN STRICT: ", tag)
  run_tag(tag,
          window_kb = cfg$strict_window_kb,
          min_overlap_gene = cfg$strict_min_overlap_gene,
          min_overlap_locus = cfg$strict_min_overlap_locus,
          suffix = "strict")
  
  if (isTRUE(cfg$run_relaxed_always)) {
    message(">>> RUN RELAXED: ", tag)
    run_tag(tag,
            window_kb = cfg$relaxed_window_kb,
            min_overlap_gene = cfg$relaxed_min_overlap_gene,
            min_overlap_locus = cfg$relaxed_min_overlap_locus,
            suffix = "relaxed")
  }
}