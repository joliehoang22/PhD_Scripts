library(dplyr)
library(data.table)

## -------------------------
## 0. Read BrainMeta (reference)
## -------------------------
brainmeta_hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz")
setDT(brainmeta_hg38)

## If your BrainMeta file uses snp_key_hg38, keep that name and use it in merge
## (Probe = gene ID column in BrainMeta)
## Ensure BrainMeta has at least: Probe, snp_key_hg38, b, SE, p

## -------------------------
## 1. Helper: preprocess a QTLtools nominal file (LIV / PM)
## -------------------------
prep_qtl <- function(path) {
  df <- fread(path, header = FALSE, sep = " ", data.table = FALSE)
  colnames(df) <- c(
    "gene_id", "gene_chr", "gene_start", "gene_end",
    "strand", "n_variants", "distance", "variant_id",
    "var_chr", "var_start", "var_end",
    "pval", "beta", "top_flag"
  )
  
  # strip gene version
  df$gene_id <- sub("\\..*$", "", df$gene_id)
  
  # extract chr:pos from variant_id (e.g., chr1:12345[A/G])
  df$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", df$variant_id))
  df$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", df$variant_id))
  df$snp_key     <- paste0("chr", df$variant_chr, ":", df$variant_pos)
  
  # z and SE from p and beta (two-sided p)
  df$z <- ifelse(
    df$pval > 0 & df$pval < 1,
    qnorm(1 - df$pval / 2),
    ifelse(df$pval == 0, 38, 1e-6)
  )
  df$se <- abs(df$beta) / abs(df$z)
  
  df
}

## -------------------------
## 2. rb function (unchanged)
## -------------------------
calcu_cor_true <- function(b1, se1, b2, se2, theta) {
  idx = which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
  if (length(idx) > 0) {
    b1    = b1[-idx];    se1 = se1[-idx]
    b2    = b2[-idx];    se2 = se2[-idx]
    theta = theta[-idx]
  }
  
  var_b1 = var(b1, na.rm = TRUE) - mean(se1^2, na.rm = TRUE)
  var_b2 = var(b2, na.rm = TRUE) - mean(se2^2, na.rm = TRUE)
  
  if (var_b1 < 0) var_b1 = var(b1, na.rm = TRUE)
  if (var_b2 < 0) var_b2 = var(b2, na.rm = TRUE)
  
  cov_b1_b2 = cov(b1, b2, use = "complete.obs") -
    mean(theta, na.rm = TRUE) * sqrt(mean(se1^2, na.rm = TRUE) * mean(se2^2, na.rm = TRUE))
  
  r = cov_b1_b2 / sqrt(var_b1 * var_b2)
  
  # jackknife
  r_jack = c()
  n = length(b1)
  for (k in 1:n) {
    b1_jack = b1[-k]; se1_jack = se1[-k]
    b2_jack = b2[-k]; se2_jack = se2[-k]
    
    var_b1_jack = var(b1_jack, na.rm = TRUE) - mean(se1_jack^2, na.rm = TRUE)
    var_b2_jack = var(b2_jack, na.rm = TRUE) - mean(se2_jack^2, na.rm = TRUE)
    
    if (var_b1_jack < 0) var_b1_jack = var(b1_jack, na.rm = TRUE)
    if (var_b2_jack < 0) var_b2_jack = var(b2_jack, na.rm = TRUE)
    
    theta_jack = theta[-k]
    
    cov_e1_jack_e2_jack = mean(theta_jack, na.rm = TRUE) *
      sqrt(mean(se1_jack^2, na.rm = TRUE) * mean(se2_jack^2, na.rm = TRUE))
    
    cov_b1_b2_jack = cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack
    
    r_tmp = cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack = c(r_jack, r_tmp)
  }
  
  r_mean = mean(r_jack, na.rm = TRUE)
  idx_na = which(is.na(r_jack))
  
  if (length(idx_na) > 0) {
    se_r = sqrt((n - 1) / n * sum((r_jack[-idx_na] - r_mean)^2))
  } else {
    se_r = sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }
  
  cbind(r, se_r)
}

## -------------------------
## 3. Compute rb for one tissue using BrainMeta as reference
## -------------------------
compute_rb_for_tissue <- function(brainmeta, tissue_df, p_thresh = 5e-8, theta = 0) {
  # 1) merge on SNP + gene
  dat <- merge(
    brainmeta,
    tissue_df,
    by.x = c("snp_key_hg38", "Probe"),  # BrainMeta columns
    by.y = c("snp_key",      "gene_id"),
    all  = FALSE
  )
  
  if (nrow(dat) == 0) {
    return(list(
      rb      = NA_real_,
      se_rb   = NA_real_,
      n_pairs = 0,
      n_ref   = 0,
      n_total = 0,
      note    = "No overlapping SNP-gene pairs"
    ))
  }
  
  # 2) require non-zero SEs in both
  dat <- dat[dat$SE != 0 & dat$se != 0 & !is.na(dat$SE) & !is.na(dat$se), ]
  if (nrow(dat) == 0) {
    return(list(
      rb      = NA_real_,
      se_rb   = NA_real_,
      n_pairs = 0,
      n_ref   = 0,
      n_total = 0,
      note    = "No non-zero SEs after merge"
    ))
  }
  
  # 3) BrainMeta-only discovery: p_brainmeta < p_thresh
  ref_idx <- which(dat$p < p_thresh & !is.na(dat$p))
  if (length(ref_idx) < 2) {
    return(list(
      rb      = NA_real_,
      se_rb   = NA_real_,
      n_pairs = length(ref_idx),
      n_ref   = length(ref_idx),
      n_total = nrow(dat),
      note    = "Too few BrainMeta discoveries"
    ))
  }
  
  dat_ref <- dat[ref_idx, ]
  
  # 4) compute rb using BrainMeta b/SE and tissue beta/se
  theta_vec <- rep(theta, nrow(dat_ref))
  
  rb_res <- calcu_cor_true(
    b1    = dat_ref$b,     # BrainMeta effect
    se1   = dat_ref$SE,    # BrainMeta SE
    b2    = dat_ref$beta,  # tissue effect
    se2   = dat_ref$se,    # tissue SE
    theta = theta_vec
  )
  
  list(
    rb      = as.numeric(rb_res[1]),
    se_rb   = as.numeric(rb_res[2]),
    n_pairs = nrow(dat_ref),
    n_ref   = length(ref_idx),
    n_total = nrow(dat),
    note    = NA_character_
  )
}

## -------------------------
## 4. nsv handling via command line args (FIXED)
## -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  # Each argument is treated as an nsv value
  nsv_vals <- as.integer(args)
} else {
  # default: full range if not specified
  nsv_vals <- 1:30
}

base_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA"

results_list <- list()

for (nsv in nsv_vals) {
  message("Processing nsv = ", nsv)
  
  ## paths
  liv_path <- file.path(base_dir, sprintf("LIV_nominal_nsv%d.txt", nsv))
  pm_path  <- file.path(base_dir, sprintf("PM_nominal_nsv%d.txt",  nsv))
  
  ## read + preprocess LIV/PM
  liv_df <- prep_qtl(liv_path)
  pm_df  <- prep_qtl(pm_path)
  
  ## compute rb for each
  rb_liv <- compute_rb_for_tissue(brainmeta_hg38, liv_df, theta = 0)
  rb_pm  <- compute_rb_for_tissue(brainmeta_hg38, pm_df,  theta = 0)
  
  ## store as rows
  results_list[[length(results_list) + 1]] <- data.frame(
    tissue  = "LIV",
    nsv     = nsv,
    rb      = rb_liv$rb,
    se_rb   = rb_liv$se_rb,
    n_pairs = rb_liv$n_pairs,
    n_ref   = rb_liv$n_ref,
    n_total = rb_liv$n_total,
    note    = rb_liv$note,
    stringsAsFactors = FALSE
  )
  
  results_list[[length(results_list) + 1]] <- data.frame(
    tissue  = "PM",
    nsv     = nsv,
    rb      = rb_pm$rb,
    se_rb   = rb_pm$se_rb,
    n_pairs = rb_pm$n_pairs,
    n_ref   = rb_pm$n_ref,
    n_total = rb_pm$n_total,
    note    = rb_pm$note,
    stringsAsFactors = FALSE
  )
}

rb_summary <- bind_rows(results_list)

out_path <- sprintf(
  "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv%dto%d_20251211.csv",
  min(nsv_vals), max(nsv_vals)
)

write.csv(rb_summary, out_path, row.names = FALSE)
message("Saved results to: ", out_path)


##############################
##############################
##############################
# Pick one NSV file manually
liv_df <- prep_qtl("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/LIV_nominal_nsv1.txt")
pm_df  <- prep_qtl("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/PM_nominal_nsv1.txt")

rb_liv <- compute_rb_for_tissue(brainmeta_hg38, liv_df)
rb_pm  <- compute_rb_for_tissue(brainmeta_hg38, pm_df)

rb_liv
rb_pm


dat_test <- merge(
  brainmeta_hg38,
  liv_df,
  by.x = c("snp_key_hg38", "Probe"),
  by.y = c("snp_key", "gene_id")
)
head(dat_test)

nrow(dat_test)


library(parallel)

compute_rb_parallel <- function(brainmeta, df, p_thresh = 5e-8, theta = 0, ncores = 30) {

  dat <- merge(
    brainmeta,
    df,
    by.x = c("snp_key_hg38", "Probe"),
    by.y = c("snp_key", "gene_id"),
    all = FALSE
  )

  dat <- dat[dat$SE != 0 & dat$se != 0, ]
  ref_idx <- which(dat$p < p_thresh)

  dat_ref <- dat[ref_idx, ]
  theta_vec <- rep(theta, nrow(dat_ref))

  b1 <- dat_ref$b
  se1 <- dat_ref$SE
  b2 <- dat_ref$beta
  se2 <- dat_ref$se

  ## ------------------------------
  ## Parallel jackknife
  ## ------------------------------
  n <- length(b1)

  r_jack <- mclapply(
    1:n,
    function(k) {
      b1_j <- b1[-k]
      se1_j <- se1[-k]
      b2_j <- b2[-k]
      se2_j <- se2[-k]
      th_j <- theta_vec[-k]

      var_b1_j <- var(b1_j) - mean(se1_j^2)
      var_b2_j <- var(b2_j) - mean(se2_j^2)
      if (var_b1_j < 0) var_b1_j <- var(b1_j)
      if (var_b2_j < 0) var_b2_j <- var(b2_j)

      cov_j <- cov(b1_j, b2_j) - mean(th_j) * sqrt(mean(se1_j^2) * mean(se2_j^2))
      cov_j / sqrt(var_b1_j * var_b2_j)
    },
    mc.cores = ncores
  )

  r_jack <- unlist(r_jack)
  r_mean <- mean(r_jack, na.rm = TRUE)
  se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2))

  list(rb = r_mean, se_rb = se_r, n_pairs = length(b1))
}

rb_liv <- compute_rb_parallel(brainmeta_hg38, liv_df, ncores = 30)
rb_pm  <- compute_rb_parallel(brainmeta_hg38, pm_df,  ncores = 30)

rb_liv
rb_pm