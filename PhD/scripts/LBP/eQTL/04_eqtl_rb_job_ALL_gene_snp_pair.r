###############################################################
# Script Name: 04_eqtl_rb_job_ALL_gene_snp_pair.r
# Author: Jolie Hoang
# Date: 2025-12-02
#
# Purpose:
#   This script computes rb (correlation of true effect sizes)
#   between BrainMeta cis-eQTL effects and LIV/PM QTLtools
#   nominal effect sizes across a range of surrogate variable
#   (nsv) settings.
#
#   It implements two analysis modes:
#
#     1) mode = "bm_sig"
#        - Uses only SNP–gene pairs that are genome-wide 
#          significant in BrainMeta (p < 5e-8).
#        - Extracts their beta/SE estimates from LIV and PM
#          regardless of significance in those tissues.
#        - This avoids LIV/PM low power issues and reflects
#          the BrainMeta "discovery → replication" design.
#
#     2) mode = "bm_all"
#        - Uses all overlapping SNP–gene pairs between
#          BrainMeta and LIV/PM (no significance filtering).
#        - Provides a more exploratory genome-wide rb estimate.
#
#
# Inputs:
#   <mode>       - "bm_sig" or "bm_all"
#   <nsv_start>  - beginning of the nsv range (integer)
#   <nsv_end>    - end of the nsv range (integer)
#   <out_csv>    - output file path for rb summary table
#
# Dependencies:
#   - BrainMeta summary file in hg38 coordinate system:
#       brainmeta_hg38.txt.gz
#   - QTLtools nominal output for LIV and PM at each nsv:
#       LIV_nominal_nsvX.txt, PM_nominal_nsvX.txt
#   - R/4.2.0 + data.table, dplyr
#
# Output:
#   CSV file summarizing rb estimates for both LIV and PM:
#     columns: mode, tissue, nsv, rb, se_rb, n_pairs, note
# example titles: rb_task2_bm_sig_nsv1to3_20251202.csv OR rb_task2_bm_all_nsv1to3_20251202.csv
#
# Usage Example:
#   Rscript run_rb_task2.R bm_sig 1 3 rb_task2_bm_sig_nsv1to3.csv
#
# This script is intended to be called via bsub jobs that loop
# over nsv chunks (e.g., 1–3, 4–6, …, 28–30).
###############################################################
library(data.table)
library(dplyr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_rb_task2.R <mode> <nsv_start> <nsv_end> <out_csv>\n  mode = bm_sig or bm_all")
}

mode      <- args[1]                # "bm_sig" or "bm_all"
nsv_start <- as.integer(args[2])
nsv_end   <- as.integer(args[3])
out_csv   <- args[4]

message("Mode: ", mode)
message("nsv range: ", nsv_start, " to ", nsv_end)
message("Output: ", out_csv)

p_thresh <- 5e-8
base_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA"

#### 1. Load BrainMeta (hg38) ####
brainmeta_hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz")
setDT(brainmeta_hg38)
brainmeta_hg38[, Probe := sub("\\..*$", "", Probe)]  # strip version

brainmeta_sig <- brainmeta_hg38[p < p_thresh]

if (mode == "bm_sig") {
  brainmeta_use <- brainmeta_sig
} else if (mode == "bm_all") {
  brainmeta_use <- brainmeta_hg38
} else {
  stop("mode must be 'bm_sig' or 'bm_all'")
}

#### 2. QTL preprocessing ####
prep_qtl <- function(path) {
  df <- fread(path, header = FALSE, sep = " ", data.table = FALSE)
  colnames(df) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                    "strand", "n_variants", "distance", "variant_id",
                    "var_chr", "var_start", "var_end",
                    "pval", "beta", "top_flag")
  
  df$gene_id <- sub("\\..*$", "", df$gene_id)
  df$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", df$variant_id))
  df$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", df$variant_id))
  df$snp_key     <- paste0("chr", df$variant_chr, ":", df$variant_pos)
  
  df$z <- ifelse(
    df$pval > 0 & df$pval < 1,
    qnorm(1 - df$pval / 2),
    ifelse(df$pval == 0, 38, 1e-6)
  )
  df$se <- abs(df$beta) / abs(df$z)
  
  setDT(df)
  df
}

#### 3. rb internals ####
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
  
  n <- length(b1)
  r_jack <- numeric(n)
  for (k in seq_len(n)) {
    b1_j  <- b1[-k]; se1_j <- se1[-k]
    b2_j  <- b2[-k]; se2_j <- se2[-k]
    theta_j <- theta[-k]
    
    var_b1_j = var(b1_j, na.rm = TRUE) - mean(se1_j^2, na.rm = TRUE)
    var_b2_j = var(b2_j, na.rm = TRUE) - mean(se2_j^2, na.rm = TRUE)
    if (var_b1_j < 0) var_b1_j = var(b1_j, na.rm = TRUE)
    if (var_b2_j < 0) var_b2_j = var(b2_j, na.rm = TRUE)
    
    cov_e1e2_j = mean(theta_j, na.rm = TRUE) *
      sqrt(mean(se1_j^2, na.rm = TRUE) * mean(se2_j^2, na.rm = TRUE))
    
    cov_b1b2_j = cov(b1_j, b2_j, use = "complete.obs") - cov_e1e2_j
    
    r_jack[k] = cov_b1b2_j / sqrt(var_b1_j * var_b2_j)
  }
  
  r_mean <- mean(r_jack, na.rm = TRUE)
  se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2, na.rm = TRUE))
  
  cbind(r, se_r)
}

# BrainMeta-driven rb: BM filter depends on 'mode'
compute_rb_for_tissue_bm <- function(brainmeta_use, tissue_df, mode, p_thresh = 5e-8, theta = 0) {
  dat <- merge(
    brainmeta_use[, .(snp_key, Probe, b, SE, p)],
    tissue_df[, .(snp_key, gene_id, beta, se, pval)],
    by.x = c("snp_key", "Probe"),
    by.y = c("snp_key", "gene_id"),
    all  = FALSE
  )
  
  dat <- dat[SE != 0 & se != 0]
  
  if (mode == "bm_sig") {
    dat_use <- dat[p < p_thresh]
  } else if (mode == "bm_all") {
    dat_use <- dat
  } else {
    stop("mode must be 'bm_sig' or 'bm_all'")
  }
  
  if (nrow(dat_use) < 10) {
    return(list(
      rb = NA_real_,
      se_rb = NA_real_,
      n_pairs = nrow(dat_use),
      note = "Too few pairs after filtering"
    ))
  }
  
  theta_vec <- rep(theta, nrow(dat_use))
  
  rb_res <- calcu_cor_true(
    b1    = dat_use$b,
    se1   = dat_use$SE,
    b2    = dat_use$beta,
    se2   = dat_use$se,
    theta = theta_vec
  )
  
  list(
    rb      = as.numeric(rb_res[1]),
    se_rb   = as.numeric(rb_res[2]),
    n_pairs = nrow(dat_use),
    note    = NA_character_
  )
}

#### 4. Loop over nsv, do LIV and PM in one shot ####
nsv_vals <- nsv_start:nsv_end

res_list <- mclapply(
  nsv_vals,
  function(nsv) {
    message("Processing nsv = ", nsv)

    ## file paths
    liv_path <- file.path(base_dir, sprintf("LIV_nominal_nsv%d.txt", nsv))
    pm_path  <- file.path(base_dir, sprintf("PM_nominal_nsv%d.txt",  nsv))

    ## preprocess LIV + PM
    liv_df <- prep_qtl(liv_path)
    pm_df  <- prep_qtl(pm_path)

    ## compute rb using the BrainMeta selection for this mode
    rb_liv <- compute_rb_for_tissue_bm(brainmeta_use, liv_df, mode = mode, theta = 0)
    rb_pm  <- compute_rb_for_tissue_bm(brainmeta_use, pm_df,  mode = mode, theta = 0)

    ## return two rows (LIV + PM)
    data.frame(
      mode    = c(mode, mode),
      tissue  = c("LIV", "PM"),
      nsv     = c(nsv, nsv),
      rb      = c(rb_liv$rb, rb_pm$rb),
      se_rb   = c(rb_liv$se_rb, rb_pm$se_rb),
      n_pairs = c(rb_liv$n_pairs, rb_pm$n_pairs),
      note    = c(rb_liv$note, rb_pm$note),
      stringsAsFactors = FALSE
    )
  },
  mc.cores = 30   
)

rb_summary <- bind_rows(res_list)
write.csv(rb_summary, out_csv, row.names = FALSE)


