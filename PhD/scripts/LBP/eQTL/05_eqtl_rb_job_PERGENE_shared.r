###############################################################################
# Script: 05_eqtl_rb_job_PERGENE_shared.r
#
# Purpose:
#   This script computes the cross-tissue effect-size correlation (rb) between
#   BrainMeta (discovery) cis-eQTL effect sizes and eQTL effects from LIV and PM
#   samples across a range of surrogate variable (nsv) models.
#
#   This corresponds to Task 3 in the analysis plan:
#     → Use only SNP–gene pairs that are shared (overlapping) between BrainMeta
#       and the target tissue (LIV or PM).
#     → For each gene, retain exactly ONE SNP–gene pair by selecting the SNP with
#       the smallest BrainMeta p-value within the shared set ("lead SNP per gene").
#     → Compute rb(BrainMeta, LIV) and rb(BrainMeta, PM) using these lead pairs.
#
#   Two analysis modes are supported:
#       1. mode = "bm_sig"
#          - Filter to BrainMeta genome-wide significant eQTLs (p < 5e-8)
#            *before* selecting the lead SNP per gene.
#
#       2. mode = "bm_all"
#          - Use all overlapping SNP–gene pairs between BrainMeta and LIV/PM,
#            without any p-value restriction, prior to selecting the lead SNP.
#
#   For each nsv in the specified range, the script:
#       (a) Loads the QTLtools nominal results for LIV and PM
#       (b) Harmonizes SNP identifiers
#       (c) Merges with BrainMeta (hg38)
#       (d) Filters depending on analysis mode
#       (e) Selects one SNP–gene pair per gene (lead SNP)
#       (f) Computes rb and jackknife SE
#       (g) Saves results as a CSV
#
# Inputs:
#   Command-line arguments:
#       <mode>       = "bm_sig" or "bm_all"
#       <nsv_start>  = first nsv value in chunk
#       <nsv_end>    = last nsv value in chunk
#       <out_csv>    = full path to output CSV
#
# Parallelization:
#   - Uses mclapply() with 30 cores to process each nsv value in parallel.
#
# Output:
#   - A CSV with rb results for each tissue (LIV, PM) and each nsv in the range.
#   - rb_task3_pergene_bm_sig_nsv1to3_20251202.csv OR rb_task3_pergene_bm_all_nsv1to3_20251202.csv
#
###############################################################################

library(data.table)
library(dplyr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 05_eqtl_rb_job_PERGENE_shared.r <mode> <nsv_start> <nsv_end> <out_csv>\n  mode = bm_sig or bm_all")
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

## use all cores data.table can
setDTthreads(30)

#### 1. Load BrainMeta (hg38) ####
brainmeta_hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz")
setDT(brainmeta_hg38)
brainmeta_hg38[, Probe := sub("\\..*$", "", Probe)]  # strip ENS version

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

## One SNP–gene per gene, based on shared BrainMeta–tissue overlap
## mode = "bm_sig": restrict to BrainMeta p < p_thresh before picking lead SNP
## mode = "bm_all": use all overlapping pairs, still picking lead SNP per gene
compute_rb_pergene_shared <- function(brainmeta, tissue_df, mode, p_thresh = 5e-8, theta = 0) {
  dat_full <- merge(
    brainmeta[, .(snp_key, Probe, b, SE, p)],
    tissue_df[, .(snp_key, gene_id, beta, se, pval)],
    by.x = c("snp_key", "Probe"),
    by.y = c("snp_key", "gene_id"),
    all  = FALSE
  )
  
  dat_full <- dat_full[SE != 0 & se != 0]
  if (nrow(dat_full) == 0) {
    return(list(
      rb = NA_real_,
      se_rb = NA_real_,
      n_pairs = 0L,
      note = "No overlap after merge"
    ))
  }
  
  if (mode == "bm_sig") {
    dat_use <- dat_full[p < p_thresh]
  } else if (mode == "bm_all") {
    dat_use <- dat_full
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
  
  ## One SNP–gene pair per gene: pick SNP with smallest BrainMeta p per Probe
  dat_lead <- dat_use[ , .SD[which.min(p)], by = Probe]
  
  theta_vec <- rep(theta, nrow(dat_lead))
  rb_res <- calcu_cor_true(
    b1    = dat_lead$b,
    se1   = dat_lead$SE,
    b2    = dat_lead$beta,
    se2   = dat_lead$se,
    theta = theta_vec
  )
  
  list(
    rb      = as.numeric(rb_res[1]),
    se_rb   = as.numeric(rb_res[2]),
    n_pairs = nrow(dat_lead),
    note    = NA_character_
  )
}

#### 4. Parallel loop over nsv (LIV + PM inside each worker) ####
nsv_vals <- nsv_start:nsv_end

res_list <- mclapply(
  nsv_vals,
  function(nsv) {
    message("Processing nsv = ", nsv, " (mode = ", mode, ")")
    
    ## paths
    liv_path <- file.path(base_dir, sprintf("LIV_nominal_nsv%d.txt", nsv))
    pm_path  <- file.path(base_dir, sprintf("PM_nominal_nsv%d.txt",  nsv))
    
    liv_df <- prep_qtl(liv_path)
    pm_df  <- prep_qtl(pm_path)
    
    rb_liv <- compute_rb_pergene_shared(brainmeta_hg38, liv_df, mode = mode, theta = 0)
    rb_pm  <- compute_rb_pergene_shared(brainmeta_hg38, pm_df,  mode = mode, theta = 0)
    
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
message("Done. Written: ", out_csv)
