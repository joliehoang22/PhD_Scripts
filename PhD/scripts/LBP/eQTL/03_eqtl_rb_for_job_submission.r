### code for job submission - this is for lead snp-gene pair per gene in the shared dataset (of brainmeta and liv/pm), lead here is the most significant one
#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_rb_lead_merged.R <tissue> <nsv_start> <nsv_end> <out_csv>\n  tissue = LIV or PM")
}

tissue   <- args[1]                # "LIV" or "PM"
nsv_start <- as.integer(args[2])
nsv_end   <- as.integer(args[3])
out_csv   <- args[4]

message("Tissue: ", tissue)
message("nsv range: ", nsv_start, " to ", nsv_end)
message("Output: ", out_csv)

p_thresh <- 5e-8
base_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA"

#### 1. Load BrainMeta (hg38) ####
brainmeta_hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz")
setDT(brainmeta_hg38)

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

compute_rb_for_tissue <- function(brainmeta, tissue_df, p_thresh = 5e-8, theta = 0) { ##CHANGED THIS
  # 1) merge on SNP + gene
  dat_full <- merge(
    brainmeta[, .(snp_key, Probe, b, SE, p)],
    tissue_df[, .(snp_key, gene_id, beta, se, pval)],
    by.x = c("snp_key", "Probe"),   # BrainMeta
    by.y = c("snp_key", "gene_id"), # LIV/PM
    all  = FALSE
  )
  
  # 2) require non-zero SEs
  dat_full <- dat_full[SE != 0 & se != 0]
  if (nrow(dat_full) == 0) {
    return(list(
      rb       = NA_real_,
      se_rb    = NA_real_,
      n_pairs  = 0L,
      n_shared = 0L,
      n_uniq1  = 0L,
      n_uniq2  = 0L,
      note     = "No overlap after merge"
    ))
  }
  
  # 3) one SNP–gene pair per gene:
  #    lead SNP based on lowest BrainMeta p among overlapping pairs
  dat <- dat_full[ , .SD[which.min(p)], by = Probe]
  
  # 4) GW-significant indices on this lead-only set
  idx1 <- which(dat$p    < p_thresh)   # BrainMeta p
  idx2 <- which(dat$pval < p_thresh)   # tissue p
  
  shared <- intersect(idx1, idx2)
  uniq1  <- setdiff(idx1, shared)
  uniq2  <- setdiff(idx2, shared)
  
  # Use all pairs that are significant in at least one dataset
  use_idx <- union(idx1, idx2)
  dat_use <- dat[use_idx]
  
  if (nrow(dat_use) < 10) {
    return(list(
      rb       = NA_real_,
      se_rb    = NA_real_,
      n_pairs  = nrow(dat_use),
      n_shared = length(shared),
      n_uniq1  = length(uniq1),
      n_uniq2  = length(uniq2),
      note     = "Too few pairs after filtering"
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
    rb       = as.numeric(rb_res[1]),
    se_rb    = as.numeric(rb_res[2]),
    n_pairs  = nrow(dat_use),
    n_shared = length(shared),
    n_uniq1  = length(uniq1),
    n_uniq2  = length(uniq2),
    note     = NA_character_
  )
}

#### 4. Loop over nsv with per-gene lead SNP based on shared merge ####
nsv_vals <- nsv_start:nsv_end
results_list <- list()

for (nsv in nsv_vals) {
  message("Processing tissue = ", tissue, ", nsv = ", nsv)
  
  qtl_path <- file.path(base_dir, sprintf("%s_nominal_nsv%d.txt", tissue, nsv))
  qtl_df   <- prep_qtl(qtl_path)
  
  rb_res <- compute_rb_for_tissue(brainmeta_hg38, qtl_df, theta = 0)
  
  results_list[[length(results_list) + 1]] <- data.frame(
    tissue  = tissue,
    nsv     = nsv,
    rb      = rb_res$rb,
    se_rb   = rb_res$se_rb,
    n_pairs = rb_res$n_pairs,
    n_shared = rb_res$n_shared,
    n_uniq1  = rb_res$n_uniq1,
    n_uniq2  = rb_res$n_uniq2,
    note     = rb_res$note,
    stringsAsFactors = FALSE
  )
}

rb_summary <- bind_rows(results_list)
write.csv(rb_summary, out_csv, row.names = FALSE)
message("Done. Written: ", out_csv)
