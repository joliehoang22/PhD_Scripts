##### code derived from blood_brain_null_form0.R | OG_null_with_empirical_p.R (this is not as comprehensive)
## code derived from LBP_blood_brain_QC_20250603.rmd
rm(list=ls())
library(data.table)
library(ggplot2)
library(readxl)
library(biomaRt)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(matrixStats)
library(purrr)
library(caret)
library(tidyr)

set.seed(2025)
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

# Number of permutations
n_perm <- 1000

#################################
############# WITH ID ###########
#################################
blood_form5_with_residID <-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
rownames(blood_form5_with_residID) <- blood_form5_with_residID$V1; blood_form5_with_residID$V1 <- NULL

brain_full_with_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full_with_residID) <- brain_full_with_residID$V1; brain_full_with_residID$V1 <- NULL

#################################
########## NO RESID ID ##########
#################################
blood_form5_no_residID <-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL

brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1; brain_full_no_residID$V1 <- NULL

##pick one
blood_form <- blood_form5_with_residID
blood_form <- blood_form5_no_residID

##pick one - this is for the naming 
#these are with base_brain
blood_id <- "blood_form5_with_ID" 
blood_id <- "blood_form5_no_residID" 
#these are with full_brain
blood_id <- "blood_form5_with_ID_full_brain" 
blood_id <- "blood_form5_no_residID_full_brain" 

##pick one
brain_form <- brain_full_with_residID #
brain_form <- brain_full_no_residID #
brain_form <- v_brain$E #

output_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2"
null_output_path <- file.path(output_dir, paste0("null_summaries_", blood_id, ".txt"))
summary_output_path <- file.path(output_dir, paste0("comparison_of_null_and_real_data_summaries_", blood_id, ".txt"))

# Number of permutations
n_perm <- 1000
summary_stats <- c(
  "min", "mean", "median", "p80", "p90", "p95", 
  "p99.5", "p99.95", "p99.995", "p99.9995", "p99.99995", "p99.999995", 
  "p99.9999995", "max"
)

null_summaries <- matrix(NA, nrow = n_perm, ncol = length(summary_stats))
colnames(null_summaries) <- summary_stats

message("Starting null permutation procedure...")
for (i in 1:n_perm) {
  if (i %% 50 == 0) message("Completed ", i, " of ", n_perm, " permutations...")
  
  v_brain_null <- brain_form[, sample(ncol(brain_form))]
  null_cor <- cor(t(blood_form), t(v_brain_null), method = "spearman")
  abs_cor <- abs(null_cor)
  
  null_summaries[i, ] <- c(
    min = min(abs_cor),
    mean = mean(abs_cor),
    median = median(abs_cor),
    p80 = quantile(abs_cor, 0.80),
    p90 = quantile(abs_cor, 0.90),
    p95 = quantile(abs_cor, 0.95),
    p99.5 = quantile(abs_cor, 0.995),
    p99.95 = quantile(abs_cor, 0.9995), #22,472,918
    p99.995 = quantile(abs_cor, 0.99995), #2,247,291
    p99.9995 = quantile(abs_cor, 0.999995), #224,729
    p99.99995 = quantile(abs_cor, 0.9999995), #22,472
    p99.999995 = quantile(abs_cor, 0.99999995), #2247
    p99.9999995 = quantile(abs_cor, 0.999999995), #224 (eight 9s)
    max = max(abs_cor)
  )
}

null_summaries_df <- as.data.frame(null_summaries)
write.table(null_summaries_df, file = null_output_path, sep = "\t", row.names = FALSE, quote = FALSE)

message("Permutation complete. Next step: compute empirical p-values.")

summarize_cor_matrix <- function(cor_mat, null_summaries_df) {
  cor_mat <- abs(as.matrix(cor_mat))
  stat_names <- colnames(null_summaries_df)  # ensure consistency
  # compute observed values matching those stats
  get_stat <- function(stat) {
    if (stat == "min")    return(min(cor_mat, na.rm = TRUE))
    if (stat == "mean")   return(mean(cor_mat, na.rm = TRUE))
    if (stat == "median") return(median(cor_mat, na.rm = TRUE))
    if (stat == "max")    return(max(cor_mat, na.rm = TRUE))
    if (grepl("^p", stat)) {
      p <- as.numeric(sub("^p", "", stat))
      # p like 99.9995 -> convert to proportion
      q <- p/100
      return(quantile(cor_mat, q, na.rm = TRUE))
    }
    stop("Unknown stat: ", stat)
  }
  observed <- vapply(stat_names, get_stat, numeric(1))
  # one-sided upper tail for "large-is-extreme" stats, lower for min
  empirical_p <- mapply(function(stat, val) {
    if (stat == "min") {
      mean(null_summaries_df[[stat]] <= val)
    } else {
      mean(null_summaries_df[[stat]] >= val)
    }
  }, stat_names, observed)
  data.frame(stat = stat_names, observed_value = observed, empirical_p_value = empirical_p, row.names = NULL)
}

# Now pass the actual matrix
cor_mat <- cor(t(blood_form), t(brain_form), method = "spearman")
cor_mat <- abs(as.matrix(cor_mat))

result_df <- summarize_cor_matrix(cor_mat, null_summaries_df)
write.table(result_df, file = summary_output_path, sep = "\t", row.names = FALSE, quote = FALSE)

