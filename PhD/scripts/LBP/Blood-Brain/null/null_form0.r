######################################################
######################################################
###### This script is for the creation of null models 
######################################################
######################################################
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

##### DYNAMIC CODE,
## Change 3 things: CHANGE 1, 2, 3

###### BLOOD #####
##################
#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE) ##CHANGE 1
#rownames(blood_form3) <- blood_form3$V1
#blood_form3$V1 <- NULL
#blood_form0[1:5,1:5]
#dim(blood_form0) #21046   225

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_voom_20250522.RData")
blood_samples_to_be_removed <- c("LBPSEMA4BLOOD795", "LBPSEMA4BLOOD142", "LBPSEMA4BLOOD049", "LBPSEMA4BLOOD755",
                                  "LBPSEMA4BLOOD335","LBPSEMA4BLOOD431","LBPSEMA4BLOOD174","LBPSEMA4BLOOD738")
v_blood$E <- v_blood$E[, !colnames(v_blood$E) %in% blood_samples_to_be_removed]
blood_metadata <- blood_metadata[!rownames(blood_metadata) %in% blood_samples_to_be_removed, ]
blood_only_metadata <- blood_metadata[, grepl("_blood$", names(blood_metadata))]
blood_only_metadata <- blood_only_metadata[!rownames(blood_only_metadata) %in% blood_samples_to_be_removed, ]
blood_only_metadata$IID_ISMMS <- blood_metadata$IID_ISMMS[match(rownames(blood_only_metadata), rownames(blood_metadata))]

#blood_form = v_blood$E
blood_form = v_blood$E ##CHANGE 2

#cor_matrix <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",data.table=FALSE) ##CHANGE 3
#cor_matrix <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_brain_baseline_spearman_cor_matrix.txt",data.table=FALSE)
cor_matrix <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_brain_spearman_cor_matrix.txt",data.table=FALSE)
cor_matrix[1:5,1:5]
row.names(cor_matrix) <- cor_matrix$V1
cor_matrix$V1 <- NULL

##################
blood_id <- "blood_form0" ##CHANGE 4 - this is the name of the outputs
#blood_id <- "blood_form3"
output_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2"
null_output_path <- file.path(output_dir, paste0("null_summaries_", blood_id, ".txt"))
summary_output_path <- file.path(output_dir, paste0("comparison_of_null_and_real_data_summaries_", blood_id, ".txt"))

base_null <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form0.txt", data.table=FALSE)
form3_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form3_no_residID.txt", data.table=FALSE)
form3_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form3_with_residID.txt", data.table=FALSE)

base_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form0.txt", data.table=FALSE)
form3_no_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_no_residID.txt", data.table=FALSE)
form3_with_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_with_residID.txt", data.table=FALSE)

form3_no_residID_summary
#           stat observed_value empirical_p_value
#11        p99.5     0.20201327             0.063
#12       p99.95     0.24839444             0.046
#13      p99.995     0.29069225             0.018
#14     p99.9995     0.43836985             0.000

form3_with_residID_summary
#           stat observed_value empirical_p_value
#13      p99.995     0.26650653             0.273
#14     p99.9995     0.31362213             0.060
#15    p99.99995     0.41299296             0.000

#null_summaries_base <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form0.txt
#", data.table=FALSE)
#/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form0.txt

###### BRAIN #####
##################
#load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_voom_20250522.RData")
## Remove outliers
brain_samples_to_be_removed <- c("LBPSEMA4BRAIN364", "LBPSEMA4BRAIN318", "LBPSEMA4BRAIN564", "LBPSEMA4BRAIN170",
                                  "LBPSEMA4BRAIN703", "LBPSEMA4BRAIN341","LBPSEMA4BRAIN391","LBPSEMA4BRAIN017")

v_brain$E <- v_brain$E[, !colnames(v_brain$E) %in% brain_samples_to_be_removed]
#dim(v_brain$E) #21356   225
brain_metadata <- brain_metadata[!rownames(brain_metadata) %in% brain_samples_to_be_removed, ]
#dim(brain_metadata) #225 339
brain_only_metadata <- brain_metadata[, grepl("_brain$", names(brain_metadata))]
brain_only_metadata <- brain_only_metadata[!rownames(brain_only_metadata) %in% brain_samples_to_be_removed, ]
#dim(brain_only_metadata) #225 169
brain_only_metadata$IID_ISMMS <- brain_metadata$IID_ISMMS[match(rownames(brain_only_metadata), rownames(brain_metadata))]
#dim(brain_only_metadata) #225 170

# Number of permutations
n_perm <- 1000
summary_stats <- c(
  "min", "p25", "mean", "median", "p80", "p90", "p92.5", "p95", "p97.5", "p98.5",
  "p99.5", "p99.95", "p99.995", "p99.9995", "p99.99995", "p99.999995", 
  "p99.9999995", "p99.99999995", "max"
)

null_summaries <- matrix(NA, nrow = n_perm, ncol = length(summary_stats))
colnames(null_summaries) <- summary_stats

message("Starting null permutation procedure...")
for (i in 1:n_perm) {
  if (i %% 50 == 0) message("Completed ", i, " of ", n_perm, " permutations...")
  
  v_brain_null <- v_brain$E[, sample(ncol(v_brain$E))]
  null_cor <- cor(t(blood_form), t(v_brain_null), method = "spearman")
  abs_cor <- abs(null_cor)
  
  null_summaries[i, ] <- c(
    min = min(abs_cor),
    p25 = quantile(abs_cor, 0.25),
    mean = mean(abs_cor),
    median = median(abs_cor),
    p80 = quantile(abs_cor, 0.80),
    p90 = quantile(abs_cor, 0.90),
    p92.5 = quantile(abs_cor, 0.925),
    p95 = quantile(abs_cor, 0.95),
    p97.5 = quantile(abs_cor, 0.975),
    p98.5 = quantile(abs_cor, 0.985),
    p99.5 = quantile(abs_cor, 0.995),
    p99.95 = quantile(abs_cor, 0.9995), #22,472,918
    p99.995 = quantile(abs_cor, 0.99995), #2,247,291
    p99.9995 = quantile(abs_cor, 0.999995), #224,729
    p99.99995 = quantile(abs_cor, 0.9999995), #22,472
    p99.999995 = quantile(abs_cor, 0.99999995), #2247
    p99.9999995 = quantile(abs_cor, 0.999999995), #224 (eight 9s)
    p99.99999995 = quantile(abs_cor, 0.9999999995), #24 (nine 9s)
    max = max(abs_cor)
  )
}

null_summaries_df <- as.data.frame(null_summaries)
write.table(null_summaries_df, file = null_output_path, sep = "\t", row.names = FALSE, quote = FALSE)

message("Permutation complete. Next step: compute empirical p-values.")

summarize_cor_matrix <- function(cor_mat, null_summaries_df) {
  cor_mat <- abs(as.matrix(cor_mat))

  stat_names <- c(
    "min", "p25", "mean", "median", "p80", "p90", "p92.5", "p95", "p97.5", "p98.5",
    "p99.5", "p99.95", "p99.995", "p99.9995", "p99.99995", "p99.999995",
    "p99.9999995", "p99.99999995", "max"
  )
  probs <- c(
    NA, 0.25, NA, NA, 0.80, 0.90, 0.925, 0.95, 0.975, 0.985,
    0.995, 0.9995, 0.99995, 0.999995, 0.9999995, 0.99999995,
    0.999999995, 0.9999999995, NA
  )

  values <- mapply(function(stat, prob) { #applies a function in parallel over multiple input arguments
    if (stat == "min") return(min(cor_mat))
    if (stat == "mean") return(mean(cor_mat))
    if (stat == "median") return(median(cor_mat))
    if (stat == "max") return(max(cor_mat))
    quantile(cor_mat, prob)
  }, stat_names, probs)

  empirical_p <- mapply(function(stat, val) {
    mean(null_summaries_df[[stat]] >= val)
  }, stat_names, values)

  data.frame(
    stat = stat_names,
    observed_value = as.numeric(values),
    empirical_p_value = empirical_p
  )
}

# Now pass the actual matrix
result_df <- summarize_cor_matrix(cor_matrix, null_summaries_df)
write.table(result_df, file = summary_output_path, sep = "\t", row.names = FALSE, quote = FALSE)

