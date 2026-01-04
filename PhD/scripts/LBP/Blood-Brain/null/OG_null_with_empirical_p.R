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

# Number of permutations
n_perm <- 1000

# Initialize storage
null_summaries <- matrix(NA, nrow = n_perm, ncol = 11) #used 7 col for the original null_summaries
#colnames(null_summaries) <- c("mean", "median", "p80", "p90", "p92.5", "p95", "p97.5")
colnames(null_summaries) <- c("mean", "median", "p80", "p90", "p92.5", "p95", "p97.5","p98.5","p99.5","p99.95","max")
message("Starting null permutation procedure...")
for (i in 1:n_perm) {
  if (i %% 50 == 0) message("Completed ", i, " of ", n_perm, " permutations...")
  
  # Shuffle brain sample columns to break subject alignment
  v_brain_null <- v_brain$E[, sample(ncol(v_brain$E))]

  # Compute Spearman correlation between blood and shuffled brain expression
  null_cor <- cor(t(v_blood$E), t(v_brain_null), method = "spearman")
  abs_cor <- abs(null_cor)
  
  # Store summary statistics
  null_summaries[i, ] <- c(
    mean = mean(abs_cor),
    median = median(abs_cor),
    p80 = quantile(abs_cor, 0.80),
    p90 = quantile(abs_cor, 0.90),
    p92.5 = quantile(abs_cor, 0.925),
    p95 = quantile(abs_cor, 0.95),
    p97.5 = quantile(abs_cor, 0.975),
    p98.5 = quantile(abs_cor, 0.985),
    p99.5 = quantile(abs_cor, 0.995),
    p99.95 = quantile(abs_cor, 0.9995),
    max = max(abs_cor)
  )
}
null_summaries_df <- as.data.frame(null_summaries)
write.table(null_summaries_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
message("Permutation complete. Next step: compute empirical p-values.")

##this has the original null summaries, colnames(null_summaries) <- c("mean", "median", "p80", "p90", "p92.5", "p95", "p97.5")
#write.table(null_summaries_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries.txt", sep = "\t", row.names = FALSE, quote = FALSE)

null_summaries_df <- read.table(
  file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries2.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)
dim(null_summaries_df)#1000   11
colnames(null_summaries_df)
# [1] "mean"   "median" "p80"    "p90"    "p92.5"  "p95"    "p97.5"  "p98.5" 
# [9] "p99.5"  "p99.95" "max" 

range(null_summaries_df$mean) #0.04415266 0.07808573
range(null_summaries_df$median) #0.03685841 0.06817109
range(null_summaries_df$p80) #0.07082596 0.12703329
range(null_summaries_df$p90) #0.0916930 0.1584893
range(null_summaries_df$p92.5) #0.0995965 0.1693279
range(null_summaries_df$p95) #0.1101496 0.1859777
range(null_summaries_df$p97.5) #0.1268911 0.2107859
range(null_summaries_df$p98.5) #0.1383765 0.2255847
range(null_summaries_df$p99.5) #0.1611820 0.2511968
range(null_summaries_df$p99.95) #0.2027223 0.2914412
range(null_summaries_df$max) #0.3408239 0.4281216

plot <- ggplot(null_summaries_df, aes(x = p99.95)) +
  geom_histogram(binwidth = 0.002, fill = "steelblue", color = "white") +
  geom_vline(xintercept = blood_by_brain_cor_99.95, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(
    breaks = seq(0, 0.3, by = 0.01),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "Distribution of Values",
    x = "Value",
    y = "Count"
  ) +
  coord_cartesian(xlim = c(0, 0.3)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.line.x = element_line(color = "black")
  )
ggsave("/hpc/users/hoangd02/www/plots/lbp/null_distribution_p99.95.pdf", plot = plot, width = 16, height = 8) 

### LOAD IN THE REAL STATS: BASELINE BLOOD AND BASELINE BRAIN 
blood_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_brain_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_by_brain_cor) <- blood_by_brain_cor$V1
blood_by_brain_cor$V1 <- NULL
blood_by_brain_cor <- abs(blood_by_brain_cor) 
blood_by_brain_cor_mean <- mean(as.matrix(blood_by_brain_cor)); blood_by_brain_cor_mean #0.05351357
blood_by_brain_cor_median <- median(as.matrix(blood_by_brain_cor)); blood_by_brain_cor_median  #0.04445627
blood_by_brain_cor_80 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.8); blood_by_brain_cor_80 #0.08596214
blood_by_brain_cor_90 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.9); blood_by_brain_cor_90 #0.1116572
blood_by_brain_cor_92.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.925); blood_by_brain_cor_92.5 #0.1213765
blood_by_brain_cor_95 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.95); blood_by_brain_cor_95 #0.1343279
blood_by_brain_cor_97.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.975); blood_by_brain_cor_97.5 #0.1547531 
blood_by_brain_cor_98.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.985); blood_by_brain_cor_98.5 #0.1686314
blood_by_brain_cor_99.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.995); blood_by_brain_cor_99.5 #0.1957751
blood_by_brain_cor_99.95 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.9995); blood_by_brain_cor_99.95 #0.2440944
blood_by_brain_cor_max <- max(as.matrix(blood_by_brain_cor)); blood_by_brain_cor_max #0.822524

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH original BASELINE BLOOD AND BASELINE BRAIN 
empirical_p_mean <- mean(null_summaries_df$mean >= blood_by_brain_cor_mean); empirical_p_mean #0.427
#In 42.7% of permuted datasets, the mean absolute correlation was at least as high as in the real paired dataset, so not significant.
empirical_p_median <- mean(null_summaries_df$median >= blood_by_brain_cor_median); empirical_p_median #0.533
empirical_p_80 <- mean(null_summaries_df$p80 >= blood_by_brain_cor_80); empirical_p_80 #0.43
empirical_p_90 <- mean(null_summaries_df$p90 >= blood_by_brain_cor_90); empirical_p_90 #0.337
empirical_p_92.5 <- mean(null_summaries_df$p92.5 >= blood_by_brain_cor_92.5); empirical_p_92.5 #0.306
empirical_p_95 <- mean(null_summaries_df$p95 >= blood_by_brain_cor_95); empirical_p_95 #0.265
empirical_p_97.5 <- mean(null_summaries_df$p97.5 >= blood_by_brain_cor_97.5); empirical_p_97.5 #0.219
empirical_p_98.5 <- mean(null_summaries_df$p98.5 >= blood_by_brain_cor_98.5); empirical_p_98.5 #0.191
empirical_p_99.5 <- mean(null_summaries_df$p99.5 >= blood_by_brain_cor_99.5); empirical_p_99.5 #0.147
empirical_p_99.95 <- mean(null_summaries_df$p99.95 >= blood_by_brain_cor_99.95); empirical_p_99.95 #0.078
empirical_p_max <- mean(null_summaries_df$max >= blood_by_brain_cor_max); empirical_p_max #0

### LOAD IN THE REAL STATS: k5_combo90 -- BE CAREFUL WITH THE NAME HERE 
blood_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc/blood_k5_combo90_brain_baseline_spearman_cor_matrix.txt",data.table=FALSE)
dim(blood_by_brain_cor) #21046 * 21356 = 449,458,376
blood_by_brain_cor[1:5,1:5]
row.names(blood_by_brain_cor) <- blood_by_brain_cor$V1
blood_by_brain_cor$V1 <- NULL
blood_by_brain_cor <- abs(blood_by_brain_cor)
blood_by_brain_cor_mean <- mean(as.matrix(blood_by_brain_cor)); blood_by_brain_cor_mean #0.06131221
blood_by_brain_cor_median <- median(as.matrix(blood_by_brain_cor)); blood_by_brain_cor_median #0.05226507
blood_by_brain_cor_80 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.8); blood_by_brain_cor_80 #0.09864096
blood_by_brain_cor_90 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.9); blood_by_brain_cor_90 #0.1258407
blood_by_brain_cor_92.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.925); blood_by_brain_cor_92.5 #0.1358534
blood_by_brain_cor_95 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.95); blood_by_brain_cor_95 #0.1489855
blood_by_brain_cor_97.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.975); blood_by_brain_cor_97.5 #0.1693015
blood_by_brain_cor_99.5 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.995); blood_by_brain_cor_99.5 #0.2094638
blood_by_brain_cor_99.95 <- quantile(as.matrix(blood_by_brain_cor), probs = 0.9995); blood_by_brain_cor_99.95 #0.25685 

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH k5_combo90 
empirical_p_mean <- mean(null_summaries_df$mean >= blood_by_brain_cor_mean); empirical_p_mean #0.063
empirical_p_median <- mean(null_summaries_df$median >= blood_by_brain_cor_median); empirical_p_median #0.08
empirical_p_80 <- mean(null_summaries_df$p80 >= blood_by_brain_cor_80); empirical_p_80 #0.061
empirical_p_90 <- mean(null_summaries_df$p90 >= blood_by_brain_cor_90); empirical_p_90 #0.054
empirical_p_92.5 <- mean(null_summaries_df$p92.5 >= blood_by_brain_cor_92.5); empirical_p_92.5 #0.053
empirical_p_95 <- mean(null_summaries_df$p95 >= blood_by_brain_cor_95); empirical_p_95 #0.05
empirical_p_97.5 <- mean(null_summaries_df$p97.5 >= blood_by_brain_cor_97.5); empirical_p_97.5 #0.044

#######TRY TO COMPARE THIS WITH DATA THAT HAS BEEN "CLEANED": FORM 3 WITH RESID ID | BRAIN BASE
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_mean <- mean(as.matrix(blood_form3_by_brain_cor)); blood_form3_mean #0.05487351
blood_form3_median <- median(as.matrix(blood_form3_by_brain_cor)); blood_form3_median  #0.04690055
blood_form3_80 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.8); blood_form3_80 #0.0881279
blood_form3_90 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9); blood_form3_90 #0.1122777
blood_form3_92.5 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.925); blood_form3_92.5 #0.1212158
blood_form3_95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.95); blood_form3_95 #0.1329857
blood_form3_97.5 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.975); blood_form3_97.5 #0.1513148
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2300379

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH blood_form3_brain_baseline_spearman_cor_matrix
empirical_p_mean <- mean(null_summaries_df$mean >= blood_form3_mean); empirical_p_mean #0.321
empirical_p_median <- mean(null_summaries_df$median >= blood_form3_median); empirical_p_median #0.325
empirical_p_80 <- mean(null_summaries_df$p80 >= blood_form3_80); empirical_p_80 #0.321
empirical_p_90 <- mean(null_summaries_df$p90 >= blood_form3_90); empirical_p_90 #0.319 
empirical_p_92.5 <- mean(null_summaries_df$p92.5 >= blood_form3_92.5); empirical_p_92.5 #0.311
empirical_p_95 <- mean(null_summaries_df$p95 >= blood_form3_95); empirical_p_95 #0.303
empirical_p_97.5 <- mean(null_summaries_df$p97.5 >= blood_form3_97.5); empirical_p_97.5 #0.299

####### I RAN THIS MORE EXTENSIVELY ELSE WHERE -- HERE'S A SUMMARY OF EMPIRICAL P VALUE BEYOND THE 97.5TH PERCENTILES
#code is in LBP_blood_brain_QC_20250603.rmd
form3_no_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_no_residID.txt", data.table=FALSE)
form3_with_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_with_residID.txt", data.table=FALSE)

####### FORM 5 WITH RESID ID BRAIN BASE
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_mean <- mean(as.matrix(blood_form5_by_brain_cor)); blood_form5_mean #0.05458049
blood_form5_median <- median(as.matrix(blood_form5_by_brain_cor)); blood_form5_median  #0.04686789
blood_form5_99.5 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.995); blood_form5_99.5 #0.1845849
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.226457

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH blood_form5_brain_baseline_spearman_cor_matrix
empirical_p_mean <- mean(null_summaries_df$mean >= blood_form5_mean); empirical_p_mean #0.338
empirical_p_median <- mean(null_summaries_df$median >= blood_form5_median); empirical_p_median #0.326
empirical_p_99.5 <- mean(null_summaries_df$p99.5 >= blood_form5_99.5); empirical_p_99.5 #0.358
empirical_p_99.95 <- mean(null_summaries_df$p99.95 >= blood_form5_99.95); empirical_p_99.95 #0.359

####### FORM 5 WITH RESID ID BRAIN FULL
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_mean <- mean(as.matrix(blood_form5_by_brain_cor)); blood_form5_mean #0.05381891
blood_form5_median <- median(as.matrix(blood_form5_by_brain_cor)); blood_form5_median  #0.04566582
blood_form5_97.5 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.975); blood_form5_97.5 #0.1501275

blood_form5_99.5 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.995); blood_form5_99.5 #0.1869606
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2304362

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH blood_form5_brain_baseline_spearman_cor_matrix
empirical_p_mean <- mean(null_summaries_df$mean >= blood_form5_mean); empirical_p_mean #0.4
empirical_p_median <- mean(null_summaries_df$median >= blood_form5_median); empirical_p_median #0.423
empirical_p_99.5 <- mean(null_summaries_df$p99.5 >= blood_form5_99.5); empirical_p_99.5 #0.297
empirical_p_99.95 <- mean(null_summaries_df$p99.95 >= blood_form5_99.95); empirical_p_99.95 #0.259

####### FORM 5 NO RESID ID BRAIN BASE
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_mean <- mean(as.matrix(blood_form5_by_brain_cor)); blood_form5_mean #0.05794256
blood_form5_median <- median(as.matrix(blood_form5_by_brain_cor)); blood_form5_median  #0.04984197
blood_form5_99.5 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.995); blood_form5_99.5 #0.194803
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2384292

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH blood_form5_brain_baseline_spearman_cor_matrix
empirical_p_mean <- mean(null_summaries_df$mean >= blood_form5_mean); empirical_p_mean #0.164
empirical_p_median <- mean(null_summaries_df$median >= blood_form5_median); empirical_p_median #0.155
empirical_p_99.5 <- mean(null_summaries_df$p99.5 >= blood_form5_99.5); empirical_p_99.5 #0.16
empirical_p_99.95 <- mean(null_summaries_df$p99.95 >= blood_form5_99.95); empirical_p_99.95 #0.128

####### FORM 5 NO RESID ID BRAIN FULL
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_mean <- mean(as.matrix(blood_form5_by_brain_cor)); blood_form5_mean #0.05463909
blood_form5_median <- median(as.matrix(blood_form5_by_brain_cor)); blood_form5_median  #0.04636009
blood_form5_99.5 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.995); blood_form5_99.5 #0.1896629 
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2337137

#How often would I observe a statistic as extreme or more extreme than my real value just by chance? WITH blood_form5_brain_baseline_spearman_cor_matrix
empirical_p_mean <- mean(null_summaries_df$mean >= blood_form5_mean); empirical_p_mean #0.334
empirical_p_median <- mean(null_summaries_df$median >= blood_form5_median); empirical_p_median #0.363
empirical_p_99.5 <- mean(null_summaries_df$p99.5 >= blood_form5_99.5); empirical_p_99.5 #0.235
empirical_p_99.95 <- mean(null_summaries_df$p99.95 >= blood_form5_99.95); empirical_p_99.95 #0.192
