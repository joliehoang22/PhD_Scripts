library(data.table)
library(readxl)
library(edgeR)
library(limma)
library(matrixStats)
library(dplyr)
library(purrr)
library(variancePartition)
library(BiocParallel)
library(tidyr)

####### PREVIOUS RELEVANT CODE
##code derived from blood-brain LBP_clean.R
# load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline.RData")

# ##code derived from LBP_clean_blood_brain.ipynb and LBP_blood_brain_QC_20250603.rmd
# load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_voom_20250522.RData")

##code derived from LBP_clean_blood_brain.ipynb and LBP_blood_brain_QC_20250603.rmd
#load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")
##this version is blood-brain_sample_baseline_with_voom_20250522.RData + removing samples 
## this doesn't residualized out of covars for v_brain or v_blood yet; v_blood <- voomWithDreamWeights(dge_filtered_blood, formula = ~1, data = blood_metadata)

##code derived from 01_two_pairs_exploration.R (not in the preprocessing folder yet)
#load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
## this is similar to the blood-brain_sample_baseline_with_dictionary_20250729.RData but with 2 additional columns on surgeryDate and timepoint (0 and 1)
dim(v_blood$E) #21046   225
dim(v_brain$E) #21356   225
dim(blood_only_metadata) #225 171
dim(brain_only_metadata) #225 171

############## BLOOD FORMULA ##############
###########################################
##code derived from LBP_clean_blood_brain.ipynb
# form3 <- ~ (1|IID_ISMMS) + mymet_rin_blood + STAR_Insertion_average_length_blood
# fit_blood <- dream(v_blood$E, form3, blood_metadata)
# resid_expr_blood_form3 <- residuals(fit_blood)

# write.csv(resid_expr_blood_form3, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3.csv", row.names = TRUE)
# write.table(resid_expr_blood_form3, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3.txt", sep = "\t", quote = FALSE, col.names = NA)

############## BRAIN FORMULA ##############
###########################################
##from LBP_blood_brain_QC_20250603.rmd
# form_full <- ~ (1|IID_ISMMS) + (1|mymet_depletionbatch_brain) + 
# mymet_rin_brain + (1|mymet_bank_brain) + (1|mymet_sex_brain) + RNASeqMetrics_INTRONIC_BASES_brain + ODC_brain + 
# (1|mymet_tissue_brain)

# fit_brain<- dream(v_brain$E, form_full, brain_metadata)
# resid_expr_brain_form_full <- residuals(fit_brain)
# write.table(resid_expr_brain_form_full, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", sep = "\t", quote = FALSE, col.names = NA)

#### NEW BRAIN FORM WITH + GABA_brain
########################################## with resid ID 
form_full <- ~ (1|IID_ISMMS) + (1|mymet_depletionbatch_brain) + mymet_rin_brain + (1|mymet_bank_brain) + (1|mymet_sex_brain) +
RNASeqMetrics_INTRONIC_BASES_brain + ODC_brain + (1|mymet_tissue_brain) + GABA_brain

fit_brain<- dream(v_brain$E, form_full, brain_only_metadata)
resid_expr_brain_form_full <- residuals(fit_brain)
write.table(resid_expr_brain_form_full, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_with_GABA.txt", sep = "\t", quote = FALSE, col.names = NA)

########################################## no resid ID 
form_full <- ~ (1|mymet_depletionbatch_brain) + mymet_rin_brain + (1|mymet_bank_brain) + (1|mymet_sex_brain) +
RNASeqMetrics_INTRONIC_BASES_brain + ODC_brain + (1|mymet_tissue_brain) + GABA_brain

fit_brain<- dream(v_brain$E, form_full, brain_only_metadata)
resid_expr_brain_form_full <- residuals(fit_brain)
write.table(resid_expr_brain_form_full, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_with_GABA_no_indivdualID.txt", sep = "\t", quote = FALSE, col.names = NA)

brain_full <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_with_GABA.txt",data.table=FALSE)
brain_full_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_with_GABA_no_indivdualID.txt",data.table=FALSE)

###### CHOOSE ONE
rownames(brain_full) <- brain_full$V1
brain_full$V1 <- NULL
dim(brain_full) #21356   225
brain_full [1:3,1:3]

rownames(brain_full_no_residID) <- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225
brain_full_no_residID[1:3,1:3]

#########FORMULAS with residID 
blood_form3 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3.txt",data.table=FALSE)
blood_form5 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)

#########FORMULAS with NO residID 
blood_form3_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
blood_form5_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)

###### CHOOSE ONE
#blood_form3 
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#blood_form5 
rownames(blood_form5) <- blood_form5$V1
blood_form5$V1 <- NULL
dim(blood_form5) #21046   225

#blood_form3_no_residID
rownames(blood_form3_no_residID) <- blood_form3_no_residID$V1
blood_form3_no_residID$V1 <- NULL
dim(blood_form3_no_residID) #21046   225

#blood_form5_no_residID
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1
blood_form5_no_residID$V1 <- NULL
dim(blood_form5_no_residID) #21046   225

############## CORRELATIONS ##############
###########################################
### BLOOD FORM 5 AND BASELINE BRAIN 

cor_matrix <- cor(t(blood_form3_no_residID), t(brain_full_no_residID), method = "spearman")
write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cor_matrix <- cor(t(blood_form5_no_residID), t(brain_full_no_residID), method = "spearman")
write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cor_matrix <- cor(t(blood_form3), t(brain_full), method = "spearman")
write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_with_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cor_matrix <- cor(t(blood_form5), t(brain_full), method = "spearman")
write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_with_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

#write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cor_matrix_mean <- mean(as.matrix(cor_matrix)); cor_matrix_mean #
cor_matrix_97.5 <- quantile(as.matrix(cor_matrix), probs = 0.975); cor_matrix_97.5 #

######################## FOR DR. RHODES ########################
## extracting top correlations 
#cor_mat <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_with_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt",data.table=FALSE)
#cor_mat <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_with_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt",data.table=FALSE)
#cor_mat <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt",data.table=FALSE)
cor_mat <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_withGABA_spearman_cor_matrix.txt",data.table=FALSE)

dim(cor_mat) #21046 x 21356

# assume first column is blood gene IDs and column names (from 2nd col onward) are brain gene IDs
blood_genes <- cor_mat$V1
brain_genes <- colnames(cor_mat)[-1]

# keep only genes present in both (intersection)
common_genes <- intersect(blood_genes, brain_genes)

length(common_genes) #17533

# subset to the diagonal elements (same gene in both)
diag_vals <- sapply(common_genes, function(g)
  cor_mat[cor_mat$V1 == g, g])

# convert to a data frame
diag_df <- data.frame(Gene = common_genes, Correlation = diag_vals)
summary(diag_df$Correlation)
#######blood_form5_with_indivdualID_brain_full
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.27440 -0.02320  0.02429  0.02809  0.07407  0.43031

#######blood_form3_with_indivdualID_brain_full
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.26154 -0.02185  0.02613  0.02950  0.07726  0.47186

#######blood_form3_no_indivdualID_brain_full
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.46782 -0.01222  0.03945  0.05092  0.09569  0.87493 

#######blood_form5_no_indivdualID_brain_full
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.48897 -0.01222  0.03882  0.05111  0.09472  0.87562

# filter for correlations > 0.3
diag_filtered <- subset(diag_df, Correlation > 0.3)

# sort descending
diag_filtered <- diag_filtered[order(-diag_filtered$Correlation), ]
rownames(diag_filtered) <- seq_len(nrow(diag_filtered))
length(diag_filtered$Correlation) 
#54 genes blood_form5_with_indivdualID_brain_full
#53 genes blood_form3_with_indivdualID_brain_full
#428 genes blood_form3_no_indivdualID_brain_full
#468 genes blood_form5_no_indivdualID_brain_full

head(diag_filtered)
#                 Gene Correlation
# 1  ENSG00000238083.7   0.4303076
# 2 ENSG00000176681.14   0.4170702
# 3  ENSG00000285534.1   0.4092183
# 4  ENSG00000223496.3   0.4070333
# 5  ENSG00000260671.2   0.4039265
# 6  ENSG00000247498.9   0.3985598

write.csv(diag_filtered, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_top_correlated_same_genes.csv", row.names = FALSE)
#write.csv(diag_filtered, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_top_correlated_same_genes.csv", row.names = FALSE)
#write.csv(diag_filtered, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_with_indivdualID_brain_full_top_correlated_same_genes.csv", row.names = FALSE)
#write.csv(diag_filtered, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_with_indivdualID_brain_full_top_correlated_same_genes.csv", row.names = FALSE)

###comparing the results, seems like the top correlated gene pairs are actually pretty robust to each other!! 
### 

test <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_top_correlated_same_genes.csv")

######## EXPLORE SAMPLE-WISE SPEARMAN CORRELATIONS ACROSS GENES BTW T0 AND T1 AGAIN 
##code derived from 02_building_foundations.R

########## DEFINE NEW BLOOD AND BRAIN --- 
##original -- no residualizign out GABA 
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full) <- brain_full$V1
brain_full$V1 <- NULL

blood <- blood_form5_no_residID  #blood_form5, blood_form3, blood_form3_no_residID or blood_form5_no_residID
brain <- brain_full_no_residID #brain_full or brain_full_no_residID

blood <- blood_form5  #blood_form5, blood_form3, blood_form3_no_residID or blood_form5_no_residID
brain <- brain_full #brain_full or brain_full_no_residID

#load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData") #already loaded up there
dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)

two_pairs <- filter(dictionary, number_of_pair_brain == "2_pairs") #146, correct
dim(two_pairs) #146   8

# --- Helper: correlate one blood–brain column pair on common genes ---
get_corr <- function(blood_id, brain_id, blood_mat, brain_mat) {
  common <- intersect(rownames(blood_mat), rownames(brain_mat))
  if (length(common) < 3) return(NA_real_)
  x <- blood_mat[common, blood_id, drop = TRUE]
  y <- brain_mat[common, brain_id, drop = TRUE]
  suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

# --- 1) Keep only rows that exist in matrices ---
two_pairs_clean <- two_pairs %>%
  filter(SAMPLE_ISMMS_blood %in% colnames(blood),
         SAMPLE_ISMMS_brain %in% colnames(brain)) %>%
  filter(timepoint %in% c(0, 1))

# --- 2) Same-timepoint correlations (T0_only, T1_only) ---
within_same <- two_pairs_clean %>%
  rowwise() %>%
  mutate(
    corr = get_corr(SAMPLE_ISMMS_blood, SAMPLE_ISMMS_brain, blood, brain),
    type = paste0("T", timepoint, "_only")
  ) %>%
  ungroup() %>%
  select(IID_ISMMS, type, corr) %>%
  distinct()

# --- 3) Cross-timepoint correlations (within person; avoid duplicates) ---
within_cross <- two_pairs_clean %>%
  inner_join(two_pairs_clean, by = "IID_ISMMS", suffix = c("_b", "_br")) %>%
  filter(timepoint_b != timepoint_br) %>%
  rowwise() %>%
  mutate(corr = get_corr(SAMPLE_ISMMS_blood_b, SAMPLE_ISMMS_brain_br, blood, brain)) %>%
  ungroup() %>%
  distinct(IID_ISMMS, corr) %>%
  mutate(type = "Cross_timepoint")

# --- 4) Long table (no pivot_wider) ---
cor_long <- bind_rows(within_same, within_cross)

table(cor_long$type)
# Cross_timepoint         T0_only         T1_only 
#             146              73              73 

head(cor_long)
#   IID_ISMMS type        corr
#   <chr>     <chr>      <dbl>
# 1 PT-0018   T1_only  0.0608 
# 2 PT-0018   T0_only  0.00202
# 3 PT-0021   T0_only -0.00861
# 4 PT-0021   T1_only -0.0489 
# 5 PT-0022   T1_only -0.0336 
# 6 PT-0022   T0_only  0.151

###making a new cor_diff df to see which individuals have elevated vs deflated corr
### POSITIVE diff = higher T1 than T0 | T1 > T0
### NEGATIVE diff = lower T1 than T0
cor_diff <- cor_long %>%
  mutate(corr = as.numeric(corr)) %>%
  pivot_wider(
    names_from = type, values_from = corr,
    values_fn = ~ mean(as.numeric(.x), na.rm = TRUE)
  ) %>%
  mutate(diff = T1_only - T0_only)

head(cor_diff)

#########correlating T1_only, T0_only, and diff to all metadata variables
# join metrics to numeric metadata across subjects (BRAIN METADATA has both blood and brain)
joined <- cor_diff %>%
  select(IID_ISMMS, T1_only, T0_only, diff) %>%
  inner_join(brain_metadata %>% select(IID_ISMMS, where(is.numeric)),
             by = "IID_ISMMS")

# join metrics to numeric metadata across subjects BLOOD ONLY
# joined <- cor_diff %>%
#   select(IID_ISMMS, T1_only, T0_only, diff) %>%
#   inner_join(blood_only_metadata %>% select(IID_ISMMS, where(is.numeric)),
#              by = "IID_ISMMS")

# correlate each numeric metadata column with T0_only, T1_only, and diff
features <- setdiff(names(joined), c("IID_ISMMS","T1_only","T0_only","diff"))

cors_mat <- cor(
  joined[, features],
  joined[, c("T0_only","T1_only","diff")],
  use = "pairwise.complete.obs"  # handles NAs
)

# tidy df with feature names
cors_df <- data.frame(feature = rownames(cors_mat), cors_mat, row.names = NULL)
cors_df  # columns are feature, T0_only, T1_only, diff

###sort T0 from highest to lowest
head(cors_df[order(cors_df$diff, decreasing = TRUE), ],4) ##BRAIN

################### ORIGINAL ################### WITHOUT RESIDUALIZING OUT FOR GABA
#################################################
#                                  feature     T0_only    T1_only      diff
# 138                           GABA_brain -0.19204561 0.14635978 0.2358734
# 137                            GLU_brain -0.22371681 0.05102085 0.2019587
# 23  STAR_Uniquely_mapped_reads_pct_brain -0.08396409 0.13119338 0.1440828
# 115            mymet_rna_conc_ngul_brain -0.12840468 0.07401698 0.1431786



######### blood_form5_no_residID and brain_full_no_residID
#                                feature    T0_only    T1_only      diff
# 137                          GLU_brain -0.2026156 0.03948529 0.1769886
# 138                         GABA_brain -0.1078835 0.12647566 0.1598646
# 6     STAR_Average_mapped_length_brain -0.1012157 0.09863871 0.1375249
# 16  STAR_Number_of_splices_AT_AC_brain -0.0151611 0.18982555 0.1295401

######### blood_form5 and brain_full
#                                                        feature     T0_only
# 14                   STAR_Number_of_reads_unmapped_other_brain -0.07335767
# 230                InsertSizeMetrics_WIDTH_OF_95_PERCENT_blood -0.11329755
# 38            AlignmentSummaryMetrics_PF_INDEL_RATE_PAIR_brain -0.04701712
# 67  AlignmentSummaryMetrics_PF_INDEL_RATE_SECOND_OF_PAIR_brain -0.06200817
#        T1_only      diff
# 14  0.10709493 0.1339540
# 230 0.04355454 0.1238762
# 38  0.11111556 0.1147914
# 67  0.08637667 0.1104028

######### blood_form3_no_residID and brain_full_no_residID
#                                         feature    T0_only    T1_only      diff
# 137                                   GLU_brain -0.1980611 0.05536045 0.1844502
# 138                                  GABA_brain -0.1375432 0.11962649 0.1801114
# 6              STAR_Average_mapped_length_brain -0.1006578 0.08696854 0.1314404
# 229 InsertSizeMetrics_WIDTH_OF_90_PERCENT_blood -0.1033049 0.08092115 0.1295564

######### blood_form3 and brain_full
#                                         feature     T0_only    T1_only
# 138                                  GABA_brain -0.11786835 0.08839575
# 6              STAR_Average_mapped_length_brain -0.11479006 0.08162971
# 137                                   GLU_brain -0.15624805 0.02422643
# 231 InsertSizeMetrics_WIDTH_OF_99_PERCENT_blood -0.01001655 0.16681825
#          diff
# 138 0.1565056
# 6   0.1492606
# 137 0.1416186
# 231 0.1261528

cor(brain_metadata[,c("GABA_brain","GLU_brain","ODC_brain")])
           GABA_brain  GLU_brain  ODC_brain
GABA_brain   1.000000  0.4299940 -0.4496800
GLU_brain    0.429994  1.0000000 -0.9882269
ODC_brain   -0.449680 -0.9882269  1.0000000

#### Boxplot and matching points between T0 and T1
library(dplyr)
library(tidyr)
library(ggplot2)

# Keep only T0/T1, drop NA, and require both timepoints per subject
df_pair <- cor_long %>%
  filter(type %in% c("T0_only", "T1_only")) %>%            # exclude cross-timepoint or others
  filter(!is.na(corr)) %>%
  group_by(IID_ISMMS) %>%
  filter(n_distinct(type) == 2) %>%                        # keep only IDs with both T0 & T1
  ungroup() %>%
  mutate(type = factor(type, levels = c("T0_only", "T1_only"),
                       labels = c("T0", "T1")))

# reshape the data to wide format: one row per IID, columns for T0 and T1
df_wide <- df_pair %>%
  pivot_wider(names_from = type, values_from = corr)
# run paired Wilcoxon signed-rank test
wilcox.test(df_wide$T0, df_wide$T1, paired = TRUE)

# compute Cohen's d for paired samples
d <- (mean(df_wide$T1 - df_wide$T0)) / sd(df_wide$T1 - df_wide$T0); d

# Add annotation bar coordinates
annot_df <- data.frame(
  xmin = 1,   # T0
  xmax = 2,   # T1
  y = 0.45,   # height of the bar (adjust depending on your data range)
  #label = "d = -0.092, p = 0.571" #form5 addGABA with_residID
  #label = "d = -0.158, p = 0.145" #form5 addGABA no_residID
  label =  "d = -0.093, p = 0.320" #form5_no_residID (original)
)

p <- ggplot(df_pair, aes(x = type, y = corr)) +
  # Boxplots (one per timepoint)
  geom_boxplot(aes(group = type),
               width = 0.35, outlier.shape = NA,
               fill = "grey90", color = "black") +
  # Connecting lines between T0 and T1 per subject
  geom_line(aes(group = IID_ISMMS),
            alpha = 0.25, linewidth = 0.6, color = "gray50") +
  # Points with transparency (easier to see overlap)
  geom_point(aes(group = IID_ISMMS),
             size = 2.8, alpha = 0.5, color = "black") +
  # Add bar and label above
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmin, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmax, xend = xmax, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_text(data = annot_df,
            aes(x = (xmin + xmax)/2, y = y + 0.02, label = label),
            inherit.aes = FALSE, size = 5.5) +
  labs(x = NULL, y = "Spearman correlation") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 10)
  )
ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_formula5_no_residID_original.pdf", p, width = 8, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_formula5_addGABA_no_residID.pdf", p, width = 8, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_formula5_addGABA_with_residID.pdf", p, width = 8, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_matching_points_wcpg.pdf", p, width = 8, height = 5)
