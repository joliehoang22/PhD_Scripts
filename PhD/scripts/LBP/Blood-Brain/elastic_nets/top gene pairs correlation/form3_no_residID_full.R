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
library(furrr)      
library(future)    
library(purrr)
library(glmnet)   
set.seed(2025)

source("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/scripts/elastic_net_function.R")

################################### no_indivdualID_brain_full ###############################
##### discovered to be data leakage on 8-7-2025
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_form3_by_brain_cor) #21046 21357
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
summary_stats <- summary(as.vector(as.matrix(blood_form3_by_brain_cor))); summary_stats
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02266 0.04790 0.05637 0.08143 0.87792 
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2414433 | strongest correlations (top 0.05%, 224k gene-pairs) are ≥ 0.2414433
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2862253

# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= blood_form3_99.95, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)

dim(top_gene_pairs) #224734      3
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
head(top_gene_pairs)
length(unique(top_gene_pairs$brain_gene)) #16789
length(unique(top_gene_pairs$blood_gene)) #19683

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#168716 ENSG00000226259.10 ENSG00000226259.10   0.8779214
#169105  ENSG00000226752.9  ENSG00000226752.9   0.8499747
#180166  ENSG00000241945.8  ENSG00000241945.8   0.8441066
#85321  ENSG00000145736.14 ENSG00000145736.14   0.8389223
#178047  ENSG00000237541.3  ENSG00000237541.3   0.8304720
#207509  ENSG00000274602.5  ENSG00000274602.5   0.8243647

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #621
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224113

# Print the result
cat("Same gene pairs:", n_same, "\nDifferent gene pairs:", n_different, "\n")
#Same gene pairs: 507 
#Different gene pairs: 224227 

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

blood_expr <- blood_form3
brain_expr <- brain_full

# Train/test split
set.seed(2025)
n <- nrow(X_all)  # number of samples
train_idx <- sample(1:n, size = floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

# Brain genes to evaluate
brain_genes_to_test <- unique(top_gene_pairs$brain_gene)

results_df <- predict_brain_gene_expression(
  X_all = t(blood_expr),
  Y_all = t(brain_expr),
  top_gene_pairs = top_gene_pairs,
  brain_genes = brain_genes_to_test,
  train_idx = train_idx,
  test_idx = test_idx
)

write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

model <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID.txt", data.table=FALSE)
dim(model) #12732     6

summary(model$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.09999 0.14069 0.15902 0.19461 0.88891 

summary(model$cor_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.1843  0.2465  0.3362  0.3384  0.4293  0.9661      24

summary(model$r2_test) #data leakage
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.97452  0.02229  0.08584  0.08718  0.14924  0.92274 


################################################################################################################################
############################################ TRY 2 WITH nestedCV function ######################################################
##Improvements: 
#Feature selection on training data only | CV with single train/test split approach

# Full expression matrices (genes x samples)
blood_expr<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_expr) <- blood_expr$V1; blood_expr$V1 <- NULL

brain_expr<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
rownames(brain_expr) <- brain_expr$V1; brain_expr$V1 <- NULL

# Transpose to samples x genes (needed for modeling)
X_all <- t(blood_expr)
Y_all <- t(brain_expr)

set.seed(2025)
n <- nrow(X_all)
train_idx <- sample(1:n, size = floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

brain_genes_to_test <- colnames(Y_all)  # or a subset if testing

results_df <- predict_brain_gene_expression(
  X_all = X_all,
  Y_all = Y_all,
  brain_genes = brain_genes_to_test,
  train_idx = train_idx,
  test_idx = test_idx,
  top_percent = 0.05  # or try 0.01, 0.1, etc.
)

write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID_nestedCV.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

model <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID_nestedCV.txt", data.table=FALSE)
dim(model) #21356     6
summary(model$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.6310  0.7489  0.7147  0.8435  0.9985 
summary(model$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-3.9328 -0.4089 -0.2463 -0.2777 -0.1108  0.8909 
summary(model$cor_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.5217 -0.0643  0.0392  0.0443  0.1438  0.9095     416

################################################################################################################################
############################################ TRY 3 WITH multiple_nestedCV function #############################################
##Improvements: 
#Feature selection on training data only | CV with multiple train/test split approach
# Full expression matrices (genes x samples)
blood_expr<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_expr) <- blood_expr$V1; blood_expr$V1 <- NULL

brain_expr<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
rownames(brain_expr) <- brain_expr$V1; brain_expr$V1 <- NULL

results <- nested_cv_brain_prediction(
    X_all = t(blood_expr),  # samples x blood_genes
    Y_all = t(brain_expr),  # samples x brain_genes  
    brain_genes = colnames(t(brain_expr)),  
    n_outer_folds = 5,
    n_inner_folds = 5,
    top_percent = 0.05,
    alpha_values = c(0, 0.25, 0.5, 0.75, 1),
    seed = 2025
    #parallel = TRUE,  # set to FALSE if you don't want parallel processing
    #n_cores = 20
    )

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID_nestedCV_multiple.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 
# # View summary results
# head(results$summary_stats)
# 
# # View detailed results for specific brain gene
# results$detailed_results[results$detailed_results$brain_gene == "GENE_NAME", ]

##for results look at nested_CV_investigation.R

