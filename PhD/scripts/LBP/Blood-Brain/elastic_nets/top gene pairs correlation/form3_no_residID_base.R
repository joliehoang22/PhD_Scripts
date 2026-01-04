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

source("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/scripts/elastic_net_function.R")
set.seed(2025)

blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_form3_by_brain_cor) #21046 21357
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
summary_stats <- summary(as.vector(as.matrix(blood_form3_by_brain_cor))); summary_stats
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02399 0.05052 0.05901 0.08524 0.83440 
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2484 | strongest correlations (top 0.05%, 224k gene-pairs) are ≥ 0.2484
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2906922 

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
length(unique(top_gene_pairs$brain_gene)) #15695
length(unique(top_gene_pairs$blood_gene)) #15366

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#167582 ENSG00000226259.10 ENSG00000226259.10   0.8344016
#181336  ENSG00000241945.8  ENSG00000241945.8   0.8261620
#210686  ENSG00000274602.5  ENSG00000274602.5   0.8222672
#208259 ENSG00000233327.10  ENSG00000273018.6   0.8079362
#208262  ENSG00000273018.6  ENSG00000273018.6   0.8060124
#175112 ENSG00000233327.10 ENSG00000233327.10   0.8029794

set.seed(2025)
brain_genes_to_test <- unique(top_gene_pairs$brain_gene)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

blood_expr <- blood_form3
brain_expr <- v_brain$E

# Transpose expression matrices: samples × genes
X_all <- t(blood_expr)   # blood gene expression (samples x blood genes)
Y_all <- t(brain_expr)   # brain gene expression (samples x brain genes)

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

# Save to file
write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_base_genes_elastic_net_blood_form3_no_indivdualID.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################################################################
############################################ TRY 2 WITH nestedCV function ######################################################
##Improvements: 
#Feature selection on training data only | CV with single train/test split approach
source("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/scripts/elastic_net_funtion_nestedCV_single_test_and_train.R")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_expr <- blood_form3
brain_expr <- v_brain$E

# Transpose expression matrices: samples × genes
X_all <- t(blood_expr)   # blood gene expression (samples x blood genes)
Y_all <- t(brain_expr)   # brain gene expression (samples x brain genes)

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

write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_base_genes_elastic_net_blood_form3_no_indivdualID_nestedCV.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

model <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_base_genes_elastic_net_blood_form3_no_indivdualID_nestedCV.txt", data.table=FALSE)
dim(model) #21356     6
summary(model$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.6103  0.7069  0.6941  0.8005  0.9985
summary(model$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-2.5912 -0.4086 -0.2656 -0.2879 -0.1410  0.9679
summary(model$cor_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.55020 -0.06350  0.03123  0.03843  0.13070  0.86983      238

################################################################################################################################
############################################ TRY 3 WITH multiple_nestedCV function #############################################
##Improvements: 
#Feature selection on training data only | CV with multiple train/test split approach
source("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/scripts/elastic_net_function_nestedCV_multiple_test_and_train.R")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_expr <- blood_form3
brain_expr <- v_brain$E

results <- nested_cv_brain_prediction(
    X_all = t(blood_expr),  # samples x blood_genes
    Y_all = t(brain_expr),  # samples x brain_genes  
    brain_genes = colnames(t(brain_expr)),  # test subset first
    n_outer_folds = 5,
    n_inner_folds = 5,
    top_percent = 0.05,
    alpha_values = c(0, 0.25, 0.5, 0.75, 1),
    seed = 2025
    #parallel = TRUE, 
    #n_cores = 20
    )

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_base_genes_elastic_net_blood_form3_no_indivdualID_nestedCV_multiple.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

##for results look at nested_CV_investigation.R
