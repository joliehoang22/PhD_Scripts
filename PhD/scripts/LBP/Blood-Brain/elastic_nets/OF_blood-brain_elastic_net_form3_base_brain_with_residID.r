########################################################################
#################### ELASTIC NET for form 3 base brain #################
########################################################################
# Load required packages
library(glmnet)   # for elastic net
library(furrr)    # for parallelized map
library(tibble)   # for tibble output
library(purrr)    # for map functions
library(data.table)
library(ggplot2)
library(readxl)
library(biomaRt)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(matrixStats)
library(caret)
library(tidyverse)

#cor_mat <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",data.table=FALSE)

set.seed(2025)
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_voom_20250522.RData")
dim(v_blood$E) #21046   233
dim(v_brain$E) #21356   233

## Remove outliers
blood_samples_to_be_removed <- c("LBPSEMA4BLOOD795", "LBPSEMA4BLOOD142", "LBPSEMA4BLOOD049", "LBPSEMA4BLOOD755",
                                  "LBPSEMA4BLOOD335","LBPSEMA4BLOOD431","LBPSEMA4BLOOD174","LBPSEMA4BLOOD738")

brain_samples_to_be_removed <- c("LBPSEMA4BRAIN364", "LBPSEMA4BRAIN318", "LBPSEMA4BRAIN564", "LBPSEMA4BRAIN170",
                                  "LBPSEMA4BRAIN703", "LBPSEMA4BRAIN341","LBPSEMA4BRAIN391","LBPSEMA4BRAIN017")


v_blood$E <- v_blood$E[, !colnames(v_blood$E) %in% blood_samples_to_be_removed]
dim(v_blood$E) #21046   225
blood_metadata <- blood_metadata[!rownames(blood_metadata) %in% blood_samples_to_be_removed, ]
dim(blood_metadata) #225 339
blood_only_metadata <- blood_metadata[, grepl("_blood$", names(blood_metadata))]
blood_only_metadata <- blood_only_metadata[!rownames(blood_only_metadata) %in% blood_samples_to_be_removed, ]
dim(blood_only_metadata) #225 169
blood_only_metadata$IID_ISMMS <- blood_metadata$IID_ISMMS[match(rownames(blood_only_metadata), rownames(blood_metadata))]
dim(blood_only_metadata) #225 170

#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5) <- blood_form5$V1
blood_form5$V1 <- NULL

v_brain$E <- v_brain$E[, !colnames(v_brain$E) %in% brain_samples_to_be_removed]
dim(v_brain$E) #21356   225
brain_metadata <- brain_metadata[!rownames(brain_metadata) %in% brain_samples_to_be_removed, ]
dim(brain_metadata) #225 339
brain_only_metadata <- brain_metadata[, grepl("_brain$", names(brain_metadata))]
brain_only_metadata <- brain_only_metadata[!rownames(brain_only_metadata) %in% brain_samples_to_be_removed, ]
dim(brain_only_metadata) #225 169
brain_only_metadata$IID_ISMMS <- brain_metadata$IID_ISMMS[match(rownames(brain_only_metadata), rownames(brain_metadata))]
dim(brain_only_metadata) #225 170

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

pca_blood <- prcomp(t(blood_form5), center = TRUE, scale. = TRUE)

# Transpose gene expression matrices so samples are rows
#blood_expr_t <- pca_blood$x    # 225 samples x 21046 blood genes
blood_expr_t <- t(v_blood$E)   # 225 samples x 21046 blood genes
brain_expr_t <- t(brain_full)   # 225 samples x 21356 brain genes

# Set up parallel backend
plan(multisession, workers = 20)  # Adjust based on your machine

# Run elastic net models in parallel
#i=10 this is just to test
model_results <- future_map_dfr(1:ncol(brain_expr_t), function(i) { #applies a function to each brain gene and row-binds the results into a single data frame
  if (i %% 100 == 0) message("Processing brain gene ", i, " of ", ncol(brain_expr_t))
  
  y <- brain_expr_t[, i] #columns=brain gene, rows=samples
  x <- blood_expr_t #use all blood genes
  
  # Split 80% train, 20% test 
  n <- nrow(x) #225
  train_idx <- sample(1:n, size = floor(0.8 * n)) #indices of training samples | randomly selects 1:n indice, floor to round (0.8*n) down
  test_idx <- setdiff(1:n, train_idx)
  
  x_train <- x[train_idx, ]
  y_train <- y[train_idx]
  
  x_test <- x[test_idx, ]
  y_test <- y[test_idx]
  
  # Cross-validated elastic net on training set
  cv_fit <- cv.glmnet(x_train, y_train, alpha = 0.5, nfolds = 5) #can add type.measure="mse" | Find the lambda that minimizes the CV error
  #cv_fit$lambda      cv_fit$cvup        cv_fit$call        cv_fit$lambda.min
  #cv_fit$cvm         cv_fit$cvlo        cv_fit$name        cv_fit$lambda.1se
  #cv_fit$cvsd        cv_fit$nzero       cv_fit$glmnet.fit  cv_fit$index
  
  best_lambda <- cv_fit$lambda.min #selecting the best lambda that gives the lowest CV error
  # Predictions
  preds_train <- predict(cv_fit, s = best_lambda, newx = x_train)
  preds_test <- predict(cv_fit, s = best_lambda, newx = x_test)
  
  # R² calculations
  r2_train <- if (sd(y_train) == 0 || sd(preds_train) == 0) NA else cor(y_train, preds_train)^2 #if either vector of sd of actual and predicted training values have 0 variance,return NAs, otherwise calculate R^2
  r2_test <- if (sd(y_test) == 0 || sd(preds_test) == 0) NA else cor(y_test, preds_test)^2
  
  tibble( #df but compact version
    brain_gene = colnames(brain_expr_t)[i],
    best_lambda = best_lambda,
    r2_train = r2_train,
    r2_test = r2_test
  )
}, .options = furrr_options(seed = TRUE)) #tells furrr to use reproducible random sampling, this is important because I'm running code in parallel,  in parallel, and I used function like sample(): I want the randomness to behave as if it were done sequentially, with a fixed seed."

write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_base_blood_PC_full_brain_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_base_blood_PC_base_brain_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_base_brain_no_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_full_brain_no_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_full_brain_with_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_base_brain_with_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_full_brain_with_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_full_brain_no_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_base_brain_no_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_base_brain_with_residID_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(model_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_brain_predicted_from_blood.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Save full model objects in case you want to predict later
model_objects <- future_map(1:ncol(brain_expr_t), function(i) {
  if (i %% 100 == 0) message("Fitting model for brain gene ", i, " of ", ncol(brain_expr_t))
  
  y <- brain_expr_t[, i]
  x <- blood_expr_t
  
  # Return the fitted cross-validated elastic net model
  cv.glmnet(x, y, alpha = 0.5, nfolds = 5)
}, .options = furrr_options(seed = TRUE))

# Save all models in a single RDS file
saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_base_blood_PC_objects_all_full_brain_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_base_blood_PC_objects_all_base_brain_genes.rds")

#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form5_PC_objects_all_base_brain_no_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form5_PC_objects_all_full_brain_no_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form5_PC_objects_all_full_brain_with_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form5_PC_objects_all_base_brain_with_residID_genes.rds")

#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_PC_objects_all_full_brain_no_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_PC_objects_all_base_brain_no_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_PC_objects_all_base_brain_with_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_objects_all_base_brain_with_residID_genes.rds")
#saveRDS(model_objects, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_objects_all_brain_genes.rds")

library(data.table)
############### PCs ################### -- this runs really fast, like less than 10 mins
#######################################
################ BASE BLOOD AND BASE BRAIN ############### 
## this is all genes results not PC -- change names please
elastic_PC_base_blood_base_brain <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_base_blood_PC_base_brain_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_blood_base_brain$r2_train)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.063   0.157   0.253   0.310   0.411   0.998   11509
summary(elastic_PC_base_blood_base_brain$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.011   0.038   0.035   0.970   11509

elastic_PC_base_blood_full_brain <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_base_blood_PC_full_brain_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_blood_full_brain$r2_train)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.056   0.125   0.187   0.223   0.281   0.888   13783
summary(elastic_PC_base_blood_full_brain$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.011   0.026   0.033   0.513   13783

################ FORM 3 ############### 
elastic_PC_base_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_base_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_with_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.038   0.076   0.137   0.204   0.269   0.981   14177 
summary(elastic_PC_base_with_residID$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.011   0.023   0.031   0.348   14177 

elastic_PC_base_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_base_brain_no_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_no_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.039   0.091   0.175   0.240   0.330   0.986   12656 
summary(elastic_PC_base_no_residID$r2_test)
#.  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.011   0.026   0.032   0.810   12656

elastic_PC_full_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_full_brain_no_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_full_no_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.039   0.096   0.182   0.246   0.345   0.993   11778
summary(elastic_PC_full_no_residID$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.000   0.002   0.010   0.022   0.029   0.325   11778 

elastic_PC_full_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_PC_full_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_full_with_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.038   0.093   0.178   0.241   0.336   0.979   12039
summary(elastic_PC_full_with_residID$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.011   0.024   0.031   0.449   12039 

################ FORM 5 ############### 
elastic_PC_base_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_base_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_with_residID$r2_train)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.040   0.085   0.159   0.227   0.312   0.999   12571
summary(elastic_PC_base_with_residID$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.011   0.023   0.030   0.552   12571 

elastic_PC_base_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_base_brain_no_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_base_no_residID$r2_train)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.040   0.089   0.172   0.236   0.328   0.992   12629
summary(elastic_PC_base_no_residID$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.010   0.023   0.029   0.670   12629 
  
elastic_PC_full_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_full_brain_no_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_full_no_residID$r2_train)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.040   0.092   0.172   0.236   0.329   0.988   11877 
summary(elastic_PC_full_no_residID$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.011   0.023   0.031   0.639   11877 

elastic_PC_full_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_PC_full_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_PC_full_with_residID$r2_train)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.039   0.092   0.177   0.242   0.340   0.999   11740
summary(elastic_PC_full_with_residID$r2_test)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.010   0.023   0.031   0.377   11740


############### ALL GENES ###################
#############################################
elastic_base_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_base_brain_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_base_no_residID$r2_train)
#.  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.062   0.158   0.262   0.322   0.433   1.000   11760
summary(elastic_base_no_residID$r2_test)
#.   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.012   0.038   0.038   0.878   11760

elastic_base_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_base_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_base_with_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.058   0.144   0.230   0.272   0.356   0.990   13093 
summary(elastic_base_with_residID$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.002   0.011   0.027   0.032   0.866   13093 

elastic_full_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_full_brain_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_full_no_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.061   0.161   0.264   0.321   0.431   0.999   11733
summary(elastic_full_no_residID$r2_test) ##############highest
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.013   0.040   0.040   0.898   11733 

elastic_full_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_full_brain_with_residID_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_full_with_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.060   0.130   0.204   0.248   0.321   0.997   14149 
summary(elastic_full_with_residID$r2_test)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.011   0.027   0.034   0.538   14149 

model_objects <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_objects_all_base_brain_genes.rds")
str(model_objects)       # Overview of the structure
length(model_objects)    # Number of elements if it's a list
model_objects[[1]]       # Access the first model or object

# Load the models
model_objects <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_model_form3_objects_all_base_brain_genes.rds")
# Predict brain gene 42 using its fitted model
cv_fit_42 <- model_objects[[42]]
preds <- predict(cv_fit_42, s = cv_fit_42$lambda.min, newx = new_blood_expr)  # make sure dimensions match

