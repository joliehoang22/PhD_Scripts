library(energy)    # for dcor()
library(furrr)     # for future_map()
library(tibble)    # for cleaner output
library(future)
library(data.table)

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

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
#blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
#rownames(blood_form5) <- blood_form5$V1
#blood_form5$V1 <- NULL

v_brain$E <- v_brain$E[, !colnames(v_brain$E) %in% brain_samples_to_be_removed]
dim(v_brain$E) #21356   225
brain_metadata <- brain_metadata[!rownames(brain_metadata) %in% brain_samples_to_be_removed, ]
dim(brain_metadata) #225 339
brain_only_metadata <- brain_metadata[, grepl("_brain$", names(brain_metadata))]
brain_only_metadata <- brain_only_metadata[!rownames(brain_only_metadata) %in% brain_samples_to_be_removed, ]
dim(brain_only_metadata) #225 169
brain_only_metadata$IID_ISMMS <- brain_metadata$IID_ISMMS[match(rownames(brain_only_metadata), rownames(brain_metadata))]
dim(brain_only_metadata) #225 170

#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

# Transpose gene expression matrices so samples are rows
blood_expr_t <- t(blood_form3)   # 225 samples x 21046 blood genes
brain_expr_t <- t(v_brain$E)  # 225 samples x 21356 brain genes

# ========================
# 2. Subset to 5x5 Genes for Testing
# ========================
#blood_test <- blood_expr_t[, 1:5]   # First 5 blood genes
#brain_test <- brain_expr_t[, 1:5]   # First 5 brain genes

# Convert to data frames (columns = genes, rows = samples)
blood_df <- as.data.frame(blood_expr_t)
brain_df <- as.data.frame(brain_expr_t)

# ========================
# 3. Define Correlation Functions
# ========================
get_dcor <- function(x, y) {
  tryCatch(dcor(x, y), error = function(e) NA)
}

#get_signed_dcor <- function(x, y) {
#  dcor_val <- tryCatch(dcor(x, y), error = function(e) NA)
#  pearson_val <- tryCatch(cor(x, y), error = function(e) NA)
#  return(sign(pearson_val) * dcor_val)
#}

# ========================
# 4. Parallel Setup
# ========================
#plan(multisession, workers = 30)  # Use more workers for full runs

# Row-wise calculation functions: 1 brain gene vs all blood genes
get_dcor_row <- function(brain_gene_vec) {
  apply(blood_df, 2, function(blood_gene_vec) get_dcor(blood_gene_vec, brain_gene_vec))
}

#get_signed_dcor_row <- function(brain_gene_vec) {
#  apply(blood_df, 2, function(blood_gene_vec) get_signed_dcor(blood_gene_vec, brain_gene_vec))
#}

# ========================
# 5. Compute Distance Correlation Matrices
# ========================
dcor_matrix <- future_map_dfr(brain_df, get_dcor_row, .progress = TRUE)
#signed_dcor_matrix <- future_map_dfr(brain_df, get_signed_dcor_row, .progress = TRUE)

# Convert to data.frame before setting row/col names
dcor_matrix <- as.data.frame(dcor_matrix)
#signed_dcor_matrix <- as.data.frame(signed_dcor_matrix)

rownames(dcor_matrix) <- colnames(brain_expr_t)
colnames(dcor_matrix) <- colnames(blood_expr_t)

#rownames(signed_dcor_matrix) <- colnames(brain_expr_t)
#colnames(signed_dcor_matrix) <- colnames(blood_expr_t)

# ========================
# 6. Save Results
# ========================
write.table(dcor_matrix,file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_full_brain_no_residID_dcor_matrix.txt",sep = "\t", row.names = TRUE, quote = FALSE)

#write.table(dcor_matrix,file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_full_brain_with_residID_dcor_matrix.txt",sep = "\t", row.names = TRUE, quote = FALSE)
#write.table(dcor_matrix,file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_base_brain_no_residID_dcor_matrix.txt",sep = "\t", row.names = TRUE, quote = FALSE)
#write.table(dcor_matrix,file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_base_brain_with_residID_dcor_matrix.txt",sep = "\t", row.names = TRUE, quote = FALSE)

#write.table(signed_dcor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_base_brain_with_residID_signed_dcor_matrix.txt",sep = "\t", row.names = TRUE, quote = FALSE)
# ========================
# 7. Print Result
# ========================
#print("Distance Correlation Matrix:")
#print(round(dcor_matrix, 3))

#print("Signed Distance Correlation Matrix:")
#print(round(signed_dcor_matrix, 3))

dcor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/form3_base_brain_with_residID_dcor_matrix.txt", data.table = FALSE)
row.names(dcor) <- dcor$V1
dcor$V1 <- NULL
dcor[1:5,1:5]
dcor <- t(dcor)
dim(dcor) #21046 21356

blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_by_brain_cor[1:5,1:5]
blood_form3_by_brain_cor <- as.matrix(blood_form3_by_brain_cor)

identical(rownames(dcor), rownames(blood_form3_by_brain_cor)) #TRUE 
identical(colnames(dcor), colnames(blood_form3_by_brain_cor)) #TRUE 

cor_matrix <- cor(as.vector(dcor), as.vector(blood_form3_by_brain_cor), method = "spearman")
#0.7024809

#Running inside a screen session: 86036.im_working_here

library(xgboost)
library(dplyr)

# Select 5 random brain genes to test
brain_expr <- t(v_brain$E) 
X <- t(blood_form3) 
#test_brain_genes <- sample(colnames(brain_expr), 10); test_brain_genes

#top_n <- 100  # Number of features to keep

# Make sure sample IDs are aligned across X and brain_expr
#results <- list()

#X <- as.matrix(v_brain$E)  # predictors: blood gene expression

#for (gene in test_brain_genes) {
  cat("Processing gene:", gene, "\n")
  
  y <- brain_expr[, gene]
  
  # Remove samples with NA
  complete_idx <- complete.cases(X, y)
  X_complete <- X[complete_idx, ]
  y_complete <- y[complete_idx]
  
  n <- nrow(X_complete)
  set.seed(123)
  train_idx <- sample(1:n, size = floor(0.8 * n))  # 80% training
  
  X_train <- X_complete[train_idx, ]
  y_train <- y_complete[train_idx]
  
  X_test <- X_complete[-train_idx, ]
  y_test <- y_complete[-train_idx]
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test)
  
  # 5-fold cross-validation to tune nrounds (early stopping for safety)
  cv <- xgb.cv(
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    verbose = 0,
    objective = "reg:squarederror",
    metrics = "rmse"
  )
  
  best_nrounds <- cv$best_iteration
  
  # Train final model
  model <- xgboost(
    data = dtrain,
    nrounds = best_nrounds,
    objective = "reg:squarederror",
    verbose = 0
  )
  
  # Predictions
  preds_train <- predict(model, dtrain)
  preds_test <- predict(model, dtest)
  
  # R² calculation
  r2_train <- cor(preds_train, y_train)^2
  r2_test <- cor(preds_test, y_test)^2
  
  # Save results
  results[[gene]] <- list(
    r2_train = r2_train,
    r2_test = r2_test,
    model = model,
    best_nrounds = best_nrounds
  )
}

# Print R² for each gene
#cat("\nR² summary:\n")
#sapply(results, function(x) c(train = round(x$r2_train, 3), test = round(x$r2_test, 3)))

# Print test R² for each gene
#sapply(results, function(x) round(x$r2, 3))

#ENSG00000225648.5 ENSG00000068137.15 ENSG00000260563.3 ENSG00000012983.11
#train             0.998              0.998             0.733                  1
#test              0.031              0.017             0.001                  0
#ENSG00000155858.6
#train              1.00
#test               0.02

###################################
###################################
###################################
# Set dCor threshold
plan(multisession, workers = 30) 

dcor_threshold <- 0.125 #this is likely at 3rd quadrant
#summary(dcor)
#dcor_sub <- dcor[,1:5]
# Initialize results list
results <- list()
cat("XGBoost modeling started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#for (gene in test_brain_genes) {
for (gene in colnames(brain_expr)) {
  cat("Processing gene:", gene, "\n")
  
  y <- brain_expr[, gene]
  
  # 1. Get blood genes with dCor > threshold
  if (!gene %in% colnames(dcor)) next
  gene_dcor <- dcor[, gene]
  high_dcor_genes <- names(gene_dcor[gene_dcor > dcor_threshold])
  
  # 2. Filter to genes that exist in X
  selected_genes <- intersect(high_dcor_genes, colnames(X))
  if (length(selected_genes) < 5) {
    warning(paste("Too few features for", gene, "- skipping"))
    next
  }
  
  X_gene <- X[, selected_genes, drop = FALSE]
  
  # 3. Remove samples with NA
  complete_idx <- complete.cases(X_gene, y)
  X_complete <- X_gene[complete_idx, ]
  y_complete <- y[complete_idx]
  
  n <- nrow(X_complete)
  if (n < 20) {
    warning(paste("Too few samples for", gene, "- skipping"))
    next
  }
  
  # 4. Train/test split
  set.seed(123)
  train_idx <- sample(1:n, floor(0.8 * n))
  X_train <- X_complete[train_idx, ]
  y_train <- y_complete[train_idx]
  X_test <- X_complete[-train_idx, ]
  y_test <- y_complete[-train_idx]
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test)
  
  # 5. CV for tuning nrounds
  params <- list(
    objective = "reg:squarederror",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.7,
    colsample_bytree = 0.7,
    lambda = 1,
    alpha = 1
  )
  
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    verbose = 0
  )
  
  best_nrounds <- cv$best_iteration
  
  # 6. Final model
  model <- xgboost(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  # 7. Predict
  preds_train <- predict(model, dtrain)
  preds_test <- predict(model, dtest)
  
  r2_train <- cor(preds_train, y_train)^2
  r2_test <- cor(preds_test, y_test)^2
  
  results[[gene]] <- list(
    r2_train = r2_train,
    r2_test = r2_test,
    n_features = length(selected_genes),
    model = model
  )
}

# Summary output
#cat("\nR² summary for top 10 test genes:\n")
r2_summary <- sapply(results, function(x) c(train = round(x$r2_train, 3), test = round(x$r2_test, 3), features = x$n_features))

r2_df <- as.data.frame(t(r2_summary))
r2_df$gene <- rownames(r2_df)  # Add gene ID as a column

# Save as tab-delimited text
write.table(r2_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/xgboost_r2_summary_all_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("XGBoost modeling finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


#features = The number of blood genes used as input to predict each brain gene.

##################### CODE WITH PARALLELIZATION #####################
library(xgboost)
library(dplyr)

# Select 5 random brain genes to test
brain_expr <- t(v_brain$E) 
X <- t(blood_form3) 

plan(multisession, workers = 30) 
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB

dcor_threshold <- 0.125

# Log start time
cat("XGBoost modeling started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Set seed for reproducibility
set.seed(123)

# Parallelized loop over genes
results <- future_map(colnames(brain_expr), function(gene) {
  cat("Processing gene:", gene, "\n")
  
  y <- brain_expr[, gene]
  
  # 1. Get top 100 blood genes with highest dCor
  if (!gene %in% colnames(dcor)) return(NULL)
  gene_dcor <- dcor[, gene]
  top_dcor_genes <- names(sort(gene_dcor, decreasing = TRUE))[1:100]
  
  # 2. Filter to genes that exist in X
  selected_genes <- intersect(top_dcor_genes, colnames(X))
  if (length(selected_genes) < 5) {
    warning(paste("Too few features for", gene, "- skipping"))
    return(NULL)
  }
  
  X_gene <- X[, selected_genes, drop = FALSE]
  
  # 3. Remove samples with NA
  complete_idx <- complete.cases(X_gene, y)
  X_complete <- X_gene[complete_idx, ]
  y_complete <- y[complete_idx]
  
  n <- nrow(X_complete)
  if (n < 20) {
    warning(paste("Too few samples for", gene, "- skipping"))
    return(NULL)
  }
  
  # 4. Train/test split
  train_idx <- sample(1:n, floor(0.8 * n))
  X_train <- X_complete[train_idx, ]
  y_train <- y_complete[train_idx]
  X_test  <- X_complete[-train_idx, ]
  y_test  <- y_complete[-train_idx]
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test)
  
  # 5. CV to tune nrounds
  params <- list(
    objective = "reg:squarederror",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.7,
    colsample_bytree = 0.7,
    lambda = 1,
    alpha = 1
  )
  
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    verbose = 0
  )
  
  best_nrounds <- cv$best_iteration
  
  # 6. Final model
  model <- xgboost(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  # 7. Predict
  preds_train <- predict(model, dtrain)
  preds_test  <- predict(model, dtest)
  
  r2_train <- if (sd(preds_train) == 0 || sd(y_train) == 0) NA else cor(preds_train, y_train)^2
  r2_test  <- if (sd(preds_test) == 0 || sd(y_test) == 0) NA else cor(preds_test, y_test)^2
  
  list(
    gene = gene,
    r2_train = r2_train,
    r2_test = r2_test,
    n_features = length(selected_genes),
    model = model
  )
}, .options = furrr_options(seed = TRUE, globals = TRUE))

#Error in future::getGlobalsAndPackages(fn, envir = env_fn, globals = TRUE) : 
#  The total size of the 4 globals exported for future expression (‘function (gene); 
#{; cat("Processing gene:", gene, "\n"); y <- brain_expr[, gene]; if (!gene %in% colnames(dcor)); 
#return(NULL); ...; list(gene = gene, r2_train = r2_train, r2_test = r2_test,; 
#n_features = length(selected_genes), model = model); }’) is 3.42 GiB. 
#This exceeds the maximum allowed size 500.00 MiB per by R option "future.globals.maxSize". 
#This limit is set to protect against transfering too large objects to parallel workers by mistake, 
#which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") 
#for further explainations and how to adjust or remove this threshold 
#The three largest globals are ‘dcor’ (3.35 GiB of class ‘numeric’), ‘brain_expr’ (37.19 MiB of class ‘numeric’) 
#and ‘X’ (36.65 MiB of class ‘numeric’)

# Convert to summary dataframe
r2_summary <- map_dfr(results, function(x) {
  if (is.null(x)) return(NULL)
  tibble(
    gene = x$gene,
    r2_train = round(x$r2_train, 3),
    r2_test = round(x$r2_test, 3),
    n_features = x$n_features
  )
})

# Save summary results
write.table(
  r2_summary,
  file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/xgboost_r2_summary_all_genes2.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Save model objects
saveRDS(
  results,
  file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/xgboost_model_objects_all_genes2.rds"
)

# Log end time
cat("XGBoost modeling finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

############# ############# ############# ############# #############
############# TESTING ONE GENE AT A TIME  ############# #############
############# ############# ############# ############# #############
# Choose one gene to test
gene <- colnames(brain_expr)[1]
cat("Processing gene:", gene, "\n")

# Extract expression for the gene
y <- brain_expr[, gene]

dcor_threshold <- 0.125

# Step 1: Get blood genes with distance correlation > threshold
if (!gene %in% colnames(dcor)) stop("Gene not found in dcor")
gene_dcor <- dcor[, gene]
high_dcor_genes <- names(gene_dcor[gene_dcor > dcor_threshold])
cat("Number of high-dCor genes:", length(high_dcor_genes), "\n")

# Step 2: Filter to genes that are in X
selected_genes <- intersect(high_dcor_genes, colnames(X))
if (length(selected_genes) < 5) stop("Too few features for this gene")

X_gene <- X[, selected_genes, drop = FALSE]

# Step 3: Remove samples with missing data
complete_idx <- complete.cases(X_gene, y)
X_complete <- X_gene[complete_idx, ]
y_complete <- y[complete_idx]

n <- nrow(X_complete)
if (n < 20) stop("Too few samples after removing NA")

# Step 4: Split into train and test sets (80/20)
set.seed(123)
train_idx <- sample(1:n, floor(0.8 * n))
X_train <- X_complete[train_idx, ]
y_train <- y_complete[train_idx]
X_test  <- X_complete[-train_idx, ]
y_test  <- y_complete[-train_idx]

# Step 5: Prepare data for xgboost
dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test)

# Step 6: Cross-validation to choose best nrounds
params <- list(
  objective = "reg:squarederror",
  max_depth = 3,
  eta = 0.1,
  subsample = 0.7,
  colsample_bytree = 0.7,
  lambda = 1,
  alpha = 1
)

cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 1000,
  nfold = 5,
  early_stopping_rounds = 10,
  verbose = 1
)

best_nrounds <- cv$best_iteration
cat("Best nrounds from CV:", best_nrounds, "\n")

#Multiple eval metrics are present. Will use test_rmse for early stopping.
#Will train until test_rmse hasn't improved in 10 rounds.

##this is with dcor threshold 0.125
#Stopping. Best iteration:
#[89]	train-rmse:0.066535+0.000820	test-rmse:0.315296+0.029833

#this is with top 100 genes
#[65]	train-rmse:0.146873+0.004808	test-rmse:0.293073+0.043839

# Step 7: Train final model
model <- xgboost(
  params = params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

# Step 8: Predict and calculate R²
preds_train <- predict(model, dtrain)
preds_test <- predict(model, dtest)

r2_train <- if (sd(preds_train) == 0 || sd(y_train) == 0) NA else cor(preds_train, y_train)^2
r2_test  <- if (sd(preds_test) == 0 || sd(y_test) == 0) NA else cor(preds_test, y_test)^2

# Output
cat("R² (train):", r2_train, "\n")
cat("R² (test):", r2_test, "\n")
cat("Number of features used:", length(selected_genes), "\n")

## for one gene, ENSG00000001497.16
#[65]	train-rmse:0.146873+0.004808	test-rmse:0.293073+0.043839
#Best nrounds from CV: 65 
#R² (train): 0.8105018 
#R² (test): 0.1303048 
#Number of features used: 100

## for another gene, ENSG00000001497.16
#[63]	train-rmse:0.099462+0.002524	test-rmse:0.425381+0.063182
#Best nrounds from CV: 63 
#R² (train): 0.9616129 
#R² (test): 0.04784178 
#Number of features used: 6748
































