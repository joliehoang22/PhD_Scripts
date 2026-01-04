library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(variancePartition)
library(matrixStats)
library(purrr)
library(tidyr)
library(furrr)      
library(future)    
library(purrr)
library(glmnet)   
library(future.apply)
library(nestedcv, lib.loc="/sc/arion/projects/mscic1/results/jolie/Rlibrary")
set.seed(2025)

## Job submussions
#path=/sc/arion/projects/mscic1/results/jolie/LBP/scripts/
#cd $path
#ml R
#bsub -q premium -P acc_mscic1 -n 20 -W 144:00 -R rusage[mem=4000] -R span[hosts=1] -R himem -o %J.stdout -eo %J.stderr Rscript ${path}nestedCV_job_form3_no_residID_base_brain.R

#Job <199223798> is submitted to queue <premium>.

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

# Transpose data so samples are rows and genes are columns (required format)
brain_t <- t(v_brain$E)  # Now 225 samples x 21356 genes
blood_t <- t(blood_form3) # Now 225 samples x 21046 genes

# Get all brain genes and their variances
all_brain_genes <- colnames(brain_t)
gene_vars <- apply(brain_t, 2, var)

# Initialize results storage
results <- list()
failed_genes <- character()

cat("Starting prediction for", length(all_brain_genes), "brain genes...\n")
t0 <- Sys.time()

# Loop through all brain genes
for (i in seq_along(all_brain_genes)) {
  target_gene <- all_brain_genes[i]
  
  # Progress update every 100 genes
  if (i %% 100 == 0) {
    cat("Processing gene", i, "of", length(all_brain_genes), ":", target_gene, "\n")
  }
  
  # Skip genes with very low variance to avoid degenerate fits
  if (gene_vars[target_gene] < 1e-6) {
    failed_genes <- c(failed_genes, target_gene)
    next
  }
  
  tryCatch({
    set.seed(123)
    res <- nestcv.glmnet(
      y = as.numeric(brain_t[, target_gene]),
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      alphaSet = seq(0.1, 1, 0.1),
      cv.cores = 20, ##change this too for //
      filterFUN = correl_filter,
      filter_options = list(
        nfilter = min(500, ncol(blood_t)),
        method = "spearman"
      )
    )
    
    # Store results
    results[[target_gene]] <- list(
      gene = target_gene,
      cv_performance = res$summary,
      coefficients = res$final_fit$beta,
      intercept = res$final_fit$a0,
      lambda = res$final_fit$lambda
    )
    
  }, error = function(e) {
    cat("Failed for gene:", target_gene, "- Error:", e$message, "\n")
    failed_genes <<- c(failed_genes, target_gene)
  })
}

t1 <- Sys.time()
cat("Total runtime:", round(difftime(t1, t0, units = "hours"), 2), "hours\n")

# Summary
successful_genes <- length(results)
cat("Successfully processed:", successful_genes, "genes\n")
cat("Failed genes:", length(failed_genes), "\n")

# Save results
saveRDS(results, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_no_residID_results_20250818.rds")
saveRDS(failed_genes, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_no_residID_failed_genes_20250818.rds")

# Extract performance metrics for all genes
performance_df <- do.call(rbind, lapply(names(results), function(gene) {
  data.frame(
    gene = gene,
    cor = results[[gene]]$cv_performance["cor"],
    rmse = results[[gene]]$cv_performance["rmse"],
    mae = results[[gene]]$cv_performance["mae"]
  )
}))

write.csv(performance_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_no_residID_prediction_performance.csv", row.names = FALSE)

cat("Results saved to:\n")
cat("- brain_gene_prediction_results.rds (full results)\n")
cat("- gene_prediction_performance.csv (performance metrics)\n")
cat("- failed_genes.rds (genes that failed)\n")

nestedcv_brain_base_blood_no_residID_failed_genes_20250818.rds
nestedcv_brain_base_blood_no_residID_results_20250818.rds

nestedcv_brain_base_blood_with_residID_failed_genes_20250818.rds
nestedcv_brain_base_blood_with_residID_results_20250818.rds

nestedcv_brain_full_blood_with_residID_failed_genes_20250818.rds
nestedcv_brain_full_blood_with_residID_results_20250818.rds

##this is probably brain_full_no_residID hereee
nestedcv_brain_gene_prediction_results_20250818.rds
nestedcv_failed_genes_20250818.rds

nestedcv_try_20250818.rds

fail <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_no_residID_failed_genes_20250818.rds")

data <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_no_residID_results_20250818.rds")
head(data, 3)

performance_df <- data.frame(
  gene = names(data),
  RMSE = sapply(data, function(x) x$cv_performance["RMSE"]),
  R_squared = sapply(data, function(x) x$cv_performance["R.squared"]),
  Pearson_r2 = sapply(data, function(x) x$cv_performance["Pearson.r^2"]),
  MAE = sapply(data, function(x) x$cv_performance["MAE"]),
  stringsAsFactors = FALSE
)
rownames(performance_df) <- NULL
dim(performance_df) #21356     5
head(performance_df)

summary(performance_df$R_squared)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.6192 -0.2350 -0.1708 -0.1684 -0.1091  0.9437 
summary(performance_df$RMSE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1194  0.3650  0.5606  0.6877  0.8793  3.6805 
summary(performance_df$Pearson_r2)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.001143 0.005051 0.016098 0.015154 0.945887 
summary(performance_df$MAE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.09256 0.28905 0.44095 0.54124 0.69270 3.11827 

good <- performance_df[performance_df$R_squared > 0.01, ]
dim(good) #655   5
summary(good)
#      gene                RMSE          R_squared         Pearson_r2     
#  Length:655         Min.   :0.1269   Min.   :0.01017   Min.   :0.04307  
#  Class :character   1st Qu.:0.3487   1st Qu.:0.02791   1st Qu.:0.07210  
#  Mode  :character   Median :0.5551   Median :0.06162   Median :0.09813  
#                     Mean   :0.6708   Mean   :0.19002   Mean   :0.21705  
#                     3rd Qu.:0.9050   3rd Qu.:0.26082   3rd Qu.:0.26523  
#                     Max.   :3.6805   Max.   :0.94366   Max.   :0.94589  
#       MAE        
#  Min.   :0.1032  
#  1st Qu.:0.2729  
#  Median :0.4298  
#  Mean   :0.5204  
#  3rd Qu.:0.6886  
#  Max.   :3.1183

## MAYBE THIS IS NOT WRONG - THERE IS A SMALL SUBSET WITH GENUINE PREDICTIVE RELATIONSHIP BUT MOST GENES DO NOT HAVE ONE?

data <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_with_residID_results_20250818.rds")
head(data, 3)
performance_df <- data.frame(
  gene = names(data),
  RMSE = sapply(data, function(x) x$cv_performance["RMSE"]),
  R_squared = sapply(data, function(x) x$cv_performance["R.squared"]),
  Pearson_r2 = sapply(data, function(x) x$cv_performance["Pearson.r^2"]),
  MAE = sapply(data, function(x) x$cv_performance["MAE"]),
  stringsAsFactors = FALSE
)
rownames(performance_df) <- NULL
dim(performance_df) #21356     5
head(performance_df)

summary(performance_df$R_squared)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.6544 -0.2829 -0.2151 -0.2192 -0.1533  0.7852 
summary(performance_df$RMSE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1200  0.3741  0.5724  0.7037  0.8962  4.0065 
summary(performance_df$Pearson_r2)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0007627 0.0033814 0.0088273 0.0099823 0.7945640 
summary(performance_df$MAE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09219 0.29704 0.45221 0.55623 0.71172 3.60164 

good <- performance_df[performance_df$R_squared > 0.01, ]
dim(good) #129   5
summary(good)
#      gene                RMSE          R_squared         Pearson_r2     
#  Length:129         Min.   :0.1697   Min.   :0.01084   Min.   :0.04412  
#  Class :character   1st Qu.:0.5410   1st Qu.:0.03407   1st Qu.:0.06973  
#  Mode  :character   Median :0.9104   Median :0.07976   Median :0.10935  
#                     Mean   :1.2256   Mean   :0.22229   Mean   :0.24703  
#                     3rd Qu.:1.7998   3rd Qu.:0.33653   3rd Qu.:0.33948  
#                     Max.   :3.2483   Max.   :0.78516   Max.   :0.79456  
#       MAE        
#  Min.   :0.1297  
#  1st Qu.:0.4230  
#  Median :0.7176  
#  Mean   :0.9570  
#  3rd Qu.:1.4461  
#  Max.   :2.6823

data <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_full_blood_with_residID_results_20250818.rds")
head(data, 3)

performance_df <- data.frame(
  gene = names(data),
  RMSE = sapply(data, function(x) x$cv_performance["RMSE"]),
  R_squared = sapply(data, function(x) x$cv_performance["R.squared"]),
  Pearson_r2 = sapply(data, function(x) x$cv_performance["Pearson.r^2"]),
  MAE = sapply(data, function(x) x$cv_performance["MAE"]),
  stringsAsFactors = FALSE
)
rownames(performance_df) <- NULL
dim(performance_df) #21356     5
head(performance_df)

summary(performance_df$R_squared)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.7984 -0.3164 -0.2421 -0.2477 -0.1740  0.3830 
summary(performance_df$RMSE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.08824 0.23789 0.36026 0.43613 0.55513 3.95084 
summary(performance_df$Pearson_r2)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000000 0.0009189 0.0043355 0.0096377 0.0122903 0.3852548
summary(performance_df$MAE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06781 0.18623 0.27988 0.33447 0.42385 3.45440 

good <- performance_df[performance_df$R_squared > 0.01, ]
dim(good) #108   5
summary(good)
#     gene                RMSE          R_squared         Pearson_r2     
# Length:108         Min.   :0.1109   Min.   :0.01013   Min.   :0.03361  
# Class :character   1st Qu.:0.3678   1st Qu.:0.02567   1st Qu.:0.06635  
# Mode  :character   Median :0.5729   Median :0.05589   Median :0.09096  
#                    Mean   :0.8002   Mean   :0.08828   Mean   :0.11536  
#                    3rd Qu.:1.0157   3rd Qu.:0.12398   3rd Qu.:0.14348  
#                    Max.   :3.0137   Max.   :0.38296   Max.   :0.38525  
#      MAE         
# Min.   :0.08378  
# 1st Qu.:0.28846  
# Median :0.44265  
# Mean   :0.61727  
# 3rd Qu.:0.75317  
# Max.   :2.57974

data <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_results_20250818.rds")
head(data, 3)

performance_df <- data.frame(
  gene = names(data),
  RMSE = sapply(data, function(x) x$cv_performance["RMSE"]),
  R_squared = sapply(data, function(x) x$cv_performance["R.squared"]),
  Pearson_r2 = sapply(data, function(x) x$cv_performance["Pearson.r^2"]),
  MAE = sapply(data, function(x) x$cv_performance["MAE"]),
  stringsAsFactors = FALSE
)
rownames(performance_df) <- NULL
dim(performance_df) #21356     5
head(performance_df)

summary(performance_df$R_squared)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.6131 -0.2643 -0.1978 -0.1947 -0.1324  0.8559
summary(performance_df$RMSE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.08519 0.23151 0.35166 0.42510 0.54541 3.59100
summary(performance_df$Pearson_r2)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.000962 0.004316 0.014262 0.013065 0.857572
summary(performance_df$MAE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06311 0.18156 0.27257 0.32559 0.41553 3.10795

### SOME DX 
# Load your data
data <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_results_20250818.rds")

# Function to check if coefficients and intercept are not NULL
has_model_params <- function(gene_result) {
  return(!is.null(gene_result$coefficients) && !is.null(gene_result$intercept))
}

# Find genes with non-NULL coefficients and intercepts
genes_with_models <- sapply(data, has_model_params)
valid_gene_names <- names(genes_with_models)[genes_with_models]

## ALL 0!! BAD BAD BAD.
# Simple test run with 2 brain genes
library(nestedcv, lib.loc="/sc/arion/projects/mscic1/results/jolie/Rlibrary")
library(glmnet)
library(data.table)

# Load data
blood_form3 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt", data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

# Transpose data
brain_t <- t(v_brain$E)
blood_t <- t(blood_form3)

# Get first 2 brain genes for testing
test_genes <- colnames(brain_t)[1:2]
cat("Testing genes:", test_genes, "\n")

# Simple extraction function
extract_coefficients <- function(nestcv_result) {
  if (is.null(nestcv_result$final_fit)) {
    return(list(coefficients = NULL, intercept = NULL))
  }
  
  # Try the most common method: cv.glmnet with lambda.min
  tryCatch({
    coefs <- coef(nestcv_result$final_fit, s = "lambda.min")
    intercept <- as.numeric(coefs[1, 1])
    coefficients <- as.matrix(coefs[-1, 1, drop = FALSE])
    return(list(coefficients = coefficients, intercept = intercept))
  }, error = function(e) {
    return(list(coefficients = NULL, intercept = NULL))
  })
}

# Test run with 2 genes
results <- list()

for (i in 1:2) {
  target_gene <- test_genes[i]
  cat("\nProcessing gene", i, ":", target_gene, "\n")
  
  tryCatch({
    set.seed(123)
    res <- nestcv.glmnet(
      y = as.numeric(brain_t[, target_gene]),
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      alphaSet = seq(0.1, 1, 0.1),
      cv.cores = 1,  # Use 1 core for testing
      filterFUN = correl_filter,
      filter_options = list(
        nfilter = min(1000, ncol(blood_t)),
        method = "spearman"
      )
    )
    
    # Debug: Check what we got
    cat("  Final fit class:", class(res$final_fit), "\n")
    cat("  Final fit is NULL:", is.null(res$final_fit), "\n")
    
    # Extract coefficients
    model_info <- extract_coefficients(res)
    
    # Store results
    results[[target_gene]] <- list(
      gene = target_gene,
      cv_performance = res$summary,
      coefficients = model_info$coefficients,
      intercept = model_info$intercept
    )
    
    # Print results
    cat("  Performance:\n")
    print(res$summary)
    cat("  Coefficients extracted:", !is.null(model_info$coefficients), "\n")
    cat("  Intercept extracted:", !is.null(model_info$intercept), "\n")
    
    if (!is.null(model_info$coefficients)) {
      cat("  Non-zero coefficients:", sum(model_info$coefficients != 0), "\n")
      cat("  Intercept value:", model_info$intercept, "\n")
    }
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
}

# Summary
cat("\n=== SUMMARY ===\n")
for (i in 1:length(results)) {
  gene <- names(results)[i]
  has_coef <- !is.null(results[[gene]]$coefficients)
  has_int <- !is.null(results[[gene]]$intercept)
  cat("Gene", i, "(", gene, "): Coefficients =", has_coef, ", Intercept =", has_int, "\n")
}

      RMSE     R.squared   Pearson.r^2           MAE   
   0.308733     -0.270371      0.003579      0.243593


## start
# Compare improved modeling approaches on 500 genes
library(nestedcv)
library(glmnet)
library(data.table)

# Load data
blood_form3 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt", data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

brain_t <- t(v_brain$E)
blood_t <- t(blood_form3)

# Get first 500 brain genes for testing
test_genes <- colnames(brain_t)[1:500]
cat("Testing", length(test_genes), "brain genes\n")

# Extract coefficients function
extract_coefficients <- function(nestcv_result) {
  if (is.null(nestcv_result$final_fit)) {
    return(list(coefficients = NULL, intercept = NULL))
  }
  
  tryCatch({
    coefs <- coef(nestcv_result$final_fit, s = "lambda.min")
    intercept <- as.numeric(coefs[1, 1])
    coefficients <- as.matrix(coefs[-1, 1, drop = FALSE])
    return(list(coefficients = coefficients, intercept = intercept))
  }, error = function(e) {
    return(list(coefficients = NULL, intercept = NULL))
  })
}

# Custom variance-based filter function
variance_filter <- function(y, x, nfilter = 1000) {
  # Remove low-variance features first (bottom 10%)
  var_threshold <- quantile(apply(x, 2, var), 0.1)
  high_var_features <- which(apply(x, 2, var) > var_threshold)
  
  # Then select by correlation among high-variance features
  if (length(high_var_features) > nfilter) {
    x_filtered <- x[, high_var_features]
    cors <- abs(cor(y, x_filtered, use = "complete.obs"))
    top_indices <- order(cors, decreasing = TRUE)[1:nfilter]
    selected_features <- high_var_features[top_indices]
  } else {
    selected_features <- high_var_features
  }
  
  return(selected_features)
}

# =============================================================================
# APPROACH 1: LESS AGGRESSIVE REGULARIZATION + 1000 FEATURES
# =============================================================================
cat("Starting APPROACH 1: Less Aggressive Regularization + 1000 Features\n")

results_approach1 <- list()
failed_genes_approach1 <- character()
start_time <- Sys.time()

for (i in seq_along(test_genes)) {
  target_gene <- test_genes[i]
  
  if (i %% 50 == 0) {
    cat("Approach 1 - Processing gene", i, "of", length(test_genes), ":", target_gene, "\n")
  }
  
  tryCatch({
    set.seed(123)
    res <- nestcv.glmnet(
      y = as.numeric(brain_t[, target_gene]),
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      # LESS AGGRESSIVE: More Ridge, less LASSO
      alphaSet = seq(0.01, 0.5, 0.1),
      cv.cores = 20,
      filterFUN = correl_filter,
      filter_options = list(
        nfilter = min(1000, ncol(blood_t)),  # 1000 features instead of 500
        method = "spearman"
      )
    )
    
    model_info <- extract_coefficients(res)
    
    results_approach1[[target_gene]] <- list(
      gene = target_gene,
      cv_performance = res$summary,
      coefficients = model_info$coefficients,
      intercept = model_info$intercept,
      approach = "less_aggressive_1000_features"
    )
    
  }, error = function(e) {
    cat("Approach 1 failed for gene:", target_gene, "- Error:", e$message, "\n")
    failed_genes_approach1 <- c(failed_genes_approach1, target_gene)
  })
}

approach1_time <- Sys.time() - start_time
cat("Approach 1 completed in", round(approach1_time, 2), "minutes\n")

# =============================================================================
# APPROACH 2: CUSTOM VARIANCE-BASED FEATURE SELECTION
# =============================================================================
cat("Starting APPROACH 2: Variance-Based Feature Selection\n")

results_approach2 <- list()
failed_genes_approach2 <- character()
start_time <- Sys.time()

for (i in seq_along(test_genes)) {
  target_gene <- test_genes[i]
  
  if (i %% 50 == 0) {
    cat("Approach 2 - Processing gene", i, "of", length(test_genes), ":", target_gene, "\n")
  }
  
  tryCatch({
    set.seed(123)
    res <- nestcv.glmnet(
      y = as.numeric(brain_t[, target_gene]),
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      # Regular alpha range
      alphaSet = seq(0.1, 1, 0.1),
      cv.cores = 20,
      # CUSTOM VARIANCE-BASED FILTERING
      filterFUN = variance_filter,
      filter_options = list(nfilter = 1000)
    )
    
    model_info <- extract_coefficients(res)
    
    results_approach2[[target_gene]] <- list(
      gene = target_gene,
      cv_performance = res$summary,
      coefficients = model_info$coefficients,
      intercept = model_info$intercept,
      approach = "variance_based_selection"
    )
    
  }, error = function(e) {
    cat("Approach 2 failed for gene:", target_gene, "- Error:", e$message, "\n")
    failed_genes_approach2 <- c(failed_genes_approach2, target_gene)
  })
}

approach2_time <- Sys.time() - start_time
cat("Approach 2 completed in", round(approach2_time, 2), "minutes\n")

# =============================================================================
# CREATE PERFORMANCE COMPARISON DATAFRAMES
# =============================================================================
cat("Creating performance comparison dataframes...\n")

# Function to extract performance metrics
extract_performance <- function(results_list, approach_name) {
  performance_df <- data.frame(
    gene = character(),
    approach = character(),
    RMSE = numeric(),
    R_squared = numeric(),
    Pearson_r2 = numeric(),
    MAE = numeric(),
    n_coefficients = integer(),
    stringsAsFactors = FALSE
  )
  
  for (gene in names(results_list)) {
    gene_data <- results_list[[gene]]
    perf <- gene_data$cv_performance
    
    n_coeff <- ifelse(!is.null(gene_data$coefficients), 
                      sum(gene_data$coefficients != 0, na.rm = TRUE), 
                      0)
    
    performance_df <- rbind(performance_df, data.frame(
      gene = gene,
      approach = approach_name,
      RMSE = as.numeric(perf["RMSE"]),
      R_squared = as.numeric(perf["R.squared"]),
      Pearson_r2 = as.numeric(perf["Pearson.r^2"]),
      MAE = as.numeric(perf["MAE"]),
      n_coefficients = n_coeff,
      stringsAsFactors = FALSE
    ))
  }
  
  return(performance_df)
}

# Extract performance for both approaches
perf_approach1 <- extract_performance(results_approach1, "less_aggressive_1000_features")
perf_approach2 <- extract_performance(results_approach2, "variance_based_selection")

# Combine performance data
combined_performance <- rbind(perf_approach1, perf_approach2)

# =============================================================================
# SAVE RESULTS
# =============================================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/"

# Save detailed results
saveRDS(results_approach1, paste0(save_dir, "approach1_less_aggressive_1000feat_results_", timestamp, ".rds"))
saveRDS(results_approach2, paste0(save_dir, "approach2_variance_based_results_", timestamp, ".rds"))

# Save performance comparison
saveRDS(combined_performance, paste0(save_dir, "performance_comparison_500genes_", timestamp, ".rds"))

# Save failed genes
saveRDS(failed_genes_approach1, paste0(save_dir, "approach1_failed_genes_", timestamp, ".rds"))
saveRDS(failed_genes_approach2, paste0(save_dir, "approach2_failed_genes_", timestamp, ".rds"))

# =============================================================================
# PERFORMANCE SUMMARY AND COMPARISON
# =============================================================================
cat("\n=== PERFORMANCE SUMMARY ===\n")

# Summary statistics for each approach
cat("APPROACH 1 (Less Aggressive + 1000 Features):\n")
approach1_summary <- perf_approach1[, c("R_squared", "Pearson_r2", "RMSE", "MAE")]
print(summary(approach1_summary))

cat("\nAPPROACH 2 (Variance-Based Selection):\n")
approach2_summary <- perf_approach2[, c("R_squared", "Pearson_r2", "RMSE", "MAE")]
print(summary(approach2_summary))

# Count genes with positive R-squared
pos_r2_approach1 <- sum(perf_approach1$R_squared > 0, na.rm = TRUE)
pos_r2_approach2 <- sum(perf_approach2$R_squared > 0, na.rm = TRUE)

cat("\nGenes with positive R-squared:\n")
cat("Approach 1:", pos_r2_approach1, "out of", nrow(perf_approach1), "\n")
cat("Approach 2:", pos_r2_approach2, "out of", nrow(perf_approach2), "\n")

# Best performing genes from each approach
cat("\nTop 5 genes by R-squared - Approach 1:\n")
top5_approach1 <- perf_approach1[order(perf_approach1$R_squared, decreasing = TRUE)[1:5], ]
print(top5_approach1[, c("gene", "R_squared", "Pearson_r2", "RMSE", "n_coefficients")])

cat("\nTop 5 genes by R-squared - Approach 2:\n")
top5_approach2 <- perf_approach2[order(perf_approach2$R_squared, decreasing = TRUE)[1:5], ]
print(top5_approach2[, c("gene", "R_squared", "Pearson_r2", "RMSE")])

# Save summary
summary_stats <- list(
  approach1_summary = summary(approach1_summary),
  approach2_summary = summary(approach2_summary),
  positive_r2_counts = c(approach1 = pos_r2_approach1, approach2 = pos_r2_approach2),
  top5_approach1 = top5_approach1,
  top5_approach2 = top5_approach2,
  failed_counts = c(approach1 = length(failed_genes_approach1), approach2 = length(failed_genes_approach2))
)

saveRDS(summary_stats, paste0(save_dir, "summary_statistics_", timestamp, ".rds"))

cat("\nAll results saved to:", save_dir, "\n")
cat("Files saved with timestamp:", timestamp, "\n")
cat("\nSuccessfully processed:\n")
cat("Approach 1:", length(results_approach1), "genes\n")
cat("Approach 2:", length(results_approach2), "genes\n")

# Load your results (replace TIMESTAMP with actual timestamp)
results <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/approach1_less_aggressive_1000feat_results_20250821_2016.rds")

# Get the top gene
top_gene <- "ENSG00000012817.15"
gene_data <- results[[top_gene]]

# See the coefficients
coeffs <- gene_data$coefficients
nonzero_coeffs <- coeffs[coeffs != 0, , drop = FALSE]

print(paste("Gene:", top_gene))
print(paste("R-squared:", 0.917))
print(paste("Non-zero coefficients:", length(nonzero_coeffs)))

# Show the actual coefficients (blood gene names and their values)
print("Top coefficients:")
print(head(nonzero_coeffs[order(abs(nonzero_coeffs), decreasing = TRUE), , drop = FALSE], 10))

                           s1
ENSG00000126012.11 -0.4133055
ENSG00000005889.15 -0.3707701
ENSG00000183943.5  -0.3222381
ENSG00000100297.16  0.2188596