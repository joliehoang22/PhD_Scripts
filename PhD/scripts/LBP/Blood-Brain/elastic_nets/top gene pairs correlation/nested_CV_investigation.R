# 
# # View summary results
# head(results$summary_stats)
# 
# # View detailed results for specific brain gene
# results$detailed_results[results$detailed_results$brain_gene == "GENE_NAME", ]

library(data.table)
test<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_with_indivdualID_nestedCV_multiple_job.txt", data.table=FALSE)
test2<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID_nestedCV_multiple.txt", data.table=FALSE)
#outputs of _job is the identical to that without _job

dim(test) #106,780     25
identical(test,test2) #TRUE
colnames(test)
 [1] "detailed_results.brain_gene"         
 [2] "detailed_results.n_predictors"       
 [3] "detailed_results.r2_train"           
 [4] "detailed_results.r2_test"            
 [5] "detailed_results.cor_test_spearman"  
 [6] "detailed_results.cor_test_pearson"   
 [7] "detailed_results.best_alpha"         
 [8] "detailed_results.lambda_min"         
 [9] "detailed_results.best_inner_cv_score"
[10] "detailed_results.outer_fold"         
[11] "summary_stats.brain_gene"            
[12] "summary_stats.n_folds"               
[13] "summary_stats.mean_r2_test"          
[14] "summary_stats.sd_r2_test"            
[15] "summary_stats.mean_cor_spearman"     
[16] "summary_stats.sd_cor_spearman"       
[17] "summary_stats.mean_cor_pearson"      
[18] "summary_stats.sd_cor_pearson"        
[19] "summary_stats.mean_n_predictors"     
[20] "summary_stats.most_common_alpha"     
[21] "parameters.n_outer_folds"            
[22] "parameters.n_inner_folds"            
[23] "parameters.top_percent"              
[24] "parameters.alpha_values"             
[25] "parameters.seed"

summary(test)
summary(test2)


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

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

# Transpose data so samples are rows and genes are columns (required format)
brain_t <- t(brain_full)  # Now 225 samples x 21356 genes
blood_t <- t(blood_form3) # Now 225 samples x 21046 genes

# Find common genes between brain and blood data
common_genes <- intersect(colnames(brain_t), colnames(blood_t))
cat("Number of common genes:", length(common_genes), "\n")

# Filter to common genes only
brain_common <- brain_t[, common_genes]
blood_common <- blood_t[, common_genes]

###### TEST for ONE GENE
# pick a gene with high variance to avoid degenerate fits
gene_vars <- apply(brain_t, 2, var)
target_gene <- names(sort(gene_vars, decreasing = TRUE))[1]
cat("Target brain gene:", target_gene, "\n")

set.seed(123)
res <- nestcv.glmnet(
  y = as.numeric(brain_t[, target_gene]),  # response vector
  x = blood_t,                             # predictors (all blood genes)
  family = "gaussian",
  outer_method = "cv",
  n_outer_folds = 5,
  n_inner_folds = 5,
  alphaSet = seq(0.1, 1, 0.1),             # elastic-net alpha grid
  cv.cores = 1,                            # can increase later
  # feature selection *inside* each outer training fold
  filterFUN = correl_filter,
  filter_options = list(
    nfilter = min(500, ncol(blood_t)),     # cap features per fold
    method  = "spearman"                   
  )
)
res
#Result:
#       RMSE     R.squared   Pearson.r^2           MAE   
#     1.7189        0.7655        0.7685        1.4219

predSummary(res$output)  
#       RMSE     R.squared   Pearson.r^2           MAE   
#     1.7189        0.7655        0.7685        1.4219 

head(res$output)
#                      testy      predy
#LBPSEMA4BLOOD326 -1.7697217 -1.3046834
#LBPSEMA4BLOOD635 -2.8900203 -0.5652066
#LBPSEMA4BLOOD760 -2.6702488 -3.5006184
#LBPSEMA4BLOOD306 -0.3987434 -0.9250373
#LBPSEMA4BLOOD160  4.3159093  5.1423924
#LBPSEMA4BLOOD151  4.7894369  3.2197577
#Row name (e.g., LBPSEMA4BLOOD326): the sample ID (same for blood/brain after your alignment).
#testy: the observed (held-out) brain expression for the target gene in that sample.
#predy: the model’s prediction of that brain expression from blood genes, trained only on the corresponding outer-fold training data (so this is out-of-sample).

##this is fine | res$output uses the row names of the samples you fed into the model. 
#If your row names were blood IDs (e.g., LBPSEMA4BLOOD326), that’s what you’ll see—these are just sample keys. 
#It doesn’t affect correctness as long as brain_t and blood_t had identical row names (same subjects/timepoints).


################## SCALE UP ##################
###### TEST for ALL BRAIN GENES ##############
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

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

# Transpose data so samples are rows and genes are columns (required format)
brain_t <- t(brain_full)  # Now 225 samples x 21356 genes
blood_t <- t(blood_form3) # Now 225 samples x 21046 genes

plan(multisession, workers = 25)

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
      cv.cores = 25, ##change this too for //
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
saveRDS(results, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_results_20250818.rds")
saveRDS(failed_genes, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_failed_genes_20250818.rds")

# Extract performance metrics for all genes
performance_df <- do.call(rbind, lapply(names(results), function(gene) {
  data.frame(
    gene = gene,
    cor = results[[gene]]$cv_performance["cor"],
    rmse = results[[gene]]$cv_performance["rmse"],
    mae = results[[gene]]$cv_performance["mae"]
  )
}))

write.csv(performance_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_performance.csv", row.names = FALSE)

cat("Results saved to:\n")
cat("- brain_gene_prediction_results.rds (full results)\n")
cat("- gene_prediction_performance.csv (performance metrics)\n")
cat("- failed_genes.rds (genes that failed)\n")

## blood_form3_no_ID vs brain_full
saveRDS(results, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_results_20250818.rds")
saveRDS(failed_genes, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_failed_genes_20250818.rds")
write.csv(performance_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_gene_prediction_performance.csv", row.names = FALSE)







######################
target_genes <- colnames(brain_t)
n_targets <- length(target_genes)
cat("Number of target brain genes:", n_targets, "\n")

plan(multisession, workers = 25)

t0 <- Sys.time()

summary_all <- future_lapply(seq_along(target_genes), function(i) {
    library(tibble)
  g <- target_genes[i]
  if (i %% 100 == 0) message("Finished ", i, "/", n_targets, " (", round(100*i/n_targets,1), "%)")

  fit <- tryCatch(
    nestcv.glmnet(
      y = as.numeric(brain_t[, g]),
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      alphaSet = seq(0.1, 1, 0.1),
      cv.cores = 1,
      filterFUN = correl_filter,
      filter_options = list(
        nfilter = min(500, ncol(blood_t)),
        method  = "spearman"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$output) || nrow(fit$output) == 0) {
    return(tibble(gene = g, RMSE = NA_real_, R2 = NA_real_))
  }

  preds <- fit$output
  testy <- as.numeric(preds$testy)
  predy <- as.numeric(preds$predy)

  RMSE <- sqrt(mean((testy - predy)^2))
  R2   <- {
    RSS <- sum((testy - predy)^2)
    TSS <- sum((testy - mean(testy))^2)
    1 - RSS/TSS
  }

  tibble(gene = g, RMSE = RMSE, R2 = R2)
}, future.seed = TRUE) |>
  bind_rows()

t1 <- Sys.time()
cat("Total runtime:", round(difftime(t1, t0, units = "mins"), 2), "minutes\n")
summary_all %>% arrange(desc(R2)) %>% slice_head(n = 20)

saveRDS(summary_all, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_try_20250818.rds")


######################################### TRY 2
# ---- packages ----
library(data.table)
library(nestedcv)
library(glmnet)
library(future)
library(future.apply)
library(dplyr)
library(tibble)
library(purrr)

# ---- load & prepare (your code) ----
blood_form3 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt", data.table = FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL

brain_full <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table = FALSE)
rownames(brain_full) <- brain_full$V1
brain_full$V1 <- NULL

brain_t <- t(brain_full)   # 225 x 21356
blood_t <- t(blood_form3)  # 225 x 21046

# ---- targets & parallel ----
target_genes <- colnames(brain_t)
n_targets <- length(target_genes)
cat("Number of target brain genes:", n_targets, "\n")

plan(multisession, workers = 25)

# ---- timer ----
t0 <- Sys.time()

# (tiny helper just for coefficients; keeps things simple)
.extract_coefs <- function(fit) {
  # Try the most common places nestedcv stores the final cv.glmnet
  if (!is.null(fit$bestModel) && inherits(fit$bestModel, "cv.glmnet")) {
    cf <- as.matrix(coef(fit$bestModel, s = "lambda.min"))[, 1]
    return(cf)
  }
  if (!is.null(fit$modlist) && length(fit$modlist) > 0 &&
      inherits(fit$modlist[[1]], "cv.glmnet")) {
    cf <- as.matrix(coef(fit$modlist[[1]], s = "lambda.min"))[, 1]
    return(cf)
  }
  NULL
}

# ---- run all genes in parallel ----
set.seed(123)
results <- future_lapply(seq_along(target_genes), function(i) {
  g <- target_genes[i]
  y <- as.numeric(brain_t[, g])

  # progress
  if (i %% 100 == 0) message("Finished ", i, "/", n_targets, " (", round(100*i/n_targets, 1), "%)")

  # fit (same settings as your single-gene code)
  fit <- tryCatch(
    nestcv.glmnet(
      y = y,
      x = blood_t,
      family = "gaussian",
      outer_method = "cv",
      n_outer_folds = 5,
      n_inner_folds = 5,
      alphaSet = seq(0.1, 1, 0.1),
      cv.cores = 1,
      filterFUN = correl_filter,
      filter_options = list(
        nfilter = min(500, ncol(blood_t)),
        method  = "spearman"
      )
    ),
    error = function(e) NULL
  )

  # if fit failed or no predictions, return NA metrics
  if (is.null(fit) || is.null(fit$output) || nrow(fit$output) == 0) {
    return(list(
      summary = tibble(gene = g, RMSE = NA_real_, R2 = NA_real_,
                       Pearson_r2 = NA_real_, MAE = NA_real_),
      coefs = NULL
    ))
  }

  preds <- fit$output  # has columns: testy, predy
  testy <- as.numeric(preds$testy)
  predy <- as.numeric(preds$predy)

  # metrics to match your single-gene readout
  RMSE <- sqrt(mean((testy - predy)^2))
  RSS  <- sum((testy - predy)^2)
  TSS  <- sum((testy - mean(testy))^2)
  R2   <- 1 - RSS/TSS
  Pearson_r2 <- suppressWarnings(cor(testy, predy, method = "pearson")^2)
  MAE  <- mean(abs(testy - predy))

  # coefficients (if available)
  cf <- tryCatch(.extract_coefs(fit), error = function(e) NULL)
  cf_tbl <- if (!is.null(cf)) {
    tibble(gene = g, feature = names(cf), coef = as.numeric(cf))
  } else NULL

  list(
    summary = tibble(gene = g, RMSE = RMSE, R2 = R2,
                     Pearson_r2 = Pearson_r2, MAE = MAE),
    coefs = cf_tbl
  )
}, future.seed = TRUE)

# ---- collect & report ----
summary_all <- bind_rows(lapply(results, `[[`, "summary"))
coefs_all   <- bind_rows(Filter(Negate(is.null), lapply(results, `[[`, "coefs")))

t1 <- Sys.time()
cat("Total runtime:", round(difftime(t1, t0, units = "mins"), 2), "minutes\n")
cat("Models with NA metrics:", sum(is.na(summary_all$RMSE)), "of", n_targets, "\n")

# peek best genes
summary_all %>% arrange(desc(R2)) %>% slice_head(n = 20)

saveRDS(summary_all, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_try_20250818.rds")
saveRDS(coefs_all, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_coefs_all_20250818.rds")

# ---- save ----
# saveRDS(summary_all, "brain_from_blood_ALL_brain_genes_summary.rds")
# saveRDS(coefs_all,   "brain_from_blood_ALL_brain_genes_coefs.rds")





####################################################
################## ARCHIVED ########################
####################################################

print(ms)
cat("Mean R2:", mean(ms$R2, na.rm = TRUE),
    "  Mean RMSE:", mean(ms$RMSE, na.rm = TRUE), "\n")

# Custom filter function for correlation-based feature selection
# This will be applied within each training fold
correlation_filter <- function(y, x, top_n = 100) {
  # y: target brain gene expression (vector)
  # x: blood gene expression matrix (samples x genes)
  
  # Calculate correlations between target brain gene and all blood genes
  correlations <- cor(y, x, method = "spearman")
  
  # Get absolute correlations and rank them
  abs_cors <- abs(correlations[1, ])
  
  # Select top_n genes with highest absolute correlation
  top_genes <- order(abs_cors, decreasing = TRUE)[1:min(top_n, length(abs_cors))]
  
  return(top_genes)
}

# Function to predict a single brain gene using nested CV
predict_brain_gene <- function(brain_gene_name, brain_data, blood_data, 
                              top_features = 100, outer_folds = 5, inner_folds = 5) {
  
  cat("Predicting brain gene:", brain_gene_name, "\n")
  
  # Extract target brain gene expression
  y <- brain_data[, brain_gene_name]
  
  # Use blood data as predictors
  x <- blood_data
  
  # Set up nested CV with correlation-based feature selection
  set.seed(123)  # For reproducibility
  
  result <- nestcv.glmnet(
    y = y,
    x = x,
    family = "gaussian",  # For continuous gene expression
    alphaSet = seq(0.1, 1, 0.1),  # Test different alpha values for elastic net
    cv.cores = 1,  # Adjust based on your system
    outer_method = "cv",  # Cross-validation for outer loop
    n_outer_folds = outer_folds,
    n_inner_folds = inner_folds,
    filterFUN = correlation_filter,  # Our custom filter function
    filter_options = list(top_n = top_features),  # Parameters for filter
    finalCV = TRUE,  # Perform final cross-validation
    verbose = TRUE
  )
  
  return(result)
}

# Example: Test prediction for one brain gene first
# Let's pick a brain gene that exists in our common genes
brain_genes_available <- intersect(colnames(brain_common), colnames(brain_common))
test_gene <- brain_genes_available[1]  # Pick first available gene
test_gene <- colnames(brain_t)[1:1]

cat("Testing with brain gene:", test_gene, "\n")

# Run nested CV for one gene (test run)
result_test <- predict_brain_gene(
  brain_gene_name = test_gene,
  brain_data = brain_common,
  blood_data = blood_common,
  top_features = 50,  # Start with fewer features for testing
  outer_folds = 3,    # Fewer folds for testing
  inner_folds = 3
)

# Examine results
cat("Test Results Summary:\n")
summary(result_test)

# Plot results
plot(result_test)

# Get performance metrics
cat("R-squared:", result_test$summary$r_squared, "\n")
cat("RMSE:", result_test$summary$rmse, "\n")

# Function to predict multiple brain genes
predict_multiple_genes <- function(brain_gene_list, brain_data, blood_data, 
                                  top_features = 100, outer_folds = 5, inner_folds = 5) {
  
  results_list <- list()
  performance_df <- data.frame(
    gene = character(),
    r_squared = numeric(),
    rmse = numeric(),
    n_features_selected = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(brain_gene_list)) {
    gene_name <- brain_gene_list[i]
    cat("Processing gene", i, "of", length(brain_gene_list), ":", gene_name, "\n")
    
    tryCatch({
      result <- predict_brain_gene(
        brain_gene_name = gene_name,
        brain_data = brain_data,
        blood_data = blood_data,
        top_features = top_features,
        outer_folds = outer_folds,
        inner_folds = inner_folds
      )
      
      results_list[[gene_name]] <- result
      
      # Extract performance metrics
      performance_df <- rbind(performance_df, data.frame(
        gene = gene_name,
        r_squared = result$summary$r_squared,
        rmse = result$summary$rmse,
        n_features_selected = length(result$final_coef) - 1,  # Exclude intercept
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("Error processing gene", gene_name, ":", e$message, "\n")
    })
  }
  
  return(list(results = results_list, performance = performance_df))
}

# Example: Predict a small subset of brain genes
# Select a few genes for initial testing
test_genes <- brain_genes_available[1:3]  # Test with 3 genes first

cat("Running nested CV for", length(test_genes), "brain genes...\n")

# Run prediction for multiple genes
multi_results <- predict_multiple_genes(
  brain_gene_list = test_genes,
  brain_data = brain_common,
  blood_data = blood_common,
  top_features = 50,
  outer_folds = 3,
  inner_folds = 3
)

# View performance results
print(multi_results$performance)

# Plot performance distribution
if (nrow(multi_results$performance) > 1) {
  ggplot(multi_results$performance, aes(x = r_squared)) +
    geom_histogram(bins = 10, fill = "lightblue", color = "black") +
    labs(title = "Distribution of R-squared values",
         x = "R-squared", y = "Count") +
    theme_minimal()
}

# Function to run full analysis (uncomment when ready for full run)
# run_full_analysis <- function() {
#   # Select top variable brain genes (optional)
#   brain_var <- apply(brain_common, 2, var)
#   top_variable_genes <- names(sort(brain_var, decreasing = TRUE))[1:100]
#   
#   cat("Running full analysis for", length(top_variable_genes), "brain genes...\n")
#   
#   full_results <- predict_multiple_genes(
#     brain_gene_list = top_variable_genes,
#     brain_data = brain_common,
#     blood_data = blood_common,
#     top_features = 100,
#     outer_folds = 5,
#     inner_folds = 5
#   )
#   
#   return(full_results)
# }

# Instructions for next steps:
cat("\n=== NEXT STEPS ===\n")
cat("1. Review the test results above\n")
cat("2. If satisfied, uncomment and run run_full_analysis() for all genes\n")
cat("3. Adjust parameters (top_features, folds) based on your computational resources\n")
cat("4. Consider filtering to most variable brain genes first to reduce computation\n")

# Key points about this implementation:
cat("\n=== KEY FEATURES ===\n")
cat("✓ Feature selection done only within training folds (no data leakage)\n")
cat("✓ Proper nested cross-validation structure\n")
cat("✓ Correlation-based feature selection for each brain gene\n")
cat("✓ Elastic net regularization with alpha tuning\n")
cat("✓ Performance metrics tracking\n")



########################################
# brain_full: 21356 genes x 225 samples
# blood_form3: 21046 genes x 225 samples

# First, let's align the genes between brain and blood data
# Find common genes
common_genes <- intersect(rownames(brain_full), rownames(blood_form3))
cat("Number of common genes:", length(common_genes), "\n") #17533 

# Subset to common genes
brain_common <- brain_full[common_genes, ]
blood_common <- blood_form3[common_genes, ]

# Transpose data so samples are rows and genes are columns
# This is the format expected by nestedcv
brain_t <- t(brain_common)  # 225 samples x genes
blood_t <- t(blood_common)  # 225 samples x genes

# ------------------------------------------------------------------------------
# predict_single_brain_gene()
# ------------------------------------------------------------------------------
# Purpose:
#   Predict expression of ONE brain gene from ALL blood genes using elastic net
#   with nested cross-validation. Feature selection (Spearman correlation) is
#   applied INSIDE the outer CV folds to avoid data leakage.
#
# Inputs:
#   gene      : character, the column name of the target brain gene in brain_t
#   brain_t   : matrix/data.frame (n_samples x n_brain_genes)
#   blood_t   : matrix/data.frame (n_samples x n_blood_genes)
#   nfilter   : integer, max number of blood genes to keep per fold (correlation-based)
#   cv_folds  : integer, K for both outer and inner K-fold CV
#   seed      : integer, RNG seed for reproducibility
#
# Returns:
#   A list with:
#     $fit     : nestedcv model object (contains folds, tuning, final fit, preds)
#     $metrics : data.frame with outer-fold metrics (R2, RMSE, etc.)
# ------------------------------------------------------------------------------

predict_single_brain_gene <- function(gene, brain_t, blood_t,
                                      nfilter = 500, cv_folds = 5, seed = 123) {
  set.seed(seed)

  # response as numeric vector
  y <- as.numeric(brain_t[, gene])
  x <- blood_t

  # nested CV with per-fold correlation filter (no leakage)
  fit <- nestcv.glmnet(
    y = y,
    x = x,
    family = "gaussian",
    outer_method = "cv",
    n_outer_folds = cv_folds,
    n_inner_folds = cv_folds,
    alphaSet = seq(0.1, 1, 0.1),
    cv.cores = 1,  # bump later
    filterFUN = correl_filter,
    filter_options = list(nfilter = min(nfilter, ncol(x)), method = "spearman")
  )

  # returns per-outer-fold performance 
  mets <- predSummary(fit)  # columns include RMSE, R2, etc.

  list(fit = fit, metrics = mets)
}

predict_multiple_brain_genes <- function(genes, brain_t, blood_t,
                                         nfilter = 500, cv_folds = 5, seed = 123) {
  res_list <- lapply(genes, function(g) {
    predict_single_brain_gene(g, brain_t, blood_t, nfilter, cv_folds, seed)
  })
  names(res_list) <- genes

  # compact summary (mean over outer folds)
  summary_df <- bind_rows(lapply(seq_along(res_list), function(i) {
    data.frame(
      gene = names(res_list)[i],
      R2_mean = mean(res_list[[i]]$metrics$R2, na.rm = TRUE),
      RMSE_mean = mean(res_list[[i]]$metrics$RMSE, na.rm = TRUE)
    )
  }))

  list(results = res_list, summary = summary_df)
}

# ---------- example run ----------
# pick a few genes to test first
test_genes <- colnames(brain_t)[1:5]

out <- predict_multiple_brain_genes(
  genes   = test_genes,
  brain_t = brain_t,
  blood_t = blood_t,
  nfilter = 500,
  cv_folds = 5,
  seed    = 123
)

print(out$summary)

# access a fitted model & its metrics, e.g. first gene:
# out$results[[1]]$fit
# out$results[[1]]$metrics


# Brain Gene Expression Prediction using Blood Gene Expression
# with Nested Cross-Validation

# Load required libraries
library(nestedcv)
library(glmnet)
library(caret)
library(dplyr)
library(corrr)

# =====================================================================
# DATA PREPARATION
# =====================================================================

# Assuming you have:

# =====================================================================
# FEATURE SELECTION FUNCTION
# =====================================================================

# Function to select top correlated blood genes for each brain gene
select_top_correlated_features <- function(brain_gene_expr, blood_data, top_n = 1000) {
  # Calculate correlations between the target brain gene and all blood genes
  correlations <- cor(brain_gene_expr, blood_data, method = "spearman")
  
  # Get indices of top correlated genes
  top_indices <- order(abs(correlations), decreasing = TRUE)[1:min(top_n, ncol(blood_data))]
  
  # Return selected blood gene expressions
  return(blood_data[, top_indices, drop = FALSE])
}

# =====================================================================
# SINGLE GENE PREDICTION FUNCTION
# =====================================================================

predict_single_brain_gene <- function(brain_gene_name, brain_data, blood_data, 
                                     top_n_features = 1000, cv_folds = 5) {
  
  cat("Predicting brain gene:", brain_gene_name, "\n")
  
  # Extract target brain gene expression
  y <- brain_data[, brain_gene_name]
  
  # Select top correlated blood genes as features
  X_selected <- select_top_correlated_features(y, blood_data, top_n_features)
  
  cat("Using", ncol(X_selected), "blood genes as features\n")
  
  # Perform nested cross-validation with glmnet
  set.seed(123)  # For reproducibility
  
  # Configure nested CV
  nested_result <- nestcv.glmnet(
    y = y,
    x = X_selected,
    family = "gaussian",
    cv.cores = 1,  # Adjust based on your system
    outer_method = "cv",
    n_outer_folds = cv_folds,
    n_inner_folds = cv_folds,
    alphaSet = seq(0.1, 1, 0.1),  # Test different alpha values for elastic net
    filterFUN = NULL,  # We've already done feature selection
    filter_options = NULL
  )
  
  return(nested_result)
}

# =====================================================================
# BATCH PREDICTION FUNCTION
# =====================================================================

predict_multiple_brain_genes <- function(brain_gene_list, brain_data, blood_data,
                                        top_n_features = 1000, cv_folds = 5) {
  
  results <- list()
  
  for (gene in brain_gene_list) {
    if (gene %in% colnames(brain_data)) {
      results[[gene]] <- predict_single_brain_gene(
        gene, brain_data, blood_data, top_n_features, cv_folds
      )
    } else {
      cat("Warning: Gene", gene, "not found in brain data\n")
    }
  }
  
  return(results)
}

# =====================================================================
# RESULTS ANALYSIS FUNCTION
# =====================================================================

analyze_prediction_results <- function(nested_results) {
  
  results_summary <- data.frame(
    gene = character(),
    outer_r2 = numeric(),
    outer_rmse = numeric(),
    best_alpha = numeric(),
    best_lambda = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (gene_name in names(nested_results)) {
    result <- nested_results[[gene_name]]
    
    # Extract performance metrics
    outer_r2 <- mean(result$outer_result$cvm)  # Mean CV R²
    outer_rmse <- sqrt(mean((result$outer_result$fit.preval - result$outer_result$response)^2))
    
    # Get best hyperparameters (from final model)
    best_alpha <- result$final_result$alpha
    best_lambda <- result$final_result$lambda.min
    
    results_summary <- rbind(results_summary, data.frame(
      gene = gene_name,
      outer_r2 = outer_r2,
      outer_rmse = outer_rmse,
      best_alpha = best_alpha,
      best_lambda = best_lambda
    ))
  }
  
  return(results_summary)
}

# =====================================================================
# EXAMPLE USAGE
# =====================================================================

# Test with a small subset of brain genes first (e.g., first 5 genes)
test_genes <- colnames(brain_t)[1:5]

cat("Testing with genes:", paste(test_genes, collapse = ", "), "\n")

# Run nested CV prediction
test_results <- predict_multiple_brain_genes(
  brain_gene_list = test_genes,
  brain_data = brain_t,
  blood_data = blood_t,
  top_n_features = 500,  # Start with 500 features
  cv_folds = 5  # 5-fold CV for faster testing
)

# Analyze results
summary_results <- analyze_prediction_results(test_results)
print(summary_results)

# =====================================================================
# VISUALIZATION FUNCTION
# =====================================================================

plot_prediction_performance <- function(summary_results) {
  library(ggplot2)
  
  # Plot R² values
  p1 <- ggplot(summary_results, aes(x = reorder(gene, outer_r2), y = outer_r2)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Prediction Performance (R²)", 
         x = "Brain Gene", y = "Cross-Validation R²") +
    theme_minimal()
  
  # Plot best alpha values
  p2 <- ggplot(summary_results, aes(x = reorder(gene, best_alpha), y = best_alpha)) +
    geom_col(fill = "coral", alpha = 0.7) +
    coord_flip() +
    labs(title = "Optimal Alpha Values", 
         x = "Brain Gene", y = "Best Alpha (Elastic Net)") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

# Plot results
if (nrow(summary_results) > 0) {
  plot_prediction_performance(summary_results)
}

# =====================================================================
# FULL ANALYSIS (RUN AFTER TESTING)
# =====================================================================

# After testing, you can run the full analysis:
# all_brain_genes <- colnames(brain_t)[1:100]  # Start with first 100 genes
# 
# full_results <- predict_multiple_brain_genes(
#   brain_gene_list = all_brain_genes,
#   brain_data = brain_t,
#   blood_data = blood_t,
#   top_n_features = 1000,
#   cv_folds = 10
# )
# 
# full_summary <- analyze_prediction_results(full_results)
# write.csv(full_summary, "brain_gene_prediction_results.csv", row.names = FALSE)

# =====================================================================
# ADDITIONAL UTILITIES
# =====================================================================

# Function to get feature importance for a specific prediction
get_feature_importance <- function(nested_result, gene_name) {
  final_model <- nested_result$final_result
  coefficients <- coef(final_model, s = "lambda.min")
  non_zero_coefs <- coefficients[coefficients[,1] != 0, , drop = FALSE]
  
  importance_df <- data.frame(
    feature = rownames(non_zero_coefs),
    coefficient = non_zero_coefs[,1],
    abs_coefficient = abs(non_zero_coefs[,1])
  )
  
  importance_df <- importance_df[order(importance_df$abs_coefficient, decreasing = TRUE), ]
  
  return(importance_df)
}

# Example usage for feature importance:
# if (length(test_results) > 0) {
#   importance <- get_feature_importance(test_results[[1]], names(test_results)[1])
#   print(head(importance, 10))
# }