#code derived from form3_with_residID_full.R
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

## Job submussions
#path=/sc/arion/projects/mscic1/results/jolie/LBP/scripts/
#cd $path
#ml R

#bsub -q premium -P acc_mscic1 -n 30 -W 144:00 -R rusage[mem=4000] -R span[hosts=1] -o %J.stdout -eo %J.stderr Rscript ${path}	null_swap.R
#bjobs

predict_brain_gene_expression_single_fold <- function(X_all, Y_all, brain_gene, 
                                                     train_idx, test_idx, 
                                                     top_percent = 0.05,
                                                     alpha_values = c(0, 0.25, 0.5, 0.75, 1),
                                                     n_inner_folds = 5) {
  library(glmnet)
  library(caret)
  
  # Skip if brain gene is missing from Y matrix
  if (!(brain_gene %in% colnames(Y_all))) return(NULL)
  
  # Get the brain expression vector
  y <- Y_all[, brain_gene]
  y_train <- y[train_idx]
  y_test <- y[test_idx]
  
  # Feature selection using only training data
  X_train_full <- X_all[train_idx, , drop = FALSE]
  
  # Calculate correlations on training data only
  cor_vec <- apply(X_train_full, 2, function(x) {
    cor(x, y_train, method = "spearman", use = "pairwise.complete.obs")
  })
  
  # Remove any NAs
  cor_vec <- cor_vec[!is.na(cor_vec)]
  
  # Select top predictors based only on training data correlations
  top_n <- ceiling(length(cor_vec) * top_percent)
  top_blood_genes <- names(sort(abs(cor_vec), decreasing = TRUE))[1:top_n]
  
  # Sanity checks
  if (length(top_blood_genes) < 2) return(NULL)
  if (!all(top_blood_genes %in% colnames(X_all))) return(NULL)
  
  # Subset data to selected features
  X_train <- X_all[train_idx, top_blood_genes, drop = FALSE]
  X_test <- X_all[test_idx, top_blood_genes, drop = FALSE]
  
  # Inner CV for hyperparameter tuning (alpha selection)
  # Create inner CV folds
  inner_folds <- createFolds(y_train, k = n_inner_folds, list = TRUE)
  
  best_alpha <- NULL
  best_inner_score <- -Inf
  
  for (alpha in alpha_values) {
    inner_scores <- numeric(n_inner_folds)
    
    for (inner_fold in 1:n_inner_folds) {
      inner_test_idx <- inner_folds[[inner_fold]]
      inner_train_idx <- setdiff(1:length(y_train), inner_test_idx)
      
      X_inner_train <- X_train[inner_train_idx, , drop = FALSE]
      y_inner_train <- y_train[inner_train_idx]
      X_inner_test <- X_train[inner_test_idx, , drop = FALSE]
      y_inner_test <- y_train[inner_test_idx]
      
      # Fit model with current alpha
      inner_model <- cv.glmnet(X_inner_train, y_inner_train, alpha = alpha, 
                              nfolds = 5, type.measure = "mse")
      
      # Predict on inner test set
      inner_pred <- predict(inner_model, newx = X_inner_test, s = "lambda.min")
      
      # Calculate R² for this inner fold
      inner_scores[inner_fold] <- 1 - sum((y_inner_test - inner_pred)^2) / 
                                     sum((y_inner_test - mean(y_inner_test))^2)
    }
    
    mean_inner_score <- mean(inner_scores, na.rm = TRUE)
    
    if (mean_inner_score > best_inner_score) {
      best_inner_score <- mean_inner_score
      best_alpha <- alpha
    }
  }
  
  # Train final model with best alpha on full training data
  final_model <- cv.glmnet(X_train, y_train, alpha = best_alpha, 
                          nfolds = n_inner_folds, type.measure = "mse")
  
  # Predict on training and test sets
  y_pred_train <- predict(final_model, newx = X_train, s = "lambda.min")
  y_pred_test <- predict(final_model, newx = X_test, s = "lambda.min")
  
  # Compute metrics
  r2_train <- 1 - sum((y_train - y_pred_train)^2) / sum((y_train - mean(y_train))^2)
  r2_test <- 1 - sum((y_test - y_pred_test)^2) / sum((y_test - mean(y_test))^2)
  cor_test_spearman <- cor(as.vector(y_pred_test), y_test, method = "spearman", use = "complete.obs")
  cor_test_pearson <- cor(as.vector(y_pred_test), y_test, method = "pearson", use = "complete.obs")
  
  # Return results
  return(list(
    brain_gene = brain_gene,
    n_predictors = length(top_blood_genes),
    r2_train = r2_train,
    r2_test = r2_test,
    cor_test_spearman = cor_test_spearman,
    cor_test_pearson = cor_test_pearson,
    best_alpha = best_alpha,
    lambda_min = final_model$lambda.min,
    best_inner_cv_score = best_inner_score
  ))
}

# Main nested CV function
nested_cv_brain_prediction <- function(X_all, Y_all, brain_genes, 
                                      n_outer_folds = 5,
                                      n_inner_folds = 5,
                                      top_percent = 0.05,
                                      alpha_values = c(0, 0.25, 0.5, 0.75, 1),
                                      seed = 123,
                                      parallel = FALSE,
                                      n_cores = NULL) {
  
  library(glmnet)
  library(caret)
  library(dplyr)
  
  if (parallel) {
    library(parallel)
    library(foreach)
    library(doParallel)
  }
  
  set.seed(seed)
  
  # Create outer CV folds
  n_samples <- nrow(X_all)
  sample_indices <- 1:n_samples
  outer_folds <- createFolds(sample_indices, k = n_outer_folds, list = TRUE)
  
  cat("Starting nested CV with", n_outer_folds, "outer folds and", n_inner_folds, "inner folds\n")
  cat("Testing", length(brain_genes), "brain genes\n")
  cat("Alpha values to test:", paste(alpha_values, collapse = ", "), "\n\n")
  
  all_results <- list()
  
  # Set up parallel processing if requested
  if (parallel) {
    if (is.null(n_cores)) n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  for (outer_fold in 1:n_outer_folds) {
    cat("Processing outer fold", outer_fold, "of", n_outer_folds, "\n")
    
    test_idx <- outer_folds[[outer_fold]]
    train_idx <- setdiff(sample_indices, test_idx)
    
    if (parallel) {
      # Parallel processing across brain genes
      fold_results <- foreach(brain_gene = brain_genes, 
                             .packages = c("glmnet", "caret"),
                             .combine = 'c') %dopar% {
        result <- predict_brain_gene_expression_single_fold(
          X_all, Y_all, brain_gene, train_idx, test_idx, 
          top_percent, alpha_values, n_inner_folds
        )
        list(result)
      }
      # Remove NULL results
      fold_results <- fold_results[!sapply(fold_results, is.null)]
    } else {
      # Sequential processing
      fold_results <- list()
      for (i in seq_along(brain_genes)) {
        if (i %% 100 == 0) cat("  Brain gene", i, "of", length(brain_genes), "\n")
        
        result <- predict_brain_gene_expression_single_fold(
          X_all, Y_all, brain_genes[i], train_idx, test_idx, 
          top_percent, alpha_values, n_inner_folds
        )
        
        if (!is.null(result)) {
          fold_results[[length(fold_results) + 1]] <- result
        }
      }
    }
    
    # Add fold information to results
    for (i in seq_along(fold_results)) {
      fold_results[[i]]$outer_fold <- outer_fold
    }
    
    all_results <- c(all_results, fold_results)
    cat("Completed outer fold", outer_fold, "- processed", length(fold_results), "brain genes\n\n")
  }
  
  # Convert results to dataframe
  if (length(all_results) > 0) {
    results_df <- do.call(rbind, lapply(all_results, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }))
    rownames(results_df) <- NULL
    
    # Calculate summary statistics across folds for each brain gene
    summary_stats <- results_df %>%
      group_by(brain_gene) %>%
      summarise(
        n_folds = n(),
        mean_r2_test = mean(r2_test, na.rm = TRUE),
        sd_r2_test = sd(r2_test, na.rm = TRUE),
        mean_cor_spearman = mean(cor_test_spearman, na.rm = TRUE),
        sd_cor_spearman = sd(cor_test_spearman, na.rm = TRUE),
        mean_cor_pearson = mean(cor_test_pearson, na.rm = TRUE),
        sd_cor_pearson = sd(cor_test_pearson, na.rm = TRUE),
        mean_n_predictors = mean(n_predictors, na.rm = TRUE),
        most_common_alpha = names(sort(table(best_alpha), decreasing = TRUE))[1],
        .groups = 'drop'
      ) %>%
      arrange(desc(mean_r2_test))
    
    cat("Nested CV completed successfully!\n")
    cat("Processed", nrow(summary_stats), "brain genes across", n_outer_folds, "outer folds\n")
    
    return(list(
      detailed_results = results_df,
      summary_stats = summary_stats,
      parameters = list(
        n_outer_folds = n_outer_folds,
        n_inner_folds = n_inner_folds,
        top_percent = top_percent,
        alpha_values = alpha_values,
        seed = seed
      )
    ))
  } else {
    warning("No valid results obtained!")
    return(NULL)
  }
}

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_expr <- blood_form3
brain_expr <-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_expr)<- brain_expr$V1; brain_expr$V1 <- NULL

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

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_with_indivdualID_nestedCV_multiple_job.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
