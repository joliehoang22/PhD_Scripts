## This code builds an elastic net 
## This code is stored in source /sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/scripts/

## main application are found in form3_no_resiID_base.R, etc.

predict_brain_gene_expression <- function(X_all, Y_all, top_gene_pairs, brain_genes, train_idx, test_idx) {
  library(glmnet)
  library(dplyr)

  results <- list()

  for (brain_gene in brain_genes) {
    
    # Get top-correlated blood genes for this brain gene
    top_blood_genes <- top_gene_pairs %>%
      filter(brain_gene == !!brain_gene) %>%
      pull(blood_gene)
    
    # Skip if fewer than 2 predictors (not enough for modeling)
    if (length(top_blood_genes) < 2) next
    
    # Check if gene exists in expression matrices
    if (!(brain_gene %in% colnames(Y_all))) next
    if (!all(top_blood_genes %in% colnames(X_all))) next
    
    # Subset X and y
    X <- X_all[, top_blood_genes, drop = FALSE]
    y <- Y_all[, brain_gene]
    
    # Training data
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    
    # Test data
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    
    # Train Elastic Net
    model <- cv.glmnet(X_train, y_train, alpha = 0.5)
    
    # Predict on training set
    y_pred_train <- predict(model, newx = X_train, s = "lambda.min")
    r2_train <- 1 - sum((y_train - y_pred_train)^2) / sum((y_train - mean(y_train))^2)
    
    # Predict on test set
    y_pred_test <- predict(model, newx = X_test, s = "lambda.min")
    cor_test <- cor(y_pred_test, y_test)
    r2_test <- 1 - sum((y_test - y_pred_test)^2) / sum((y_test - mean(y_test))^2)
    
    # Save results
    results[[brain_gene]] <- list(
      brain_gene = brain_gene,
      n_predictors = length(top_blood_genes),
      r2_train = r2_train,
      cor_test = cor_test,
      r2_test = r2_test,
      lambda_min = model$lambda.min
    )
  }

  # Return results as a dataframe
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  rownames(results_df) <- NULL
  results_df <- results_df[order(-results_df$r2_test), ] #order by highest to lowest r2_test
  return(results_df)
}
