## improved version of elastic_net_function_data_leakage.R
##Feature selection on training data only | CV with single train/test split approach

predict_brain_gene_expression <- function(X_all, Y_all, brain_genes, train_idx, test_idx, top_percent = 0.05) {
  library(glmnet)
  library(dplyr)

  results <- list()

  for (brain_gene in brain_genes) {
    
    # Skip if brain gene is missing from Y matrix
    if (!(brain_gene %in% colnames(Y_all))) next
    
    # Get the brain expression vector
    y <- Y_all[, brain_gene]

    # === [CHANGE: use only training data for feature selection] ===
    y_train <- y[train_idx]
    X_train_full <- X_all[train_idx, , drop = FALSE]

    # === [CHANGE: use Spearman correlation for robustness] ===
    cor_vec <- apply(X_train_full, 2, function(x) cor(x, y_train, method = "spearman", use = "pairwise.complete.obs"))

    # === [CHANGE: select top predictors based only on training data correlations] ===
    top_n <- ceiling(length(cor_vec) * top_percent)
    top_blood_genes <- names(sort(abs(cor_vec), decreasing = TRUE))[1:top_n]

    # Sanity checks
    if (length(top_blood_genes) < 2) next
    if (!all(top_blood_genes %in% colnames(X_all))) next

    # Subset both training and test data to selected features
    X <- X_all[, top_blood_genes, drop = FALSE]
    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]

    # === [NO CHANGE: fit model only on training data] ===
    model <- cv.glmnet(X_train, y_train, alpha = 0.5)

    # Predict on training and test
    y_pred_train <- predict(model, newx = X_train, s = "lambda.min")
    y_pred_test <- predict(model, newx = X_test, s = "lambda.min")

    # Compute R²
    r2_train <- 1 - sum((y_train - y_pred_train)^2) / sum((y_train - mean(y_train))^2)
    r2_test <- 1 - sum((y_test - y_pred_test)^2) / sum((y_test - mean(y_test))^2)

    # === [OPTIONAL: report Spearman correlation for test too] ===
    cor_test <- cor(y_pred_test, y_test, method = "spearman")

    # Save results
    results[[brain_gene]] <- list(
      brain_gene = brain_gene,
      n_predictors = length(top_blood_genes),
      r2_train = r2_train,
      r2_test = r2_test,
      cor_test = cor_test,
      lambda_min = model$lambda.min
    )
  }

  # Return as dataframe
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  rownames(results_df) <- NULL
  results_df <- results_df[order(-results_df$r2_test), ]
  return(results_df)
}
