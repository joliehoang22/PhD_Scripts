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
path=/sc/arion/projects/mscic1/results/jolie/LBP/scripts/
cd $path
ml R
bsub -q premium -P acc_mscic1 -n 20 -W 144:00 -R rusage[mem=4000] -R span[hosts=1] -R himem -o %J.stdout -eo %J.stderr Rscript ${path}nestedCV_job_form3_with_residID_base_brain.R

#Job <199223693> is submitted to queue <premium>.

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
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
saveRDS(results, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_with_residID_results_20250818.rds")
saveRDS(failed_genes, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_with_residID_failed_genes_20250818.rds")

# Extract performance metrics for all genes
performance_df <- do.call(rbind, lapply(names(results), function(gene) {
  data.frame(
    gene = gene,
    cor = results[[gene]]$cv_performance["cor"],
    rmse = results[[gene]]$cv_performance["rmse"],
    mae = results[[gene]]$cv_performance["mae"]
  )
}))

write.csv(performance_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/nestedcv_brain_base_blood_with_residID_prediction_performance.csv", row.names = FALSE)

cat("Results saved to:\n")
cat("- brain_gene_prediction_results.rds (full results)\n")
cat("- gene_prediction_performance.csv (performance metrics)\n")
cat("- failed_genes.rds (genes that failed)\n")
