## code derived from LBP_blood_brain_QC_20250729.R

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
library(furrr)      # for future_map
library(future)     # for plan()
library(purrr)      # for map-style helpers

set.seed(2025)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")
table(brain_metadata$number_of_pair_brain) ##
# 1_pair 2_pairs 
#     79     146

#1 pair = 79 samples / people
#2 pairs = 146 samples / 2 = 73 people

#swap: 0%, 10%, 25%, 50%, 75%, 100%
#For each swap, compute the blood-brain correlation, and 
#repeat the process 1000 times per swap level to get a distribution of expected correlations.

#To simulate what happens when you increasingly disrupt the true pairing — 
#if correlation drops as swaps increase, that supports the importance of specific sample-level pairing (not just person-level identity)

#null summary code: blood_brain_null_form0.r

#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE) ##CHANGE 1
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL

blood_form = blood_form3 ##CHANGE 1

# Swap levels in percentages
swap_levels <- c(0, 0.1, 0.25, 0.5, 0.75, 1.0)
swap_levels <- 1
n_perm <- 1000

# Summary stats to compute
summary_stats <- c(
  "min", "p25", "mean", "median", "p80", "p90", "p92.5", "p95", "p97.5", "p98.5",
  "p99.5", "p99.95", "p99.995", "p99.9995", "p99.99995", "p99.999995",
  "p99.9999995", "p99.99999995", "max"
)
n_stats <- length(summary_stats) #19 summary stats

two_pair_inds <- unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])
n_two_pair <- length(two_pair_inds) #73 people

# Preallocate list to store summaries
null_results <- list()

for (swap_frac in swap_levels) {
  message("Starting swap level: ", swap_frac * 100, "%")

  swap_summaries <- matrix(NA, nrow = n_perm, ncol = n_stats)
  colnames(swap_summaries) <- summary_stats

  for (i in 1:n_perm) {
    if (i %% 50 == 0) message("  Permutation ", i, " of ", n_perm)

    # Select fraction to swap
    n_swap <- floor(swap_frac * n_two_pair)
    swap_inds <- sample(two_pair_inds, n_swap)

    # Make a copy of brain expression matrix
    v_brain_perm <- v_brain$E

    # Swap the brain samples within selected individuals
    for (id in swap_inds) {
      brain_cols <- which(brain_metadata$IID_ISMMS == id)
      if (length(brain_cols) == 2) {
        v_brain_perm[, brain_cols] <- v_brain_perm[, rev(brain_cols)]
      }
    }

    # Compute Spearman correlation with blood
    null_cor <- cor(t(blood_form), t(v_brain_perm), method = "spearman")
    abs_cor <- abs(null_cor)

    # Compute all summary stats
    swap_summaries[i, ] <- c(
      min = min(abs_cor),
      p25 = quantile(abs_cor, 0.25),
      mean = mean(abs_cor),
      median = median(abs_cor),
      p80 = quantile(abs_cor, 0.80),
      p90 = quantile(abs_cor, 0.90),
      p92.5 = quantile(abs_cor, 0.925),
      p95 = quantile(abs_cor, 0.95),
      p97.5 = quantile(abs_cor, 0.975),
      p98.5 = quantile(abs_cor, 0.985),
      p99.5 = quantile(abs_cor, 0.995),
      p99.95 = quantile(abs_cor, 0.9995),
      p99.995 = quantile(abs_cor, 0.99995),
      p99.9995 = quantile(abs_cor, 0.999995),
      p99.99995 = quantile(abs_cor, 0.9999995),
      p99.999995 = quantile(abs_cor, 0.99999995),
      p99.9999995 = quantile(abs_cor, 0.999999995),
      p99.99999995 = quantile(abs_cor, 0.9999999995),
      max = max(abs_cor)
    )
  }

  # Save to list
  null_results[[paste0("swap_", swap_frac * 100, "pct")]] <- as.data.frame(swap_summaries)
}

#############################
###### PARALLELIZATION ######
#############################
# Set up parallel processing with 30 workers
plan(multisession, workers = 20)

# List to store results for each swap level
null_results <- list()

# Loop through each swap level (0%, 10%, 25%, 50%, 75%, 100%)
for (swap_frac in swap_levels) {
  swap_label <- paste0("swap_", swap_frac * 100, "pct")
  message("\n========== Starting ", swap_label, " ==========")

  # Optional: Track total time for each swap level
  start_time <- Sys.time()

  # Run 1000 permutations in parallel
  swap_summaries <- future_map(
    1:n_perm,
    function(i) {
      # Print progress every 50 permutations (safe in parallel)
      if (i %% 50 == 0) cat("[", swap_label, "] Completed permutation", i, "of", n_perm, "\n")

      # Randomly select which 2-pair individuals to swap
      n_swap <- floor(swap_frac * n_two_pair)
      swap_inds <- sample(two_pair_inds, n_swap)

      # Copy brain matrix and apply within-individual swaps
      v_brain_perm <- v_brain$E
      for (id in swap_inds) {
        brain_cols <- which(brain_metadata$IID_ISMMS == id)
        if (length(brain_cols) == 2) {
          v_brain_perm[, brain_cols] <- v_brain_perm[, rev(brain_cols)]
        }
      }

      # Compute Spearman correlation between blood and permuted brain data
      null_cor <- cor(t(blood_form), t(v_brain_perm), method = "spearman")
      abs_cor <- abs(null_cor)

      # Return vector of summary statistics
      c(
        min = min(abs_cor),
        p25 = quantile(abs_cor, 0.25),
        mean = mean(abs_cor),
        median = median(abs_cor),
        p80 = quantile(abs_cor, 0.80),
        p90 = quantile(abs_cor, 0.90),
        p92.5 = quantile(abs_cor, 0.925),
        p95 = quantile(abs_cor, 0.95),
        p97.5 = quantile(abs_cor, 0.975),
        p98.5 = quantile(abs_cor, 0.985),
        p99.5 = quantile(abs_cor, 0.995),
        p99.95 = quantile(abs_cor, 0.9995),
        p99.995 = quantile(abs_cor, 0.99995),
        p99.9995 = quantile(abs_cor, 0.999995),
        p99.99995 = quantile(abs_cor, 0.9999995),
        p99.999995 = quantile(abs_cor, 0.99999995),
        p99.9999995 = quantile(abs_cor, 0.999999995),
        p99.99999995 = quantile(abs_cor, 0.9999999995),
        max = max(abs_cor)
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  # Store results for this swap level
  null_results[[swap_label]] <- as.data.frame(do.call(rbind, swap_summaries))

  # Report time taken
  end_time <- Sys.time()
  message("Finished ", swap_label, " | Elapsed time: ", round(difftime(end_time, start_time, units = "mins"), 2), " mins")
}
save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results_swap100.RData")
#this is NO residID
save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")

#this is with residID
save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")

#start: 4:02pm 7/29/25

#start2: 5:38pm 7/29/25  - running the null shuffling with residID 
#Each swaps took 164.76 mins (2.7 hours) and i have 6 swaps --> 16.5 hours --> end at 10:08am the next day

#start3: 7:11pm 7/30/25 - running the null shuffling with NO residID 

