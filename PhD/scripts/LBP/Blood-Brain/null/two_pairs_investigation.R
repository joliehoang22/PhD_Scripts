## this code explores the discrepancy of two-pairs in the dictionary RData vs the original RData
rm(list=ls())
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
### PROBLEM
##code derived from examining_swap_results.R
#listing rows that are NAs
two_pairs[is.na(two_pairs$LR_brain_extraction_diff_days), ]
#IID_ISMMS SAMPLE_ISMMS_blood SAMPLE_ISMMS_brain number_of_pair_brain
#9     PT-0024   LBPSEMA4BLOOD728   LBPSEMA4BRAIN532              2_pairs
#14    PT-0028   LBPSEMA4BLOOD067   LBPSEMA4BRAIN783              2_pairs
#15    PT-0029   LBPSEMA4BLOOD657   LBPSEMA4BRAIN765              2_pairs
#100   PT-0113   LBPSEMA4BLOOD296   LBPSEMA4BRAIN713              2_pairs
#103   PT-0123   LBPSEMA4BLOOD307   LBPSEMA4BRAIN747              2_pairs
#106   PT-0126   LBPSEMA4BLOOD356   LBPSEMA4BRAIN460              2_pairs
#133   PT-0162   LBPSEMA4BLOOD281   LBPSEMA4BRAIN124              2_pairs
#    mymet_tissue_brain mymet_tissue_blood mymet_extractiondate_brain
#9              R_Brain            R_Blood                 2018-08-20
#14             R_Brain            R_Blood                 2018-08-22
#15             R_Brain            R_Blood                 2018-08-22
#100            R_Brain            R_Blood                 2018-09-19
#103            R_Brain            R_Blood                 2018-08-22
#106            L_Brain            L_Blood                 2018-06-08
#133            R_Brain            R_Blood                 2018-11-08
#    LR_brain_extraction_diff_days
#9                              NA
#14                             NA
#15                             NA
#100                            NA
#103                            NA
#106                            NA
#133                            NA

### SOLUTION
##go to the root 
metadata <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")          
dim(metadata) ##[1] 530 170  
sub<-filter(metadata, IID_ISMMS == "PT-0024"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0028"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0029"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0113"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0123"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0126"); sub #2-pairs XX
sub<-filter(metadata, IID_ISMMS == "PT-0162"); sub #2-pairs XX

blood_samples_to_be_removed <- c("LBPSEMA4BLOOD795", "LBPSEMA4BLOOD142", "LBPSEMA4BLOOD049", "LBPSEMA4BLOOD755",
                                  "LBPSEMA4BLOOD335","LBPSEMA4BLOOD431","LBPSEMA4BLOOD174","LBPSEMA4BLOOD738")

brain_samples_to_be_removed <- c("LBPSEMA4BRAIN364", "LBPSEMA4BRAIN318", "LBPSEMA4BRAIN564", "LBPSEMA4BRAIN170",
                                  "LBPSEMA4BRAIN703", "LBPSEMA4BRAIN341","LBPSEMA4BRAIN391","LBPSEMA4BRAIN017")

sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD795"); sub$IID_ISMMS #PT-0113 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD142"); sub$IID_ISMMS #PT-0028 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD049"); sub$IID_ISMMS #PT-0126 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD755"); sub$IID_ISMMS #PT-0162 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD335"); sub$IID_ISMMS #PT-0179 this is originally two_brain_one_blood_sample --> 1 pair
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD431"); sub$IID_ISMMS #PT-0024 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD174"); sub$IID_ISMMS #PT-0029 XX
sub<-filter(metadata, SAMPLE_ISMMS == "LBPSEMA4BLOOD738"); sub$IID_ISMMS #PT-0123 XX

## SO IT'S SAFE TO CHANGE THOSE IDs ABOVE FROM 2-PAIRS TO 1-PAIR 
## THEY WERE ORIGINALLY 2 PAIRS BUT 1 PAIR GOT DROPPED BECAUSE OF DATA QUALITY CONTROL, SO NOW THEY ARE 1-PAIR

#### OLD R DATA
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_voom_20250522.RData")
dim(brain_metadata)#233 339
test <- filter(brain_metadata, IID_ISMMS == "PT-0024") #1 pair
test$SAMPLE_ISMMS_blood #LBPSEMA4BLOOD326
test$SAMPLE_ISMMS_brain #LBPSEMA4BRAIN034

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")
head(dictionary)
sub <- filter(dictionary, SAMPLE_ISMMS_blood == "LBPSEMA4BLOOD738"); sub #0, this means that the samples above were already dropped in the creation of RData
table(brain_metadata$number_of_pair_brain) 
# 1_pair 2_pairs 
#     72     153

brain_metadata$number_of_pair_brain[
  brain_metadata$IID_ISMMS %in% c("PT-0024", "PT-0028", "PT-0029", 
                                  "PT-0113", "PT-0123", "PT-0126", "PT-0162")] <- "1_pair"

table(brain_metadata$number_of_pair_brain)
# 1_pair 2_pairs 
#     79     146
##EXPECTED  -- good this is correct

#1 pair = 79 samples / people
#2 pairs = 146 samples / 2 = 73 people
# total people = 79 + 73 = 152 people

# BEFORE
sub <- filter(dictionary, IID_ISMMS == "PT-0024"); sub

dictionary$number_of_pair_brain[
  dictionary$IID_ISMMS %in% c("PT-0024", "PT-0028", "PT-0029", 
                                  "PT-0113", "PT-0123", "PT-0126", "PT-0162")] <- "1_pair"

##SANITY CHECK - AFTER
sub <- filter(dictionary, IID_ISMMS == "PT-0024"); sub

#save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")
##SANITY CHECK - AFTER
sub <- filter(dictionary, IID_ISMMS == "PT-0024"); sub
table(brain_metadata$number_of_pair_brain) 
# 1_pair 2_pairs 
#     79     146

## RERUN THE WITHIN-INDIVIDUALS SHUFFLING TO SEE IF THE RESULTS HOLD 
# code derived from null_swap.R

## THERE ARE 2 BLOOD FORMS !! REMEMBER TO CHANGE
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE) ##CHANGE 1
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_form = blood_form3 

# Swap levels in percentages
swap_levels <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
n_perm <- 10000 #this was 1000, changed to 10000

# Summary stats to compute
summary_stats <- c(
  "min", "p25", "mean", "median", "p80", "p90", "p92.5", "p95", "p97.5", "p98.5",
  "p99.5", "p99.95", "p99.995", "p99.9995", "p99.99995", "p99.999995",
  "p99.9999995", "p99.99999995", "max"
)
n_stats <- length(summary_stats) #19 summary stats

two_pair_inds <- unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])
n_two_pair <- length(two_pair_inds); n_two_pair #73 people

# Set up parallel processing with 30 workers
plan(multisession, workers = 30)

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

save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_with_indivdualID_null_results_10000perm.RData")

save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_no_indivdualID_null_results_10000perm.RData")

### OLD RESULTS WHERE THE NUMBER OF 2 PAIRS WERE WRONG !
#save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results_swap100.RData")
#save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")
#save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")
