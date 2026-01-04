---
title: "Untitled"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## GitHub Documents

This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.

## Question of the Day: Does Pairing within the Same Person Matter? 7-29-2025
Shuffling sample lables within individuals to see whether pairing of blood-brain matters.
Creating null summaries with varying degree of swaps 

```{r cars}
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
table(brain_metadata$number_of_pair_brain)
# 1_pair 2_pairs 
#     72     153

#153 × 2 = 306 samples are candidates for shuffling
#swap: 0%, 10%, 25%, 50%, 75%, 100%
#For each swap, compute the blood-brain correlation, and 
#repeat the process 1000 times per swap level to get a distribution of expected correlations.

#To simulate what happens when you increasingly disrupt the true pairing — 
#if correlation drops as swaps increase, that supports the importance of specific sample-level pairing (not just person-level identity).

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

# Assume you already have these sample sets prepared:
# data.frame with columns: IID_ISMMS, sample_id_blood, sample_id_brain
two_pair_inds <- unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])
n_two_pair <- length(two_pair_inds) #80
#length(unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])) #80

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

########this is what results will look like 
names(null_results)
# [1] "swap_0pct"  "swap_10pct"  "swap_25pct"  "swap_50pct"  "swap_75pct"  "swap_100pct"

dim(null_results[["swap_0pct"]])
# [1] 1000   19

head(null_results[["swap_0pct"]])
#         min      p25     mean   median      p80 ...
# 1   0.12032  0.25143  0.31244  0.30321  0.38765 ...
# 2   0.11840  0.24873  0.31123  0.30400  0.38621 ...
# ...

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
  message("✅ Finished ", swap_label, " | Elapsed time: ", round(difftime(end_time, start_time, units = "mins"), 2), " mins")
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

```

## Examining the null results: 07-30-2025 

You can also embed plots, for example:

```{r pressure, echo=FALSE}
library(data.table)
library(ggplot2)
library(readxl)
library(biomaRt)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(matrixStats)
library(tibble)
library(purrr)
library(caret)
library(tidyr)
library(furrr)      # for future_map
library(future)     # for plan()
library(purrr)      # for map-style helpers
set.seed(2025)
##first thing first just check the dictionary to make sure my pairing is correct
#did this in LBP_blood_brain_QC_20250603.rmd and LBP_clean_blood_brain.ipynb
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")
## swapping 0%, 10% (8 pairs), 25% (20 pairs), 50% (40 pairs), 75% (60 pairs), and 100% (80 pairs)
names(null_results)
#[1] "swap_0pct"   "swap_10pct"  "swap_25pct"  "swap_50pct"  "swap_75pct" 
#[6] "swap_100pct"

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")

dim(null_results[["swap_0pct"]])
# [1] 1000   19

swap_0pct <- null_results[["swap_0pct"]]
swap_10pct <- null_results[["swap_10pct"]]
swap_25pct <- null_results[["swap_25pct"]]
swap_50pct <- null_results[["swap_50pct"]]
swap_75pct <- null_results[["swap_75pct"]]
swap_100pct <- null_results[["swap_100pct"]]

#summary of each column of each swap 
summary_stats_0pct <- sapply(head(swap_0pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t0 <- t(summary_stats_0pct)  #transpose to make it more readable

summary_stats_10pct <- sapply(head(swap_10pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t10 <- t(summary_stats_10pct)

summary_stats_25pct <- sapply(head(swap_25pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t25 <- t(summary_stats_25pct)

summary_stats_50pct <- sapply(head(swap_50pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t50 <- t(summary_stats_50pct)

summary_stats_75pct <- sapply(head(swap_75pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t75 <- t(summary_stats_75pct)

summary_stats_100pct <- sapply(head(swap_100pct), function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
})
t100 <- t(summary_stats_100pct)

df <- bind_rows(
  as.data.frame(t0) %>% rownames_to_column("stat") %>% mutate(swap_pct = 0),
  as.data.frame(t10) %>% rownames_to_column("stat") %>% mutate(swap_pct = 10),
  as.data.frame(t25) %>% rownames_to_column("stat") %>% mutate(swap_pct = 25),
  as.data.frame(t50) %>% rownames_to_column("stat") %>% mutate(swap_pct = 50),
  as.data.frame(t75) %>% rownames_to_column("stat") %>% mutate(swap_pct = 75),
  as.data.frame(t100) %>% rownames_to_column("stat") %>% mutate(swap_pct = 100)
)

df_long <- df %>%
  pivot_longer(cols = c("mean", "median", "min", "max"), names_to = "measure", values_to = "value")

plot<- ggplot(df_long, aes(x = swap_pct, y = value, color = stat, group = stat)) +
  geom_line() +
  facet_wrap(~ measure) +
  labs(x = "Swap Percentage", y = "Value", color = "Summary Statistic") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/shuffling_summary_no_residID.pdf", plot = plot, width = 16, height = 8) 

#ggsave("/hpc/users/hoangd02/www/plots/lbp/shuffling_summary_with_residID.pdf", plot = plot, width = 16, height = 8) 

```

## 7-31-2025: Checking RedCap to remove people who did not concent for research + looking at some imaging
Shuffling sample lables within individuals to see whether pairing of blood-brain matters.
Creating null summaries with varying degree of swaps 

```{r cars}
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
set.seed(2025)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")
head(dictionary)

table(dictionary$IID_ISMMS)

no_go <- c("PT-0003","PT-0007","PT-0008","PT-0014","PT-0050","PT-0066",
            "PT-0067", "PT-0075","PT-0077","PT-0107","PT-0129","PT-0137",
            "PT-0151","PT-0168","PT-0176","PT-0216","PT-0259","PT-0417","PT-0428")

#good to go! no one has that ID 

##A bit of imaging -- do NOT change anything in the path! copy and paste brains in your folder if needed
#Path for Raw nifty
/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output

#Path for Freesurfer outputs
/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output_fs2

imaging_map <- fread("/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/mapping_files/scan_id_date_path_map.csv")

length(unique(imaging_map$subject_id)) #179

subset(imaging_map, subject_id == "PT-0018")
#FSPGR and FSPGR_contrast

subset(imaging_map, subject_id == "PT-0019")
#FSPGR and FSPGR_contrast


##potential confounders: look for scanner differences 
#T1-weighted: 
#GE: FSPGR | FSPGR_contrast = after a contrast agent (typically Gadolinium) used to enhanced signals in tumors, inflammation, etc.
#Siemens: MPRAGE | Ax T1 MPRAGE = Axial slice only?

# Extract unique subject IDs from the dictionary
IID_ISMMS <- unique(dictionary$IID_ISMMS)
# Subset imaging_map based on matching subject_id
imaging_map_subset <- imaging_map[subject_id %in% IID_ISMMS]

dim(imaging_map_subset) #309x9

table(imaging_map_subset$SeriesDescription2)
#   Ax T1 MPRAGE           FSPGR  FSPGR_contrast          MPRAGE MPRAGE_contrast 
#              1              56              39              95             118 

subset(imaging_map_subset, SeriesDescription2 == "Ax T1 MPRAGE") #PT-0064
subset(imaging_map_subset, subject_id == "PT-0064") 
#Ax T1 MPRAGE, FSPGR, and MPRAGE_contrast

##there are scans pre-surgery and post-surgery (which has the electrode) - subset to only scans pre-surgery?
#make two columns to compare acquisition data to surgery dates first 
imaging_map_subset$surgery1_timing <- ifelse(
  imaging_map_subset$AcquisitionDate < imaging_map_subset$surgery1,
  "pre-surgery", "post-surgery"
)

imaging_map_subset$surgery2_timing <- ifelse(
  imaging_map_subset$AcquisitionDate < imaging_map_subset$surgery2,
  "pre-surgery", "post-surgery"
)
#is this actually meaningful? surgery2 is always after surgery1 and pre-surgery of surgery2 doesn't mean true pre-surgery (aka before surgery1)

presurgery_scans <- filter(imaging_map_subset, surgery1_timing == "pre-surgery")
dim(presurgery_scans) #256  11

table(presurgery_scans$surgery1_timing) #TRUE
table(presurgery_scans$surgery2_timing) #TRUE

length(unique(presurgery_scans$subject_id)) #151 people, so some have multiple scans | also i think for dictionary i have 152 people?
table(presurgery_scans$subject_id)
table(dictionary$IID_ISMMS)

setdiff(unique(dictionary$IID_ISMMS), unique(presurgery_scans$subject_id)) 
##PT-0025 is in dictionary but not presurgery_scans

table(presurgery_scans$SeriesDescription2) #total: 256 scans ~ 1 week? double check with someone
#          FSPGR  FSPGR_contrast          MPRAGE MPRAGE_contrast 
#              4              39              95             118

length(unique(presurgery_scans$scan_id)) #256, good this match the dimension of the df so this means that all pt has unique scan number 
head(presurgery_scans)
#      V1 subject_id scan_id AcquisitionDate ContentTime SeriesDescription2
#1:    18    PT-0019  scan31      2014-03-28      140642     FSPGR_contrast
presurgery_scans$V1<-NULL
write.csv(presurgery_scans, "/sc/arion/projects/mscic1/results/jolie/LBP/imaging/presurgery_scans_subset_with_blood_brain_ge.csv")

test<-read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/imaging/presurgery_scans_subset_with_blood_brain_ge.csv")
test$X <- NULL

```

## Continue to examine the null results: 08-05-2025 and 08-06-2025
To avoid redundancy, see code above, otherwise new code will be posted below.
Also explored:
- Any missingness in anesthesia and medications?

```{r}
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
set.seed(2025)
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

#Reminder 
#this is NO residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")
#this is with residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")

##looked at the graphs above- 
https://hoangd02.u.hpc.mssm.edu/plots/lbp/shuffling_summary_no_residID.pdf
https://hoangd02.u.hpc.mssm.edu/plots/lbp/shuffling_summary_with_residID.pdf

#how far apart are the surgeries
dim(blood_metadata) #225 340
dim(brain_metadata) #225 340
dim(blood_only_metadata) #225 171
dim(brain_only_metadata) #225 171

head(dictionary)
#adding mymet_extractiondate_brain and mymet_tissue_brain
library(dplyr)

# Step 1: Join brain extraction dates
dictionary <- dictionary %>%
  left_join(
    brain_metadata %>%
      select(IID_ISMMS, mymet_tissue = mymet_tissue_brain, mymet_extractiondate_brain),
    by = c("IID_ISMMS", "mymet_tissue_brain" = "mymet_tissue")
  )

# Step 2: Join blood extraction dates
dictionary <- dictionary %>%
  left_join(
    brain_metadata %>%
      select(IID_ISMMS, mymet_tissue = mymet_tissue_blood, mymet_extractiondate_blood),
    by = c("IID_ISMMS", "mymet_tissue_blood" = "mymet_tissue")
  )
tail(dictionary)
# not all blood and brain are taken on the same day it seems 

# Make sure both date columns are of Date type
dictionary$mymet_extractiondate_brain <- as.Date(dictionary$mymet_extractiondate_brain)
dictionary$mymet_extractiondate_blood <- as.Date(dictionary$mymet_extractiondate_blood)

# Calculate date difference (positive = blood after brain)
dictionary$extraction_date_diff_days <- as.numeric(
  dictionary$mymet_extractiondate_blood - dictionary$mymet_extractiondate_brain
)

#some blood are taken before and after brain, only 2 blood and brain are taken the same day? 
table(dictionary$extraction_date_diff_days)
#-103 -100  -83  -80  -60  -57  -35  -33  -30  -29  -26  -23  -22  -19  -16   -6 
#   1    2    1    2    1    1    2    1    1    1    2    2    1    1    1    1 
#  -3   -2    0    1    4    6    7   11   17   20   21   24   26   27   28   29 
#   1    2    2    3    2    1    2    2    1    4    4    1    1    4    3    6 
#  34   35   42   43   44   46   48   49   51   55   56   62   63   65   69   70 
#   1    3    1    3    3    2    2    2    1    2    2    7    3    1    3    1 
#  71   72   75   76   77   78   79   81   83   85   88   90   95   97   98   99 
#   2    1    3    1    3    3    1    4    1    4    3    2    3    4    1    3 
# 101  102  103  104  105  108  110  111  113  114  117  118  122  124  126  130 
#   1    1    5    2    3    2    3    2    3    2    3    1    2    3    4    1 
# 131  138  139  144  145  146  147  151  153  154  159  160  161  162  167  169 
#   1    3    1    2    1    2    1    3    1    3    1    3    1    1    3    3 
# 172  174  176  181  183  186  187  190  194  211  221  228  229  526  546  595 
#   1    1    2    1    2    3    1    1    1    1    1    1    2    1    1    1 
summary(dictionary$extraction_date_diff_days)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -103.0    29.0    81.0    84.4   124.0   595.0 

##calculating the difference between L_Brain and R_Brain - subset for only 2 pairs 
two_pairs <- filter(dictionary, number_of_pair_brain == "2_pairs")

# Initialize new column with NA
two_pairs$LR_brain_extraction_diff_days <- NA_real_

# Get unique patient IDs
unique_ids <- unique(two_pairs$IID_ISMMS)

# Loop over each patient
for (id in unique_ids) {
  # Subset rows for this patient
  rows <- which(two_pairs$IID_ISMMS == id)
  
  if (length(rows) == 2) {
    tissues <- two_pairs$mymet_tissue_brain[rows]
    dates <- two_pairs$mymet_extractiondate_brain[rows]
    
    # Check if both L_Brain and R_Brain are present
    if (all(c("L_Brain", "R_Brain") %in% tissues)) {
      # Get dates
      date_L <- dates[tissues == "L_Brain"]
      date_R <- dates[tissues == "R_Brain"]
      
      # Calculate difference Right minus Left 
      diff_days <- as.numeric(date_R - date_L)
      
      # Assign to both rows
      two_pairs$LR_brain_extraction_diff_days[rows] <- diff_days
    }
  }
}

table(two_pairs$LR_brain_extraction_diff_days)
#-170 -165 -147 -134 -106 -103  -93  -79  -78  -76  -75  -73  -72  -60  -51  -50 
#   2    2    2    2    2    2    2    4    2    2    2    2    2    2    2    2 
# -44  -41  -28  -21  -15  -13   -6   -4    0    1    4    6    8   10   13   14 
#   4    2    2    2    2    2    4    2   10    2    2    2    2    2    2    2 
#  22   29   33   36   40   41   43   55   57   61   68   72   75   84   98   99 
#   6    2    6    2    2    2    2    2    2    2    2    2    2    2    2    2 
# 101  104  117  120  124  128  129  147  170  173  179  186  445  484 
#   2    2    2    2    2    2    2    2    2    2    2    2    2    2

summary(two_pairs$LR_brain_extraction_diff_days)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-170.00  -44.00   10.00   25.38   72.00  484.00       7 

#484 days is almost 16 months!

## Missingness in anesthesia and medications
dop_daily_dose <- readRDS("/sc/arion/projects/psychgen/lbp/data/emr/lbp_lel2021_clinical/dopamine_daily_dosing.RDS")
dim(dop_daily_dose) #166 7
length(unique(dop_daily_dose$iid)) #108

dop_daily_dose <- dop_daily_dose[!duplicated(iid)] #getting rid of duplicates, keeping only the first row for each unique iids (the values duplicate for duplicated IDs)
dim(dop_daily_dose) #108   7

dop_daily_dose$IID_ISMMS <- dop_daily_dose$iid
dop_daily_dose_sub <- dop_daily_dose[,c("IID_ISMMS", "carbidopa","levodopa","entacapone")]
head(dop_daily_dose_sub)
#   IID_ISMMS carbidopa levodopa entacapone
#      <char>     <num>    <num>      <num>
#1:   PT-0002     25.00      100          0
#2:   PT-0010    312.50     1250          0
#3:   PT-0012    187.50      750       1000

blood_only_metadata <- merge(blood_only_metadata, dop_daily_dose_sub, by = "IID_ISMMS", all.x = TRUE)
dim(blood_only_metadata) #225 173

summary(blood_only_metadata$carbidopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   25.0   100.0   175.0   194.6   250.0   793.8      84 
summary(blood_only_metadata$levodopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    100     400     750     832    1050    3175      84 
summary(blood_only_metadata$entacapone)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.0     0.0     0.0   106.4     0.0  6000.0      84 

brain_only_metadata <- merge(brain_only_metadata, dop_daily_dose_sub, by = "IID_ISMMS", all.x = TRUE)
dim(brain_only_metadata) #225 173

summary(brain_only_metadata$carbidopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   25.0   100.0   175.0   194.6   250.0   793.8      84 
summary(brain_only_metadata$levodopa)
summary(brain_only_metadata$entacapone)
## same as above - for blood (the samples are paired so it should be the same!)

##need to take a closer look at the NAs..... and not sure if merge IID_ISMMS is the best way

####################################### TEST TEST TEST 
##99.95th percentile of form3_no_residID_summary
#aka top 0.05 correlations pair: 224,734 gene-pairs:
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_form3_by_brain_cor) #21046 21357
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
summary_stats <- summary(as.vector(as.matrix(blood_form3_by_brain_cor))); summary_stats
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02266 0.04790 0.05637 0.08143 0.87792 
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2414433 | strongest correlations (top 0.05%, 224k gene-pairs) are ≥ 0.2414433
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2862253

# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= blood_form3_99.95, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
dim(top_gene_pairs) #224734      3
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
head(top_gene_pairs)
length(unique(top_gene_pairs$brain_gene)) #16789
length(unique(top_gene_pairs$blood_gene)) #19683

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#168716 ENSG00000226259.10 ENSG00000226259.10   0.8779214
#169105  ENSG00000226752.9  ENSG00000226752.9   0.8499747
#180166  ENSG00000241945.8  ENSG00000241945.8   0.8441066
#85321  ENSG00000145736.14 ENSG00000145736.14   0.8389223
#178047  ENSG00000237541.3  ENSG00000237541.3   0.8304720
#207509  ENSG00000274602.5  ENSG00000274602.5   0.8243647

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE)

# Total number of pairs
n_total <- nrow(top_gene_pairs)

# Number of different gene pairs
n_different <- n_total - n_same

# Print the result
cat("Same gene pairs:", n_same, "\nDifferent gene pairs:", n_different, "\n")
#Same gene pairs: 507 
#Different gene pairs: 224227 

##Building Elastic net using only top genes 
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
#brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL

blood_expr <- blood_form3
brain_expr <- brain_full

# Transpose to samples x genes
blood_mat <- t(blood_expr)
brain_mat <- t(brain_expr)

# Predict brain gene ENSG00000226259.10
target_gene <- "ENSG00000226259.10"

sub<-subset(top_gene_pairs, brain_gene == "ENSG00000226259.10"); sub

# Get ALL blood genes highly correlated with this brain gene
predictors <- top_gene_pairs %>%
  filter(brain_gene == target_gene) %>%
  pull(blood_gene)

predictors
#[1] "ENSG00000226259.10" "ENSG00000145736.14" "ENSG00000230847.4" 
#[4] "ENSG00000215630.6"  "ENSG00000235558.3"  "ENSG00000050165.17"
#[7] "ENSG00000104722.14"

#pick top 1, 3, or 5: 
predictors <- top_gene_pairs %>%
  filter(brain_gene == target_gene) %>%
  top_n(5, wt = correlation) %>%
  pull(blood_gene)

# Subset predictor matrix
X <- blood_mat[, predictors, drop = FALSE]
y <- brain_mat[, target_gene]

# Use Elastic Net
library(glmnet)
model <- cv.glmnet(X, y, alpha = 0.5)  # alpha=0.5 = elastic net | 10 CV by default
pred <- predict(model, newx = X, s = "lambda.min")
cor(pred, y)  

# For glmnet predictions
r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
r2  # 0.09070857

##########################################################################################
########################## HELD OUT TEST SET FOR ONE BRAIN GENE ##########################
library(glmnet)

# Transpose matrices to samples × genes
X_all <- t(blood_expr)
Y_all <- t(brain_expr)

# Set your target brain gene
target_gene <- "ENSG00000226259.10"

# Step 1: Get top-correlated blood genes for this brain gene
top_blood_genes <- top_gene_pairs %>%
  filter(brain_gene == target_gene) %>%
  pull(blood_gene)

# Subset predictors (X) and response (y)
X <- X_all[, top_blood_genes, drop = FALSE]
y <- Y_all[, target_gene]

# Step 2: Train/test split
set.seed(2025)
n <- nrow(X)
train_idx <- sample(1:n, size = floor(0.8 * n))  # 80% train
test_idx <- setdiff(1:n, train_idx)             # 20% test

X_train <- X[train_idx, ]
y_train <- y[train_idx]

X_test <- X[test_idx, ]
y_test <- y[test_idx]

# Step 3: Fit model on training set
model <- cv.glmnet(X_train, y_train, alpha = 0.5)

# Step 4: Predict on test set
y_pred <- predict(model, newx = X_test, s = "lambda.min")

# Step 5: Evaluate on test set
cor_test <- cor(y_pred, y_test)
r2_test <- 1 - sum((y_test - y_pred)^2) / sum((y_test - mean(y_test))^2)

cat("Test set Pearson correlation:", round(cor_test, 3), "\n")
cat("Test set R²:", round(r2_test, 3), "\n")

## For 1 brain gene: ENSG00000226259.10
#Test set Pearson correlation: 0.304  
#Test set R²: 0.083 

##########################################################################################
########################## HELD OUT TEST SET FOR MANY BRAIN GENES #########################
library(glmnet)
library(dplyr)

################################### no_indivdualID_brain_full ###############################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_form3_by_brain_cor) #21046 21357
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
summary_stats <- summary(as.vector(as.matrix(blood_form3_by_brain_cor))); summary_stats
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02266 0.04790 0.05637 0.08143 0.87792 
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2414433 | strongest correlations (top 0.05%, 224k gene-pairs) are ≥ 0.2414433
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2862253

# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= blood_form3_99.95, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
dim(top_gene_pairs) #224734      3
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
head(top_gene_pairs)
length(unique(top_gene_pairs$brain_gene)) #16789
length(unique(top_gene_pairs$blood_gene)) #19683

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#168716 ENSG00000226259.10 ENSG00000226259.10   0.8779214
#169105  ENSG00000226752.9  ENSG00000226752.9   0.8499747
#180166  ENSG00000241945.8  ENSG00000241945.8   0.8441066
#85321  ENSG00000145736.14 ENSG00000145736.14   0.8389223
#178047  ENSG00000237541.3  ENSG00000237541.3   0.8304720
#207509  ENSG00000274602.5  ENSG00000274602.5   0.8243647

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #621
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224113

# Print the result
cat("Same gene pairs:", n_same, "\nDifferent gene pairs:", n_different, "\n")
#Same gene pairs: 507 
#Different gene pairs: 224227 
brain_genes_to_test <- unique(top_gene_pairs$brain_gene)

results_df <- predict_brain_gene_expression(
  X_all = t(blood_expr),
  Y_all = t(brain_expr),
  top_gene_pairs = top_gene_pairs,
  brain_genes = brain_genes_to_test,
  train_idx = train_idx,
  test_idx = test_idx
)

write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_full_genes_elastic_net_blood_form3_no_indivdualID.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

###################################################################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_form3_by_brain_cor) #21046 21357
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1
blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
summary_stats <- summary(as.vector(as.matrix(blood_form3_by_brain_cor))); summary_stats
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02399 0.05052 0.05901 0.08524 0.83440 
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2484 | strongest correlations (top 0.05%, 224k gene-pairs) are ≥ 0.2484
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2906922 

# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= blood_form3_99.95, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)

dim(top_gene_pairs) #224734      3
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
head(top_gene_pairs)
length(unique(top_gene_pairs$brain_gene)) #15695
length(unique(top_gene_pairs$blood_gene)) #15366

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#167582 ENSG00000226259.10 ENSG00000226259.10   0.8344016
#181336  ENSG00000241945.8  ENSG00000241945.8   0.8261620
#210686  ENSG00000274602.5  ENSG00000274602.5   0.8222672
#208259 ENSG00000233327.10  ENSG00000273018.6   0.8079362
#208262  ENSG00000273018.6  ENSG00000273018.6   0.8060124
#175112 ENSG00000233327.10 ENSG00000233327.10   0.8029794

set.seed(2025)
brain_genes_to_test <- unique(top_gene_pairs$brain_gene)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
dim(blood_form3) #21046   225

blood_expr <- blood_form3
brain_expr <- v_brain$E

# Transpose expression matrices: samples × genes
X_all <- t(blood_expr)   # blood gene expression (samples x blood genes)
Y_all <- t(brain_expr)   # brain gene expression (samples x brain genes)

# Train/test split
set.seed(2025)
n <- nrow(X_all)  # number of samples
train_idx <- sample(1:n, size = floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

# Brain genes to evaluate
brain_genes_to_test <- unique(top_gene_pairs$brain_gene)

results_df <- predict_brain_gene_expression(
  X_all = t(blood_expr),
  Y_all = t(brain_expr),
  top_gene_pairs = top_gene_pairs,
  brain_genes = brain_genes_to_test,
  train_idx = train_idx,
  test_idx = test_idx
)

# Save to file
write.table(results_df, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_brain_base_genes_elastic_net_blood_form3_no_indivdualID.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

```


## Predict Gene Expression Function - 08-06-2025

```{r}
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

```



