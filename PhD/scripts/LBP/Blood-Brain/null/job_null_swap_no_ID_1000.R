#code derived from two_pairs_investigation.R 
#this code is too computational expensive to run concurrently with other scripts so 
#imma submit a job for it 

## Job submussions
path=/sc/arion/projects/mscic1/results/jolie/LBP/scripts/
cd $path
ml R

bsub -q premium -P acc_mscic1 -n 30 -W 144:00 -R rusage[mem=20000] -R span[hosts=1] -R himem -o %J.stdout -eo %J.stderr Rscript ${path}job_null_swap_no_ID_1000.R
bsub -q premium -P acc_mscic1 -n 30 -W 144:00 -R rusage[mem=20000] -R span[hosts=1] -R himem -o %J.stdout -eo %J.stderr Rscript ${path}job_null_swap_with_ID_1000.R

##resubmitted on Aug 18th 2025
#Job <199127166> is submitted to queue <premium>.
#Job <199127168> is submitted to queue <premium>.
#--
#Job <198844848> is submitted to queue <premium> null_swap_with_ID.R | out of memory
#Job <198844898> is submitted to queue <premium>. null_swap_no_ID.R | out of memory
#--
#Job <198849035> is submitted to queue <premium>. | out of memory
#Job <198849036> is submitted to queue <premium>. | out of memory
#--
#Job <198874416> is submitted to queue <premium>.
#Job <198874417> is submitted to queue <premium>.
#--
#Job <198883732> is submitted to queue <premium>. | requested himem | might ran out of time with 144 hrs (6 days)
#Job <198883735> is submitted to queue <premium>. | requested himem
#bjobs

#this code runs the null_swap with 10000 permutations for blood form3 with IDs
#this script is in /sc/arion/projects/mscic1/results/jolie/LBP/scripts
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
library(tibble)
library(furrr)      # for future_map
library(future)     # for plan()
library(purrr)      # for map-style helpers

set.seed(2025)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")
table(brain_metadata$number_of_pair_brain) ##

#blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE) ##CHANGE 
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_form = blood_form3 

# Swap levels in percentages
swap_levels <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
n_perm <- 1000 #this was 1000, changed to 10000, change back to 1000

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
save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_no_indivdualID_null_results_1000perm.RData")

#save(null_results, file = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_no_indivdualID_null_results_10000perm.RData")

###### CHANGE BETWEEN THESE !!!
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_no_indivdualID_null_results_1000perm.RData")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_with_indivdualID_null_results_1000perm.RData")

##code derived from examining_swap_results.R
names(null_results)

swap_0pct <- null_results[["swap_0pct"]]
swap_10pct <- null_results[["swap_10pct"]]
swap_25pct <- null_results[["swap_25pct"]]
swap_50pct <- null_results[["swap_50pct"]]
swap_75pct <- null_results[["swap_75pct"]]
swap_100pct <- null_results[["swap_100pct"]]
dim(null_results[["swap_0pct"]])
# [1] 1000   19

##all rows are identical
all_rows_identical <- length(unique(apply(swap_0pct, 1, paste, collapse = "_"))) == 1
all_rows_identical #TRUE 

all_rows_identical <- length(unique(apply(swap_50pct, 1, paste, collapse = "_"))) == 1
all_rows_identical #FALSE

all_rows_identical <- length(unique(apply(swap_100pct, 1, paste, collapse = "_"))) == 1
all_rows_identical #TRUE 

## JUST COMPARING THE 0% AND 100%
swap_0pct <- swap_0pct[1, ]
colnames(swap_0pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(swap_0pct)) #renaming
#  min        p25       mean     median        p80       p90     p92.5       p95
#1   0 0.02399389 0.05900994 0.05051833 0.09466814 0.1205099 0.1300895 0.1427349
#      p97.5     p98.5     p99.5    p99.95   p99.995  p99.9995 p99.99995
#1 0.1624989 0.1758765 0.2020133 0.2483944 0.2906922 0.4383699  0.746393
#  p99.999995 p99.9999995 p99.99999995       max
#1  0.7723494   0.8187232    0.8325499 0.8344016

swap_100pct <- swap_100pct[1, ]
colnames(swap_100pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(swap_100pct)) #renaming
#  min        p25       mean     median       p80       p90     p92.5       p95
#1   0 0.02661926 0.06485842 0.05609566 0.1046745 0.1321018 0.1419574 0.1546555
#      p97.5     p98.5     p99.5    p99.95   p99.995 p99.9995 p99.99995
#1 0.1738243 0.1863675 0.2100801 0.2508944 0.2881974 0.434954 0.7336171
#  p99.999995 p99.9999995 p99.99999995       max
#1    0.76765    0.822049    0.8332175 0.8344617

#p99.95 ~ 200k gene-pairs
#p99.995 ~ 20k gene-pairs - this is where 0 swaps start to get higher - PAIRING STARTS TO MATTER WHEN IT COMES TO THE TOP GENES
#p99.9995 ~ 2k gene-pairs
#p99.99995 ~ 200 gene-pairs
#p99.999995 ~ 20 gene-pairs

swap_0pct$swap_pct   <- 0
swap_100pct$swap_pct <- 100

df <- bind_rows(swap_0pct, swap_100pct)

# choose which measures to plot (edit this vector as you like)
keep_measures <- c("min","median", "p25", "p99.95","p99.995","p99.9995","p99.99995","p99.999995","max")

# long format for plotting
df_long <- df %>%
  select(swap_pct, all_of(keep_measures)) %>%
  pivot_longer(cols = -swap_pct, names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = keep_measures))

plot <- ggplot(df_long, aes(x = swap_pct, y = value, color = measure, group = measure)) +
  geom_line() +  # Remove position_dodge - we want lines to connect 0 to 100
  geom_point(size = 2) +  # Remove position_dodge - we want points at exactly 0 and 100
  geom_text(aes(label = sprintf("%.3f", value)),
            vjust = -0.6, size = 3, show.legend = FALSE) +
  labs(x = "Swap Percentage", y = "Value", color = "Measure") +
  scale_x_continuous(breaks = c(0, 100), limits = c(-5, 105)) +  # Slightly wider limits for text
  coord_cartesian(ylim = c(0, max(df_long$value) * 1.1), clip = "off") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/shuffling_summary_with_residID_0and100swaps.pdf", plot = plot, width = 7, height = 6) 
#ggsave("/hpc/users/hoangd02/www/plots/shuffling_summary_no_residID_0and100swaps.pdf", plot = plot, width = 7, height = 6) 

####### for swap_10pct, swap_25pct, etc. take the average of all 1000 rows for each column
avg_swap_10pct <- data.frame(t(colMeans(swap_10pct)))
avg_swap_25pct <- data.frame(t(colMeans(swap_25pct)))
avg_swap_50pct <- data.frame(t(colMeans(swap_50pct)))
avg_swap_75pct <- data.frame(t(colMeans(swap_75pct)))

colnames(avg_swap_10pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(avg_swap_10pct)) #renaming
colnames(avg_swap_25pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(avg_swap_25pct)) #renaming
colnames(avg_swap_50pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(avg_swap_50pct)) #renaming
colnames(avg_swap_75pct) <- sub("^(p\\d+(?:\\.\\d+)?)\\..*$", "\\1", colnames(avg_swap_75pct)) #renaming

avg_swap_10pct <- avg_swap_10pct %>%
  rename(
    p25 = p25.25,
    p80 = p80.80,
    p90 = p90.90,
    p95 = p95.95
  )
avg_swap_25pct <- avg_swap_25pct %>%
  rename(
    p25 = p25.25,
    p80 = p80.80,
    p90 = p90.90,
    p95 = p95.95
  )
avg_swap_50pct <- avg_swap_50pct %>%
  rename(
    p25 = p25.25,
    p80 = p80.80,
    p90 = p90.90,
    p95 = p95.95
  )
avg_swap_75pct <- avg_swap_75pct %>%
  rename(
    p25 = p25.25,
    p80 = p80.80,
    p90 = p90.90,
    p95 = p95.95
  )
swap_0pct$swap_pct   <- 0
avg_swap_10pct$swap_pct  <- 10
avg_swap_25pct$swap_pct  <- 25
avg_swap_50pct$swap_pct  <- 50
avg_swap_75pct$swap_pct  <- 75
swap_100pct$swap_pct <- 100

df <- bind_rows(swap_0pct, avg_swap_25pct, avg_swap_50pct, avg_swap_75pct, swap_100pct) #didnt add 10%

# choose which measures to plot (edit this vector as you like)
keep_measures <- c("min","median", "p25", "p99.95","p99.995","p99.9995","p99.99995","p99.999995","max")

# long format for plotting
df_long <- df %>%
  select(swap_pct, all_of(keep_measures)) %>%
  pivot_longer(cols = -swap_pct, names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = keep_measures))

plot <- ggplot(df_long, aes(x = swap_pct, y = value, color = measure, group = measure)) +
  geom_line() +  # Remove position_dodge - we want lines to connect 0 to 100
  geom_point(size = 2) +  # Remove position_dodge - we want points at exactly 0 and 100
  geom_text(aes(label = sprintf("%.3f", value)),
            vjust = -0.6, size = 3, show.legend = FALSE) +
  labs(x = "Swap Percentage", y = "Value", color = "Measure") +
  #scale_x_continuous(breaks = c(0, 100), limits = c(-5, 105)) +  # Slightly wider limits for text
  coord_cartesian(ylim = c(0, max(df_long$value) * 1.1), clip = "off") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/shuffling_summary_with_residID_all_swaps.pdf", plot = plot, width = 7, height = 6) 
#ggsave("/hpc/users/hoangd02/www/plots/shuffling_summary_no_residID_all_swaps.pdf", plot = plot, width = 7, height = 6) 

################## ARCHIVED ################## 
#This approach is not good anymore cuz all permutation yielded the same result for swap_0pct and swap_100pct, but not the between percentiles, 
#which makes sense cuz there are different ways to swap 10% etc. but 1 way to swap 0% and 100%

#summary of each column of each swap 
summary_stats_0pct <- sapply(head(swap_0pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
})
t0 <- t(summary_stats_0pct)  #transpose to make it more readable

summary_stats_10pct <- sapply(head(swap_10pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
})
t10 <- t(summary_stats_10pct)

summary_stats_25pct <- sapply(head(swap_25pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
})
t25 <- t(summary_stats_25pct)

summary_stats_50pct <- sapply(head(swap_50pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
})
t50 <- t(summary_stats_50pct)

summary_stats_75pct <- sapply(head(swap_75pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
})
t75 <- t(summary_stats_75pct)

summary_stats_100pct <- sapply(head(swap_100pct), function(x) {
  c(mean = mean(x), median = median(x), p25 = quantile(x, 0.25), p99.95 = quantile(x, 0.9995))
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
  pivot_longer(cols = c("mean", "median", "p25.25%", "p99.95.99.95%"), names_to = "measure", values_to = "value")

plot<- ggplot(df_long, aes(x = swap_pct, y = value, color = stat, group = stat)) +
  geom_line() +
  facet_wrap(~ measure) +
  labs(x = "Swap Percentage", y = "Value", color = "Summary Statistic") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/shuffling_summary_no_residID_1000perm_try2.pdf", plot = plot, width = 16, height = 8) 

library(dplyr)
sva <- read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim.csv")
sva_oracle <- filter(sva, method == "oracle"); dim(sva_oracle)
sva_none <- filter(sva, method == "none"); dim(sva_none)
sva_known_nsv <- filter(sva, method == "known_nsv"); dim(sva_known_nsv)
sva_be <- filter(sva, method == "be"); dim(sva_be)
sva_leek <- filter(sva, method == "leek"); dim(sva_leek)

write.csv(sva_oracle, "/sc/arion/projects/mscic1/results/jolie/sva_oracle.csv")
write.csv(sva_none, "/sc/arion/projects/mscic1/results/jolie/sva_none.csv")
write.csv(sva_known_nsv, "/sc/arion/projects/mscic1/results/jolie/sva_known_nsv.csv")
write.csv(sva_be, "/sc/arion/projects/mscic1/results/jolie/sva_be.csv")
write.csv(sva_leek, "/sc/arion/projects/mscic1/results/jolie/sva_leek.csv")