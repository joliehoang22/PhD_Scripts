## code derived from LBP_blood_brain_QC_20250729.R
library(data.table)
library(ggplot2)
library(readxl)
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

## swapping 0%, 10% (8 pairs), 25% (20 pairs), 50% (40 pairs), 75% (60 pairs), and 100% (80 pairs)

#this is NO residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")
names(null_results)
#[1] "swap_0pct"   "swap_10pct"  "swap_25pct"  "swap_50pct"  "swap_75pct" 
#[6] "swap_100pct"

swap_0pct <- null_results[["swap_0pct"]]
swap_10pct <- null_results[["swap_10pct"]]
swap_25pct <- null_results[["swap_25pct"]]
swap_50pct <- null_results[["swap_50pct"]]
swap_75pct <- null_results[["swap_75pct"]]
dim(null_results[["swap_0pct"]])
# [1] 1000   19

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results_swap100.RData") #227.73 mins = 4 hrs
swap_100pct <- null_results[["swap_100pct"]]

#this is with residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")
swap_0pct <- null_results[["swap_0pct"]]
swap_10pct <- null_results[["swap_10pct"]]
swap_25pct <- null_results[["swap_25pct"]]
swap_50pct <- null_results[["swap_50pct"]]
swap_75pct <- null_results[["swap_75pct"]]
swap_100pct <- null_results[["swap_100pct"]]

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

ggsave("/hpc/users/hoangd02/www/plots/lbp/shuffling_summary_no_residID.pdf", plot = plot, width = 16, height = 8) 

#ggsave("/hpc/users/hoangd02/www/plots/lbp/shuffling_summary_no_residID.pdf", plot = plot, width = 16, height = 8) 
#ggsave("/hpc/users/hoangd02/www/plots/lbp/shuffling_summary_with_residID.pdf", plot = plot, width = 16, height = 8) 

###################
###################
###################
#this is NO residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise__no_indivdualID_null_results.RData")
#this is with residID
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/shuffling_samples_stepwise_null_results.RData")

##looked at the graphs above- 
https://hoangd02.u.hpc.mssm.edu/plots/lbp/shuffling_summary_no_residID.pdf
https://hoangd02.u.hpc.mssm.edu/plots/lbp/shuffling_summary_with_residID.pdf

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")

#how far apart are the surgeries
dim(blood_metadata) #225 340
dim(brain_metadata) #225 340
dim(blood_only_metadata) #225 171
dim(brain_only_metadata) #225 171

head(dictionary)
test<-filter(dictionary, IID_ISMMS =="PT-0024")
#adding mymet_extractiondate_brain and mymet_tissue_brain to dictionary

dictionary <- dictionary %>%
  left_join(
    brain_metadata %>%
      select(IID_ISMMS, mymet_tissue = mymet_tissue_brain, mymet_extractiondate_brain),
    by = c("IID_ISMMS", "mymet_tissue_brain" = "mymet_tissue")
  )

dictionary <- dictionary %>%
  left_join(
    brain_metadata %>%
      select(IID_ISMMS, mymet_tissue = mymet_tissue_blood, mymet_extractiondate_blood),
    by = c("IID_ISMMS", "mymet_tissue_blood" = "mymet_tissue")
  )

tail(dictionary)
# not all blood and brain are taken on the same day it seems hmmm

#make sure both date columns are of Date type
dictionary$mymet_extractiondate_brain <- as.Date(dictionary$mymet_extractiondate_brain)
dictionary$mymet_extractiondate_blood <- as.Date(dictionary$mymet_extractiondate_blood)

#calculate date difference (positive = blood after brain)
dictionary$extraction_date_diff_days <- as.numeric(
  dictionary$mymet_extractiondate_blood - dictionary$mymet_extractiondate_brain)

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
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
#-170.00  -44.00   10.00   25.38   72.00  484.00       

#484 days is almost 16 months!

### how far apart are the blood samples?
two_pairs$LR_blood_extraction_diff_days <- NA_real_

for (id in unique_ids) {
  rows <- which(two_pairs$IID_ISMMS == id)
  
  if (length(rows) == 2) {
    tissues <- two_pairs$mymet_tissue_blood[rows]          # <-- change if needed
    dates   <- two_pairs$mymet_extractiondate_blood[rows]  # <-- change if needed
    
    if (all(c("L_Blood", "R_Blood") %in% tissues) && all(!is.na(dates))) { # <-- adjust labels
      date_L <- dates[tissues == "L_Blood"]
      date_R <- dates[tissues == "R_Blood"]
      diff_days <- as.numeric(date_R - date_L)
      two_pairs$LR_blood_extraction_diff_days[rows] <- diff_days
    }
  }
}

table(two_pairs$LR_blood_extraction_diff_days)
#-114 -111 -104  -97  -91  -84  -77  -61  -54  -43  -23  -20  -17   -7    0    7 
#   2    2    2    4    2    6    2   14    2    8    4    2    2    6   38    6 
#  14   20   23   43   50   61   77   84   91   97  104  111 
#   4    4    2    2    2    8    4    6    2    2    4    4 
summary(two_pairs$LR_blood_extraction_diff_days)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
#-114.00  -43.00    0.00   -3.11   20.00  111.00       

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/matching_pairs.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/matching_pairs.RData")




blood_samples_to_be_removed <- c("LBPSEMA4BLOOD795", "LBPSEMA4BLOOD142", "LBPSEMA4BLOOD049", "LBPSEMA4BLOOD755",
                                  "LBPSEMA4BLOOD335","LBPSEMA4BLOOD431","LBPSEMA4BLOOD174","LBPSEMA4BLOOD738")

brain_samples_to_be_removed <- c("LBPSEMA4BRAIN364", "LBPSEMA4BRAIN318", "LBPSEMA4BRAIN564", "LBPSEMA4BRAIN170",
                                  "LBPSEMA4BRAIN703", "LBPSEMA4BRAIN341","LBPSEMA4BRAIN391","LBPSEMA4BRAIN017")

