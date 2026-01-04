#######
#######
#this code runs the null_swap with 10000 permutations for blood form5 no IDs
#this code is similar to job_null_swap_no_ID.R (which ran blood form3 no IDs)
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
library(ggsignif)
library(purrr)
library(caret)
library(tidyr)
library(tibble)
library(furrr)      # for future_map
library(future)     # for plan()
library(purrr)      # for map-style helpers

set.seed(2025)

#load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
head(dictionary)

#################################
############# NO ID #############
#################################
blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL

brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1; brain_full_no_residID$V1 <- NULL

blood_form = blood_form5_no_residID
brain_form = brain_full_no_residID

#################################
############# WITH ID ###########
#################################
# blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
# rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL

# brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
# row.names(brain_full) <- brain_full$V1; brain_full$V1 <- NULL

# blood_form = blood_form5
# brain_form = brain_full

# Identify the 2-pair individuals
two_pair_inds <- unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])
length(two_pair_inds) #73

# Deterministic 0% / 100% swap on brain columns
permute_brain_once <- function(brain_mat, brain_metadata, two_pair_inds, swap_frac) {
  stopifnot(swap_frac %in% c(0, 1))
  out <- brain_mat
  if (swap_frac == 1) {
    for (id in two_pair_inds) {
      idx <- which(brain_metadata$IID_ISMMS == id)
      if (length(idx) == 2) out[, idx] <- out[, rev(idx)]
    }
  }
  out
}

# Where to save
out_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_one_swap <- function(swap_frac) {
  cat("\n=== Computing full matrices for swap", swap_frac*100, "% ===\n")

  brain_once <- permute_brain_once(brain_form, brain_metadata, two_pair_inds, swap_frac)

  # Spearman correlations across samples (rows after transpose):
  # X = t(blood_form):   samples x blood_genes
  # Y = t(brain_once):   samples x brain_genes
  # cor(X, Y) returns (blood_genes x brain_genes) matrix
  cor_mat <- cor(t(blood_form), t(brain_once), method = "spearman")
  abs_mat <- abs(cor_mat)

  # Save (plain .rds)
  saveRDS(cor_mat, file.path(out_dir, sprintf("spearman_swap_%dpct_nonabs_with_residID.rds", as.integer(swap_frac*100))), compress = FALSE)
  saveRDS(abs_mat,  file.path(out_dir, sprintf("spearman_swap_%dpct_abs_with_residID.rds", as.integer(swap_frac*100))), compress = FALSE)

  rm(cor_mat, abs_mat); gc()
  cat("Saved swap", swap_frac*100, "% matrices.\n")
}

# Run once for 0% and 100%
save_one_swap(0)
save_one_swap(1)

##outputs 
spearman_swap_0pct_nonabs.rds
spearman_swap_0pct_abs.rds
spearman_swap_100pct_nonabs.rds
spearman_swap_100pct_abs.rds

############## ############## ############## ################### 
############### blood form5 no residIDs ############## #########
############## ############## ############## ################### 
no_swap_nonabs  <- readRDS(file.path(out_dir, "spearman_swap_0pct_nonabs_no_residID.rds"))
no_swap_abs <- readRDS(file.path(out_dir, "spearman_swap_0pct_abs_no_residID.rds"))
dim(no_swap_nonabs)   # 21046 x 21356
dim(no_swap_abs) #21046 21356
range(no_swap_nonabs)#-0.7925811  0.8780499
range(no_swap_abs) #0.0000000 0.8780499

all_swap_nonabs <- readRDS(file.path(out_dir, "spearman_swap_100pct_nonabs_no_residID.rds"))
all_swap_abs <- readRDS(file.path(out_dir, "spearman_swap_100pct_abs_no_residID.rds"))
dim(all_swap_nonabs)   # 21046 x 21356
dim(all_swap_abs) #21046 21356
range(all_swap_nonabs)#-0.7891351  0.8627823
range(all_swap_abs) #0.0000000 0.8627823

############## ############## ############## ################### 
############### blood form5 with residIDs ############## #######
############## ############## ############## ################### 
no_swap_nonabs  <- readRDS(file.path(out_dir, "spearman_swap_0pct_nonabs_with_residID.rds"))
no_swap_abs <- readRDS(file.path(out_dir, "spearman_swap_0pct_abs_with_residID.rds"))
dim(no_swap_nonabs)   # 21046 x 21356
dim(no_swap_abs) #21046 21356
range(no_swap_nonabs)#-0.4045438  0.4450948
range(no_swap_abs) #0.0000000 0.4450948

all_swap_nonabs <- readRDS(file.path(out_dir, "spearman_swap_100pct_nonabs_with_residID.rds"))
all_swap_abs <- readRDS(file.path(out_dir, "spearman_swap_100pct_abs_with_residID.rds"))
dim(all_swap_nonabs)   # 21046 x 21356
dim(all_swap_abs) #21046 21356
range(all_swap_nonabs)#-0.4947335  0.4527149
range(all_swap_abs) #0.0000000 0.4947335

##lets just focus on the absolute spearman correlation for now 
no_swap_abs and all_swap_abs

# ---- 1) Load matrices ----
##no residID ABSOLUTE VALUE
no_swap_abs  <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_0pct_abs_no_residID.rds")
all_swap_abs <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_100pct_abs_no_residID.rds")

#with residID ABSOLUTE VALUE
no_swap_abs  <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_0pct_abs_with_residID.rds")
all_swap_abs <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_100pct_abs_with_residID.rds")

########################################################################################################################################################################
##no residID NOT ABSOLUTE VALUE
no_swap_abs  <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_0pct_nonabs_no_residID.rds")
all_swap_abs <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_100pct_nonabs_no_residID.rds")

#with residID NOT ABSOLUTE VALUE
no_swap_abs  <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_0pct_nonabs_with_residID.rds")
all_swap_abs <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_100pct_nonabs_with_residID.rds")


# 2) Percentiles 
probs <- c(`90th`=0.90, `99th`=0.99, `99.95th`=0.9995, `99.995th`=0.99995)

# 3) Cutoffs for each matrix
cuts_no_swap  <- quantile(as.vector(no_swap_abs),  probs = probs, names = TRUE)
cuts_all_swap <- quantile(as.vector(all_swap_abs), probs = probs, names = TRUE)

##########################################
############### no residID ############### ABSOLUTE VALUE
cuts_no_swap
#       90%       99%    99.95%   99.995% 
# 0.1124294 0.1744648 0.2337137 0.2721019 
cuts_all_swap
#       90%       99%    99.95%   99.995% 
# 0.1106658 0.1722429 0.2311083 0.2688971 

############################################
############### with residID ############### ABSOLUTE VALUE
cuts_no_swap
#       90%       99%    99.95%   99.995% 
# 0.1107238 0.1719648 0.2304362 0.2677602 

cuts_all_swap
#       90%       99%    99.95%   99.995% 
# 0.1078529 0.1681553 0.2263264 0.2631753 


##########################################
############### no residID ############### NOT ABSOLUTE VALUE
cuts_no_swap
#        90%        99%     99.95%    99.995% 
# 0.08747261 0.15849558 0.22313949 0.26367678 
cuts_all_swap
#        90%        99%     99.95%    99.995% 
# 0.08778445 0.15765487 0.22056679 0.25985154 

# 4) Create “top-X%” matrices (values below cutoff → NA)
make_top_mats <- function(mat, cuts_named) {
  out <- lapply(names(cuts_named), function(nm) {
    thr <- cuts_named[[nm]]
    m   <- mat
    m[m < thr] <- NA_real_
    m
  })
  names(out) <- names(cuts_named)
  out
}

top_no_swap  <- make_top_mats(no_swap_abs,  cuts_no_swap)
top_all_swap <- make_top_mats(all_swap_abs, cuts_all_swap)

# access a specific matrix, e.g., top 99.95th for no-swap:
top_all_swap[["90%"]]       # 90th percentile matrix
top_all_swap[["99%"]]       # 99th
top_all_swap[["99.95%"]]    # 99.95th
top_all_swap[["99.995%"]]   # 99.995th

top_no_swap[["99.95%"]]    # 99.95th
top_no_swap[["99.995%"]]   # 99.995th

# 3) Quick sanity checks

##########################################
############### no residID ###############
# how many gene pairs survive at each threshold?
sapply(top_no_swap, function(m) sum(!is.na(m)))
#      90%      99%   99.95%  99.995% 
# 44946659  4494697   224737    22473 
#3 extra pairs that happened to have the same value as the threshold at 99.95%  

sapply(top_all_swap,  function(m) sum(!is.na(m)))
#      90%      99%   99.95%  99.995% 
# 44945943  4494617   224734    22473 

############################################
############### with residID ###############
# sapply(top_no_swap, function(m) sum(!is.na(m)))
#      90%      99%   99.95%  99.995% 
# 44946884  4494757   224735    22475 

sapply(top_all_swap,  function(m) sum(!is.na(m)))
#      90%      99%   99.95%  99.995% 
# 44946899  4494680   224735    22474 

# --- pick the percentile groups you want to compare ---
wanted <- c("99.95%", "99.995%")   # add "90%", "99%" if you’ve created those too

# --- helper to extract non-NA top values into a long table ---
extract_top_vals <- function(top_list, label, wanted) {
  rbindlist(lapply(wanted, function(pct) {
    v <- as.vector(top_list[[pct]])
    v <- v[!is.na(v)]
    data.table(value = v, swap = label, percentile = pct)
  }))
}

# Build tidy data: combine No swap vs All swap
dt <- rbindlist(list(
  extract_top_vals(top_no_swap,  "No swap (0%)",  wanted),
  extract_top_vals(top_all_swap, "All swap (100%)", wanted)
))

table(dt$swap, dt$percentile)    
############### no residID ###############        
  #                 99.95% 99.995%
  # No swap (0%)    224737   22473
  # All swap (100%) 224734   22473

############### with residID ###############
  #                 99.95% 99.995%
  # All swap (100%) 224735   22474
  # No swap (0%)    224735   22475

# order facet and legend nicely
# dt[, percentile := factor(percentile, levels = wanted)]
# dt[, swap := factor(swap, levels = c("No swap (0%)", "All swap (100%)"))]


# # --- Violin + boxplots, faceted by percentile ---
# k_map <- c("99.95%" = "224,734 gene pairs",
#            "99.995%" = "22,473 gene pairs")  # adjust if you include other percentiles

# lab_fun <- labeller(
#   percentile = function(x) paste0(x, " (top ", k_map[x], ")")
# )

# p <- ggplot(dt, aes(x = swap, y = value, fill = swap)) +
#   geom_violin(trim = TRUE, alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.size = 0.3) +
#   facet_wrap(~ percentile, nrow = 1, scales = "fixed", labeller = lab_fun) +
#   labs(
#     x = NULL, y = "Absolute Spearman correlation of all-gene concordance",
#     title = "Top blood brain gene pair correlations by percentile (no_residID)",
#   ) +
#   scale_y_continuous(limits = c(0.2, 1)) +   # same y-axis in all panels
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")

# ggsave("/hpc/users/hoangd02/www/plots/lbp/top_pairs_violin_box_by_percentile.pdf",
#       p, width = 9, height = 8)

########### T-TEST, WILCOXON'S TEST AND PLOTTING 

# dt must have: value, swap ("No swap (0%)"/"All swap (100%)"), percentile (e.g., "99.95%","99.995%")
# If dt is a data.frame, convert for convenience:
dt <- as.data.table(dt)

# --- per-percentile tests ---
p_fmt <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-16, "<1e-16", sprintf("p=%.1e", p)))

# per-percentile: t-test + Cohen's d + y-position for the bar
tests <- dt[, {
  x <- value[swap == "No swap (0%)"]
  y <- value[swap == "All swap (100%)"]

  # t-test p-value
  p_t <- tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
  p_w <- tryCatch(wilcox.test(value ~ swap, exact = FALSE)$p.value, error = function(e) NA_real_)

  # Cohen's d (pooled SD)
  n1 <- length(x); n2 <- length(y)
  s1 <- sd(x);     s2 <- sd(y)
  m1 <- mean(x);   m2 <- mean(y)
  sp <- sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  d  <- (m1 - m2) / sp

  y_max <- max(value, na.rm = TRUE)

  .(p_t = p_t, d = d, p_w = p_w, y_max = y_max)
}, by = percentile]

# positions for two lines (cap so they stay on the 0–1 axis)
tests[, `:=`(
  xmin = 1, xmax = 2,
  y_t  = pmin(y_max + 0.06, 0.99),
  ann_t = paste0("t-test ", p_fmt(p_t), "; d=", sprintf("%.3f", d)),
  ann_w = paste0("Wilcoxon ", p_fmt(p_w))
)]


######## 
######## Or do it individually
# Example for 99.95% only
dt_9995 <- subset(dt, percentile == "99.95%")

ttest <- t.test(value ~ swap, data = dt_9995); ttest
############### no residID ###############
# data:  value by swap
# t = 39.904, df = 449433, p-value < 2.2e-16
# alternative hypothesis: true difference in means between group No swap (0%) and group All swap (100%) is not equal to 0
# 95 percent confidence interval:
#  0.002757300 0.003042155
# sample estimates:
#    mean in group No swap (0%) mean in group All swap (100%) 
#                     0.2514516                     0.2485518 

############### with residID ###############
#  mean in group No swap (0%) mean in group All swap (100%) 
#                   0.2468403                     0.2425398 


wilcoxon <- wilcox.test(value ~ swap, data = dt_9995, exact = FALSE); wilcoxon
# W = 2.9049e+10, p-value < 2.2e-16

# Example for 99.995%
dt_99995 <- subset(dt, percentile == "99.995%")

t.test(value ~ swap, data = dt_99995)
wilcox.test(value ~ swap, data = dt_99995, exact = FALSE)

######## 
######## 
k_map <- c("90%"="44,945,838", "99%"="4,494,584", "99.95%"="224,734", "99.995%"="22,473")
lab_fun <- labeller(percentile = function(x) paste0(x, " (top ", k_map[x], " gene pairs)"))
dt$swap <- factor(dt$swap, levels = c("No swap (0%)", "All swap (100%)"))


with_residID <- ggplot(dt, aes(x = swap, y = value, fill = swap)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  facet_wrap(~ percentile, nrow = 1, scales = "fixed", labeller = lab_fun) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL, y = "Absolute Spearman correlation of all-gene concordance",
    title = "Top blood brain gene pair correlations by percentile (with_residID)",
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  # t-test line
  geom_signif(
    data = tests,
    aes(xmin = xmin, xmax = xmax, annotations = ann_t, y_position = y_t),
    manual = TRUE, tip_length = 0.01, textsize = 3,
    inherit.aes = FALSE)
  # ) +
  # # Wilcoxon line (offset above)
  # geom_signif(
  #   data = tests,
  #   aes(xmin = xmin, xmax = xmax, annotations = ann_w, y_position = y_w),
  #   manual = TRUE, tip_length = 0.01, textsize = 3,
  #   inherit.aes = FALSE
  # )

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_pairs_violin_box_by_percentile_with_residID.pdf",
      with_residID, width = 9, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/top_pairs_violin_box_by_percentile_with_residID.pdf

no_residID <- ggplot(dt, aes(x = swap, y = value, fill = swap)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  facet_wrap(~ percentile, nrow = 1, scales = "fixed", labeller = lab_fun) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL, y = "Absolute Spearman correlation of all-gene concordance",
    title = "Top blood brain gene pair correlations by percentile (no_residID)",
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  # t-test line
  geom_signif(
    data = tests,
    aes(xmin = xmin, xmax = xmax, annotations = ann_t, y_position = y_t),
    manual = TRUE, tip_length = 0.01, textsize = 3,
    inherit.aes = FALSE)
  # ) +
  # # Wilcoxon line (offset above)
  # geom_signif(
  #   data = tests,
  #   aes(xmin = xmin, xmax = xmax, annotations = ann_w, y_position = y_w),
  #   manual = TRUE, tip_length = 0.01, textsize = 3,
  #   inherit.aes = FALSE
  # )

ggsave("/hpc/users/hoangd02/www/plots/lbp/test_no_residID.pdf",
      no_residID, width = 9, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/top_pairs_violin_box_by_percentile.pdf

library(patchwork)
combo <- no_residID + with_residID

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_pairs_violin_box_by_percentile_non_absolute_values.pdf",
      combo, width =18, height = 8)

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_pairs_violin_box_by_percentile.pdf",
      combo, width =18, height = 8)
# Mean difference: practically tiny.
# Cohen’s d (0.27): standardized effect size says it’s small but non-negligible.
# t-test: with n this big, any consistent difference looks “statistically significant.”

#For interpretation, I’d rely more on Cohen’s d or even visual inspection of distributions than on the p-value here.

### hmmmmmmmmmmmm
dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)
#   0   1 
# 152  73 

## 2 timepoints
table(dictionary$number_of_pair_brain)
#  1_pair 2_pairs 
#      79     146 

two_pair <- filter(dictionary, number_of_pair_brain == "2_pairs")
table(two_pair$timepoint)
#  0  1 
# 73 73 

#################################
############# NO ID #############
#################################
#blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
#rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL

# brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
# row.names(brain_full_no_residID)<- brain_full_no_residID$V1; brain_full_no_residID$V1 <- NULL

# blood_form = blood_form5_no_residID
# brain_form = brain_full_no_residID

#################################
############# WITH ID ###########
#################################
blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full) <- brain_full$V1; brain_full$V1 <- NULL

blood_form = blood_form5
brain_form = brain_full

# Assuming:
# blood_form5 = expression matrix (21046 x 225)
# metadata = dataframe with SAMPLE_ISMMS_blood and timepoint columns

# Make sure colnames(blood_form5) match SAMPLE_ISMMS_blood values
head(colnames(blood_form5))
head(two_pair$SAMPLE_ISMMS_blood)

# Subset for timepoint == 1
samples_t1 <- two_pair$SAMPLE_ISMMS_blood[two_pair$timepoint == 1]
blood_t1 <- blood_form5[, colnames(blood_form5) %in% samples_t1, drop = FALSE]

# Subset for timepoint == 0
samples_t0 <- two_pair$SAMPLE_ISMMS_blood[two_pair$timepoint == 0]
blood_t0 <- blood_form5[, colnames(blood_form5) %in% samples_t0, drop = FALSE]

# Check dimensions
dim(blood_t1) #21046    73
dim(blood_t0) #21046    73

# Subset for timepoint == 1
samples_brain_t1 <- two_pair$SAMPLE_ISMMS_brain[two_pair$timepoint == 1]
brain_t1 <- brain_full[, colnames(brain_full) %in% samples_brain_t1, drop = FALSE]

# Subset for timepoint == 0
samples_brain_t0 <- two_pair$SAMPLE_ISMMS_brain[two_pair$timepoint == 0]
brain_t0 <- brain_full[, colnames(brain_full) %in% samples_brain_t0, drop = FALSE]

# Check dimensions
dim(brain_t1) #21356    73
dim(brain_t0) #21356    73

timepoint0_cor <- cor(t(blood_t0), t(brain_t0), method = "spearman")
timepoint0_cor_abs  <- abs(timepoint0_cor)
dim(timepoint0_cor_abs) #21046 21356

timepoint1_cor <- cor(t(blood_t1), t(brain_t1), method = "spearman")
timepoint1_cor_abs  <- abs(timepoint1_cor)
dim(timepoint1_cor_abs) #21046 21356


quantile_tp0 <- quantile(as.vector(timepoint0_cor_abs), probs = c(0.9, 0.95, 0.9995, 0.99995)); quantile_tp0
#       90%       95%    99.95%   99.995% 
# 0.1927064 0.2288350 0.3953783 0.4540602 
quantile_tp1 <- quantile(as.vector(timepoint1_cor_abs), probs = c(0.9, 0.95, 0.9995, 0.99995)); quantile_tp1
#       90%       95%    99.95%   99.995% 
# 0.2011292 0.2384611 0.4092929 0.4693632 

#### Overlap of top correlated pairs
# Flatten matrices into vectors
vals0 <- as.vector(timepoint0_cor_abs)
vals1 <- as.vector(timepoint1_cor_abs)

# Find thresholds
cut0 <- quantile(vals0, 0.9995)
cut1 <- quantile(vals1, 0.9995)

# Get index sets
top0 <- which(timepoint0_cor_abs >= cut0, arr.ind = TRUE)
top1 <- which(timepoint1_cor_abs >= cut1, arr.ind = TRUE)

# Overlap
overlap <- intersect(paste(top0[,1], top0[,2]), paste(top1[,1], top1[,2]))
length(overlap) #144 gene-pairs remain highly correlated across time

#########plotting the distributions 
# Flatten correlations into vectors
vals_tp0 <- as.vector(timepoint0_cor_abs)
vals_tp1 <- as.vector(timepoint1_cor_abs)

# Get thresholds for top 0.05% (99.95 percentile)
thr_tp0 <- quantile(vals_tp0, probs = 0.9995, na.rm = TRUE)
thr_tp1 <- quantile(vals_tp1, probs = 0.9995, na.rm = TRUE)

##subset to top pairs
top_tp0 <- vals_tp0[vals_tp0 >= thr_tp0]
top_tp1 <- vals_tp1[vals_tp1 >= thr_tp1]

length(top_tp0)  # number of pairs at timepoint 0 - 224788
length(top_tp1)  # number of pairs at timepoint 1 - 224860

library(ggplot2)

df_plot <- data.frame(
  correlation = c(top_tp0, top_tp1),
  timepoint = rep(c("Timepoint 0", "Timepoint 1"),
                  c(length(top_tp0), length(top_tp1)))
)

p <- ggplot(df_plot, aes(x = timepoint, y = correlation, fill = timepoint)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  labs(y = "Absolute Spearman correlation",
       title = "Top 0.05 percent all-gene concordance") +
  theme_minimal(base_size = 14)

ggsave("/hpc/users/hoangd02/www/plots/lbp/test_timepoint.pdf",
      p, width =8, height = 8)

# Welch t-test
t_test <- t.test(top_tp0, top_tp1)
t_test
#         Welch Two Sample t-test
# data:  top_tp0 and top_tp1
# t = -199.29, df = 449510, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.01469036 -0.01440422
# sample estimates:
# mean of x mean of y 
# 0.4212328 0.4357801 

# Cohen's d
library(effsize)
cohen_d <- cohen.d(top_tp0, top_tp1)
cohen_d
# d estimate: -0.5943946 (medium)
# 95 percent confidence interval:
#      lower      upper 
# -0.6003681 -0.5884211 






# Swap levels in percentages
#swap_levels <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
swap_levels <- c(0, 1)
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

# Loop through each swap level (0%,100%)
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

