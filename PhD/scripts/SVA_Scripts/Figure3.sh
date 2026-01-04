###################################################################################
######## Lollipop Plots from status_var_and _rates_lollipop plot.sh ##################
###################################################################################
#https://hoangd02.u.hpc.mssm.edu/plots/sva/gandal_sim_results/FPR_by_cor_combined_simData1_and_simData2_status_sd0.1.pdf

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrastr)
library(rlang)
library(stringr)
library(fs)
library(withr)
library(assertthat)
library(plyr)
library(dplyr)
library(S4Vectors)
library(variancePartition)
library(glmnet)
library(BiocParallel)
library(rctutils)
library(magrittr)
library(limma)
library(data.table)
library(sva)

res <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed_corrected.rds") 
dim(res) #1,318,660      18
res<-as.data.frame(res)

# this is for grid summary
# oracle <-filter(res, method == "oracle")
# breaks <- c(0, 0.5, 1, 1.5, 2, 3.2)  # Choose your own cut points
# labels <- c("Very Low", "Low", "Medium", "High", "Very High")

# # Cut sv_sd into categories
# oracle$sv_sd_bucket <- cut(oracle$param.sv_sd, breaks = breaks, labels = labels, include.lowest = TRUE)

# # Count frequency in each category
# table(oracle$sv_sd_bucket)
#  Very Low       Low    Medium      High Very High 
#     97509     86223     50000     16250     13750

library(dplyr)
res <- res %>%
  select(-Status_mean, -Status_sd, -Status_median, -Feature)

res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP),   # False Negative Rate: FN / (FN + TP) truth: all positive
    FDR = FP / (TP + FP),   # False Discovery Rate
    Precision = TP / (TP + FP),
    param.total_sv_sd = sqrt(param.sv_sd^2 * param.n_sv),
    param.n_batches = 15,
    param.batch_sd = 0.3,
    ## For status and batch, we multiply by (n-1)/n where n is the
    ## number of levels in the variable
    param.mean_batch_var = param.batch_sd^2 * (param.n_batches - 1) / param.n_batches,
    param.mean_sv_var = param.total_sv_sd^2 * param.fraction_sv,
    param.mean_status_var = param.status_sd^2 * param.fraction_degs * 1/2,
    param.mean_resid_var = param.resid_sd^2,
    param.total_mean_var = param.mean_batch_var + param.mean_sv_var + param.mean_status_var + param.mean_resid_var,
    param.batch_var_fraction = param.mean_batch_var / param.total_mean_var,
    param.sv_var_fraction = param.mean_sv_var / param.total_mean_var,
    param.status_var_fraction = param.mean_status_var / param.total_mean_var,
    param.resid_var_fraction = param.mean_resid_var / param.total_mean_var
  )

res$param.total_sv_sd_rounded <- round(res$param.total_sv_sd, 1)
table(res$param.total_sv_sd_rounded)
#     0   0.5     1   1.5     2   2.5     3   3.5     4   4.5     5   5.5     6 
# 68755 62505 62510 62505 62505 62510 62505 62505 62500 62500 62500 62500 62500 
#   6.5     7   7.5     8   8.5     9   9.5    10 
# 62500 62505 62500 62500 62500 62500 62500 62355 

##SUBSETTING FOR FPR >=0.05 - P-value = 0.05 means we're okay with a 5% FPR per test 
sub <- filter(res, method %in% c("oracle", "none", "known_nsv"))
dim(sub) #791196     34

sub <- sub %>% mutate(FPR_level = case_when(
    FPR >= 0.05 ~ "high",
    FPR < 0.05  ~ "low")) %>% mutate(TPR_level = case_when(
    TPR >= 0.50 ~ "high",
    TPR < 0.50  ~ "low")) %>% mutate(Precision_level = case_when(
    Precision >= 0.75 ~ "high",
    Precision < 0.75  ~ "low")) %>% mutate(FDR_level = case_when(
    FDR >= 0.05 ~ "high", FDR < 0.05  ~ "low"))

table(sub$FDR_level) #total: 301174+ 441943 = 743117
#   high    low 
# 301174 441943
#301174/791196 = 0.3806566

table(sub$FPR_level) #total: 157286+ 633910 = 791196 / 791,196 
#   high    low 
# 157286 633910 
#157286/791196 = 0.1987952

table(sub$TPR_level) #total: 143304+568764 = 712068
#   high    low 
# 568764 143304
#568764/791196 = 0.7188661

##we have all FPR info; 712068 TPR info (79128 NAs); 743117 Precision info (48079 NAs)
known <- filter(sub, method == "known_nsv")

table(known$FPR_level) 
#   high    low 
# 157286 106446
sum(known$FPR_level=="high")/length(known$FPR_level) #0.5963857
 
table(known$FDR_level) 
#   high    low 
# 195727  63015  #195727/(195727+63015) #0.7564562

sum(known$FDR_level=="high")/length(known$FDR_level) ## NAs?

table(known$TPR_level)
#   high    low 
# 234477   2879 | 234477/(234477 + 2879 ) #0.9878705

#table(sub$FPR[sub$method=="oracle"]<0.025)

sub$nsv[sub$nsv == 1] <- 0
#sub <- sub[sub$param.sv_sd != 0, ]

# Step 1: Separate 0 values and non-zero values
sub_zero <- sub[sub$param.status_var_fraction == 0, ]
sub_nonzero <- sub[sub$param.status_var_fraction != 0, ]

# Step 2: Create 9 quantile bins for non-zero values
q_breaks <- unique(quantile(
  sub_nonzero$param.status_var_fraction,
  probs = seq(0, 1, length.out = 11),  # 10 bins = 11 cut points
  na.rm = TRUE
))

# Label format: (a%, b%]
percent_labels <- paste0(
  "(", round(100 * head(q_breaks, -1), 3), "%, ",
      round(100 * tail(q_breaks, -1), 3), "%]"
)

# Step 3: Bin the non-zero subset
sub_nonzero$status_var_percent_bin <- cut(
  sub_nonzero$param.status_var_fraction,
  breaks = q_breaks,
  include.lowest = TRUE,
  labels = percent_labels
)

# Step 4: Assign "0%" label to zero values
sub_zero$status_var_percent_bin <- "0%"

# Step 5: Combine the two subsets
sub_binned <- rbind(sub_zero, sub_nonzero)

bin_levels <- c("0%", percent_labels)

sub_binned$status_var_percent_bin <- factor(sub_binned$status_var_percent_bin, levels = bin_levels)

sub_binned$method <- factor(sub_binned$method, levels = c("oracle", "none", "known_nsv"))
dim(sub_binned) #791196     38

#### ONLY FOCUSING ON THE NON-ZERO NOW 
sub_binned <- filter(sub_binned, status_var_percent_bin != "0%")
dim(sub_binned) #712068     38

#sub_binned$method <- factor(sub_binned$method, levels = c("oracle","none","known_nsv"))

sub_binned <- sub_binned %>%
  mutate(method = case_when(
    method == "oracle" ~ "Oracle",
    method == "none" ~ "No SVs",
    method == "known_nsv" ~ "Known Number of SVs",
    ))

sub_binned$method <- factor(sub_binned$method, levels = c("Oracle","No SVs","Known Number of SVs"))

#looking to see values of TPR, FPR, FDR at real signal ranges 
table(sub_binned$status_var_percent_bin)
 #               0%      (0%, 0.003%]  (0.003%, 0.009%]  (0.009%, 0.018%] 
 #                0             71208             71208             71205 
 # (0.018%, 0.032%]  (0.032%, 0.058%]  (0.058%, 0.106%]  (0.106%, 0.206%] 
 #            71208             71205             71214             71199 
 # (0.206%, 0.463%]  (0.463%, 1.429%] (1.429%, 39.936%] 
 #            71217             71202             71202 

lower_range <- filter(sub_binned, status_var_percent_bin == "(0.032%, 0.058%]")
lower_range_known <- filter(lower_range, method == "Known Number of SVs")
lower_range_known_TPR <- filter(lower_range_known, TPR_level == "high")
summary(lower_range_known_TPR$TPR)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5001  0.8644  0.9255  0.8934  0.9588  0.9945 

higher_range <- filter(sub_binned, status_var_percent_bin == "(0.463%, 1.429%]")


### LOW FPR/FDR = CONTROLLED --- GOOD
### HIGH FPR/FDR = UNCONTROLLED -- BAD

### FDR HIGH
fdr_summary <- sub_binned %>%
  group_by(status_var_percent_bin, method) %>%
  dplyr::summarise( 
    total = n(),
    n_high = sum(FDR_level == "high", na.rm = TRUE),
    frac_high = n_high / total
  ) %>%
  ungroup()

df_fdr <- as.data.frame(fdr_summary)
#    status_var_percent_bin              method total n_high    frac_high
# 3            (0%, 0.003%] Known Number of SVs 23736  23543 0.99186889
# 15       (0.032%, 0.058%] Known Number of SVs 23735  20798 0.87625869
# 27       (0.463%, 1.429%] Known Number of SVs 23734   9274 0.39074745
# 30      (1.429%, 39.936%] Known Number of SVs 23734   2460 0.10364877


x_levels <- unique(fdr_summary$status_var_percent_bin)
fdr_summary$status_var_percent_bin <- factor(fdr_summary$status_var_percent_bin, levels = x_levels)
dodge_width <- 0.6
methods <- unique(fdr_summary$method)

# Lollipop plot for High FDR control 
high_fdr <- ggplot(fdr_summary, aes(x = status_var_percent_bin, y = frac_high, color = method)) +
  geom_segment(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods), 
        xend = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods),
        y = 0, yend = frac_high),
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods)),
    size = 3
  ) +
  scale_x_continuous(
    breaks = 1:length(x_levels),
    labels = x_levels
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(values = c("Oracle" = "#d8a97d", "No SVs" = "#3466ff", "Known Number of SVs" = "#83376c"), labels = c("Oracle" = "Oracle", "No SVs" = "No SVs", "Known Number of SVs" = "Known Number\nof SVs")) +
  labs(
    x = "Total Variation Explained by Phenotype of Interest (Deciles)",
    y = "Simulations with Uncontrolled FDR",
    color = "Method",
    #title = " "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

##HIGH FPR, When it is NOT controlled for
fpr_summary <- sub_binned %>%
  group_by(status_var_percent_bin, method) %>%
  dplyr::summarise( 
    total = n(),
    n_high = sum(FPR_level == "high", na.rm = TRUE),
    frac_high = n_high / total
  ) %>%
  ungroup()

df_fpr <- as.data.frame(fpr_summary)
#   status_var_percent_bin              method total n_high  frac_high
#3            (0%, 0.003%] Known Number of SVs 23736  21952 0.92483991
#15       (0.032%, 0.058%] Known Number of SVs 23735  16984 0.71556773
#27       (0.463%, 1.429%] Known Number of SVs 23734   6043 0.25461363
#30      (1.429%, 39.936%] Known Number of SVs 23734    924 0.03893149


x_levels <- unique(fpr_summary$status_var_percent_bin)
fpr_summary$status_var_percent_bin <- factor(fpr_summary$status_var_percent_bin, levels = x_levels)
dodge_width <- 0.6
methods <- unique(fpr_summary$method)

high_fpr <- ggplot(fpr_summary, aes(x = status_var_percent_bin, y = frac_high, color = method)) +
  geom_segment(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods), 
        xend = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods),
        y = 0, yend = frac_high),
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods)),
    size = 3
  ) +
  scale_x_continuous(
    breaks = 1:length(x_levels),
    labels = x_levels
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(values = c("Oracle" = "#d8a97d", "No SVs" = "#3466ff", "Known Number of SVs" = "#83376c"), labels = c("Oracle" = "Oracle", "No SVs" = "No SVs", "Known Number of SVs" = "Known Number\nof SVs")) +
  labs(
    x = "Total Variation Explained by Phenotype of Interest (Deciles)",
    y = "Simulations with Uncontrolled FPR",
    color = "Method",
    #title = " "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")


#### HIGH TPR
tpr_summary <- sub_binned %>%
  group_by(status_var_percent_bin, method) %>%
  dplyr::summarise( 
    total = n(),
    n_high = sum(TPR_level == "high", na.rm = TRUE),
    frac_high = n_high / total
  ) %>%
  ungroup()
  
df_tpr <- as.data.frame(tpr_summary)
#    status_var_percent_bin              method total n_high    frac_high
# 3            (0%, 0.003%] Known Number of SVs 23736  22963 0.9674334344
# 15       (0.032%, 0.058%] Known Number of SVs 23735  23488 0.9895934274
# 27       (0.463%, 1.429%] Known Number of SVs 23734  23734 1.0000000000
# 30      (1.429%, 39.936%] Known Number of SVs 23734  23734 1.0000000


x_levels <- unique(tpr_summary$status_var_percent_bin)
tpr_summary$status_var_percent_bin <- factor(tpr_summary$status_var_percent_bin, levels = x_levels)
dodge_width <- 0.6
methods <- unique(tpr_summary$method)

high_tpr <- ggplot(tpr_summary, aes(x = status_var_percent_bin, y = frac_high, color = method)) +
  geom_segment(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods), 
        xend = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods),
        y = 0, yend = frac_high),
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = as.numeric(status_var_percent_bin) + 
          (as.numeric(factor(method, levels = methods)) - 2) * dodge_width / length(methods)),
    size = 3
  ) +
  scale_x_continuous(
    breaks = 1:length(x_levels),
    labels = x_levels
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(values = c("Oracle" = "#d8a97d", "No SVs" = "#3466ff", "Known Number of SVs" = "#83376c"), labels = c("Oracle" = "Oracle", "No SVs" = "No SVs", "Known Number of SVs" = "Known Number\nof SVs")) +
  labs(
    x = "Total Variation Explained by Phenotype of Interest (Deciles)",
    y = "Simulations with High TPR",
    color = "Method",
    #title = " "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),legend.position = "none")

library(patchwork)
combo <- high_fdr / high_fpr / high_tpr 

ggsave(filename = "/hpc/users/hoangd02/www/plots/FDR_TPR_FPR_Lollipop_05132025.pdf", plot = combo, width=10, height=11)

#/hpc/users/hoangd02/www/plots/FDR_TPR_FPR_Precision_Percent_Status_Variation_Lollipop.pdf
#/hpc/users/hoangd02/www/plots/FDR_TPR_FPR_Precision_Percent_Status_Variation_Lollipop_new_deciles.pdf

##load VPA from real data ---- KEEP EVERYTHING AS PERCENTAGE
## BIPOLAR DISEASE
bp_cmc_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_variance_partition.txt",data.table=FALSE)
bp_cmc_frac_status_var_mean <-mean(bp_cmc_status_var$Dx); bp_cmc_frac_status_var_mean #0.0006640974

bp_cmc_hbcc_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_variance_partition.txt",data.table=FALSE)
bp_cmc_hbcc_frac_status_var_mean <-mean(bp_cmc_hbcc_status_var$Dx); bp_cmc_hbcc_frac_status_var_mean #0.006344017

## SCHIZOPHRENIA
sz_cmc_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/CMC/sz_mssm_variance_partition.txt",data.table=FALSE)
sz_cmc_frac_status_var_mean <-mean(sz_cmc_status_var$Dx); sz_cmc_frac_status_var_mean #0.0005545046

sz_cmc_hbcc_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/CMC/sz_hbcc_variance_partition.txt",data.table=FALSE)
sz_cmc_hbcc_frac_status_var_mean <-mean(sz_cmc_hbcc_status_var$Dx); sz_cmc_hbcc_frac_status_var_mean #0.006920495

## AD
rosmap_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/rosmap_variance_partition.txt",data.table=FALSE)
rosmap_frac_status_var_mean <-mean(rosmap_status_var$ceradsc); rosmap_frac_status_var_mean  #0.002282867

msbb_status_var<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/msbb_variance_partition.txt",data.table=FALSE)
msbb_frac_status_var_mean <-mean(msbb_status_var$PlaqueMean); msbb_frac_status_var_mean #0.005673065


### WHICH DECILE ARE THESE IN 
#0.058-0.106% BP CMC & SZ CMC
#0.206-0.463% ROSMAP
#0.463-1.429% MSBB & BP CMC-HBCC & SZ CMC-HBCC

# Get unique bin labels in order
bin_labels <- levels(sub_binned$status_var_percent_bin)
bin_labels <- bin_labels[bin_labels != "0%"]
n_bins <- length(bin_labels)

# Create x boundaries so tiles touch without white space
bin_edges <- seq(0.5, n_bins + 0.5, by = 1)

# Create rect data with xmin, xmax for each bin
bin_plot_df <- data.frame(
  bin_label = bin_labels,
  xmin = bin_edges[-length(bin_edges)],
  xmax = bin_edges[-1],
  ymin = 0,
  ymax = 0.2
)

# Make the plot
ggplot() +
  geom_rect(
    data = bin_plot_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey",
    color = NA
  ) +
  # Add vertical black lines at bin borders
  geom_vline(xintercept = bin_edges, color = "black", size = 0.3) +
  scale_x_continuous(
    breaks = (bin_plot_df$xmin + bin_plot_df$xmax) / 2,
    labels = bin_plot_df$bin_label,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal() +
  labs(x = "Percent Status Variation", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank()
  )

##ADDING REAL DATA
# Create the dataset values in percentage (%)
real_data <- data.frame(
  dataset = c("BP: CMC", "BP: CMC-HBCC", "SZ: CMC", "SZ: CMC-HBCC", "AD: ROSMAP", "AD: MSBB"),
  status_var_percent = c(
    bp_cmc_frac_status_var_mean, bp_cmc_hbcc_frac_status_var_mean,
    sz_cmc_frac_status_var_mean, sz_cmc_hbcc_frac_status_var_mean,
    rosmap_frac_status_var_mean, msbb_frac_status_var_mean
  ) * 100  
)

find_bin_label <- function(value, bin_labels) {
  for (label in bin_labels) {
    range_vals <- str_extract_all(label, "\\d+\\.\\d+")[[1]] %>% as.numeric()
    if (length(range_vals) == 2 && !any(is.na(range_vals))) {
      if (value > range_vals[1] && value <= range_vals[2]) {
        return(label)
      }
    }
  }
  return(NA)
}

real_data$bin_label <- sapply(real_data$status_var_percent, find_bin_label, bin_labels = bin_labels)

real_data <- real_data %>%
  left_join(bin_plot_df, by = "bin_label") 

# Extract bin numeric range
bounds <- str_match(real_data$bin_label, "\\((\\d+\\.\\d+)%?, (\\d+\\.\\d+)%?\\]")
real_data$bin_min <- as.numeric(bounds[, 2])
real_data$bin_max <- as.numeric(bounds[, 3])

# Compute x_pos based on actual percentage location within the bin
real_data <- real_data %>%
  mutate(
    x_pos = xmin + (status_var_percent - bin_min) / (bin_max - bin_min) * (xmax - xmin)
  )

bar<-ggplot() +
  # Grey continuous bar
  geom_rect(
    data = bin_plot_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#f2f2f2"
  ) +
  # Bin borders as segments limited to grey height
  geom_segment(
    data = data.frame(x = bin_edges),
    aes(x = x, xend = x, y = 0, yend = 0.2),
    color = "black", size = 0.3
  ) +
  # Real data vertical lines inside the bar
  geom_segment(
    data = real_data,
    aes(x = x_pos, xend = x_pos, y = 0, yend = 0.2, color = dataset),
    size = 1.3
  ) +
  scale_x_continuous(
    breaks = (bin_plot_df$xmin + bin_plot_df$xmax) / 2,
    labels = bin_plot_df$bin_label,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(
  values = c(
    # Alzheimer's Disease (AD) - RED TONES
    "AD: MSBB" = "#E37C72",    # Coral red (light)
    "AD: ROSMAP" = "#A3392F",  # Brick red (dark)

    # Bipolar Disorder (BP) - GREEN TONES
    "BP: CMC" = "#BBD8A3",     # Light green (soft)
    "BP: CMC-HBCC" = "#79AC78",# Dark green

    # Schizophrenia (SZ) - COMPLEMENTARY TONES
    "SZ: CMC" = "#C4A5D6",     # Pale indigo
    "SZ: CMC-HBCC" = "#734F96" # Classic indigo
  ),
    labels = c("AD: MSBB" = "AD: MSBB (0.57%)", 
      "AD: ROSMAP" = "AD: ROSMAP (0.23%)", 
      "BP: CMC" = "BP: CMC (0.07%)",
      "BP: CMC-HBCC" = "BP: CMC-HBCC (0.63%)", 
      "SZ: CMC" = "SZ: CMC (0.06%)", 
      "SZ: CMC-HBCC" = "SZ: CMC-HBCC (0.69%)"),
    name = "Dataset"
  ) +
  theme_minimal() +
  labs(x = "Total Variation Explained by Phenotype of Interest (Deciles)", y = NULL, color = "Dataset") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank()
  ) + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none", #change this to bottom/none
        legend.text = element_text(size = 12),      
        legend.title = element_text(size = 14, face = "bold"))

combo <- wrap_plots(
  high_fdr,
  patchwork::wrap_elements(full = bar),
  high_fpr,
  patchwork::wrap_elements(full = bar),
  high_tpr,
  patchwork::wrap_elements(full = bar),
  ncol = 1,
  heights = c(1, 0.6, 1, 0.6, 1, 0.6)  # adjust height ratios for bars vs plots
)

ggsave(filename = "/hpc/users/hoangd02/www/plots/FDR_FPR_TPR_Percent_Status_Var_08042025.pdf", plot = combo, width=9, height=16)

ggsave(filename = "/hpc/users/hoangd02/www/plots/FDR_FPR_TPR_Percent_Status_Var_05132025_FOUR.pdf", plot = combo, width=9, height=16)








































