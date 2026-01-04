################################################################################################
###### This figure include VPA, p-value distribution and boxplots/violin plots for all the rates 

###################################################################################
########################################### Rates from extract_set5.sh ############
###################################################################################
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(purrr)
library(furrr)
library(patchwork)
library(ggrastr)
library(cowplot)

#OLD: res <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed_corrected.rds") 
#NEW
res <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/all_simulations_averaged_by_seeds.rds")
res<-as.data.frame(res)

res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

dim(res)
#OLD: 1318660 x 30
#NEW: 1318750 x 19

#FPR = FP / (FP + TN) FPR = NA when FP + TN = 0
#This usually happens when the simulation generates no “negative” cases at all, so the false-positive rate is undefined.


### AVERAGE TNR, FPR, TNR, AND FNR ACROSS METHODS
res_rates <- res %>%
  mutate(across(c(TNR, TPR, FNR, FPR),
                ~ ifelse(is.nan(.), NA_real_, .)))  # drop NaN

mean <- res_rates %>%
  group_by(method) %>%
  summarise(across(c(TNR, TPR, FNR, FPR), ~ mean(., na.rm = TRUE)),
            .groups = "drop")

mean
#   method      TNR   TPR   FNR     FPR
#   <chr>     <dbl> <dbl> <dbl>   <dbl>
# 1 be        0.568 0.871 0.129 0.432  
# 2 known_nsv 0.586 0.889 0.111 0.414  
# 3 leek      0.585 0.865 0.135 0.415  
# 4 none      0.996 0.412 0.588 0.00403
# 5 oracle    0.992 0.862 0.138 0.00763

sva_methods <- res[res$method %in% c("known_nsv", "be", "leek"), ]
dim(sva_methods) #791250     19

sum(is.na(sva_methods$FPR)) #there are no NAs
# Count how many simulations meet the criterion
count_below <- sum(sva_methods$FPR < 0.05, na.rm = TRUE)

# Total number of simulations
total <- sum(!is.na(sva_methods$FPR)) #791250

# Percentage
percent_below <- (count_below / total) * 100; percent_below #39.11254


# res$Feature <- NULL
# res$true_degs<-res$param.fraction_degs * 20000
# res$all_positive=res$TP + res$FP
# res$all_negative=res$TN + res$FN
# res$accuracy=(res$TP + res$TN)/(res$TP + res$TN + res$FP + res$FN)
# res$precision=res$TP/(res$TP + res$FP) 
# res$specificity=res$TN/(res$TN + res$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
# #res$falsePositiveRate=res$FP/(res$FP + res$TN) #type 1 error = false positive #SVA commits this!
# res$recall=res$TP/(res$TP + res$FN) #true positive rate; measure of a test's ability to correctly identify true positives
# res$negativePredictiveValue=res$TN/(res$TN + res$FN)
# res$totalDEGsignal=res$param.fraction_degs*res$param.status_sd

# ## change nsv for oracle and known from 1 to 0
# res$nsv <- mapply(function(nsv_val, method_val) {
#   if (method_val %in% c("known_nsv", "oracle") && length(nsv_val) == 1 && nsv_val == 1) {
#     return(0)
#   } else {
#     return(nsv_val)
#   }
# }, res$nsv, res$method, SIMPLIFY = FALSE)


str(res)

res$seeds <- vapply(res$seeds, function(x) paste0('"', paste(x, collapse = ","), '"'), FUN.VALUE = character(1))
res$Ratio_DE_non_DE <- vapply(res$Ratio_DE_non_DE, function(x) paste0('"', paste(x, collapse = ","), '"'), FUN.VALUE = character(1))
res$nsv <- vapply(res$nsv, function(x) paste0('"', paste(x, collapse = ","), '"'), FUN.VALUE = character(1))

write.csv(res, "/sc/arion/projects/mscic1/results/jolie/sim_bug/res_91percentsimulations_supplementary_figure.csv", row.names = FALSE)

res_read <- read.csv("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_91percentsimulations_supplementary_figure.csv", stringsAsFactors = FALSE)

###only run this once
res$method <- factor(res$method, 
	levels = c("oracle", "none", "known_nsv", "be", "leek"),
	labels = c('Oracle', 'No SVs','Known Number\nof SVs','SVA "BE"','SVA "Leek"'))

custom_color <- c("Oracle" = "#d8a97d",
				  "No SVs" = "#3466ff",    
  				  "Known Number\nof SVs" = "#af6ca6",
  				  'SVA "BE"' = "#a84c4c",
   				  'SVA "Leek"' = "#b4747e")

## colored version
p1 <- ggplot(res, aes(x = method, y = TPR, fill = method)) +  
  geom_violin(aes(fill = method), trim = FALSE, color = "black") + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), 
               color = "black", outlier.shape = NA) +  
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) + 
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  labs(x = "Method", y = "True Positive Rate") + 
  scale_fill_manual(name = "DEA Approaches", values = custom_color)

## colored version
p2 <- ggplot(res, aes(x = method, y = TNR, fill = method)) +
  geom_violin(aes(fill = method), trim = FALSE, color = "black") + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), #can add alpha
               color = "black", outlier.shape = NA) + 
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  scale_y_continuous(limits = c(0, 1),
    				 breaks = seq(0, 1, by = 0.25)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(x = "Method", y = "True Negative Rate") + 
  scale_fill_manual(name = "DEA Approaches", values = custom_color)


## colored version, x axis has no legend
p3 <- ggplot(res, aes(x = method, y = FPR, fill = method)) +
  geom_violin(aes(fill = method), trim = FALSE, color = "black") + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9),
               color = "black", outlier.shape = NA) +  
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  scale_y_continuous(limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Method", y = "False Positive Rate") + 
  scale_fill_manual(name = "DEA Approaches", values = custom_color)

## colored version, no x-axis
p4 <- ggplot(res, aes(x = method, y = FNR, fill = method)) +
  geom_violin(aes(fill = method), trim = FALSE, color = "black") + 
  #geom_violin(trim = FALSE, color = "black", fill= NA) +  
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9),
               color = "black", outlier.shape = NA) + 
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
    scale_y_continuous(limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none", 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(x = "Method", y = "False Negative Rate") + 
  scale_fill_manual(name = "DEA Approaches", values = custom_color)

# Combine the plots
# combo <- p1 + p2 + p3 + p4 

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_05142025.pdf", plot = combo, width = 8, height = 5) 

# combo <- (p1 + p2 + p3 + p4 + plot_layout(guides = "collect")) /
# 		  guide_area() + plot_layout(heights = c(5, 1)) &
#   		  theme(legend.position = "top") 

combo <- (p1 + p2 + p3 + p4 + plot_layout(guides = "collect")) &
		  theme(legend.position = "top")  # shared legend at top 

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_09122025.pdf", plot = combo, width = 8, height = 5) 

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_05142025_two.pdf", plot = combo, width = 8, height = 5) 

###################################################################################
######### P VALUE DISTRIBUTION from simulation_pvalue_distributions.sh ############
###################################################################################
case0_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case0_pval_out.csv")
case1_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case1_pval_out.csv")
case2_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case2_pval_out.csv")
case3_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case3_pval_out.csv")
case4_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case4_pval_out.csv")

### KNOWN, ORACLE AND NONE
case0_sub <- filter(case0_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case1_sub <- filter(case1_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case2_sub <- filter(case2_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case3_sub <- filter(case3_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case4_sub <- filter(case4_pval_out, Method %in% c("known_nsv", "none", "oracle"))

# Filter out rows where P.Value is NA or zero to avoid log issues
filtered_data <- case0_sub %>% filter(!is.na(P.Value) & P.Value > 0)
filtered_data$Method <- factor(filtered_data$Method, 
	levels = c("oracle", "none", "known_nsv"))

# Plot histogram colored by Method
case0 <- ggplot(filtered_data, aes(x = P.Value, fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", position = "identity") +
  facet_wrap(~Method, ncol = 1) +
  xlim(0, 1) +  
  labs(x = "p-values",y = "Count") +
  theme_minimal() +
  theme(legend.position = "none",
  		strip.text = element_blank(),
    	axis.text.x = element_text(angle = 0, hjust = 1), 
    	axis.ticks.x = element_blank(),
  		plot.margin = margin(5.5, 5.5, 5.5, 5.5)) + 
  scale_fill_manual(values = c("known_nsv" = "#af6ca6", 
                               "none" = "#3466ff", 
                               "oracle" = "#d8a97d")) + 
  guides(fill = "none")

# layout <- c(
#   area(1, 2, 1, 2),  # case0 at row 1, col 2
#   area(2, 1, 2, 1),  # p1 at row 2, col 1
#   area(3, 1, 3, 1),  # p3 at row 3, col 1
#   area(2, 2, 2, 2),  # p2 at row 2, col 2
#   area(3, 2, 3, 2)   # p4 at row 3, col 2
# )

# combo <- case0 + p1 + p3 + p2 + p4 +
#   plot_layout(design = layout, guides = "collect") &
#   theme(legend.position = "bottom")
# combo

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_05152025.pdf", plot = combo, width = 8, height = 8) 

https://hoangd02.u.hpc.mssm.edu/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_05152025.pdf

###################################################################################
################ VPA from sva_sim_parallel_for_vp.R in SVA --> Scripts ############
###################################################################################
https://hoangd02.u.hpc.mssm.edu/plots/sva/sim_results_oct2024/sim_vpa/

##this is the final one 
https://hoangd02.u.hpc.mssm.edu/plots/sva/sim_results_oct2024/sim_vpa/vp_plot_row_9.pdf

#load("/sc/arion/projects/mscic1/results/jolie/simulation_vp_plots_results.RData")
vp_plots <- readRDS("/sc/arion/projects/mscic1/results/jolie/simulation_vp_plots_results_20251015.rds")

length(vp_plots) #10

##main figure: vp_plots[[9]]

main_vp_plot <- vp_plots[[9]] +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1), 
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  labs(y = "Variance explained (%)", title = NULL) +
  scale_x_discrete(labels = c(
    "Batch_ID" = "Batch ID",
    "Status" = "Status",
    "SV01" = "EIV 1",
    "SV02" = "EIV 2",
    "SV03" = "EIV 3",
    "SV04" = "EIV 4",
    "SV05" = "EIV 5",
    "SV06" = "EIV 6",
    "SV07" = "EIV 7",
    "SV08" = "EIV 8",
    "SV09" = "EIV 9",
    "SV10" = "EIV 10",
    "Residuals" = "Residuals"
  )) +
  guides(fill = "none") 

layout <- c(
  area(1, 1, 1, 1),  # main_vp_plot at row 1, col 1
  area(1, 2, 1, 2),  # case0 at row 1, col 2
  area(2, 1, 2, 1),  # p1 at row 2, col 1
  area(3, 1, 3, 1),  # p3 at row 3, col 1
  area(2, 2, 2, 2),  # p2 at row 2, col 2
  area(3, 2, 3, 2)   # p4 at row 3, col 2
)

combo <- main_vp_plot + case0_shrunk+ p1 + p3 + p2 + p4 +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = "/hpc/users/hoangd02/www/plots/Figure1_20250518_2.pdf",
  plot = combo, width = 9, height = 8) 

https://hoangd02.u.hpc.mssm.edu/plots/Figure1_20250518_2.pdf

### Check if the 5 representative cases overfit:
test <- read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_results_with_summary_vpa_stats.csv")

Case 1: n_sv =10, total_sv_sd = 5, sv_sd = 1.58, fraction_degs = 0.30, status_sd = 0.2, resid_sd = 0.2, fraction_sv = 0.1, seed = 12352; 
Case 2: n_sv =10, total_sv_sd = 7.5, sv_sd = 2.37, fraction_degs = 0.08, status_sd = 0.5, resid_sd = 0.3, fraction_sv = 0.5, seed = 12351; 
Case 3: n_sv =10, total_sv_sd = 1.5, sv_sd = 0.47, fraction_degs = 0.20, status_sd = 0.2, resid_sd = 0.4, fraction_sv = 0.1, seed = 12348; 
Case 4: n_sv =10, total_sv_sd = 2.0, sv_sd = 0.63, fraction_degs = 0.30, status_sd = 0.2, resid_sd =0.3, fraction_sv = 0.7, seed = 12348

test <- test %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

## MAIN CASE 
n_sv =10, total_sv_sd = 2.5, sv_sd = 0.79, fraction_degs = 0.04, status_sd = 0.1, resid_sd = 0.4, fraction_sv = 0.3, seed = 12345 (row 9 in 10_random_vpa_set photo) 
##
sub <- test[
  test$param.n_sv == 10 &
  test$param.fraction_degs == 0.04 &
  test$param.fraction_sv == 0.3 &
  test$param.status_sd == 0.1 &
  test$param.sv_sd == 0.790569415042095 &
  test$param.resid_sd == 0.4 &
  test$param.seed == 12345, ]
sub
## 3 cases where FPR = 0.001041667 = 0.1041667% 

### Case 1:
Case 1: n_sv =10, total_sv_sd = 5, sv_sd = 1.58, fraction_degs = 0.30, status_sd = 0.2, resid_sd = 0.2, fraction_sv = 0.1, seed = 12352; 
sub <- test[
  test$param.n_sv == 10 &
  test$param.fraction_degs == 0.30 &
  test$param.fraction_sv == 0.1 &
  test$param.status_sd == 0.2 &
  test$param.sv_sd == 1.58113883008419 &
  test$param.resid_sd == 0.2 &
  test$param.seed == 12352, ]
## 3 cases where FPR = 0.23450000 = 23% !


### Case 2:
Case 2: n_sv =10, total_sv_sd = 7.5, sv_sd = 2.37, fraction_degs = 0.08, status_sd = 0.5, resid_sd = 0.3, fraction_sv = 0.5, seed = 12351; 
sub <- test[
  test$param.n_sv == 10 &
  test$param.fraction_degs == 0.08 &
  test$param.fraction_sv == 0.5 &
  test$param.status_sd == 0.5 &
  test$param.sv_sd == 2.37170824512628 &
  test$param.resid_sd == 0.3 &
  test$param.seed == 12351, ]

## 3 cases where FPR = 0.927228261 = 93% !

### Case 3: 
Case 3: n_sv =10, total_sv_sd = 1.5, sv_sd = 0.47, fraction_degs = 0.20, status_sd = 0.2, resid_sd = 0.4, fraction_sv = 0.1, seed = 12348; 
sub <- test[
  test$param.n_sv == 10 &
  test$param.fraction_degs == 0.20 &
  test$param.fraction_sv == 0.1 &
  test$param.status_sd == 0.2 &
  test$param.sv_sd == 0.474341649025257 &
  test$param.resid_sd == 0.4 &
  test$param.seed == 12348, ]
table(sub$param.sv_sd)

## 3 cases where FPR for be = 0.0085000 = 0.85% and leek = 0.0053125 = 0.53%

### Case 4: 
Case 4: n_sv =10, total_sv_sd = 2.0, sv_sd = 0.63, fraction_degs = 0.30, status_sd = 0.2, resid_sd =0.3, fraction_sv = 0.7, seed = 12348
sub <- test[
  test$param.n_sv == 10 &
  test$param.fraction_degs == 0.30 &
  test$param.fraction_sv == 0.7 &
  test$param.status_sd == 0.2 &
  test$param.sv_sd == 0.632455532033676 &
  test$param.resid_sd == 0.3 &
  test$param.seed == 12348, ]
sub

## 3 cases where FPR = 0.026714286 = 2.671429 %







