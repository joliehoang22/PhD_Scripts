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

res <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed_corrected.rds") 
res<-as.data.frame(res)

res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP),  # False Negative Rate: FN / (FN + TP) truth: all positive
    FDR = FP / (TP + FP)   # FDR: TP / ALL POSITIVE
  )

res$Feature <- NULL
res$true_degs<-res$param.fraction_degs * 20000
res$all_positive=res$TP + res$FP
res$all_negative=res$TN + res$FN
res$accuracy=(res$TP + res$TN)/(res$TP + res$TN + res$FP + res$FN)
res$precision=res$TP/(res$TP + res$FP) 
res$specificity=res$TN/(res$TN + res$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
#res$falsePositiveRate=res$FP/(res$FP + res$TN) #type 1 error = false positive #SVA commits this!
res$recall=res$TP/(res$TP + res$FN) #true positive rate; measure of a test's ability to correctly identify true positives
res$negativePredictiveValue=res$TN/(res$TN + res$FN)
res$totalDEGsignal=res$param.fraction_degs*res$param.status_sd

res$param.total_sv_sd <- sqrt(res$param.n_sv * res$param.sv_sd^2)
res$param.total_sv_var <- res$param.total_sv_sd^2
res$param.status_var <- res$param.status_sd^2
res$param.resid_var <- res$param.resid_sd^2

## change nsv for oracle and known from 1 to 0
res$nsv <- mapply(function(nsv_val, method_val) {
  if (method_val %in% c("known_nsv", "oracle") && length(nsv_val) == 1 && nsv_val == 1) {
    return(0)
  } else {
    return(nsv_val)
  }
}, res$nsv, res$method, SIMPLIFY = FALSE)

dim(res)
#1318660 x 30@

#######################
sva_methods <- res[res$method %in% c("known_nsv", "leek", "be"), ]
dim(sva_methods) #791196     35

#how many have good FPR and FDR
sum(sva_methods$FPR < 0.05, na.rm = TRUE) #309506

#309506 / 791196 = 0.3911875

sum(sva_methods$FDR < 0.05, na.rm = TRUE) #190493

#190493 / 791196 = 0.2407659

########################

##mean of each rate per method 
library(dplyr)

res_summary <- res %>%
  group_by(method) %>%
  summarise(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_TNR = mean(TNR, na.rm = TRUE),
    mean_FNR = mean(FNR, na.rm = TRUE),
    mean_FDR = mean(FDR, na.rm = TRUE)
  )

res_summary

############divergence start
sub <- res[res$method %in% c("be", "known_nsv", "leek"), ]

ten_eiv <- filter(sub, param.n_sv == "10")
summary(ten_eiv$TP)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  0.0   721.2  1558.8  3093.1  5591.6  9942.9 

hundred_eiv <- filter(sub, param.n_sv == "100")
summary(hundred_eiv$TP)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  0.0   692.6  1519.9  2953.9  5327.8  9848.6 

##for known only
sub2 <- filter(res, method == "known_nsv")
ten_eiv <- filter(sub2, param.n_sv == "10")
summary(ten_eiv$TP)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  0.0   732.1  1569.3  3115.5  5636.2  9942.9 

hundred_eiv <- filter(sub2, param.n_sv == "100")
summary(hundred_eiv$TP)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  0.0   708.4  1542.5  3024.3  5465.3  9848.6 

##########
known <- filter(res, method == "known_nsv")

known_filtered <- known[
  known$param.resid_var == 0.25 &
  known$param.status_var == 0.25 &
  known$param.total_sv_var == 1, 
]
summary(known_filtered$FPR)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001913 0.004420 0.008226 0.014420 0.025070

############divergence ends

############### DIVERGENCE START ##############
#pull results for supp. figure 3 + supplementary materials
oracle <- filter(res, method == "oracle")
oracle_0 <- filter(oracle, param.n_sv == "0")
dim(oracle_0) #1250   30
table(oracle_0$nsv)

be <- filter(res, method == "be")
dim(be) #263732     30
# Convert list column to length of each element
be$nsv_count <- sapply(be$nsv, length)

# Total number of rows
total_rows <- nrow(be)

# Underestimation
underestimate_count <- sum(be$nsv_count < be$param.n_sv, na.rm = TRUE)
underestimate_percent <- 100 * underestimate_count / total_rows

# Overestimation
overestimate_count <- sum(be$nsv_count > be$param.n_sv, na.rm = TRUE)
overestimate_percent <- 100 * overestimate_count / total_rows

# Output
cat("Underestimation: ", underestimate_count, "(", round(underestimate_percent, 2), "% )\n")
cat("Overestimation: ", overestimate_count, "(", round(overestimate_percent, 2), "% )\n")



############### DIVERGENCE ENDS ###############

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

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions_05142025_two.pdf", plot = combo, width = 8, height = 5) 

###################################################################################
######### P VALUE DISTRIBUTION from simulation_pvalue_distributions.sh ############
###################################################################################
case0_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case0_pval_out.csv")
case1_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case1_pval_out.csv")
case2_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case2_pval_out.csv")
case3_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case3_pval_out.csv")
case4_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case4_pval_out.csv")

case1_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case1_out.csv")


### KNOWN, ORACLE AND NONE
case0_sub <- filter(case0_pval_out, Method %in% c("known_nsv", "none", "oracle")) #FPR = 0.001041667 = 0.1041667% | vp_plots[[9]]
case1_sub <- filter(case1_pval_out, Method %in% c("known_nsv", "none", "oracle")) #FPR = 0.23450000 = 23% | this one has bimodal vpa | vp_plots[[2]]
case2_sub <- filter(case2_pval_out, Method %in% c("known_nsv", "none", "oracle")) #FPR = 0.927228261 = 93% | vp_plots[[7]]
case3_sub <- filter(case3_pval_out, Method %in% c("known_nsv", "none", "oracle")) #FPR for be = 0.0085000 = 0.85% and leek = 0.0053125 = 0.53% | vp_plots[[8]]
case4_sub <- filter(case4_pval_out, Method %in% c("known_nsv", "none", "oracle")) #FPR = 0.026714286 = 2.671429% | vp_plots[[10]]

# Filter out rows where P.Value is NA or zero to avoid log issues
filtered_data <- case1_sub %>% filter(!is.na(P.Value) & P.Value > 0)
filtered_data$Method <- factor(filtered_data$Method, 
	levels = c("oracle", "none", "known_nsv"))

# Plot histogram colored by Method
case2 <- ggplot(filtered_data, aes(x = P.Value, fill = Method)) +
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

##case0: vp_plots[[9]]
#case1: vp_plots[[2]]
#case2: vp_plots[[7]]
#case3: vp_plots[[8]]
#case4: vp_plots[[10]]

case2_vpa <- vp_plots[[7]] +
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

ggsave(filename = "/hpc/users/hoangd02/www/plots/test.png",
  plot = case2_vpa, width = 9, height = 8) 

layout <- c(
  area(1, 1, 1, 1),  # main_vp_plot at row 1, col 1
  area(1, 2, 1, 2),  # case0 at row 1, col 2
  area(2, 1, 2, 1),  # p1 at row 2, col 1
  area(3, 1, 3, 1),  # p3 at row 3, col 1
  area(2, 2, 2, 2),  # p2 at row 2, col 2
  area(3, 2, 3, 2)   # p4 at row 3, col 2
)

#combo <- main_vp_plot + case0_shrunk + p1 + p3 + p2 + p4 +

combo <- main_vp_plot + case0_shrunk + p1 + p3 + p2 + p4 +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = "/hpc/users/hoangd02/www/plots/Figure1_20250610.pdf",
  plot = combo, width = 9, height = 8) 

https://hoangd02.u.hpc.mssm.edu/plots/Figure1_20250518_2.pdf

### Check if the 5 representative cases overfit:
test <- read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_results_with_summary_vpa_stats.csv")

Case 0: n_sv =10, total_sv_sd = 2.5, sv_sd = 0.79, fraction_degs = 0.04, status_sd = 0.1, resid_sd = 0.4, fraction_sv = 0.3, seed = 12345 #row 9 in 10_random_vpa_set photo
Case 1: n_sv =10, total_sv_sd = 5, sv_sd = 1.58, fraction_degs = 0.30, status_sd = 0.2, resid_sd = 0.2, fraction_sv = 0.1, seed = 12352; #row 2 in 10_random_vpa_set photo
Case 2: n_sv =10, total_sv_sd = 7.5, sv_sd = 2.37, fraction_degs = 0.08, status_sd = 0.5, resid_sd = 0.3, fraction_sv = 0.5, seed = 12351; #row 7 in 10_random_vpa_set photo
Case 3: n_sv =10, total_sv_sd = 1.5, sv_sd = 0.47, fraction_degs = 0.20, status_sd = 0.2, resid_sd = 0.4, fraction_sv = 0.1, seed = 12348; #row 8 in 10_random_vpa_set photo
Case 4: n_sv =10, total_sv_sd = 2.0, sv_sd = 0.63, fraction_degs = 0.30, status_sd = 0.2, resid_sd =0.3, fraction_sv = 0.7, seed = 12348 #row 10 in 10_random_vpa_set photo

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

all_data<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_results_with_summary_vpa_stats.csv")

test<-read.csv("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_91percentsimulations_supplementary_figure.csv")














