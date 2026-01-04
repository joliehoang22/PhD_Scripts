###################################################################################
############################# Snake plot from simData_plot.R ######################
###################################################################################
#https://hoangd02.u.hpc.mssm.edu/plots/sva/gandal_sim_results/FPR_by_cor_combined_simData1_and_simData2_status_sd0.1.pdf

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrastr)

out1 <- read.csv('/sc/arion/projects/mscic1/results/jolie/sva_sim_gandal2/gandal_sim_results_status_sd_0.1_20240911.csv') 
out2 <- read.csv('/sc/arion/projects/mscic1/results/jolie/sva_sim_gandal2/gandal_sim_results_status_sd_0.2_20240911.csv') 
out3 <- read.csv('/sc/arion/projects/mscic1/results/jolie/sva_sim_gandal2/gandal_sim_results_status_sd_0.3_20240911.csv') 
out4 <- read.csv('/sc/arion/projects/mscic1/results/jolie/sva_sim_gandal2/gandal_sim_results_status_sd_0.4_20240911.csv') 
out5 <- read.csv('/sc/arion/projects/mscic1/results/jolie/sva_sim_gandal2/gandal_sim_results_status_sd_0.5_20240911.csv') 

out1 <- out1 %>% 
  mutate(
    FPR_1 = FP_1 / (FP_1 + TN_1),
    FNR_1 = FN_1 / (FN_1 + TP_1),
    TNR_1 = TN_1 / (TN_1 + FP_1),
    TPR_1 = TP_1 / (TP_1 + FN_1))

out2 <- out2 %>% 
  mutate(
    FPR_1 = FP_1 / (FP_1 + TN_1),
    FNR_1 = FN_1 / (FN_1 + TP_1),
    TNR_1 = TN_1 / (TN_1 + FP_1),
    TPR_1 = TP_1 / (TP_1 + FN_1))

out3 <- out3 %>% 
  mutate(
    FPR_1 = FP_1 / (FP_1 + TN_1),
    FNR_1 = FN_1 / (FN_1 + TP_1),
    TNR_1 = TN_1 / (TN_1 + FP_1),
    TPR_1 = TP_1 / (TP_1 + FN_1))

out4 <- out4 %>% 
  mutate(
    FPR_1 = FP_1 / (FP_1 + TN_1),
    FNR_1 = FN_1 / (FN_1 + TP_1),
    TNR_1 = TN_1 / (TN_1 + FP_1),
    TPR_1 = TP_1 / (TP_1 + FN_1))

out5 <- out5 %>% 
  mutate(
    FPR_1 = FP_1 / (FP_1 + TN_1),
    FNR_1 = FN_1 / (FN_1 + TP_1),
    TNR_1 = TN_1 / (TN_1 + FP_1),
    TPR_1 = TP_1 / (TP_1 + FN_1))

out1 <- out1 %>% 
  mutate(
    FPR_2 = FP_2 / (FP_2 + TN_2),
    FNR_2 = FN_2 / (FN_2 + TP_2),
    TNR_2 = TN_2 / (TN_2 + FP_2),
    TPR_2 = TP_2 / (TP_2 + FN_2))

out2 <- out2 %>% 
  mutate(
    FPR_2 = FP_2 / (FP_2 + TN_2),
    FNR_2 = FN_2 / (FN_2 + TP_2),
    TNR_2 = TN_2 / (TN_2 + FP_2),
    TPR_2 = TP_2 / (TP_2 + FN_2))

out3 <- out3 %>% 
  mutate(
    FPR_2 = FP_2 / (FP_2 + TN_2),
    FNR_2 = FN_2 / (FN_2 + TP_2),
    TNR_2 = TN_2 / (TN_2 + FP_2),
    TPR_2 = TP_2 / (TP_2 + FN_2))

out4 <- out4 %>% 
  mutate(
    FPR_2 = FP_2 / (FP_2 + TN_2),
    FNR_2 = FN_2 / (FN_2 + TP_2),
    TNR_2 = TN_2 / (TN_2 + FP_2),
    TPR_2 = TP_2 / (TP_2 + FN_2))

out5 <- out5 %>% 
  mutate(
    FPR_2 = FP_2 / (FP_2 + TN_2),
    FNR_2 = FN_2 / (FN_2 + TP_2),
    TNR_2 = TN_2 / (TN_2 + FP_2),
    TPR_2 = TP_2 / (TP_2 + FN_2))

#create a new column with values of each status sd 
out1$dataset <- 'sd_0.1'
out2$dataset <- 'sd_0.2'
out3$dataset <- 'sd_0.3'
out4$dataset <- 'sd_0.4'
out5$dataset <- 'sd_0.5'

#combine all the datasets tgt!
combined <- rbind(out1, out2, out3, out4, out5)
combined$total_positive1 <- combined$TP_1 + combined$FP_1
combined$total_positive2 <- combined$TP_2 + combined$FP_2

table(combined$dataset)
# sd_0.1 sd_0.2 sd_0.3 sd_0.4 sd_0.5 
#   2600   2600   2600   2600   2600 
#TPR_vs_FPR <- ggplot(combined, aes(x = FPR_1, y = TPR_1, color = model_nsv)) +

dim(combined) #13000    34

##################################### 9/6/2025
##### ##### ##### NEW PLOT FOR NOAM ##### ##### ##### 
# plot <- ggplot(combined, aes(x = cor, y = model_nsv, color = factor(param.status_sd))) +
#   geom_point(alpha = 0.7, size = 2) +
#   scale_color_brewer(palette = "Set1", name = "Status SD") +
#   labs(
#     x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
#     y = "Model nsv",
#     title = "Correlation vs Model nsv colored by Status SD"
#   ) +
#   theme_minimal(base_size = 14)

plot <- ggplot(combined, aes(x = model_nsv, y = cor, color = factor(param.status_sd))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_viridis_d(name = "Magnitude of\nthe differences between\ncase and controls") +
  labs(
    x = "Number of SVs",
    y = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

ggsave("/hpc/users/hoangd02/www/plots/testtest.pdf", plot = plot, width = 12, height = 7)

ggsave("/hpc/users/hoangd02/www/plots/gandal_model_nsv_correlation_simData1_simData2_20250906.pdf", plot = plot, width = 12, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_test_20250906.pdf", plot = plot, width = 12, height = 5)

# Add loop_seed column that cycles from 1 to 100
combined$loop_seed <- rep(1:100, length.out = nrow(combined))

# double checking
#combined[99:102, ]
#combined[199:202, ]
#combined[599:602, ]

#subset to rows with the highest cor per seed per status_sd, dim should be 500 rows | REWRITING COMBINED
combined <- combined %>%
  group_by(loop_seed, param.status_sd) %>%
  slice_max(cor, n = 1, with_ties = FALSE) %>%
  ungroup()

dim(combined) #500  35  

## create cor_bin for ALL not for individual param.status_sd
combined <- filter(combined, param.status_sd == "0.3") ##this needs to be before cor_bin or else there will not be 10 bins for status sd of 0.3
## the problem with 0.4 is that for 1 cor_bin it's empty

########################################
########## THE CHOSEN STATUS SD IS 0.3
########################################

combined <- combined %>%
  mutate(cor_bin = cut(cor, breaks = 10)) %>%
  ungroup()

combined<-as.data.frame(combined)
table(combined$cor_bin)

sub <- filter(combined, cor_bin == "(0.464,0.472]")

#combined <- combined[order(combined$model_nsv), ]

## create cor_bin per param.status_sd to look at TPR and FPR per param.status_sd for supplementary
#combined <- combined %>%
#  group_by(param.status_sd) %>%
#  mutate(cor_bin = cut(cor, breaks = 5)) %>%
#  ungroup()
##just plotting status_sd = 0.1

#################################### plot TPR color with both simData1 and simData1 rbinded for only status_sd == 0.1 
plot_data <- combined %>%
  select(TPR_1, TPR_2, cor_bin) %>%
  pivot_longer(cols = c(TPR_1, TPR_2), 
               names_to = "Dataset", 
               values_to = "TPR") %>%
  select(cor_bin, TPR)

tpr <- ggplot(plot_data, aes(x = cor_bin, y = TPR)) +
  geom_violin(fill = "#9932CC", alpha = 0.5) +
  geom_boxplot(fill = "#9932CC", alpha = 0.5, width = 0.2, outlier.alpha = 0.5) +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2 (Deciles)"),
       y = "True Positive Rate") +
  ylim(NA,0.95) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(drop = FALSE) 

ggsave("/hpc/users/hoangd02/www/plots/gandal_test_status_sd0.3.pdf", plot = tpr, width = 12, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_test_status_sd0.4.pdf", plot = tpr, width = 12, height = 5)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_test2.pdf", plot = plot, width = 12, height = 5)
plot_data <- combined %>%
  select(FPR_1, FPR_2, cor_bin) %>%
  pivot_longer(cols = c(FPR_1, FPR_2), 
               names_to = "Dataset", 
               values_to = "FPR") %>%
  select(cor_bin, FPR)

fpr <- ggplot(plot_data, aes(x = cor_bin, y = FPR)) +
  geom_violin(fill = "#9932CC", alpha = 0.3) +
  geom_boxplot(fill = "#9932CC", alpha = 0.5, width = 0.2, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2 (Deciles)"),
       y = "False Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(drop = FALSE) 

combo <- (tpr + fpr)
ggsave("/hpc/users/hoangd02/www/plots/gandal_test_combo_status_sd0.3.pdf", plot = combo, width = 13, height = 4)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_test_combo_status_sd0.4.pdf", plot = combo, width = 15, height = 5)

##quick check, can ignore in the future
# combined <- combined %>%
#   mutate(
#     param.total_sv_sd = sqrt(param.sv_sd^2 * param.n_sv),
#     param.n_batches = 15,
#     param.batch_sd = 0.3,
#     ## For status and batch, we multiply by (n-1)/n where n is the
#     ## number of levels in the variable.
#     param.mean_batch_var = param.batch_sd^2 * (param.n_batches - 1) / param.n_batches,
#     param.mean_sv_var = param.total_sv_sd^2 * param.fraction_sv,
#     param.mean_status_var = param.status_sd^2 * param.fraction_degs * 1/2,
#     param.mean_resid_var = param.resid_sd^2,
#     param.total_mean_var = param.mean_batch_var + param.mean_sv_var + param.mean_status_var + param.mean_resid_var,
#     param.batch_var_fraction = param.mean_batch_var / param.total_mean_var,
#     param.sv_var_fraction = param.mean_sv_var / param.total_mean_var,
#     param.status_var_fraction = param.mean_status_var / param.total_mean_var,
#     param.resid_var_fraction = param.mean_resid_var / param.total_mean_var,
#     param.status_var_percent = param.status_var_fraction * 100
#   )

#################################### plot TPR color by status sd and facet wrap by simData1 and simData2
#combine all the datasets tgt!
combined <- rbind(out1, out2, out3, out4, out5)
combined$total_positive1 <- combined$TP_1 + combined$FP_1
combined$total_positive2 <- combined$TP_2 + combined$FP_2

table(combined$dataset)
# sd_0.1 sd_0.2 sd_0.3 sd_0.4 sd_0.5 
#   2600   2600   2600   2600   2600 
#TPR_vs_FPR <- ggplot(combined, aes(x = FPR_1, y = TPR_1, color = model_nsv)) +

dim(combined) #13000    34

# Add loop_seed column that cycles from 1 to 100
combined$loop_seed <- rep(1:100, length.out = nrow(combined))

# double checking
#combined[99:102, ]
#combined[199:202, ]
#combined[599:602, ]

#subset to rows with the highest cor per seed per status_sd, dim should be 500 rows | REWRITING COMBINED
combined <- combined %>%
  group_by(loop_seed, param.status_sd) %>%
  slice_max(cor, n = 1, with_ties = FALSE) %>%
  ungroup()

dim(combined) #500  35  

combined <- combined %>%
  mutate(cor_bin = cut(cor, breaks = 10)) %>%
  ungroup()

combined_long <- combined %>%
  # Create a long format with simData1 and simData2
  pivot_longer(cols = c(TPR_1, TPR_2), 
               names_to = "sim_data", 
               values_to = "TPR") %>%
  # Clean up the sim_data labels
  mutate(sim_data = case_when(
    sim_data == "TPR_1" ~ "simData1",
    sim_data == "TPR_2" ~ "simData2"
  ))

combined_long<-as.data.frame(combined_long)
head(combined_long)

# TPR color by status sd and facet wrap by simData1 and simData2
#this used to be called plot
p1<- ggplot(combined_long, aes(x = cor_bin, y = TPR, fill = factor(param.status_sd))) +
  geom_violin(
    alpha = 0.5,
    scale = "width", 
    width = 0.9,
    position = position_dodge(width = 0.95)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 0.95),
    outlier.shape = NA
  ) +
  facet_wrap(~ sim_data, scales = "free") +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
      y = "True Positive Rate", fill = "Magnitude of\nthe differences between\ncase and controls") +
  scale_fill_viridis_d() +
  ylim(NA,1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    panel.spacing = unit(0.5, "cm")
  )
#ggsave("/hpc/users/hoangd02/www/plots/p1_test.pdf", plot = p1, width = 12, height = 5)

ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_cor_bin_status_sd_simData1_simData2.pdf", plot = plot, width = 12, height = 5)

https://hoangd02.u.hpc.mssm.edu/plots/gandal_test.pdf

combined_long <- combined %>%
  # Create a long format with simData1 and simData2
  pivot_longer(cols = c(FPR_1, FPR_2), 
               names_to = "sim_data", 
               values_to = "FPR") %>%
  # Clean up the sim_data labels
  mutate(sim_data = case_when(
    sim_data == "FPR_1" ~ "simData1",
    sim_data == "FPR_2" ~ "simData2"
  ))

combined_long<-as.data.frame(combined_long)
head(combined_long)

sub <- filter(combined_long, cor_bin == "(0.464,0.472]")
mean(sub$FPR < 0.05, na.rm = TRUE) * 100 #75%
sum(sub$FPR < 0.05, na.rm = TRUE) #9 #number of rows where the false positive rate is below 0.05.

summary(sub$FPR)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00800 0.01241 0.02006 0.03707 0.04823 0.09775

# FPR color by status sd and facet wrap by simData1 and simData2
#this used to be plot
p2<- ggplot(combined_long, aes(x = cor_bin, y = FPR, fill = factor(param.status_sd))) +
  geom_violin(
    alpha = 0.5,
    scale = "width", 
    width = 0.9,
    position = position_dodge(width = 0.95)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 0.95),
    outlier.shape = NA
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ sim_data, scales = "free") +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
      y = "False Positive Rate", fill = "Magnitude of\nthe differences between\ncase and controls") +
  scale_fill_viridis_d() +
  ylim(NA,1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    panel.spacing = unit(0.5, "cm")
  )

ggsave("/hpc/users/hoangd02/www/plots/p2_test.pdf", plot = p2, width = 12, height = 5)

ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_cor_bin_status_sd_simData1_simData2.pdf", plot = plot, width = 12, height = 5)

###making extended data figure 4: 
library(patchwork)

combo <- (plot / p1 / p2) +
  plot_annotation(tag_levels = 'A')  # adds a, b, c

ggsave("/hpc/users/hoangd02/www/plots/combo_test.pdf", 
       plot = combo, width = 12, height = 15)

#############################
#############################

table(combined$param.status_sd)
# 0.1 0.2 0.3 0.4 0.5 
# 100 100 100 100 100

combined_long <- combined %>%
  select(cor_bin, TPR_1, TPR_2, FPR_1, FPR_2) %>%
  pivot_longer(
    cols = c(TPR_1, TPR_2, FPR_1, FPR_2),
    names_to = c("Metric", "SimData"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    SimData = recode(SimData, "1" = "SimData1", "2" = "SimData2"),
    Metric = recode(Metric, "TPR" = "TPR", "FPR" = "FPR") 
  )

combined_long <-as.data.frame(combined_long)
dim(combined_long)
head(combined_long)

tpr_plot <- combined_long %>%
  filter(Metric == "TPR") %>%
  ggplot(aes(x = cor_bin, y = Value, fill = SimData)) +
  geom_violin(alpha = 0.5,scale = "width",width = 0.9,position = position_dodge(width = 0.95)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.95),outlier.shape = NA) +
  scale_fill_manual(
    values = c("SimData1" = "#d633a6", "SimData2" = "#a0f600"),
    name = "Dataset") +
    labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
      y = "True Positive Rate", fill = "Dataset") +
  ylim(NA, 1) + 
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_param_status_sd01.pdf", plot = tpr_plot, width = 8, height = 5, dpi = 300)
##right


# highest<- filter(combined, cor_bin == "(0.464,0.484]") ## this is at 12 cor bins btw!!! 
# summary(highest$param.status_var_percent)
# #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.2411  0.4278  0.6668  0.5595  0.6668  0.6668 


# ##now i have param.status_var_fraction and param.status_var_percent
# table(combined$param.status_sd)

status_sd0.1<- filter(combined, param.status_sd == "0.1") ## 
summary(status_sd0.1$param.status_var_percent)
# #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #0.02685 0.02685 0.02685 0.02685 0.02685 0.02685 

status_sd0.2<- filter(combined, param.status_sd == "0.2") ## 
summary(status_sd0.2$param.status_var_percent)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.1073  0.1073  0.1073  0.1073  0.1073  0.1073

status_sd0.3<- filter(combined, param.status_sd == "0.3") ## 
summary(status_sd0.3$param.status_var_percent)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.2411  0.2411  0.2411  0.2411  0.2411  0.2411 

status_sd0.4<- filter(combined, param.status_sd == "0.4") ## 
summary(status_sd0.4$param.status_var_percent)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.4278  0.4278  0.4278  0.4278  0.4278  0.4278 

status_sd0.5<- filter(combined, param.status_sd == "0.5") ## 
summary(status_sd0.5$param.status_var_percent)
# #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.6668  0.6668  0.6668  0.6668  0.6668  0.6668 

## SUPPLEMENTARY -- 
#sub <- filter(combined, param.status_sd == "0.4")
#plot <- ggplot(combined, aes(x = model_nsv)) +
#  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
#  labs(title = "Distribution of Number of SVs Used by Model", 
#       x = "model_nsv", 
#       y = "Count") +
#  theme_minimal()

#ggsave("/hpc/users/hoangd02/www/plots/gandal_dist_of_nsv.pdf", plot = plot, width = 10, height = 6)

#################################### JUMP TO PLOTTING
####################################

#combined$method <- factor(combined$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))

#combined <- combined[order(combined$FPR_1), ]
#combined <- combined[order(combined$FPR_2), ]
combined <- rbind(out1, out2, out3, out4, out5)
sd_0.3 <- filter(combined, dataset=="sd_0.3")
#sd_0.1$seed <- rep(1:100, 26) 
sd_0.3 <- sd_0.3 %>%
  group_by(model_nsv) %>%
  mutate(seed_id = row_number()) %>%
  ungroup()

sd_0.3 <- as.data.frame(sd_0.3)
### SHOWING A RANDOM 5 SEEDS
set.seed(123)  # For reproducibility

set.seed(2025)
selected_seeds <- sample(1:100, 5) # Pick 5 random unique seeds: 31 79 51 14 67
selected_seeds #13 76 36 26 65 set.seed(2025)good
#28  5 30 37 85  set.seed(2017)
#28 87 22 88 65 2020 

sd_subset <- sd_0.3 %>% 
  filter(seed_id %in% selected_seeds)

table(sd_subset$model_nsv) 

# Find min and max values to set consistent scales
y_min <- min(c(sd_subset$FPR_1, sd_subset$FPR_2), na.rm = TRUE)
y_max <- max(c(sd_subset$FPR_1, sd_subset$FPR_2), na.rm = TRUE)

color_min <- min(c(sd_subset$TPR_1, sd_subset$TPR_2), na.rm = TRUE)
color_max <- max(c(sd_subset$TPR_1, sd_subset$TPR_2), na.rm = TRUE)

# Pick label points at max correlation for each line
# label_points <- sd_subset %>%
#   group_by(seed_id, model_nsv) %>%
#   slice_max(order_by = cor, n = 1)

# a_OG <- ggplot(sd_subset, aes(x = cor, y = FPR_1, color = TPR_1, group = seed_id)) +
#   geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
#   geom_text_repel(
#     data = label_points,
#     aes(label = model_nsv),
#     size = 5, fontface = "bold",
#     box.padding = 0.5,
#     force = 2,
#     force_pull = 0.5,
#     max.overlaps = Inf,
#     min.segment.length = 0,
#     segment.color = "gray30",
#     segment.alpha = 0.3
#   ) +
#   facet_wrap(~as.factor(seed_id), nrow = 1) +
#   scale_color_viridis_c(
#     option = "viridis",
#     name = "True Positive Rate\nof SimData1",
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.25)
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.25)
#   ) +
#   labs(
#     x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
#     y = "False Positive Rate of SimData1"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.grid.major = element_line(color = "gray90", size = 0.3),
#     panel.grid.minor = element_line(color = "gray95", size = 0.2),
#     legend.position = "right",
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 16),
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
#     strip.text = element_text(size = 16)
#   )

##start
label_points <- sd_subset %>%
  filter(model_nsv %in% seq(0, 25, 5)) %>%
  group_by(seed_id) %>%
  nest() %>%
  left_join(
    sd_subset %>%
      group_by(seed_id) %>%
      nest() %>%
      mutate(loess_fit = purrr::map(data, ~ loess(FPR_1 ~ cor, data = .x, span = 0.75))) %>%
      select(seed_id, loess_fit),
    by = "seed_id"
  ) %>%
  mutate(FPR_1_smooth = purrr::map2(data, loess_fit, ~ predict(.y, newdata = .x))) %>%
  unnest(cols = c(data, FPR_1_smooth)) %>%
  ungroup() %>%
  mutate(FPR_1_smooth_offset = FPR_1_smooth + 0.015) ## adjust this so the line/text do not overlap with loess

a <- ggplot(sd_subset, aes(x = cor, y = FPR_1, color = TPR_1, group = seed_id)) +
  geom_smooth(aes(group = seed_id), method = "loess", se = FALSE, linewidth = 1) +
  geom_text_repel(
    data = label_points,
    aes(x = cor, y = FPR_1_smooth, label = model_nsv),
    size = 5, fontface = "bold",
    box.padding = 1,
    force =10,
    force_pull = 2,
    nudge_y = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "gray30",
    segment.alpha = 0.3
  )+
  facet_wrap(~as.factor(seed_id), nrow = 1) +
  scale_color_viridis_c(
    option = "viridis",
    name = "True Positive Rate\nof SimData1",
    limits = c(color_min, color_max),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
    y = "False Positive Rate of SimData1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 16)
  )

ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a_20250821.pdf", plot = a, width = 22, height = 4)

#ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a_08042025.png", plot = a, width = 22, height = 4)

#ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a_four.png", plot = a, width = 22, height = 4)

#https://hoangd02.u.hpc.mssm.edu/plots/static_gandal_snake_plot_a_two.png

#ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a.png", plot = a_OG, width = 22, height = 4)

# Plot for SimData2
# label_points_b <- sd_subset %>%
#   group_by(seed_id, model_nsv) %>%
#   slice_max(order_by = cor, n = 1)

# b <- ggplot(sd_subset, aes(x = cor, y = FPR_2, color = TPR_2, group = seed_id)) +
#   geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
#   geom_text_repel(
#     data = label_points_b,
#     aes(label = model_nsv),
#     size = 5, fontface = "bold",
#     box.padding = 0.5,
#     force = 2,
#     force_pull = 0.5,
#     max.overlaps = Inf,
#     min.segment.length = 0,
#     segment.color = "gray30",
#     segment.alpha = 0.3
#   ) +
#   facet_wrap(~as.factor(seed_id), nrow = 1) +
#   scale_color_viridis_c(
#     option = "viridis",
#     name = "True Positive Rate\nof SimData2",
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.25)
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.25)
#   ) +
#   labs(
#     x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
#     y = "False Positive Rate of SimData2"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.grid.major = element_line(color = "gray90", size = 0.3),
#     panel.grid.minor = element_line(color = "gray95", size = 0.2),
#     legend.position = "right",
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 16),
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
#     strip.text = element_text(size = 16)
#   )

label_points_b <- sd_subset %>%
  filter(model_nsv %in% seq(0, 25, 5)) %>%
  group_by(seed_id) %>%
  nest() %>%
  left_join(
    sd_subset %>%
      group_by(seed_id) %>%
      nest() %>%
      mutate(loess_fit = purrr::map(data, ~ loess(FPR_2 ~ cor, data = .x, span = 0.75))) %>%
      select(seed_id, loess_fit),
    by = "seed_id"
  ) %>%
  mutate(FPR_2_smooth = purrr::map2(data, loess_fit, ~ predict(.y, newdata = .x))) %>%
  unnest(cols = c(data, FPR_2_smooth)) %>%
  ungroup() %>%
  mutate(FPR_2_smooth_offset = FPR_2_smooth + 0.01) ## adjust this so the line/text do not overlap with loess

b <- ggplot(sd_subset, aes(x = cor, y = FPR_2, color = TPR_2, group = seed_id)) +
  geom_smooth(aes(group = seed_id), method = "loess", se = FALSE, linewidth = 1) +
  geom_text_repel(
    data = label_points_b,
    aes(x = cor, y = FPR_2_smooth, label = model_nsv), #FPR_2_smooth_offset
    size = 5, fontface = "bold",
    box.padding = 1,
    force = 5,
    force_pull = 5,
    nudge_y = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "gray30",
    segment.alpha = 0.3,
  ) +
  facet_wrap(~as.factor(seed_id), nrow = 1) +
  scale_color_viridis_c(
    option = "viridis",
    name = "True Positive Rate\nof SimData2",
    limits = c(color_min, color_max),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_y_continuous(
    limits = c(-0.01, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
    y = "False Positive Rate of SimData2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 16)
  )
ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_b_20250821.pdf", plot = b, width = 22, height = 4) ## FINAL

#ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_b_08042025.png", plot = b, width = 22, height = 4)

#ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_b_tw0.png", plot = b, width = 22, height = 4)

# ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_b.png", plot = b_OG, width = 22, height = 4)

# ggsave(filename = "/hpc/users/hoangd02/www/plots/test.pdf",
#   plot = b, width=20, height=4) 

#ggsave("/hpc/users/hoangd02/www/plots/sva/gandal_sim_results/FPR_by_cor_color_by_TPR_0.1_simData2.pdf", plot=b, width = 20, height = 5)

#combined <- combined %>%
#  mutate(cor_bin = cut(cor, breaks = 10))

#aa <- filter(combined, param.status_sd == "0.1")
#summary(aa$cor)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1289  0.2169  0.2521  0.2549  0.2871  0.4012 

# ab <- filter(combined, param.status_sd == "0.5")
# summary(ab$cor)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3838  0.4264  0.4357  0.4369  0.4473  0.4840 


# Reshape both TPR and FPR into long format - THIS IS THE ORIGNAL no status_sd info
combined_long <- combined %>%
  select(cor_bin, TPR_1, TPR_2, FPR_1, FPR_2) %>%
  pivot_longer(
    cols = c(TPR_1, TPR_2, FPR_1, FPR_2),
    names_to = c("Metric", "SimData"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    SimData = recode(SimData, "1" = "SimData1", "2" = "SimData2"),
    Metric = recode(Metric, "TPR" = "TPR", "FPR" = "FPR") 
  )

combined_long <-as.data.frame(combined_long)
head(combined_long)
##quantifying FPR
sub <- filter(combined_long, Metric == "FPR")
dim(sub) #1000 4

#how many are < 0.05
sum(sub$Value < 0.05) #242 = 24.2%

##find median of highest cor_bin 
sub <- filter(combined_long, cor_bin == "(0.464,0.484]")
sub_fpr <- filter(sub, Metric == "FPR")
summary(sub_fpr$Value)

################### RBINDING THINGS TOGETHER ################### 
################### ################### ################### ####
# Extract SimData1 (columns ending with _1)
combined_long <- combined %>%
  # Reshape to separate the two simulation datasets
  pivot_longer(cols = c(FPR_1, TPR_1, FPR_2, TPR_2),
               names_to = c("Metric", "SimData"),
               names_pattern = "(.+)_([12])",
               values_to = "Value") %>%
  # Convert SimData numbers to descriptive names
  mutate(SimData = case_when(
    SimData == "1" ~ "SimData1",
    SimData == "2" ~ "SimData2"
  ))

head(combined_long)
combined_long <- as.data.frame(combined_long)
dim(combined_long) #2000   35
table(combined_long$SimData)
# SimData1 SimData2 
#     1000     1000 

#original data 
# Row 1: cor_bin="(0.362,0.386]", FPR_1=0.264, TPR_1=0.773, FPR_2=0.148, TPR_2=0.762

#transformed data
# Row 1a: cor_bin="(0.362,0.386]", Metric="FPR", SimData="SimData1", Value=0.264
# Row 1b: cor_bin="(0.362,0.386]", Metric="TPR", SimData="SimData1", Value=0.773  
# Row 1c: cor_bin="(0.362,0.386]", Metric="FPR", SimData="SimData2", Value=0.148
# Row 1d: cor_bin="(0.362,0.386]", Metric="TPR", SimData="SimData2", Value=0.762

fpr_plot <- combined_long %>%
  filter(Metric == "FPR") %>%
  ggplot(aes(x = cor_bin, y = Value)) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    width = 0.9,
    fill = "#9932CC"
    #fill = "#FF8C00"  #orange
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#8C1515", linewidth = 1) +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
       y = "False Positive Rate") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot_combined.pdf", plot = fpr_plot, width = 13, height = 8, dpi = 300)

tpr_plot <- combined_long %>%
  filter(Metric == "TPR") %>%
  ggplot(aes(x = cor_bin, y = Value)) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    width = 0.9,
    fill = "#9932CC"  # Dark orchid
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
       y = "True Positive Rate") +
  ylim(NA,1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

combo <- (tpr_plot + fpr_plot)
ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_tpr_plot_combined_20250820.pdf", plot = combo, width = 18, height = 6, dpi = 300)

################### ################### ################### ####
################### ################### ################### ####

########### FPR ############
fpr_plot <- combined_long %>%
  filter(Metric == "FPR") %>%
  ggplot(aes(x = cor_bin, y = Value, fill = SimData)) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    width = 0.9,
    position = position_dodge(width = 0.95)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 0.95),
    outlier.shape = NA
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#8C1515", linewidth = 1) +
  scale_fill_manual(
    values = c("SimData1" = "#d633a6", "SimData2" = "#a0f600"),
    name = "Dataset") +
    labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),  	   
    	y = "False Positive Rate", fill = "Dataset") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
    )
ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot_dodeciles.pdf", plot = fpr_plot, width = 15, height = 8, dpi = 300)

#12 cor_bins
ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot5.pdf", plot = fpr_plot, width = 15, height = 8, dpi = 300)

#10 cor_bins
ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot4.pdf", plot = fpr_plot, width = 15, height = 8, dpi = 300)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot_nsv=25.pdf", plot = fpr_plot, width = 14, height = 8, dpi = 300)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot.pdf", plot = fpr_plot, width = 10, height = 8, dpi = 300)


########### TPR ############
tpr_plot <- combined_long %>%
  filter(Metric == "TPR") %>%
  ggplot(aes(x = cor_bin, y = Value, fill = SimData)) +
  geom_violin(alpha = 0.5,scale = "width",width = 0.9,position = position_dodge(width = 0.95)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.95),outlier.shape = NA) +
  scale_fill_manual(
    values = c("SimData1" = "#d633a6", "SimData2" = "#a0f600"),
    name = "Dataset") +
    labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2"),
  		y = "True Positive Rate", fill = "Dataset") +
  ylim(NA, 1) + 
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_plot_dodeciles.pdf", plot = tpr_plot, width = 15, height = 8, dpi = 300)

#ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_plot_nsv=25.pdf", plot = tpr_plot, width = 14, height = 8, dpi = 300)

ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_plot.pdf", plot = tpr_plot, width = 10, height = 8, dpi = 300)

ggsave("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a.png", plot = a, width = 20, height = 4, dpi = 300)

library(png)
library(grid)

a_raster <- rasterGrob(readPNG("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_a.png"), interpolate = TRUE)

b_raster <- rasterGrob(readPNG("/hpc/users/hoangd02/www/plots/static_gandal_snake_plot_b.png"), interpolate = TRUE)

combo <- (tpr_plot + fpr_plot) / wrap_elements(a_raster) / wrap_elements(b_raster) 

combo <- (tpr_plot + fpr_plot) / 
         wrap_elements(a_raster) / 
         wrap_elements(b_raster) +
  plot_layout(heights = c(0.8, 1.2, 1.2))  # or tweak to c(1.2, 1, 1)

##final 
combo <- (tpr_plot + fpr_plot) / a / b +
  plot_layout(heights = c(1.3, 1, 1))

ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_combo_20250804.pdf",
  plot = combo, width = 20, height = 13) 

combo <- (tpr_plot + fpr_plot)
ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_tpr_fpr_combo_20250820.pdf",plot = combo, width = 20, height = 6) 

#ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_combo.pdf",
#  plot = combo, width = 20, height = 14) 

#ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_combo_20250524_two.pdf",
#  plot = combo, width = 20, height = 14) 

#ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_combo_20250524_three.pdf",
#  plot = combo, width = 20, height = 14) 

#ggsave(filename = "/hpc/users/hoangd02/www/plots/gandal_simData1_simData2_combo_20250524.pdf",
#  plot = combo, width = 20, height = 14) 

#https://hoangd02.u.hpc.mssm.edu/plots/gandal_simData1_simData2_combo_20250524_four.pdf

######################## FINAL ONE ##################
https://hoangd02.u.hpc.mssm.edu/plots/gandal_simData1_simData2_combo_20250804.pdf

#########################
# UPDATED AUG 19 2025
# Add param.status_sd to the pivot_longer transformation
# First, create a proper ordering factor for cor_bin within each param.status_sd

combined_long <- combined %>%
  group_by(param.status_sd) %>%
  mutate(
    # Extract the lower bound of each interval for sorting
    cor_lower = as.numeric(sub("\\((.+),.*", "\\1", cor_bin)),
    cor_bin = factor(cor_bin, levels = unique(cor_bin[order(cor_lower)]))
  ) %>%
  select(-cor_lower) %>%
  ungroup() %>%
  select(cor_bin, TPR_1, TPR_2, FPR_1, FPR_2, param.status_sd) %>%  # Keep param.status_sd
  pivot_longer(
    cols = c(TPR_1, TPR_2, FPR_1, FPR_2),
    names_to = c("Metric", "SimData"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    SimData = recode(SimData, "1" = "SimData1", "2" = "SimData2"),
    Metric = recode(Metric, "TPR" = "TPR", "FPR" = "FPR") 
  )

########### FPR ############
fpr_plot <- combined_long %>%
  filter(Metric == "FPR") %>%
  ggplot(aes(x = cor_bin, y = Value, fill = SimData)) +
  geom_violin(
    alpha = 0.5,
    scale = "width", 
    width = 0.9,
    position = position_dodge(width = 0.95)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 0.95),
    outlier.shape = NA
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#8C1515", linewidth = 1) +
  scale_fill_manual(
    values = c("SimData1" = "#d633a6", "SimData2" = "#a0f600"),
    name = "Dataset") +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2 (Quintiles)"),      
       y = "False Positive Rate", fill = "Dataset") +
  facet_wrap(~ param.status_sd, labeller = label_both, scales = "free_x") +  # Add scales = "free_x"
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Made smaller
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    panel.spacing = unit(0.5, "cm")  # Add space between panels
  )

ggsave("/hpc/users/hoangd02/www/plots/gandal_fpr_plot3.pdf", plot = fpr_plot, width = 15, height = 8, dpi = 300)

########### TPR ############ IN ONE ROW
tpr_plot <- combined_long %>%
  filter(Metric == "TPR") %>%
  ggplot(aes(x = cor_bin, y = Value, fill = SimData)) +
  geom_violin(
    alpha = 0.5,
    scale = "width", 
    width = 0.9,
    position = position_dodge(width = 0.95)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 0.95),
    outlier.shape = NA
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#8C1515", linewidth = 1) +
  scale_fill_manual(
    values = c("SimData1" = "#d633a6", "SimData2" = "#a0f600"),
    name = "Dataset") +
  labs(x = expression("Correlation between " * log[2] * "FC of SimData1 and SimData2 (Quintiles)"),      
       y = "True Positive Rate", fill = "Dataset") +
  facet_wrap(~ param.status_sd, labeller = label_both, scales = "free_x", nrow = 1) +  # Added nrow = 1
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    panel.spacing = unit(0.5, "cm")
  )

ggsave("/hpc/users/hoangd02/www/plots/gandal_tpr_plot_facetwrap_by_statussd.pdf", plot = tpr_plot, width = 16, height = 5, dpi = 300)

https://hoangd02.u.hpc.mssm.edu/plots/gandal_tpr_plot_facetwrap_by_statussd.pdf

##try plotting 2 ways
#1/ TPR color by status sd and facet wrap by simData1 and simData2
#2/ line plots color by status sd rbind simData1 and simData2
#otherwise pick one statusSD and show that 

#add the facetwrap tgt as supplementary











