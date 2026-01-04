###OCTOBER 17TH 2024
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(purrr)
library(furrr)
library(ggrastr)

res <- fread("/sc/arion/projects/mscic1/results/jolie/sva_sim/all_results_of_sva_simulation_copy.txt", data.table = FALSE)
dim(res) #13,185,631 x 13 | this should have 2887500 - 14 = 2,887,486 sets of simulations

set6 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set6.rds") #3161 sets
set7 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set7.rds") #721 sets
set8 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set8.rds") #265 sets
set9 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set9.rds") #17 sets
set10 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set10.rds") #6 sets

dim(res)
#dim: 13,185,631 x 20

#checking number of rows - this works!
13159285 + 18966 + 5047 + 2120 + 153 + 60 = 13,185,631

#turning lists into dataframes: 
set6_df <- do.call(rbind, set6) #18966 x 13
set7_df <- do.call(rbind, set7) #5047 x 13
set8_df <- do.call(rbind, set8) #2120 x 13
set9_df <- do.call(rbind, set9) #153 x 13
set10_df <- do.call(rbind, set10) #60 x 13
##the math checks out 

combo_df_to_remove <- rbind(set6_df,set7_df,set8_df,set9_df,set10_df)
dim(combo_df_to_remove) #26346 x 13

#remove matching rows
set5 <- anti_join(res, combo_df_to_remove, by = colnames(res))

dim(set5) #13,159,285 x 13
head(set5)

#REMEMBER THE OG GRID HAS 2,887,500 SETS OF SIMULATIONS BUT
#14 out of 2,887,500 (0.005%) sets of parameters in the grid, the simulations encountered an unknown error in sva for 1 out of 10 seeds

#AND THERE IS SOME ISSUES WITH SAVING OUTPUTS, SAVED IT IN A WAY WHERE WE DID NOT CONSIDER THE TOTAL_SV_SD
#total_sv_sd = seq(from = 0,to = 10, by=0.5)
#grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
#grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
#grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0

##FOR ASHG, USE SET5 which has 2631857 SETS OF SIMULATIONS (91.15%)
##save this as set5 
saveRDS(set5, "/sc/arion/projects/mscic1/results/jolie/sim_bug/set5.rds")

set5 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/set5.rds") #3161 sets
dim(set5) #13,159,285 x 13

all_data<-fread("/sc/arion/projects/mscic1/results/jolie/sva_sim_results_with_summary_vpa_stats.csv",data.table=FALSE)

all_data<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_results_with_summary_vpa_stats.csv")
dim(all_data) #43,593,283 x 20
#5 methods x 3 features = 15 rows per sim
#15 x 1sim x 4688 files = 70,320

##cleaning data for downstreat analyses
row_number <- which(all_data$Feature == 'DE' & !is.na(all_data$Ratio_DE_non_DE) & all_data$Ratio_DE_non_DE != 0)[1]

all_data[18450:18478, ]

all_data[18469:18483, ]

# Check if Ratio_DE_non_DE is always the same for DE and non_DE features
discrepancies <- all_data %>%
  filter(Feature %in% c("DE", "non_DE")) %>%
  summarise(
    ratio_de = unique(Ratio_DE_non_DE[Feature == "DE"]),
    ratio_non_de = unique(Ratio_DE_non_DE[Feature == "non_DE"]),
    discrepancy = any(!is.na(ratio_de) & !is.na(ratio_non_de) & ratio_de != ratio_non_de)
  ) %>%
  filter(discrepancy)

# Output the result
if (nrow(discrepancies) == 0) {
  cat("All entries have consistent Ratio_DE_non_DE values for DE and non_DE features.\n")
} else {
  cat("Discrepancies found in the dataset.\n")
}

##the answer is YES!!!


#######Impute Ratio_DE_non_DE values for "all" features based on the corresponding non_DE values
all_data2 <- all_data %>%
  group_by(param.n_sv, param.fraction_degs, param.fraction_sv, param.status_sd, 
           param.sv_sd, param.resid_sd, param.seed, File, Simulation) %>%
  mutate(Ratio_DE_non_DE = ifelse(is.na(Ratio_DE_non_DE) & Feature == "all", 
                                  Ratio_DE_non_DE[Feature == "non_DE"], 
                                  Ratio_DE_non_DE)) %>%
  ungroup()

all_data2 <- as.data.frame(all_data2)
#all_data2 have the imputed ratio of DE_non_DE

##so right now sub2 contains all the simulations, including the one with bugs (at least one extra row)
#so only look at set5 for now
##FOR ASHG, USE SET5 which has 2631857 SETS OF SIMULATIONS (91.15%)
sub <- all_data2[, c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed", "Status_mean", "Status_sd", "Status_median","Feature","Ratio_DE_non_DE")] #remove File and Simulation
dim(sub) #43593283 x 18

table(sub$Feature)
#     all       DE   non_DE 
#14531091 14531096 14531096 

##kinda strange that all has 5 less rows than DE and non_DE!

all_feature <- sub[sub$Feature == "all", ]
DE_feature <- sub[sub$Feature == "DE", ]
non_DE_feature <- sub[sub$Feature == "non_DE", ]


summary(sub$Ratio_DE_non_DE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04167 0.08696 0.26679 0.42857 1.00000 


summary(sub$Status_mean) ##this grabbed info of DE features somehow
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.0     0.0     0.0     0.0     0.0     0.4 1453049


summary(all_feature$Status_mean)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#2.000e-08 1.711e-04 1.416e-03 1.094e-02 8.947e-03 2.040e-01

summary(DE_feature$Status_mean)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.0     0.0     0.0     0.1     0.1     0.4 1453049 (10% missing)

summary(non_DE_feature$Status_mean)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#1.580e-08 5.071e-06 4.931e-05 1.844e-04 3.385e-04 8.815e-04 


summary(DE_feature$Ratio_DE_non_DE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04167 0.08696 0.26679 0.42857 1.00000 
hist_obj <- hist(DE_feature$Ratio_DE_non_DE)
hist_obj$breaks


summary(non_DE_feature$Ratio_DE_non_DE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04167 0.08696 0.26679 0.42857 1.00000
hist_obj <- hist(non_DE_feature$Ratio_DE_non_DE)
hist_obj$breaks


#THERE ARE DUPLICATED ROWS prob due to the total_sv_sd problem?
any(duplicated(sub)) #TRUE
any(duplicated(all_feature))

##since there are duplicates, regular merging will not work
# res3 <- merge(set5, sub2, 
#                    by = c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed"),
#                    all.x = TRUE) #all.x=TRUE ensures that all rows from set5 are kept and rows from sub2 that dont match set5 are excluded
# dim(res3)#14063980 x 17 whereas the EXPECTED IS 13,159,285 x 13

##APPROACH 1 IS TO REMOVE DUPLICATED ROWS AND APPROACH 2 IS TO MATCH THE INDICES DIRECTLY BETWEEN SET5 AND SUB2
#find matching rows in sub2 based on the parameters in set5
sub2_unique <- sub[!duplicated(sub[, c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed")]), ]
dim(sub2_unique) #old: 13185576 x 17, new: 13185581 x 18

all_feature_unique <- all_feature[!duplicated(all_feature[, c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed")]), ]
dim(all_feature_unique) #13185576 x 18
all_feature_unique<-as.data.frame(all_feature_unique)


DE_feature_unique <- DE_feature[!duplicated(DE_feature[, c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed")]), ]
dim(DE_feature_unique) #13185581 x 18

non_DE_feature_unique <- non_DE_feature[!duplicated(non_DE_feature[, c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed")]), ]
dim(non_DE_feature_unique) #13185581 x 18

#merge set5 with the unique feature_unique -- if want to merge all features (all, DE, and non-DE) use sub but idk when you'd want to use sub 
res <- merge(set5, all_feature_unique, 
                   by = c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed"),all.x = TRUE)

res_DE <- merge(set5, DE_feature_unique, 
                   by = c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed"),all.x = TRUE)

res_nonDE <- merge(set5, non_DE_feature_unique, 
                   by = c("nsv","TP","TN","FP","FN","method","param.n_sv", "param.fraction_degs", "param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", "param.seed"),all.x = TRUE)

dim(res) #13159285 x 18
dim(res_DE) #13159285 x 18
dim(res_nonDE) #13159285 x 18
dim(set5)# 13159285 x 13

#check if the number of rows matches set5
nrow(res) == nrow(set5)  #TRUE

saveRDS(res, "/sc/arion/projects/mscic1/results/jolie/sim_bug/all_feature_res_for_set5.rds")

saveRDS(res_DE, "/sc/arion/projects/mscic1/results/jolie/sim_bug/DE_feature_res_for_set5.rds")

saveRDS(res_nonDE, "/sc/arion/projects/mscic1/results/jolie/sim_bug/nonDE_feature_res_for_set5.rds")

############################# ASHG PLOTTING ############################################
########################################################################################
#ASHG START HERE 
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(purrr)
library(furrr)
library(patchwork)

#just read in res and start from here --- (1) need to look at ratio DE/ non DE 
res1 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/all_feature_res_for_set5.rds")

table(res1$param.seed) #the seeds are not all the same - maybe due to the bugs of set6, set7, etc. issue

#calculate averages grouped by the specified parameters
res1_avg <- res1 %>%
  group_by(param.n_sv, param.fraction_degs, param.fraction_sv, 
           param.status_sd, param.sv_sd, param.resid_sd, method, Feature, Ratio_DE_non_DE) %>%
  summarise(
    seeds = list(unique(param.seed)), # List all unique seeds
    nsv = mean(nsv, na.rm = TRUE),
    TP = mean(TP, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    Status_mean = mean(Status_mean, na.rm = TRUE),
    Status_sd = mean(Status_sd, na.rm = TRUE),
    Status_median = mean(Status_median, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(res1_avg, file = "/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed.rds")

##to plot, start here
res <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed.rds") 

res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

res$true_degs<-res$param.fraction_degs * 20000
res$all_positive=res$TP + res$FP
res$all_negative=res$TN + res$FN
res$accuracy=(res$TP + res$TN)/(res$TP + res$TN + res$FP + res$FN)
res$precision=res$TP/(res$TP + res$FP) 
res$specificity=res$TN/(res$TN + res$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
#res$falsePositiveRate=res$FP/(res$FP + res$TN) #type 1 error = false positive #SVA commits this!
res$recall=res$TP/(res$TP + res$FN) #true positive rate; measure of a test's ability to correctly identify true positives
res$negativePredictiveValue=res$TN/(res$TN + res$FN)
#res$falseNegativeRate=res$FN/(res$FN + res$TP) #type 2 error = false negative
res$totalDEGsignal=res$param.fraction_degs*res$param.status_sd

dim(res)
#1318660 x 31


######################################################################################
res_DE <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/DE_feature_res_for_set5.rds")

res_DE <- res_DE %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

res_DE$accuracy=(res_DE$TP + res_DE$TN)/(res_DE$TP + res_DE$TN + res_DE$FP + res_DE$FN)
res_DE$precision=res_DE$TP/(res_DE$TP + res_DE$FP) 
res_DE$specificity=res_DE$TN/(res_DE$TN + res_DE$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
#res_DE$falsePositiveRate=res_DE$FP/(res_DE$FP + res_DE$TN) #type 1 error = false positive #SVA commits this!
res_DE$recall=res_DE$TP/(res_DE$TP + res_DE$FN) #true positive rate; measure of a test's ability to correctly identify true positives
res_DE$negativePredictiveValue=res_DE$TN/(res_DE$TN + res_DE$FN)
#res_DE$falseNegativeRate=res_DE$FN/(res_DE$FN + res_DE$TP) #type 2 error = false negative
res_DE$totalDEGsignal=res_DE$param.fraction_degs*res_DE$param.status_sd

##use this approach to find identity of the simulation set!
summary(res_DE$FPR)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.001771 0.011750 0.254827 0.649348 0.997250

low_FPR <- res_DE[res_DE$FPR == 0, ] #11.9% of the time

table(low_FPR$Ratio_DE_non_DE)
#                 0 0.0204081632653061 0.0416666666666667 0.0638297872340425 
#            820063             126491             106377              95990 
#0.0869565217391304  0.111111111111111               0.25  0.428571428571429 
#             88859              83662              69452              62560 
# 0.666666666666667                  1 
#             58158              55390 

high_sn_in_low_FPR <- low_FPR[low_FPR$Ratio_DE_non_DE == 0.0204081632653061, ] #


#check their FPR


########################################################################################
res_non_DE <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/nonDE_feature_res_for_set5.rds")

res_non_DE <- res_non_DE %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

res_non_DE$accuracy=(res_non_DE$TP + res_non_DE$TN)/(res_non_DE$TP + res_non_DE$TN + res_non_DE$FP + res_non_DE$FN)
res_non_DE$precision=res_non_DE$TP/(res_non_DE$TP + res_non_DE$FP) 
res_non_DE$specificity=res_non_DE$TN/(res_non_DE$TN + res_non_DE$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
#res_non_DE$falsePositiveRate=res_non_DE$FP/(res_non_DE$FP + res_non_DE$TN) #type 1 error = false positive #SVA commits this!
res_non_DE$recall=res_non_DE$TP/(res_non_DE$TP + res_non_DE$FN) #true positive rate; measure of a test's ability to correctly identify true positives
res_non_DE$negativePredictiveValue=res_non_DE$TN/(res_non_DE$TN + res_non_DE$FN)
#res_non_DE$falseNegativeRate=res_non_DE$FN/(res_non_DE$FN + res_non_DE$TP) #type 2 error = false negative
res_non_DE$totalDEGsignal=res_non_DE$param.fraction_degs*res_non_DE$param.status_sd


###
res_grouped <- res %>%
  group_by(param.n_sv, param.fraction_degs, param.fraction_sv, param.status_sd, 
           param.sv_sd, param.resid_sd, param.seed, Feature) %>%
  summarize(
   	nsv = list(nsv),
    method = list(unique(method)),  
   	TP = list(TP),
    TN = list(TN),
    FP = list(FP),
    FN = list(FN),
    TPR = list(TPR),
    TNR = list(TNR),
    FPR = list(FPR),
    FNR = list(FNR),
    Status_mean = list(Status_mean),
    Status_sd = list(Status_sd),
    Status_median = list(Status_median),
    .groups = 'drop'
  )

# View the grouped results
aaa<-as.data.frame(res_grouped)


############################## TABLE 1 ################################################
#calculate mean rates (TPR, TNR, FPR, FNR) by method
mean_rates <- res %>%
  group_by(method) %>%
  summarize(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_TNR = mean(TNR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_FNR = mean(FNR, na.rm = TRUE)
  )
mean_rates

#chord flip - add color?
############################## FIGURE 1 ################################################
############Violin Plot showing TPR distribution ###################################
########################################################################################

#colors are accessible: https://venngage.com/tools/accessible-color-palette-generator
###only run this once
res$method <- factor(res$method, levels = c("oracle", "none", "known_nsv", "be", "leek"),labels = c('Oracle', 'No SVs','Known Number of SVs','SVA "BE"','SVA "leek"'))

# custom_colors <- c('SVA "BE"' = "#9b8bf4", 
#                    'SVA "leek"' = "#f8b8d0", 
#                    'Known Number of SVs' = "#f194b8", 
#                    'No SVs' = "#d9e4ff", 
#                    'Oracle' = "#8babf1")

###only run this once

##################TPR##########################################
#color version
# p1 <- ggplot(res, aes(x = method, y = TPR, fill = method)) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#   		legend.position = "none",
#   		axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) + 
#   labs(x = "Method", y = "True Positive Rate")
# #2: Removed 1315820 rows containing non-finite outside the scale range
# #(`stat_boxplot()`). 

#colorless version
p1 <- ggplot(res, aes(x = method, y = TPR)) +  
  geom_violin(trim = FALSE, color = "black", fill = NA) +  
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), color = "black", fill = NA) +  
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) + 
  scale_y_continuous(
    limits = c(0, 1),  # Set the range of the y-axis
    breaks = seq(0, 1, by = 0.25)  # Set the labels at 0.25 intervals
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  labs(x = "Method", y = "True Positive Rate")


ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/p1.pdf", plot = p1, width = 8, height = 5) 

##################TNR##########################################
#color version
# p2 <- ggplot(res, aes(x = method, y = TNR, fill = method)) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#   		legend.position = "none",
#   		axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#         ) +
#         scale_y_continuous(limits = c(0, 1)) + 
#   labs(x = "Method", y = "True Negative Rate")

#colorless version
p2 <- ggplot(res, aes(x = method, y = TNR)) +  # Removed 'fill = method'
  geom_violin(trim = FALSE, color = "black", fill = NA) +  # Transparent fill with black border for violin
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), color = "black", fill = NA) +  # Transparent fill with black border for boxplot
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) + 
  labs(x = "Method", y = "True Negative Rate")

##################FPR##########################################
#color version
# p3 <- ggplot(res, aes(x = method, y = FPR, fill = method)) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#   		legend.position = "none") +
#   labs(x= "Method", y = "False Positive Rate")

#colorless version
p3 <- ggplot(res, aes(x = method, y = FPR)) +  # Removed 'fill = method'
  geom_violin(trim = FALSE, color = "black", fill = NA) +  # Transparent fill with black border for violin
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), color = "black", fill = NA) +  # Transparent fill with black border for boxplot
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) + 
  scale_y_continuous(
    limits = c(0, 1),  # Set the range of the y-axis
    breaks = seq(0, 1, by = 0.25)  # Set the labels at 0.25 intervals
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Method", y = "False Positive Rate")

# ##tried histogram or ridge plot to display distribution
# library(ggridges)

# ##HISTOGRAM
# p3 <- ggplot(res, aes(x = FPR, fill = method)) +
#   geom_histogram(bins = 30, alpha = 0.7, position = "identity") + 
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   facet_wrap(~ method, scales = "free") + # Separate histograms by method
#   labs(x= "False Positive Rate", y = "Count")

# ##RIDGE PLOT
# p3 <- ggplot(res, aes(x = TNR, y = method, fill = method)) +
#   geom_density_ridges(alpha = 0.7) + # Ridge plot
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#       legend.position = "none") +
#   labs(x= "True Negative Rate", y = "Method")

# #ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/p4.png", plot = p3, width = 8, height = 5) 
# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/ridge_plot_TNR_distrubution.pdf", plot = p3, width = 8, height = 5) 

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/ridge_plot_FPR_distrubution.pdf", plot = p3, width = 8, height = 5) 

# ##VIOLIN PLOT ALONE -- 
# p3 <- ggplot(res, aes(x = method, y = FPR, fill = method)) +
#   geom_violin(trim = FALSE, alpha = 0.7) + # Violin plot
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         legend.position = "none") +
#   labs(x = "Method", y = "False Positive Rate") +
#   coord_flip() # Flip coordinates if you prefer horizontal violins

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/violin_plot_FPR_distribution.pdf", plot = p3, width = 8, height = 5)

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/violin_plot_TNR_distribution.pdf", plot = p3, width = 8, height = 5)



##################FNR##########################################
#color version
# p4 <- ggplot(res, aes(x = method, y = FNR, fill = method)) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
#   theme_bw() +
#   scale_fill_manual(values = custom_colors) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#   		legend.position = "none",
#   		axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#   		) + 
#   labs(x = "Method", y = "False Negative Rate")

#colorless version
p4 <- ggplot(res, aes(x = method, y = FNR)) +  # Removed 'fill = method'
  geom_violin(trim = FALSE, color = "black", fill = NA) +  # Transparent fill with black border for violin
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), color = "black", fill = NA) +  # Transparent fill with black border for boxplot
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 25)) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  labs(x = "Method", y = "False Negative Rate")

# Combine the plots
combo <- p1 + p2 + p3 + p4 + plot_annotation(tag_levels = 'A')

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions2.pdf", plot = combo, width = 8, height = 5) 

# Warning messages:
# 1: Removed 131880 rows containing non-finite outside the scale range
# (`stat_ydensity()`). 
# 2: Removed 131880 rows containing non-finite outside the scale range
# (`stat_boxplot()`). 
# 3: Removed 262 rows containing missing values or values outside the scale range
# (`geom_violin()`). 
# 4: Removed 131880 rows containing non-finite outside the scale range
# (`stat_ydensity()`). 
# 5: Removed 131880 rows containing non-finite outside the scale range
# (`stat_boxplot()`). 

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_violin_plots_of_rate_distributions2.png", plot = combo, width = 8, height = 5) 


summary_by_method <- by(res$TNR, res$method, summary); summary_by_method


############################## FIGURE 2A ################################################
##### Comparison of SV Detection Rates Across Simulation Methods ######################
#######################################################################################
res$method <- factor(res$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))



plot <- ggplot(res, aes(y = nsv, x = as.factor(param.n_sv))) + 
  geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs simulated") + 
  ylab("Number of SVs used") + 
  ggtitle("Comparison of SV Used by SV Simulated Across Simulation Methods") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/nsv_used_by_nsv_simulated.pdf", plot = plot, width = 8, height = 5) 

#plot all negative 
#pick 1 set of sim per SV and check TP

##Comparison of SV Simulated by FPR Across Simulation Methods
res_filtered <- res %>% 
  filter(method == "leek")

plot <- ggplot(res_filtered, aes(y = FPR, x = as.factor(param.n_sv), fill = as.factor(true_degs))) + 
#plot <- ggplot(res, aes(y = TP, x = as.factor(nsv))) + 
  #geom_violin(trim = FALSE) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs simulated") + 
  ylab("False Positive Rate") + 
  ggtitle('Comparison of SV Used by False Positive Rate in SVA "leek"') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/leek_FPR_by_sva_simulated_color_true_degs.pdf", plot = plot, width = 14, height = 5) 


library(ggplot2)
library(dplyr)

# Filter for only the "be" method
res_filtered <- res %>% 
  filter(method == "be")

# Generate the plot
plot <- ggplot(res_filtered, aes(y = TP, x = as.factor(param.n_sv), fill = as.factor(param.fraction_degs))) + 
  geom_boxplot(outlier.shape = NA) +  # Removes outliers from the boxplot
  xlab("Number of SVs Simulated") + 
  ylab("True Positives") + 
  ggtitle("Comparison of SV Used by True Positives (SVA 'BE' Method Only)") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_brewer(palette = "Set3", name = "Fraction DEGs")

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/TP_by_sva_simulated_be_no_outlier.pdf", plot = plot, width = 8, height = 5) 



#summary statistics
tp_summary <- res %>%
  group_by(method, param.n_sv) %>%
  summarise(
    Mean_TP = mean(TP),
    Variance_TP = var(TP),
    .groups = "drop"
  )

tp_summary <- as.data.frame(tp_summary)

# Plot mean TP by method and number of SVs
plot<-ggplot(tp_summary, aes(x = param.n_sv, y = Mean_TP, color = method, group = method)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Mean TP Across Methods and SVs",
    x = "Number of SVs",
    y = "Mean True Positives"
  ) +
  theme_minimal()
ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/meanTP_by_numbers_of_SV.pdf", plot = plot, width = 8, height = 5) 

summary(res$TP[res$param.n_sv == 0])
summary(res$TP[res$param.n_sv == 60])

plot<-ggplot(res, aes(x = param.fraction_degs, y = TP)) +
  geom_boxplot() +
  facet_wrap(~ method) +
  labs(title = "TP Distribution by Fraction of DE Genes", x = "Fraction of DE Genes", y = "True Positives")

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/TP_by_fraction_DEGs.pdf", plot = plot, width = 8, height = 5) 


plot<-ggplot(res, aes(x = param.status_sd, y = TP)) +
  geom_boxplot() +
  facet_wrap(~ method) +
  labs(title = "TP Distribution by Status Effect Size", x = "Status Effect Size", y = "True Positives")
ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/TP_by_status_sd.pdf", plot = plot, width = 8, height = 5) 


#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/FPR_by_sva_simulated2.pdf", plot = plot, width = 8, height = 5) 


plot <- ggplot(res, aes(y = TP, x = as.factor(nsv))) + 
  geom_violin(trim = FALSE) +
  #geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs simulated") + 
  ylab("True Positives") + 
  ggtitle("Comparison of SV Used by True Positives Across Simulation Methods") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/TP_by_sva_used_violin.pdf", plot = plot, width = 8, height = 5) 




##Comparison of SV Simulated by Total Positives Across Simulation Methods
all_pos <- ggplot(res, aes(y = all_positive, x = as.factor(param.n_sv))) + 
  geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs simulated") + 
  ylab("Total Positives") + 
  ggtitle("Comparison of SV Simulated by Total Positives Across Simulation Methods") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/all_pos_by_sva_simulated.pdf", plot = all_pos, width = 8, height = 5) 


##Comparison of SV Used by Total Positives Across Simulation Methods - BOXPLOT
all_pos <- ggplot(res, aes(y = all_positive, x = as.factor(nsv))) + 
  geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs used") + 
  ylab("Total Positives") + 
  ggtitle("Comparison of SV Used by Total Positives Across Simulation Methods") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/all_pos_by_sva_used.pdf", plot = all_pos, width = 8, height = 5) 

##Comparison of SV Used by Total Positives Across Simulation Methods - SCATTERPLOT
all_pos <- ggplot(res, aes(y = all_positive, x = as.factor(nsv))) + 
  geom_point(aes(color = method), alpha = 0.6, position = position_jitter(width = 0.2, height = 0)) + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs used") + 
  ylab("Total Positives") + 
  ggtitle("Comparison of SV Used by Total Positives Across Simulation Methods") + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )
ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/all_pos_by_sva_used_scatterplot.pdf", plot = all_pos, width = 8, height = 5) 




##Comparison of 
all_pos <- ggplot(res, aes(y = all_positive, x = as.factor(param.n_sv))) + 
  geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("Number of SVs simul") + 
  ylab("Total Positives") + 
  ggtitle("Comparison of SV Simulated by Total Positives Across Simulation Methods") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/all_pos_by_sva_simulated.pdf", plot = all_pos, width = 8, height = 5) 





plot <- ggplot(res, aes(y = FPR, x = as.factor(all_positive))) + 
  geom_boxplot() + 
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) + 
  xlab("All Positives") + 
  ylab("FPR") + 
  ggtitle("Title") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/FPR_vs_all_positive.pdf", plot = plot, width = 8, height = 5) 

############################## FIGURE 2AA ################################################
##### Comparison of SV Detection Rates Across Simulation Methods ######################
########################################################################################
##need to reorder oracle etc. 
res$method <- factor(res$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))

p <- ggplot(res, aes(x = as.factor(all_positive), y = FPR, fill = method, color = as.factor(nsv))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "All Positive (Categorical)",
    y = "False Positive Rate (FPR)",
    title = "Violin Plot of FPR by All Positive and Method"
  ) +
  facet_wrap(~method, 
             labeller = labeller(method = c(
               oracle = "Oracle", 
               none = "No SVs", 
               known_nsv = "Known Number of SVs", 
               be = 'SVA "BE"', 
               leek = 'SVA "leek"'))) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/violin_plot_optimized.pdf", 
       plot = p, width = 12, height = 6)


############################## FIGURE 2B ################################################
##### NUMBER OF SVs USED and FPR Across Simulation Methods ######################
########################################################################################
res_subs <- res[res$nsv >= 0 & res$nsv <= 20, ]
res_subs2 <- res[res$nsv >= 0 & res$nsv <= 10, ]

#subset simulated SVs to be btw 0 and 30 and all positive > 100
res_subs3 <- res[res$param.n_sv >= 0 & res$param.n_sv <= 30, ] 
res_subs3 <- res_subs3[res_subs3$all_positive >100, ]

plot <- ggplot(res_subs3, aes(y=FPR, x=as.factor(nsv), color=as.factor(param.n_sv))) + 
    geom_boxplot() + 
    facet_wrap(~method, labeller = labeller(method = c(oracle = "Oracle", none = "No SVs", known_nsv = "Known Number of SVs", be = "SVA \"be\"", leek = "SVA \"leek\"")), ncol = 3) + 
    xlab("Number of SVs used") + 
    ylab("False Positive Rate") + 
    ggtitle("Title") + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/a.pdf", plot = plot, width = 8, height = 5) 

#scatterplot of     
sub1<-res_subs3[res_subs3$method == "be",]
sub2<-res_subs3[res_subs3$method == "leek",]

dput(colnames(sub1))

col<-c("param.n_sv", "param.fraction_degs", 
"param.fraction_sv", "param.status_sd", "param.sv_sd", "param.resid_sd", 
"param.seed")
colnames(sub1)[!(colnames(sub1) %in% col)] = paste0(colnames(sub1)[!(colnames(sub1) %in% col)], "_be")

colnames(sub2)[!(colnames(sub2) %in% col)] = paste0(colnames(sub2)[!(colnames(sub2) %in% col)], "_leek")

mer<-merge(sub1, sub2, by=col)

plot <- ggplot(mer, aes(y=nsv_leek, x=nsv_be)) + 
    geom_hex(bins = 30) +
    scale_fill_viridis_c(trans="log10") +
    facet_wrap(~param.n_sv, ncol = 2) + 
    xlab("nsv be") + 
    ylab("nsv leek") + 
    coord_fixed() +
    ggtitle("Comparison of False Positive Rates by Numbers of SVs Used Across SV Methods") + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/nsv_leek_vs_be_hexplot3.pdf", plot = plot, width = 8, height = 5) 

table(table(mer$all_positive_be[mer$nsv_be>mer$param.n_sv]))

# plot <- ggplot(mer, aes(y=FPR_leek, x=nsv_be, color=FPR_be)) + 
#     geom_point() + 
#     facet_wrap(~param.n_sv, ncol = 3) + 
#     xlab("nsv be") + 
#     ylab("FPR leek") + 
#     ggtitle("Comparison of False Positive Rates by Numbers of SVs Used Across SV Methods") + 
#     theme_bw() + 
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/FPR_leek_vs_be3.pdf", plot = plot, width = 8, height = 5) 


hexbin_for_TPR_vs_FPR <- ggplot(res, aes(x = FPR, y = TPR)) +
    geom_hex(bins = 15) +  # Adjust `bins` based on desired granularity
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",
        y = "True Positive Rate",
        fill = "Density"
    ) + 
    #ggtitle("Density of True Positive Rate \nvs False Positive Rates by Method") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        known_nsv = "Known Number of SVs",
        be = "SVA 'be'",
        leek = "SVA 'leek'"
    )), ncol = 1)


############################# VPA PLOTTING #############################################
########################################################################################

###plotting code is in SVA Part 2 DEA.sh

#code to run VPA: sva_sim_parallel_for_vp.R
form <- as.formula("~ Status + Batch_ID")
remaining_svs <- setdiff(sv_columns, sv_columns[cols_to_drop])
if (length(remaining_svs) > 0) {
  form <- as.formula(paste("~ Status + Batch_ID +", paste(remaining_svs, collapse = "+")))
}
varPart <- fitExtractVarPartModel(simData$E, form, simData$targets)
vp <- sortCols(varPart)
#Check if "Status" column exists in vp
if ("Status" %in% colnames(vp)) {
  vp_res <- data.frame(Status = vp[, "Status"], feature = rownames(vp))
} else {
  stop("Status column not found in variance partitioning results")
} #just adding this here to close the loop 

#code to combine VPA results: job_submission_sva_sim_vp.sh


############################# PLOT 1 - OG ###################################################
#############################################################################################
##PLOT1: FPR by Precision across methods by Fraction DEGs (sanity check, made this before)

plot<-ggplot(res, aes(x = falsePositiveRate, y = precision, color = param.fraction_degs)) +
  geom_point() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
  	legend.position = "none")  + 
   xlab("False Positive Rate") + 
   ylab("Precision") + 
   ggtitle("Precision and False Positive Rates Across Methods by Fraction of DEGs") + facet_wrap(~method, labeller = labeller(method = c(oracle = "Oracle",be = "SVA \"be\"", leek = "SVA \"leek\"", known_nsv = "Known Number of SVs", none = "No SVs")))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Precision_FPR_plot3.png", plot = plot, width = 8, height = 5)


############################# PLOT 1 A ###################################################
########################################################################################
##PLOT1: FPR by Precision across methods by STATUS MEAN, RATIO DE TO NON DE, TOTAL DEG SIGNAL, AND RECALL
#(change color to Status_mean, Ratio_DE_non_DE, and total DEGsignal) 

#plot<-ggplot(res, aes(x = FPR, y = precision, color = Ratio_DE_non_DE)) +
#plot<-ggplot(res, aes(x = FPR, y = precision, color = TPR)) +

res <- res[sample(nrow(res)), ] #shuffle row randomly

res$method <- factor(res$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))

precision_vs_FPR <- ggplot(res, aes(x = FPR, y = precision, color = Ratio_DE_non_DE)) +
    geom_point_rast(alpha = 0.5) +  # Rasterized points
    #geom_point(alpha = 0.5) + 
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",          
        y = "Precision",           
        color = "Ratio of DEGs to non-DEGs" 
    ) + 
    #ggtitle("Precision vs False Positive Rates \nAcross Methods by Ratio of DEGs to non-DEGs") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        known_nsv = "Known Number of SVs",
        be = "SVA 'BE'",
        leek = "SVA 'leek'"
    )), ncol = 1)

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Precision_FPR_by_ratio3.png", plot = plot, width = 11, height = 5)

##PLOT FPR by TPR across methods by STATUS MEAN, RATIO DE TO NON DE, TOTAL DEG SIGNAL, AND RECALL
#(change color to Status_mean, Ratio_DE_non_DE, and total DEGsignal) 


TPR_vs_FPR <- ggplot(res, aes(x = FPR, y = TPR, color = Ratio_DE_non_DE)) +
    geom_point_rast(alpha = 0.5) +  # Rasterized points
    #geom_point(alpha = 0.5) + 
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",          
        y = "True Positive Rate",           
        color = "Ratio of DEGs to non-DEGs" 
    ) + 
    #ggtitle("True Positive Rate vs False Positive Rates \nAcross Methods by Ratio of DEGs to non-DEGs") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        known_nsv = "Known Number of SVs",
        be = "SVA 'BE'",
        leek = "SVA 'leek'"
    )), ncol = 1)

library(hexbin)

hexbin_for_TPR_vs_FPR <- ggplot(res, aes(x = FPR, y = TPR)) +
    geom_hex(bins = 15) +  # Adjust `bins` based on desired granularity
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",
        y = "True Positive Rate",
        fill = "Density"
    ) + 
    #ggtitle("Density of True Positive Rate \nvs False Positive Rates by Method") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        known_nsv = "Known Number of SVs",
        be = "SVA 'BE'",
        leek = "SVA 'leek'"
    )), ncol = 1)

#TPR_by_FPR_with_hexbin.png has bin = 30
#TPR_by_FPR_with_hexbin2.png has bin = 100
#TPR_by_FPR_with_hexbin3.png has bin = 10

library(patchwork)

combo_sim <-TPR_vs_FPR + hexbin_for_TPR_vs_FPR + precision_vs_FPR + plot_annotation(tag_levels = 'A')

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_sim_FPR_TPR_Precision_rasterized_SUPP.pdf", plot = combo_sim, width = 9, height = 11)

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_sim_plot2.png", plot = combo_sim, width = 9, height = 11)

####SUBSETTING THE PLOTS FOR ONLY ORACLE, NONE, AND BE
filtered_res <- res %>% 
  filter(method %in% c("oracle", "none", "be"))

filtered_res$method <- factor(filtered_res$method, levels = c("oracle", "none", "be"))

precision_vs_FPR <- ggplot(filtered_res, aes(x = FPR, y = precision, color = Ratio_DE_non_DE)) +
    geom_point_rast(alpha = 0.5) +  # Rasterized points
    #geom_point(alpha = 0.5) + 
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",          
        y = "Precision",           
        color = "Ratio of DEGs to non-DEGs" 
    ) + 
    #ggtitle("Precision vs False Positive Rates \nAcross Methods by Ratio of DEGs to non-DEGs") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        be = "SVA 'BE'"
    )), ncol = 1)

TPR_vs_FPR <- ggplot(filtered_res, aes(x = FPR, y = TPR, color = Ratio_DE_non_DE)) +
    geom_point_rast(alpha = 0.5) +  # Rasterized points
    #geom_point(alpha = 0.5) + 
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",          
        y = "True Positive Rate",           
        color = "Ratio of DEGs to non-DEGs" 
    ) + 
    #ggtitle("True Positive Rate vs False Positive Rates \nAcross Methods by Ratio of DEGs to non-DEGs") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        be = "SVA 'BE'"
    )), ncol = 1)

hexbin_for_TPR_vs_FPR <- ggplot(filtered_res, aes(x = FPR, y = TPR)) +
    geom_hex(bins = 15) +  # Adjust `bins` based on desired granularity
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    ) + 
    labs(
        x = "False Positive Rate",
        y = "True Positive Rate",
        fill = "Density"
    ) + 
    #ggtitle("Density of True Positive Rate \nvs False Positive Rates by Method") + 
    facet_wrap(~method, labeller = as_labeller(c(
        oracle = "Oracle",
        none = "No SVs",
        be = "SVA 'BE'"
    )), ncol = 1)

combo_sim2 <-TPR_vs_FPR + hexbin_for_TPR_vs_FPR + precision_vs_FPR + plot_annotation(tag_levels = 'A')

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_sim_FPR_TPR_Precision_rasterized.pdf", plot = combo_sim2, width = 8, height = 8)




aggregate(TPR ~ method, data = res, summary)
#      method   TPR.Min. TPR.1st Qu. TPR.Median   TPR.Mean TPR.3rd Qu.   TPR.Max.
# 1        be 0.04816667  0.84400000 0.89500000 0.87054902  0.93175000 1.00000000
# 2 known_nsv 0.00012500  0.86437500 0.92120000 0.88895268  0.95625000 1.00000000
# 3      leek 0.02000000  0.83583333 0.90700000 0.86537931  0.94675000 1.00000000
# 4      none 0.00000000  0.03000000 0.37525000 0.41187910  0.77250000 0.99500000
# 5    oracle 0.35250000  0.82500000 0.89950000 0.86164549  0.94237500 0.99500000

aggregate(FPR ~ method, data = res, summary)
#      method     FPR.Min.  FPR.1st Qu.   FPR.Median     FPR.Mean  FPR.3rd Qu.	 FPR.Max.
# 1        be 0.0000000000 0.0090217391 0.4897959184 0.4327260359 0.8250000000	0.9957500000
# 2 known_nsv 0.0000000000 0.0060000000 0.2825625000 0.4147207567 0.8707000000	0.9972500000
# 3      leek 0.0000000000 0.0049444444 0.2865625000 0.4150413730 0.8694285714	0.9972500000
# 4      none 0.0000000000 0.0000000000 0.0004081633 0.0040200999 0.0038043478	0.2003000000
# 5    oracle 0.0000000000 0.0017187500 0.0040217391 0.0076285245 0.0135000000	0.0309000000

############################# PLOT 3 ###################################################
########################################################################################
##PLOT3: FPR by Status Mean across methods by Method (saved)
#higher status_mean = larger overall average effect or signal associated with the "Status" variable (e.g., a condition or treatment effect) across the samples.


##stronger signal is associated with lower FPR
##no SVs and oracle --> FPR is always low 

#In variance partition analysis, the goal is to determine the proportion of total variance in gene expression data explained by different variables (e.g., batch, individual, condition). A higher Status_mean implies that the condition (status) contributes significantly to the variability observed in the data, potentially increasing the proportion of variance attributed to this variable in the analysis.

#SVA "be" and "leek" ---
#Both methods show a clear trend where the FPR decreases as Status_mean increases. This suggests that these methods are more prone to false positives when the condition has a low impact on the variance (low Status_mean), but they become more reliable as the condition's contribution to variance increases.

#Known Number of SVs ---
#This method also shows a decrease in FPR as Status_mean increases, but it appears to have a higher FPR in general compared to SVA "leek" and SVA "be" at lower Status_mean values.

# No SVs ---
# This method exhibits relatively low FPR values regardless of Status_mean, but the distribution is more dispersed compared to other methods. This suggests that not accounting for any surrogate variables might lead to inconsistent results.

# Oracle ---
# This method shows a consistently low FPR across all levels of Status_mean, which is expected as the "Oracle" method likely has prior knowledge or an ideal adjustment for the variance structure. It serves as a benchmark for how other methods should perform.

#The plot suggests that methods like SVA "be" and SVA "leek" are sensitive to the contribution of Status_mean to the total variance. They perform better (lower FPR) as Status_mean increases, indicating their ability to adjust appropriately when the condition has a stronger effect.

##THIS IS INFORMATIVE !!!!!
##this plot is different compared to sub (have all, non-DE and DE)
plot <- ggplot(res, aes(x = Status_mean, y = FPR, color = Ratio_DE_non_DE)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Status Mean on FPR Across Methods") +
  facet_wrap(~method)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR_with_ratio2.png", plot = plot, width = 20, height = 16)


plot <- ggplot(res_DE, aes(x = Status_mean, y = FPR, color = Ratio_DE_non_DE)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Status Mean on FPR Across Methods") +
  facet_wrap(~method)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR_with_ratio_res_DE.png", plot = plot, width = 20, height = 16)


plot <- ggplot(res_non_DE, aes(x = Status_mean, y = FPR, color = Ratio_DE_non_DE)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Status Mean on FPR Across Methods") +
  facet_wrap(~method)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR_with_ratio_res_non_DE.png", plot = plot, width = 20, height = 16)


############################# PLOT 3 ###################################################
########################################################################################


plot <- ggplot(res, aes(x = Status_mean, y = FPR, color = Ratio_DE_non_DE)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Status Mean on FPR Across Methods") +
  facet_wrap(~method)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR_with_ratio.png", plot = plot, width = 20, height = 16)







plot <- ggplot(res, aes(x = Status_mean, y = FPR, color = method)) +
  geom_point() + 
  #stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.4) +
  #geom_density_2d(bandwidth = c(0.5, 0.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Status Mean on FPR Across Methods") +
  facet_wrap(~method); plot

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR_with_plot_density.png", plot = plot, width = 20, height = 16)
#Warning message:
#Removed 55 rows containing missing values or values outside the scale range (`geom_point()`). 


###making BINS for low, medium, and high status mean based on 33rd and 67th percentiles

#define quantile-based bins for Status_mean with NA removal
res <- res %>%
  mutate(Status_bin = cut(Status_mean, 
                          breaks = quantile(Status_mean, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE),
                          labels = c("Low", "Medium", "High"),
                          include.lowest = TRUE))

label_data <- data.frame(
  Status_bin = c("Low", "Medium", "High"),
  falsePositiveRate = c(0.8, 0.8, 0.8),  # Adjust Y position if needed for better visibility
  label = c("Low", "Medium", "High")
)

plot <- ggplot(res, aes(x = Status_bin, y = falsePositiveRate, color = method)) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Jitter to show distribution within bins
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean (Binned)") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Stratified Status Mean on FPR Across Methods") +
  facet_wrap(~method) +
  geom_text(data = label_data, aes(x = Status_bin, y = falsePositiveRate, label = label), 
            size = 6, color = "black", fontface = "bold", vjust = -0.5); plot

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_bins_by_FPR.png", plot = plot, width = 20, height = 16)

###plot without labels
plot2 <- ggplot(res, aes(x = Status_bin, y = falsePositiveRate, color = method)) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Jitter to show distribution within bins
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  ) +
  xlab("Status Mean (Binned)") + 
  ylab("False Positive Rate") + 
  ggtitle("Effect of Stratified Status Mean on FPR Across Methods") +
  facet_wrap(~method); plot2

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_bins_by_FPR_without_labels.png", plot = plot, width = 20, height = 16)


##PLOT 4: Scatter Plot of Status_mean vs. TotalDEGsignal by Method, screenshotted
#This plot examines the relationship between Status_mean (indicating the effect size) and the total signal of differentially expressed genes (DEGs) across methods
plot2<-ggplot(res, aes(x = Status_mean, y = totalDEGsignal, color = method)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  xlab("Status Mean") +
  ylab("Total DEG Signal") +
  ggtitle("Effect of Status Mean on Total DEG Signal Across Methods") +
  facet_wrap(~method) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )
plot2

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/Status_mean_by_FPR.png", plot = plot, width = 20, height = 16)

##PLOT 5: Violin Plot of Precision by Method
ggplot(res, aes(x = method, y = precision, fill = method)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  xlab("Method") +
  ylab("Precision") +
  ggtitle("Distribution of Precision Across Methods") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )


##PLOT6: 
library(reshape2)

# Bin Status_mean
res$Status_mean_bin <- cut(res$Status_mean, breaks = 10)

# Create heatmap data
heatmap_data <- aggregate(falsePositiveRate ~ method + Status_mean_bin, data = res, mean)

ggplot(heatmap_data, aes(x = Status_mean_bin, y = method, fill = falsePositiveRate)) +
  geom_tile() +
  theme_bw() +
  xlab("Status Mean Bin") +
  ylab("Method") +
  ggtitle("Heatmap of FPR Across Methods and Status Mean Bins") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )



####SUBSET RES BY AVERAGING ACROSS SEEDS ###################################################
############################################################################################

param_columns <- c("param.n_sv", "param.fraction_degs", 
                   "param.fraction_sv", "param.status_sd", 
                   "param.sv_sd", "param.resid_sd")

# Group by these columns and calculate the mean across all other columns
res_avg <- res %>%
  group_by(across(all_of(param_columns))) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# Display the resulting averaged dataframe
head(res_avg)

res_avg <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/res_avg_across_seed.rds") 

##########################################################################################
#Future Plots to Explore 
g = ggplot(res[sample(1:nrow(res)),],aes(x = param.n_sv, y = FPR, color = TPR)) + geom_point()
ggsave(g, filename="~/www/plots/sva/tmp.png")

cor(res$TPR[res$method=="known_nsv"],res$FPR[res$method=="known_nsv"],method="spearman")
cor(res$TPR[res$method=="known_nsv"],res$FPR[res$method=="known_nsv"],method="spearman",use="p[1] 0.5865317ete.obs")
g = ggplot(res[sample(1:nrow(res)),],aes(x = TPR, y = FPR, color = param.n_sv)) + geom_point()

ggsave(g, filename="~/www/plots/sva/tmp2.png")

g = ggplot(res,aes(x = param.n_sv, y = FPR, color = TPR)) + geom_point()
ggsave(g, "~/www/plots/sva/tmp.png")
















