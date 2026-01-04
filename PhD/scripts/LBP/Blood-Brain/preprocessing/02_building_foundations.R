library(data.table)
library(dplyr)
library(ggplot2)
library(variancePartition)
library(tidyr)
library(stringr)
library(tibble)
library(limma)
library(glmnet)
library(reshape2)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
#head(blood_only_metadata)
head(dictionary)
#v_blood$E and v_brain$E

##some demographic data
table(dictionary$mymet_tissue_brain)
summary(blood_only_metadata$mymet_age_blood)

df_long <- reshape2::melt(v_blood$E)
summary(df_long$value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.9404  0.0213  2.1121  2.0118  4.1664 19.9860 

df_long_highest_first <- df_long %>%
  arrange(desc(value))
head(df_long_highest_first)
#                 Var1             Var2    value
# 1  ENSG00000276168.1 LBPSEMA4BLOOD710 19.98597
# 2  ENSG00000274012.1 LBPSEMA4BLOOD710 19.94300
# 3  ENSG00000276168.1 LBPSEMA4BLOOD668 19.47055
# 4  ENSG00000274012.1 LBPSEMA4BLOOD668 19.34698
# 5  ENSG00000276168.1 LBPSEMA4BLOOD756 19.32443
# 6  ENSG00000274012.1 LBPSEMA4BLOOD756 19.24718

df_long_lowest_first <- df_long %>%
  arrange(value)
head(df_long_lowest_first)
#                  Var1             Var2     value
# 1  ENSG00000001626.15 LBPSEMA4BLOOD646 -7.940354
# 2  ENSG00000002079.14 LBPSEMA4BLOOD646 -7.940354
# 3   ENSG00000003137.8 LBPSEMA4BLOOD646 -7.940354
# 4  ENSG00000005001.10 LBPSEMA4BLOOD646 -7.940354
# 5  ENSG00000005981.13 LBPSEMA4BLOOD646 -7.940354
# 6  ENSG00000006047.13 LBPSEMA4BLOOD646 -7.940354
# 7  ENSG00000006071.13 LBPSEMA4BLOOD646 -7.940354


#after melt it becomes
| Var1  | Var2    | value |
| ----- | ------- | ----- |
| GeneA | Sample1 | 2.1   |
| GeneB | Sample1 | 0.5   |
| GeneA | Sample2 | 3.5   |
| GeneB | Sample2 | 0.8   |
| GeneA | Sample3 | 1.8   |
| GeneB | Sample3 | 1.2   |

plot <- ggplot(df_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(x = "Blood gene expression", 
  y = "Density", 
  title = "Overall blood gene expression distribution across all samples") +
  theme_bw() +
  scale_x_continuous(breaks = seq(-8, 20, by = 2), limits = c(-8, 20))

ggsave("/hpc/users/hoangd02/www/plots/lbp/voom_blood_distribution.pdf", plot, width = 10, height = 6)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/voom_blood_distribution.pdf

##winsorized?

df_long <- reshape2::melt(v_brain$E)
summary(df_long$value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -7.952   0.998   3.163   2.998   4.900  18.566 

df_long_highest_first <- df_long %>%
  arrange(desc(value))
head(df_long_highest_first)
#                 Var1             Var2    value
# 1  ENSG00000274012.1 LBPSEMA4BRAIN524 18.56562
# 2  ENSG00000274012.1 LBPSEMA4BRAIN360 18.45295
# 3  ENSG00000274012.1 LBPSEMA4BRAIN302 18.24060
# 4  ENSG00000274012.1 LBPSEMA4BRAIN259 18.22504
# 5  ENSG00000274012.1 LBPSEMA4BRAIN618 18.18487

df_long_lowest_first <- df_long %>%
  arrange(value)
head(df_long_lowest_first)
#                 Var1             Var2     value
# 1  ENSG00000226278.1 LBPSEMA4BRAIN202 -7.951930
# 2  ENSG00000264943.1 LBPSEMA4BRAIN202 -7.951930
# 3  ENSG00000274655.1 LBPSEMA4BRAIN202 -7.951930
# 4  ENSG00000277577.1 LBPSEMA4BRAIN202 -7.951930
# 5  ENSG00000286179.1 LBPSEMA4BRAIN202 -7.951930
# 6  ENSG00000223779.6 LBPSEMA4BRAIN239 -7.550252
# 7  ENSG00000226278.1 LBPSEMA4BRAIN239 -7.550252
# 8  ENSG00000264943.1 LBPSEMA4BRAIN239 -7.550252

plot <- ggplot(df_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(
    x = "Brain gene expression",
    y = "Density",
    title = "Overall brain gene expression distribution across all samples"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-8, 20, by = 2), limits = c(-8, 20))

ggsave("/hpc/users/hoangd02/www/plots/lbp/voom_brain_distribution.pdf", plot, width = 10, height = 6)


## do this for full brain + blood form 3 + blood form 5

########### CREATING A TEMPORAL COMPONENT ########### ################
########### calculating time between surgeries ########### ########### 
##Timepoint values are just 0 and 1 (baseline vs follow-up) 

dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)
#   0   1 
# 152  73
# 79-73 is correct

##### creating time_between_surgeries column
library(lubridate)
dictionary <- dictionary %>%
  mutate(surgeryDate = as.Date(surgeryDate)) %>%
  group_by(IID_ISMMS) %>%
  mutate(
    time_between_surgeries = ifelse(
      n() > 1,
      as.numeric(max(surgeryDate, na.rm = TRUE) - min(surgeryDate, na.rm = TRUE)),
      0
    )
  ) %>%
  ungroup()
dictionary <- as.data.frame(dictionary)

table(dictionary$time_between_surgeries)
#  0 21 23 24 25 26 28 30 31 32 33 35 39 53 65 
# 79  2  2 10  2 30 56 10  4  8  6  6  6  2  2 

sub <- filter(dictionary, time_between_surgeries == 53); sub 
#PT-0082 is 65 days
##PT-0103 is 53 days

sub <- filter(dictionary, time_between_surgeries == 21); sub
#   IID_ISMMS SAMPLE_ISMMS_blood SAMPLE_ISMMS_brain number_of_pair_brain
# 1   PT-0041   LBPSEMA4BLOOD219   LBPSEMA4BRAIN121              2_pairs
# 2   PT-0041   LBPSEMA4BLOOD252   LBPSEMA4BRAIN266              2_pairs
#   mymet_tissue_brain mymet_tissue_blood surgeryDate timepoint
# 1            L_Brain            L_Blood  2014-10-29         0
# 2            R_Brain            R_Blood  2014-11-19         1
#   time_between_surgeries
# 1                     21
# 2                     21

##each of this can be treated as their own timepoint as well - way 2 

### way 1: compare 79 peole with baseline (0s for time_between_surgeries) and 73 people (non 0s)

### SAMPLE-LEVEL CONCORDANCE OVER TIME 
### GENExGENE CROSS-TISSUE CONCORDANCE + FUNCTIONAL ANNOTATION (GO TERMS) FOR HIGHEST COR PAIRS
### LATERALITY (which is sample-level)

############################################
############# SAMPLE-LEVEL ################# ORIGINAL 
############################################
## aka correlation of all 17k genes in blood samples to correlation of all 17k genes in brain samples
#brain full vs blood form5 no_residID
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
blood <- blood_form5_no_residID   # genes x blood_samples
brain <- brain_full_no_residID  # genes x brain_samples

# # brain full vs blood form5 with_residID
# brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
# row.names(brain_full) <- brain_full$V1
# brain_full$V1 <- NULL

# blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
# rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL
# blood <- blood_form5   # genes x blood_samples
# brain <- brain_full # genes x brain_samples

# blood <- v_blood$E   # genes x blood_samples
# brain <- v_brain$E   # genes x brain_samples

# Align genes (rows) between matrices
common_genes <- intersect(rownames(blood), rownames(brain))
length(common_genes) #17533 genes
blood <- blood[common_genes, , drop = FALSE]
brain <- brain[common_genes, , drop = FALSE]

# Sample–sample correlation matrix: rows = blood samples, cols = brain samples
#    (Spearman over genes)
cor_mat <- cor(blood, brain, method = "spearman")
dim(cor_mat) #225 225 ~~~ i think this is okay ~ and you'll subset after. i think it's better if you subset the df first though
#blood profile x = [2,5,1,3]
#brain profile y = [10,5,3,6]

# Absolute correlation
abs_cor_mat <- abs(cor_mat)

# Extract concordance for each row in dictionary by matching IDs
blood_ids <- dictionary$SAMPLE_ISMMS_blood
brain_ids <- dictionary$SAMPLE_ISMMS_brain

row_idx <- match(blood_ids, rownames(abs_cor_mat))
col_idx <- match(brain_ids, colnames(abs_cor_mat))

# Add concordance column to dictionary
dictionary$abs_concordance <- abs_cor_mat[cbind(row_idx, col_idx)]
dictionary$concordance <- cor_mat[cbind(row_idx, col_idx)]

## Creating a cross-time btw T1 and T2 for within person 
table(dictionary$number_of_pair_brain)
#  1_pair 2_pairs 
#      79     146 


##measuring concordance across timepoint would have to be within the same person
#sub <- filter(dictionary, time_between_surgeries == 21); sub

#sub <- filter(dictionary, time_between_surgeries == 65); sub
#   IID_ISMMS SAMPLE_ISMMS_blood SAMPLE_ISMMS_brain number_of_pair_brain
# 1   PT-0082   LBPSEMA4BLOOD357   LBPSEMA4BRAIN050              2_pairs
# 2   PT-0082   LBPSEMA4BLOOD475   LBPSEMA4BRAIN655              2_pairs
#   mymet_tissue_brain mymet_tissue_blood surgeryDate timepoint
# 1            R_Brain            R_Blood  2016-02-10         1
# 2            L_Brain            L_Blood  2015-12-07         0
#   time_between_surgeries concordance
# 1                     65   0.4956386
# 2                     65   0.5286571

# cor <- cor(t(v_blood$E), t(v_brain$E), method = "spearman")
# abs_cor <- abs(null_cor)
two_pairs <- filter(dictionary, number_of_pair_brain == "2_pairs") #146, correct

#Paired t-test: tests whether the mean difference (follow-up − baseline) is 0. Assumes the differences are ~normally distributed (not each group). Most powerful if that’s true; reasonably robust with n≈146.
#Wilcoxon signed-rank: tests whether the median difference is 0 (more generally, symmetry about 0). No normality assumption, uses ranks; more robust to outliers/heavy tails; slightly less power if differences are truly normal.

# Keep only IIDs that truly have both timepoints
paired_wide <- two_pairs %>%
  distinct(IID_ISMMS, timepoint, concordance) %>%   # guard against duplicates
  filter(timepoint %in% c(0, 1)) %>%
  group_by(IID_ISMMS) %>%
  filter(n_distinct(timepoint) == 2) %>%
  ungroup() %>%
  pivot_wider(names_from = timepoint, values_from = concordance, names_prefix = "tp_") %>%
  rename(baseline = tp_0, followup = tp_1)

head(paired_wide,10)
######################## blood_form5 vs full_brain no_residID 
#    IID_ISMMS followup baseline
#    <chr>        <dbl>    <dbl>
#  1 PT-0018   0.0453     0.0228
#  2 PT-0021   0.0352     0.0280
#  3 PT-0022   0.0523     0.0916
#  4 PT-0023   0.200      0.109 
#  5 PT-0026   0.000574   0.137 

######################## blood_form5 vs full_brain with_residID 
#    IID_ISMMS followup baseline
#    <chr>        <dbl>    <dbl>
#  1 PT-0018    0.0188    0.0330
#  2 PT-0021    0.121     0.0415
#  3 PT-0022    0.0445    0.0116
#  4 PT-0023    0.173     0.0366

# Paired t-test (mean difference)
ttest <- t.test(paired_wide$followup, paired_wide$baseline, paired = TRUE); ttest

######################## blood_form5 vs full_brain no_residID 
# data:  paired_wide$followup and paired_wide$baseline
# t = -2.5752, df = 72, p-value = 0.01207
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#  -0.065182352 -0.008299487
# sample estimates:
# mean difference 
#     -0.03674092 

# In a paired t-test, each subject (or sample pair) contributes its own before–after difference.
# If those differences are consistently in the same direction (most follow-up concordances < baseline concordances), 
#the standard error of the mean difference becomes small.
# Even if the absolute difference is tiny, consistency makes it statistically detectable.

# Mean difference: –0.0367 (≈ 3.7 percentage points on concordance).
# 95% CI: [–0.065, –0.008] → excludes 0, hence p ≈ 0.012.
# Interpretation: We’re 95% confident the true mean difference is a small but nonzero drop.
# Whether this is biologically meaningful is a separate (practical) question — statistics only tells you the difference is unlikely to be pure noise.

######################### blood_form5 vs full_brain with_residID 
# data:  paired_wide$followup and paired_wide$baseline
# t = -1.5949, df = 72, p-value = 0.1151
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#  -0.048813155  0.005421186
# sample estimates:
# mean difference 
#     -0.02169598 

#########
#########
#########
# Wilcoxon signed-rank (median difference), if you prefer a nonparametric test
wilcoxontest <- wilcox.test(paired_wide$followup, paired_wide$baseline, paired = TRUE, exact = FALSE); wilcoxontest
######################## blood_form5 vs full_brain no_residID 
# data:  paired_wide$followup and paired_wide$baseline
# V = 936, p-value = 0.02285
# alternative hypothesis: true location shift is not equal to 0

# p = 0.023 means the test rejects the null → there is evidence that the median concordance difference between follow-up and baseline is not 0.
# In plain words: concordance at follow-up tends to be lower than at baseline (since your paired t-test difference was negative).

# The Wilcoxon test doesn’t assume normality of differences.
# It’s more robust when differences are skewed or have outliers.
# Since both the paired t-test (p = 0.012) and Wilcoxon (p = 0.023) agree, the finding 
#(a small but consistent decrease in concordance at follow-up) is statistically reliable.

######################### blood_form5 vs full_brain with_residID 
# data:  paired_wide$followup and paired_wide$baseline
# V = 1064, p-value = 0.1159
# alternative hypothesis: true location shift is not equal to 0

library(ggsignif)

p_txt <- paste0("paired t: p = ", signif(ttest$p.value, 3),
                "  |  Wilcoxon: p = ", signif(wilcoxontest$p.value, 3))
p <- ggplot(two_pairs, aes(x = factor(timepoint), y = concordance, fill = factor(timepoint))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(
    x = "Timepoint: Baseline (0) and Follow-up (1)",
    y = "Concordance",
    title = "Sample-level blood brain GE concordance at baseline vs follow-up",
    subtitle = "Concordance = Spearman cor across 17,533 shared genes per blood brain pair"
  ) +
  #coord_cartesian(ylim = c(0, max(two_pairs$concordance, na.rm = TRUE) * 1.1)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ypos <- max(two_pairs$concordance, na.rm = TRUE) * 1.05
p_annot <- p + geom_signif(
      comparisons = list(c("0","1")),
      annotations = p_txt,
      y_position  = ypos,
      tip_length  = 0.4,
      vjust       = -1
    )

ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_baseline_follow_up_blood_form5_brain_full_with_residID_not_absolute_spearman.pdf", p_annot, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_baseline_follow_up_blood_form5_brain_full_no_residID_not_absolute_spearman.pdf", p_annot, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_baseline_follow_up_blood_form5_brain_full_with_residID.pdf", p_annot, width = 10, height = 8)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_baseline_follow_up_blood_form5_brain_full.pdf", p_annot, width = 10, height = 8)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_baseline_follow_up.pdf", p_annot, width = 10, height = 6)

table(two_pairs$timepoint)
#  0  1 
# 73 73 

######################## 
########################
##this is not using timepoint but time_between_surgeries
##for each person changing time_between_surgeries back to 0 for those with timepoint = 0 
two_pairs <- two_pairs %>%
  mutate(time_between_surgeries = if_else(timepoint == 0, 0, time_between_surgeries))

table(two_pairs$time_between_surgeries) ## this is only within two_pairs 
#  0 21 23 24 25 26 28 30 31 32 33 35 39 53 65 
# 73  1  1  5  1 15 28  5  2  4  3  3  3  1  1 

##with baseline
model <- lm(concordance ~ time_between_surgeries, data = two_pairs)
summary(model) 
######################## v_blood$E vs v_brain$E
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.091645 -0.038890 -0.004185  0.035027  0.149307 

# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.5261509  0.0059404  88.572   <2e-16 ***
# time_between_surgeries -0.0000156  0.0002802  -0.056    0.956 

######################## blood_form5 vs full_brain no_residID 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.12460 -0.06503 -0.02840  0.05153  0.31008 

# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.1247124  0.0103359  12.066  < 2e-16 ***
# time_between_surgeries -0.0012903  0.0004875  -2.647  0.00903 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.09016 on 144 degrees of freedom
# Multiple R-squared:  0.04639,	Adjusted R-squared:  0.03977 
# F-statistic: 7.005 on 1 and 144 DF,  p-value: 0.009034

######################## blood_form5 vs full_brain with_residID 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.11113 -0.06765 -0.02380  0.05183  0.28950 

# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.1118443  0.0094789  11.799   <2e-16 ***
# time_between_surgeries -0.0008301  0.0004471  -1.857   0.0654 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.08268 on 144 degrees of freedom
# Multiple R-squared:  0.02338,	Adjusted R-squared:  0.0166 
# F-statistic: 3.447 on 1 and 144 DF,  p-value: 0.0654

coefs <- coef(model)
a <- round(coefs[1], 3)
b <- round(coefs[2], 3)

r2 <- round(summary(model)$r.squared, 3)
pval <- signif(summary(model)$coefficients[2,4], 3)  # slope p-value

plot <- ggplot(two_pairs, aes(x = time_between_surgeries, y = concordance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate(
    "text",
    x = min(two_pairs$time_between_surgeries),
    y = max(two_pairs$concordance),
    hjust = 0,
    label = paste0("y = ", a, " + ", b, "x\nR² = ", r2, ", p = ", pval)
  ) +
  labs(
    x = "Days Between Surgeries",
    y = "Concordance",
    title = "Effect of Time Between Surgeries on Concordance"
  ) +
  theme_minimal(base_size = 14)

ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_blood_form5_brain_full_with_residID.pdf", plot, width = 10, height = 6)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_blood_form5_brain_full.pdf", plot, width = 10, height = 6)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up.pdf", plot, width = 10, height = 6)

##without baselines 
two_pairs_not_baseline <- two_pairs[two_pairs$time_between_surgeries != 0, ]
table(two_pairs_not_baseline$time_between_surgeries)
# 21 23 24 25 26 28 30 31 32 33 35 39 53 65 
#  1  1  5  1 15 28  5  2  4  3  3  3  1  1 

model <- lm(concordance ~ time_between_surgeries, data = two_pairs_not_baseline)
summary(model)
######################## v_blood$E vs v_brain$E
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.092837 -0.035617 -0.001371  0.027124  0.153361 

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.4993667  0.0292574   17.07   <2e-16 ***
# time_between_surgeries 0.0008586  0.0009758    0.88    0.382  

######################## blood_form5 vs full_brain no_residID 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.09260 -0.06303 -0.02803  0.03936  0.30894 

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)   
# (Intercept)             0.137962   0.048998   2.816   0.0063 **
# time_between_surgeries -0.001723   0.001634  -1.054   0.2954

######################## blood_form5 vs full_brain with_residID 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.08581 -0.05203 -0.01637  0.03752  0.20807 

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.143349   0.040274   3.559 0.000668 ***
# time_between_surgeries -0.001858   0.001343  -1.384 0.170833    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.07073 on 71 degrees of freedom
# Multiple R-squared:  0.02625,	Adjusted R-squared:  0.01254 
# F-statistic: 1.914 on 1 and 71 DF,  p-value: 0.1708

###PLOTTING
coefs <- coef(model)
a <- round(coefs[1], 3)
b <- round(coefs[2], 3)

r2 <- round(summary(model)$r.squared, 3)
pval <- signif(summary(model)$coefficients[2,4], 3)  # slope p-value

plot <- ggplot(two_pairs_not_baseline, aes(x = time_between_surgeries, y = concordance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate(
    "text",
    x = min(two_pairs_not_baseline$time_between_surgeries),
    y = max(two_pairs_not_baseline$concordance),
    hjust = 0,
    label = paste0("y = ", a, " + ", b, "x\nR² = ", r2, ", p = ", pval)
  ) +
  labs(
    x = "Days Between Surgeries",
    y = "Concordance",
    title = "Effect of Time Between Surgeries on Concordance"
  ) +
  theme_minimal(base_size = 14)

ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_no_baseline_blood_form5_brain_full_with_residID.pdf", plot, width = 10, height = 6)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_no_baseline_blood_form5_brain_full.pdf", plot, width = 10, height = 6)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_no_baseline.pdf", plot, width = 10, height = 6)


#two_pairs is the new dictionary 
## dropping the last two time 53 & 65 WITH BASELINE
two_pairs_1_month <- two_pairs[!two_pairs$time_between_surgeries %in% c(53, 65), ]
model <- lm(concordance ~ time_between_surgeries, data = two_pairs_1_month)
summary(model)
############################# brain full vs blood form5 with_residID removing IDs > 39 days
# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.1118359  0.0096632  11.573   <2e-16 ***
# time_between_surgeries -0.0008355  0.0004791  -1.744   0.0833 . 

############################# brain full vs blood form5 no_residID removing IDs > 39 days
# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.1243110  0.0105380    11.8   <2e-16 ***
# time_between_surgeries -0.0012540  0.0005224    -2.4   0.0177 * 

## dropping the last two time 53 & 65 WITH NO BASELINE removing IDs > 39 days
############################# brain full vs blood form5 no_residID
two_pairs_not_baseline <- two_pairs[two_pairs$time_between_surgeries != 0, ]
two_pairs_1_month <- two_pairs_not_baseline[!two_pairs_not_baseline$time_between_surgeries %in% c(53, 65), ]
model <- lm(concordance ~ time_between_surgeries, data = two_pairs_1_month)
summary(model)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
# (Intercept)             0.135011   0.083819   1.611    0.112
# time_between_surgeries -0.001624   0.002918  -0.556    0.580

############################# brain full vs blood form5 with_residID key:askg
two_pairs_not_baseline <- two_pairs[two_pairs$time_between_surgeries != 0, ]
two_pairs_1_month <- two_pairs_not_baseline[!two_pairs_not_baseline$time_between_surgeries %in% c(53, 65), ]
model <- lm(concordance ~ time_between_surgeries, data = two_pairs_1_month)
summary(model)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)   
# (Intercept)             0.201177   0.068266   2.947  0.00437 **
# time_between_surgeries -0.003922   0.002376  -1.650  0.10342 

coefs <- coef(model)
a <- round(coefs[1], 3)
b <- round(coefs[2], 3)

r2 <- round(summary(model)$r.squared, 3)
pval <- signif(summary(model)$coefficients[2,4], 3)  # slope p-value

plot <- ggplot(two_pairs_1_month, aes(x = time_between_surgeries, y = concordance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate(
    "text",
    x = min(two_pairs_1_month$time_between_surgeries),
    y = max(two_pairs_1_month$concordance),
    hjust = 0,
    label = paste0("y = ", a, " + ", b, "x\nR² = ", r2, ", p = ", pval)
  ) +
  labs(
    x = "Days Between Surgeries",
    y = "Concordance",
    title = "Effect of Time Between Surgeries on Concordance"
  ) +
  theme_minimal(base_size = 14)

ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_no_baseline_blood_form5_brain_full_no_residID_1_month.pdf", plot, width = 10, height = 6)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_with_baseline_blood_form5_brain_full_no_residID_1_month.pdf", plot, width = 10, height = 6)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/concordance_at_follow_up_with_baseline_blood_form5_brain_full_with_residID_1_month.pdf", plot, width = 10, height = 6)

########### COMPARING CONCORDANCE IN HEMISPHERES ########### ################
head(dictionary)

# # brain full vs blood form5 no_residID
# brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
# row.names(brain_full_no_residID)<- brain_full_no_residID$V1; brain_full_no_residID$V1 <- NULL
# dim(brain_full_no_residID) #21356   225

# blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
# rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
# blood <- blood_form5_no_residID # genes x blood_samples
# brain <- brain_full_no_residID # genes x brain_samples

# brain full vs blood form5 with_residID
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full) <- brain_full$V1; brain_full$V1 <- NULL

blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL
blood <- blood_form5   # genes x blood_samples
brain <- brain_full # genes x brain_samples

# Align genes (rows) between matrices
common_genes <- intersect(rownames(blood), rownames(brain))
length(common_genes) #17533 genes
blood <- blood[common_genes, , drop = FALSE]
brain <- brain[common_genes, , drop = FALSE]

# Sample–sample correlation matrix: rows = blood samples, cols = brain samples
#    (Spearman over genes)
cor_mat <- cor(blood, brain, method = "spearman")
#blood profile x = [2,5,1,3]
#brain profile y = [10,5,3,6]

# Absolute correlation
abs_cor_mat <- abs(cor_mat)

# Extract concordance for each row in dictionary by matching IDs
blood_ids <- dictionary$SAMPLE_ISMMS_blood
brain_ids <- dictionary$SAMPLE_ISMMS_brain

row_idx <- match(blood_ids, rownames(abs_cor_mat))
col_idx <- match(brain_ids, colnames(abs_cor_mat))

# Add concordance column to dictionary
dictionary$abs_concordance <- abs_cor_mat[cbind(row_idx, col_idx)]
dictionary$concordance <- cor_mat[cbind(row_idx, col_idx)]

##measuring concordance across timepoint would have to be within the same person
#sub <- filter(dictionary, time_between_surgeries == 21); sub
head(dictionary)

table(dictionary$mymet_tissue_brain)
# L_Blood L_Brain R_Blood R_Brain 
#       0     119       0     106 

########### NO RESID_ID
wilcoxontest <- wilcox.test(concordance ~ mymet_tissue_brain, data = dictionary); wilcoxontest
#W = 6909, p-value = 0.2172
ttest <-t.test(concordance ~ mymet_tissue_brain, data = dictionary); ttest
#t = 0.95071, df = 223, p-value = 0.3428

########### WITH RESID_ID
# Wilcoxon test (non-parametric, safer for skewed data)
wilcoxontest <- wilcox.test(concordance ~ mymet_tissue_brain, data = dictionary); wilcoxontest
#W = 6797, p-value = 0.3152
ttest <- t.test(concordance ~ mymet_tissue_brain, data = dictionary); ttest
# t = 0.89491, df = 223, p-value = 0.3718

library(ggsignif)

p_txt <- paste0("paired t: p = ", signif(ttest$p.value, 3),
                "  |  Wilcoxon: p = ", signif(wilcoxontest$p.value, 3))
p <- ggplot(dictionary, aes(x = mymet_tissue_brain, y = concordance, fill = mymet_tissue_brain)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(
    x = "Hemisphere",
    y = "Concordance",
    title = "Concordance by Hemisphere",
    subtitle = "Concordance = Spearman cor across 17,533 shared genes per blood brain pair"
  ) +
  #coord_cartesian(ylim = c(0, max(dictionary$concordance, na.rm = TRUE) * 1.1)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ypos <- max(dictionary$concordance, na.rm = TRUE) * 1.05
p_annot <- p + geom_signif(
      comparisons = list(c("L_Brain", "R_Brain")),
      annotations = p_txt,
      y_position  = ypos,
      tip_length  = 0.4,
      vjust       = -1
    )
ggsave("/hpc/users/hoangd02/www/plots/lbp/laterality_blood_form5_full_brain_with_residID_no_absolute.pdf", p_annot, width = 10, height = 8)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/laterality_blood_form5_full_brain_no_residID_no_absolute.pdf", p_annot, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/laterality_blood_form5_full_brain_with_residID.pdf",  p_annot, width = 10, height = 8)
#ggsave("/hpc/users/hoangd02/www/plots/lbp/laterality_blood_form5_full_brain_no_residID.pdf", p_annot, width = 10, height = 8)

############################################
############# GENE-LEVEL ###################
############################################
############################################
############# SAMPLE-LEVEL #################
############################################
# brain full vs blood form5 no_residID
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
blood <- blood_form5_no_residID   # genes x blood_samples
brain <- brain_full_no_residID  # genes x brain_samples

# brain full vs blood form5 with_residID
brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full) <- brain_full$V1
brain_full$V1 <- NULL

blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL
blood <- blood_form5   # genes x blood_samples
brain <- brain_full # genes x brain_samples

### ALL GENES ~~~
cor <- cor(t(blood), t(brain), method = "spearman")
dim(cor) #21046 21356
cor_mean <- mean(as.matrix(cor)); cor_mean #-0.001142294
cor_median <- median(as.matrix(cor)); cor_median #-0.001142294
cor_99.95 <- quantile(as.matrix(cor), probs = 0.9995); cor_99.95 #0.2184869
cor_99.995 <- quantile(as.matrix(cor), probs = 0.99995); cor_99.995 #0.2583218

cor_99.9995 <- quantile(as.matrix(cor), probs = 0.999995); cor_99.9995 #0.3069401

min(as.matrix(cor)) # -0.4045438
max(as.matrix(cor)) # 0.4450948

h <- hist(cor, breaks = 100, plot = FALSE)  # compute histogram counts only
df_hist <- data.frame(
  mids = h$mids,
  counts = h$counts
)

p<- ggplot(df_hist, aes(x = mids, y = counts)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(
    x = "Spearman correlation",
    y = "Frequency",
    #title = "Distribution of blood brain gene to gene correlations"
  ) +
  theme_bw()

ggsave("/hpc/users/hoangd02/www/plots/lbp/all_gene_pair_distribution_blood_form5_full_brain_no_residID.pdf",  p, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/all_gene_pair_distribution_blood_form5_full_brain_with_residID.pdf",  p, width = 10, height = 8)

### COMMON GENES ~~~
common_genes <- intersect(rownames(blood), rownames(brain))
length(common_genes) #17533 genes
blood <- blood[common_genes, , drop = FALSE]
brain <- brain[common_genes, , drop = FALSE]

cor <- cor(t(blood), t(brain), method = "spearman")
dim(cor) #17533 x 17533

#key: asbh
cor_mean <- mean(as.matrix(cor)); cor_mean #-0.001671532
cor_99.95 <- quantile(as.matrix(cor), probs = 0.9995); cor_99.95 #0.2173272
cor_99.995 <- quantile(as.matrix(cor), probs = 0.99995); cor_99.995 #0.2567991
cor_99.9995 <- quantile(as.matrix(cor), probs = 0.999995); cor_99.9995 #0.3100742

min(as.matrix(cor)) # -0.4045438
max(as.matrix(cor)) # 0.4450948

################## JUST THE DIAGONAL, A GENE TO ITSELF
gene_wise_cor <- diag(cor)
length(gene_wise_cor)  # should be 17533

df_gene_cor <- data.frame(
  gene = common_genes,
  correlation = gene_wise_cor
)

summary(df_gene_cor$correlation)
#with_residID
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.26374 -0.02554  0.02102  0.02501  0.06977  0.43333 

#no_residID
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.48069 -0.01421  0.03584  0.04831  0.08977  0.87805

gene_wise_cor_99.95 <- quantile(as.matrix(gene_wise_cor), probs = 0.9995); gene_wise_cor_99.95 #0.3689219 
gene_wise_cor_99.995 <- quantile(as.matrix(gene_wise_cor), probs = 0.99995); gene_wise_cor_99.995 #0.4174699

p <- ggplot(df_gene_cor, aes(x = correlation)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(
    x = "Spearman correlation (Blood vs. Brain, Same Gene)",
    y = "Count",
    #title = "Distribution of gene_wise blood brain correlations"
  ) +
  theme_bw()

ggsave("/hpc/users/hoangd02/www/plots/lbp/gene_to_itself_distribution_blood_form5_full_brain_no_residID.pdf",  p, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/gene_to_itself_distribution_blood_form5_full_brain_with_residID.pdf",  p, width = 10, height = 8)

##with red dotted line = mean 
mean_cor <- mean(df_gene_cor$correlation, na.rm = TRUE)
p <- ggplot(df_gene_cor, aes(x = correlation)) +
  geom_histogram(bins = 50, fill = "cornflowerblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean_cor, 
             color = "red", 
             linetype = "dotted", 
             linewidth = 1.5) +    
  annotate("text",
           x = mean_cor,
           y = Inf, vjust = -1,
           label = paste0("Mean = ", round(mean_cor, 3)),
           color = "red", size = 8, fontface = "bold") +  # increase annotation text size
  labs(
    x = "Spearman correlation (Blood vs. Brain, Same Gene)",
    y = "Count"
  ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  theme_bw(base_size = 12) +  # <-- overall base font size
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 16),  # axis labels
    axis.text  = element_text(size = 14)                 # tick labels
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/gene_to_itself_distribution_blood_form5_full_brain_no_residID_with_mean.pdf",  p, width = 6, height = 5)

################ gene-level (all genes and same-genes) and sample-level cross-time T0 and T1
################
# brain full vs blood form5 no_residID
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
blood <- blood_form5_no_residID   # genes x blood_samples
brain <- brain_full_no_residID  # genes x brain_samples

##run the temporal component of dictionary
dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)

two_pairs <- filter(dictionary, number_of_pair_brain == "2_pairs") #146, correct
dim(two_pairs) #146   8

two_pairs_timepoint0 <- filter(two_pairs, timepoint == "0")
two_pairs_timepoint1 <- filter(two_pairs, timepoint == "1")

dim(two_pairs_timepoint0) #73  8
dim(two_pairs_timepoint1) #73  8

########## TIMEPOINT 0 ##########
selected_samples_blood <- intersect(colnames(blood), two_pairs_timepoint0$SAMPLE_ISMMS_blood)
blood_subset_t0 <- blood[, selected_samples_blood, drop = FALSE]
dim(blood_subset_t0) #21046    73

selected_samples_brain <- intersect(colnames(brain), two_pairs_timepoint0$SAMPLE_ISMMS_brain)
brain_subset_t0 <- brain[, selected_samples_brain, drop = FALSE]
dim(brain_subset_t0) #21356    73

########## TIMEPOINT 1 ##########
selected_samples_blood <- intersect(colnames(blood), two_pairs_timepoint1$SAMPLE_ISMMS_blood)
blood_subset_t1 <- blood[, selected_samples_blood, drop = FALSE]

selected_samples_brain <- intersect(colnames(brain), two_pairs_timepoint1$SAMPLE_ISMMS_brain)
brain_subset_t1 <- brain[, selected_samples_brain, drop = FALSE]

########## CORRELATIONS ########## GENE-LEVEL SAME GENES
## ## ## ## ## ## ## ## T0 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
common_genes_t0 <- intersect(rownames(blood_subset_t0), rownames(brain_subset_t0))
length(common_genes_t0) #17533 genes

blood_subset_t0_common_genes <- blood_subset_t0[common_genes_t0, , drop = FALSE]
dim(blood_subset_t0_common_genes) #17533    73
brain_subset_t0_common_genes <- brain_subset_t0[common_genes_t0, , drop = FALSE]
dim(brain_subset_t0_common_genes) #17533    73

cor_common_genes_t0 <- cor(t(blood_subset_t0_common_genes), t(brain_subset_t0_common_genes), method = "spearman")
dim(cor_common_genes_t0) #17533 x 17533

## ## ## ## ## ## ## ## T1 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
common_genes_t1 <- intersect(rownames(blood_subset_t1), rownames(brain_subset_t1))
length(common_genes_t1) #17533 genes

blood_subset_t1_common_genes <- blood_subset_t1[common_genes_t1, , drop = FALSE]
dim(blood_subset_t1_common_genes) # 17533    73 perfect
brain_subset_t1_common_genes <- brain_subset_t1[common_genes_t1, , drop = FALSE]
dim(brain_subset_t1_common_genes) # 17533    73 perfect

cor_common_genes_t1 <- cor(t(blood_subset_t1_common_genes), t(brain_subset_t0_common_genes), method = "spearman")
dim(cor_common_genes_t1) #17533 x 17533

########## CORRELATIONS ########## GENE-LEVEL ALL GENES
cor_t0 <- cor(t(blood_subset_t0), t(brain_subset_t0), method = "spearman")
dim(cor_t0) #21046 x 21356
cor_t0[1:5,1:5]

cor_t1 <- cor(t(blood_subset_t1), t(brain_subset_t1), method = "spearman")
dim(cor_t1) #21046 x 21356
cor_t1[1:5,1:5]

################################### GENE-LEVEL 
###################################
# Create timepoint-aware swap function
permute_timepoints <- function(blood_t0, brain_t0, blood_t1, brain_t1, swap_frac) {
  stopifnot(swap_frac %in% c(0, 1))
  if (swap_frac == 0) {
    # Same timepoint pairing
    blood_combined <- cbind(blood_t0, blood_t1)
    brain_combined <- cbind(brain_t0, brain_t1)
  } else {
    # Cross timepoint pairing  
    blood_combined <- cbind(blood_t0, blood_t1)
    brain_combined <- cbind(brain_t1, brain_t0)  # Swapped!
  }
  
  return(list(blood = blood_combined, brain = brain_combined))
}

# Enhanced function to handle combined + separate timepoint analyses
save_timepoint_analyses <- function(blood_t0, brain_t0, blood_t1, brain_t1, 
                                   swap_frac, out_dir) {
  cat("\n=== Computing correlations for swap", swap_frac*100, "% ===\n")
  
  if (swap_frac == 0) {
    # SAME TIMEPOINT analyses
    
    # 1. Combined same timepoints
    blood_combined <- cbind(blood_t0, blood_t1)
    brain_combined <- cbind(brain_t0, brain_t1)
    cor_combined <- cor(t(blood_combined), t(brain_combined), method = "spearman")
    
    # 2. Timepoint 0 only
    cor_t0_only <- cor(t(blood_t0), t(brain_t0), method = "spearman")
    
    # 3. Timepoint 1 only  
    cor_t1_only <- cor(t(blood_t1), t(brain_t1), method = "spearman")
    
    # Save all three
    ## ALL-GENES
    # saveRDS(cor_combined, file.path(out_dir, "spearman_same_timepoint_combined_20250923.rds"))
    # saveRDS(cor_t0_only, file.path(out_dir, "spearman_same_timepoint_t0_only_20250923.rds"))
    # saveRDS(cor_t1_only, file.path(out_dir, "spearman_same_timepoint_t1_only_20250923.rds"))
    ## SAME-GENES
    saveRDS(cor_combined, file.path(out_dir, "spearman_same_timepoint_same_genes_combined_20250923.rds"))
    saveRDS(cor_t0_only, file.path(out_dir, "spearman_same_timepoint_same_genes_t0_only_20250923.rds"))
    saveRDS(cor_t1_only, file.path(out_dir, "spearman_same_timepoint_same_genes_t1_only_20250923.rds"))
    
    # Also save absolute values
    ## ALL-GENES
    # saveRDS(abs(cor_combined), file.path(out_dir, "spearman_same_timepoint_combined_abs_20250923.rds"))
    # saveRDS(abs(cor_t0_only), file.path(out_dir, "spearman_same_timepoint_t0_only_ab_20250923.rds"))
    # saveRDS(abs(cor_t1_only), file.path(out_dir, "spearman_same_timepoint_t1_only_abs_20250923.rds"))
    ## SAME-GENES
    saveRDS(abs(cor_combined), file.path(out_dir, "spearman_same_timepoint_same_genes_combined_abs_20250923.rds"))
    saveRDS(abs(cor_t0_only), file.path(out_dir, "spearman_same_timepoint_same_genes_t0_only_ab_20250923.rds"))
    saveRDS(abs(cor_t1_only), file.path(out_dir, "spearman_same_timepoint_same_genes_t1_only_abs_20250923.rds"))
    
    return(list(combined = cor_combined, t0_only = cor_t0_only, t1_only = cor_t1_only))
    
  } else {
    # CROSS TIMEPOINT analysis
    blood_combined <- cbind(blood_t0, blood_t1)
    brain_combined <- cbind(brain_t1, brain_t0)  # Swapped!
    cor_cross <- cor(t(blood_combined), t(brain_combined), method = "spearman")
    ## ALL-GENES
    # saveRDS(cor_cross, file.path(out_dir, "spearman_cross_timepoint_20250923.rds"))
    # saveRDS(abs(cor_cross), file.path(out_dir, "spearman_cross_timepoint_abs_20250923.rds"))
    ## SAME-GENES
    saveRDS(cor_cross, file.path(out_dir, "spearman_cross_timepoint_same_genes_20250923.rds"))
    saveRDS(abs(cor_cross), file.path(out_dir, "spearman_cross_timepoint_same_genes_abs_20250923.rds"))
    
    return(list(cross = cor_cross))
  }
}

# Set output directory
out_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap"

# Run analyses - GENE-WISE ALL GENES 
same_tp_results <- save_timepoint_analyses(blood_subset_t0, brain_subset_t0, 
                                          blood_subset_t1, brain_subset_t1, 
                                          swap_frac = 0, out_dir)

cross_tp_results <- save_timepoint_analyses(blood_subset_t0, brain_subset_t0, 
                                           blood_subset_t1, brain_subset_t1, 
                                           swap_frac = 1, out_dir)

##GENE-WISE SAME GENES 
same_tp_results <- save_timepoint_analyses(blood_subset_t0_common_genes, brain_subset_t0_common_genes,  
                                          blood_subset_t1_common_genes, brain_subset_t1_common_genes, 
                                          swap_frac = 0, out_dir)

cross_tp_results <- save_timepoint_analyses(blood_subset_t0_common_genes, brain_subset_t0_common_genes,  
                                          blood_subset_t1_common_genes, brain_subset_t1_common_genes, 
                                           swap_frac = 1, out_dir)
# spearman_cross_timepoint_20250923.rds
# spearman_cross_timepoint_abs_20250923.rds
# spearman_same_timepoint_combined_20250923.rds
# spearman_same_timepoint_combined_abs_20250923.rds
# spearman_same_timepoint_t0_only_20250923.rds
# spearman_same_timepoint_t0_only_ab_20250923.rds
# spearman_same_timepoint_t1_only_20250923.rds
# spearman_same_timepoint_t1_only_abs_20250923.rds

# Load all results
library(effsize)

out_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap"

#cor_combined <- readRDS(file.path(out_dir, "spearman_same_timepoint_combined_20250923.rds"))
cor_t0_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t0_only_20250923.rds"))
cor_t1_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t1_only_20250923.rds"))
cor_cross <- readRDS(file.path(out_dir, "spearman_cross_timepoint_20250923.rds"))

#dim for all is 21046 21356 ~~ it's just across different numbers of samples

# Basic dimensions and summaries
cat("Dimensions:\n")
cat("Combined:", dim(cor_combined), "\n")
cat("T0 only:", dim(cor_t0_only), "\n") 
cat("T1 only:", dim(cor_t1_only), "\n")
cat("Cross timepoint:", dim(cor_cross), "\n")

######################################################################
################################### OCTOBER 10 
# Convert matrices to vectors (this preserves the pairing)
vec_t0 <- as.vector(cor_t0_only)
vec_t1 <- as.vector(cor_t1_only)
vec_cross <- as.vector(cor_cross)

# Perform paired Wilcoxon tests
wilcox_t0_vs_cross <- wilcox.test(vec_t0, vec_cross, paired = TRUE) #479037.work_work_work
wilcox_t1_vs_cross <- wilcox.test(vec_t1, vec_cross, paired = TRUE) #1739991.wcpg_plotting
wilcox_t0_vs_t1 <- wilcox.test(vec_t0, vec_t1, paired = TRUE) #cbipm01-2

# Cohen's d for paired samples
d_t0_t1 <- cohen.d(vec_t0, vec_t1, paired = TRUE)
d_combined_t0 <- cohen.d(vec_t0, vec_cross, paired = TRUE)
d_combined_t1 <- cohen.d(vec_t1, vec_cross, paired = TRUE)

### RESULTS
wilcox_t0_vs_t1
#         Wilcoxon signed rank test with continuity correction
# data:  vec_t0 and vec_t1
# V = 4.9995e+16, p-value < 2.2e-16

d_t0_t1
# d estimate: -0.01363887 (negligible)
# 95 percent confidence interval:
#       lower       upper 
# -0.01376914 -0.01350861

#####
wilcox_t0_vs_cross

d_combined_t0
# d estimate: -0.03509636 (negligible)
# 95 percent confidence interval:
#       lower       upper 
# -0.03521082 -0.03498189 

#####
wilcox_t1_vs_cross

d_combined_t1
# d estimate: -0.01818533 (negligible)
# 95 percent confidence interval:
#       lower       upper 
# -0.01829246 -0.01807821

# First, identify which genes are in BOTH blood (rows) and brain (columns)
common_genes <- intersect(rownames(cor_t0_only), colnames(cor_t0_only))

# Check how many genes overlap
length(common_genes)

# Extract diagonal values for common genes only
diag_t0 <- numeric(length(common_genes))
diag_t1 <- numeric(length(common_genes))
diag_cross <- numeric(length(common_genes))

for (i in seq_along(common_genes)) {
  gene <- common_genes[i]
  diag_t0[i] <- cor_t0_only[gene, gene]
  diag_t1[i] <- cor_t1_only[gene, gene]
  diag_cross[i] <- cor_cross[gene, gene]
}

# Now perform paired Wilcoxon tests on diagonal elements only
wilcox_t0_vs_t1 <- wilcox.test(diag_t0, diag_t1, paired = TRUE)
wilcox_t0_vs_cross <- wilcox.test(diag_t0, diag_cross, paired = TRUE)
wilcox_t1_vs_cross <- wilcox.test(diag_t1, diag_cross, paired = TRUE)

wilcox_t0_vs_t1
#         Wilcoxon signed rank test with continuity correction
# data:  diag_t0 and diag_t1
# V = 91994222, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox_t0_vs_cross
#         Wilcoxon signed rank test with continuity correction
# data:  diag_t0 and diag_cross
# V = 82447292, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox_t1_vs_cross
#         Wilcoxon signed rank test with continuity correction
# data:  diag_t1 and diag_cross
# V = 62034029, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

cohen_d_paired <- function(x1, x2) {
  diff <- x1 - x2
  d <- mean(diff) / sd(diff)
  return(d)
}

# Calculate Cohen's d for each comparison
d_t0_vs_t1 <- cohen_d_paired(diag_t0, diag_t1);d_t0_vs_t1 #0.1752861
d_t0_vs_cross <- cohen_d_paired(diag_t0, diag_cross); d_t0_vs_cross #0.06412488
d_t1_vs_cross <- cohen_d_paired(diag_t1, diag_cross); d_t1_vs_cross #-0.1691533

######################################################################
########### END OCT 10

# Distribution summaries
#summary(as.vector(cor_combined))
summary(as.vector(cor_t0_only))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.798686 -0.079816 -0.001018 -0.001010  0.077811  0.853449 
summary(as.vector(cor_t1_only))
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.8275022 -0.0833333  0.0008639  0.0006206  0.0847526  0.8946686 
summary(as.vector(cor_cross))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.778516 -0.054106  0.002416  0.002566  0.059127  0.852934
# Create comparison dataframe
cor_comparison <- data.frame(
  #Combined = as.vector(abs(cor_combined)),
  T0_only = as.vector(abs(cor_t0_only)),
  T1_only = as.vector(abs(cor_t1_only)),
  Cross_timepoint = as.vector(abs(cor_cross))
)

# Boxplot comparison
cor_long <- melt(cor_comparison)
plot <- ggplot(cor_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Absolute Correlation Distributions", 
       x = "Analysis Type", y = "Absolute Correlation") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/AbsoluteCorrelationDistributions.pdf", plot, width = 10, height = 6)

library(ggrastr)
# Rasterize the entire plot
cor_long <- melt(cor_comparison)
plot <- ggplot(cor_long, aes(x = variable, y = value)) +
  rasterise(geom_boxplot(), dpi = 300) +  # Rasterize boxplot at 300 DPI
  labs(title = "Absolute Correlation Distributions", 
       x = "Analysis Type", y = "Absolute Correlation") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/AbsoluteCorrelationDistributions_rasterized.pdf", plot, width = 10, height = 8)
https://hoangd02.u.hpc.mssm.edu/plots/lbp/AbsoluteCorrelationDistributions_rasterized.pdf

############ NON ABSOLUTE ############
######################################
cor_comparison <- data.frame(
  #Combined = as.vector(cor_combined),
  T0_only = as.vector(cor_t0_only),
  T1_only = as.vector(cor_t1_only),
  Cross_timepoint = as.vector(cor_cross)
)

# cor_long <- melt(cor_comparison)
# plot <- ggplot(cor_long, aes(x = variable, y = value)) +
#   rasterise(geom_boxplot(), dpi = 300) +  # Rasterize boxplot at 300 DPI
#   labs(title = "Correlation Distributions", 
#        x = "Analysis Type", y = "Correlation") +
#   theme_minimal()

# ggsave("/hpc/users/hoangd02/www/plots/lbp/CorrelationDistributions_rasterized_without_combined.pdf", plot, width = 10, height = 8)

# #ggsave("/hpc/users/hoangd02/www/plots/lbp/CorrelationDistributions_rasterized.pdf", plot, width = 10, height = 8)
# https://hoangd02.u.hpc.mssm.edu/plots/lbp/CorrelationDistributions_rasterized.pdf
# ### i like this !! 

######## computing some stats 
library(ggpubr)
library(ggrastr)
library(effsize)

# Sample data for faster computation
set.seed(123)
#cor_sample <- cor_comparison[sample(nrow(cor_comparison), 1000000), ] #1 million rows for all 
#cor_long_sample <- melt(cor_sample)

#1 million rows per group 
cor_sample <- melt(cor_comparison)
cor_sample <- cor_sample %>%
  group_by(variable) %>%
  slice_sample(n = 1e7)
dim(cor_sample) #3000000       2

# Calculate Cohen's d for all comparisons
# d_combined_cross <- cohen.d(cor_sample$Combined, cor_sample$Cross_timepoint)
##for stats purposes, use the WHOLE MATRIX ~ 449M
cor_sample <- melt(cor_comparison)

library(dplyr)
library(effsize)
library(purrr)
library(tibble)

# ---- helpers ----
get_group <- function(df, g) df$value[df$variable == g]

cohen_d_pair <- function(df, g1, g2,
                         hedges.correction = FALSE, pooled = FALSE) {
  x <- get_group(df, g1); y <- get_group(df, g2)
  effsize::cohen.d(x, y, hedges.correction = hedges.correction, pooled = pooled)
}

wilcox_pair <- function(df, g1, g2, paired = FALSE) {
  x <- get_group(df, g1); y <- get_group(df, g2)
  wilcox.test(x, y, paired = paired, exact = FALSE)
}

# ---- your three comparisons (unpaired) ----
d_t0_t1        <- cohen_d_pair(cor_sample, "T0_only", "T1_only")
d_combined_t0  <- cohen_d_pair(cor_sample, "T0_only", "Cross_timepoint")
d_combined_t1  <- cohen_d_pair(cor_sample, "T1_only", "Cross_timepoint")

w_t0_vs_x  <- wilcox_pair(cor_sample, "T0_only", "Cross_timepoint", paired = FALSE)
w_t1_vs_x  <- wilcox_pair(cor_sample, "T1_only", "Cross_timepoint", paired = FALSE)
w_t0_vs_t1 <- wilcox_pair(cor_sample, "T0_only", "T1_only",        paired = FALSE)

# If you want just numeric d values:
as.numeric(d_t0_t1$estimate)
as.numeric(d_combined_t0$estimate)
as.numeric(d_combined_t1$estimate)

# ---- Optional: run many pairs + return a tidy table ----
pairs <- tribble(
  ~group1,          ~group2,
  "T0_only",        "T1_only",
  "T0_only",        "Cross_timepoint",
  "T1_only",        "Cross_timepoint"
)

results <- pairs %>%
  mutate(
    d_obj = map2(group1, group2, ~ cohen_d_pair(cor_sample, .x, .y)),
    d     = map_dbl(d_obj, ~ as.numeric(.x$estimate)),
    p     = map2_dbl(group1, group2, ~ wilcox_pair(cor_sample, .x, .y, paired = FALSE)$p.value)
  ) %>%
  select(group1, group2, d, p)

results 
## THIS IS FOR ALL-GENE PAIRS!!! not just 1 million 
#   group1  group2                d     p
#   <chr>   <chr>             <dbl> <dbl>
# 1 T0_only T1_only         -0.0133     0
# 2 T0_only Cross_timepoint -0.0427     0
# 3 T1_only Cross_timepoint -0.0233     0


################# OLD

d_t0_t1 <- cohen.d(cor_sample$T0_only, cor_sample$T1_only)
d_combined_t0 <- cohen.d(cor_sample$T0_only, cor_sample$Cross_timepoint)
d_combined_t1 <- cohen.d(cor_sample$T1_only, cor_sample$Cross_timepoint)

wilcox.test(cor_sample$T0_only, cor_sample$Cross_timepoint, paired = FALSE)
wilcox.test(cor_sample$T1_only, cor_sample$Cross_timepoint, paired = FALSE)
wilcox.test(cor_sample$T0_only, cor_sample$T1_only, paired = FALSE)

t.test(cor_sample$T0_only, cor_sample$Cross_timepoint, paired = FALSE)
t.test(cor_sample$T1_only, cor_sample$Cross_timepoint, paired = FALSE)
t.test(cor_sample$T0_only, cor_sample$T1_only, paired = FALSE)

# Calculate y positions that align with bracket levels
y_max <- max(cor_long_sample$value)
bracket_levels <- c(
  y_max * 1.15,  # Top bracket (Combined vs Cross_timepoint)
  y_max * 1.10,  # Second bracket (T0_only vs T1_only)  
  y_max * 1.05,  # Third bracket (Combined vs T0_only)
  y_max * 1.00   # Fourth bracket (Combined vs T1_only)
)

plot <- ggplot(cor_long_sample, aes(x = variable, y = value)) +
  rasterise(geom_boxplot(), dpi = 300) +
  
  # Keep your existing statistical comparisons
  stat_compare_means(
    comparisons = list(
      c("Combined", "Cross_timepoint"),  # Position 1 & 4
      c("T0_only", "T1_only"),          # Position 2 & 3
      c("Combined", "T0_only"),         # Position 1 & 2
      c("Combined", "T1_only")          # Position 1 & 3
    ),
    method = "wilcox.test",
    label = "p.format",
    step.increase = 0.05  # Reduced spacing to match your plot
  ) +
  
  # Add Cohen's d aligned with brackets
  annotate("text", x = 2.5, y = bracket_levels[1] + 0.02, 
           label = sprintf("d = %.3f", d_combined_cross$estimate), 
           size = 3, color = "blue", hjust = 0.5) +
  
  annotate("text", x = 2.5, y = bracket_levels[2] + 0.02, 
           label = sprintf("d = %.3f", d_t0_t1$estimate), 
           size = 3, color = "blue", hjust = 0.5) +
  
  annotate("text", x = 1.5, y = bracket_levels[3] + 0.02, 
           label = sprintf("d = %.3f", d_combined_t0$estimate), 
           size = 3, color = "blue", hjust = 0.5) +
  
  annotate("text", x = 2, y = bracket_levels[4] + 0.02, 
           label = sprintf("d = %.3f", d_combined_t1$estimate), 
           size = 3, color = "blue", hjust = 0.5) +
  
  # Keep your existing overall comparison
  stat_compare_means(method = "kruskal.test", 
                    label.y = max(cor_long_sample$value) * 0.1) +
  
  labs(title = "Correlation Distributions with Statistical Tests", 
       x = "Analysis Type", y = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Expand y-axis to accommodate the annotations
  coord_cartesian(ylim = c(min(cor_long_sample$value), y_max * 1.2))

ggsave("/hpc/users/hoangd02/www/plots/lbp/CorrelationDistributions_rasterized_with_wilcoxontest.pdf", plot, width = 10, height = 12)

#With your massive dataset (449M points), 
#you'll likely see mostly **** because even tiny differences 
#become "statistically significant" with huge sample sizes. 
#That's why effect sizes (Cohen's d) are more meaningful than p-values for your data - 
#they tell you if the differences are practically important, not just statistically detectable.

##try 2
## ---- settings you can tweak ----
paired_tests <- FALSE        # set TRUE if rows are matched by subject
p_adjust     <- "BH"       # use "BH" for FDR if you like
x_levels     <- c("T0_only","T1_only","Cross_timepoint")

## ---- 1) long format ----
library(dplyr)

cor_long <- cor_sample %>%
  rename(Group = variable, Correlation = value) %>%
  mutate(Group = factor(Group, levels = x_levels))

## ---- 2) define comparisons ----
comparisons <- list(
  c("T0_only", "Cross_timepoint"),
  c("T1_only", "Cross_timepoint"),
  c("T0_only", "T1_only")
)

## ---- 3) Wilcoxon p-values (paired or unpaired) ----
library(rstatix)
wilx_df <- cor_long %>%
  ungroup() %>%
  pairwise_wilcox_test(
    Correlation ~ Group,
    p.adjust.method = p_adjust,
    paired = paired_tests,
    comparisons = comparisons
  ) %>%
  select(group1, group2, p) %>%
  mutate(p_label = rstatix::p_format(p, accuracy = 0.001))
#hereee
## ---- 4) Cohen's d computed in the SAME order as wilcoxon ----
library(dplyr)
library(purrr)
library(effsize)

compute_d <- function(g1, g2) {
  x <- cor_long %>% filter(Group == g1) %>% pull(Correlation)
  y <- cor_long %>% filter(Group == g2) %>% pull(Correlation)
  effsize::cohen.d(x, y, hedges.correction = FALSE, pooled = FALSE)$estimate |> as.numeric()
}

d_df <- wilx_df %>%
  mutate(
    d = purrr::map2_dbl(group1, group2, compute_d),
    d_label = sprintf("d = %.2f", d)
  )

# 5) annotation table: labels + vertical positions
yr  <- range(cor_long$Correlation, na.rm = TRUE)
pad <- 0.12 * diff(yr)

annot_df <- d_df %>%   # <- use d_df directly (already has d & p)
  mutate(
    label      = paste0(d_label, ", p ", p_label),
    y.position = seq(yr[2] + pad, by = pad, length.out = n())
  )

# Ensure plotting order uses current factor levels
lvls <- levels(cor_long$Group)
cor_long <- cor_long %>% mutate(Group = factor(Group, levels = lvls))

# Map group names to x positions
xpos <- setNames(seq_along(lvls), lvls)

annot_plot <- annot_df %>%
  mutate(
    xmin = xpos[group1],
    xmax = xpos[group2],
    xmid = (xmin + xmax) / 2,
    y    = y.position,
    ycap = y + pad * 0.18
  )

# 6) base plot
ylim_top <- max(yr[2], max(annot_plot$ycap, na.rm = TRUE)) + pad * 0.2

p <- ggplot(cor_long, aes(x = Group, y = Correlation, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian(ylim = c(yr[1], ylim_top), clip = "off") +
  theme_minimal(base_size = 12) +
  labs(title = "Correlations by Group", x = NULL, y = "Correlation") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 30, 10, 10)
  )

# 7) brackets + labels
p <- p +
  geom_segment(data = annot_plot,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_plot,
               aes(x = xmin, xend = xmin, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_plot,
               aes(x = xmax, xend = xmax, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_label(data = annot_plot,
             aes(x = xmid, y = ycap, label = label),
             inherit.aes = FALSE, size = 3.5,
             label.size = 0, fill = "white", alpha = 0.92)

ggsave("/hpc/users/hoangd02/www/plots/lbp/timepoint_cor_1million_per_group_with_wilcoxontest.pdf", p, width = 10, height = 8)

#ggsave("/hpc/users/hoangd02/www/plots/lbp/timepoint_cor_rasterized_with_wilcoxontest.pdf", p, width = 10, height = 8)

####looking at the top percentiles
cor_t0_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t0_only_20250923.rds"))
cor_t1_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t1_only_20250923.rds"))
cor_cross <- readRDS(file.path(out_dir, "spearman_cross_timepoint_20250923.rds"))

probs <- c(`99.95th` = 0.9995, `99.995th` = 0.99995)

## Helper: thresholds on full (rectangular) matrix
get_thresholds <- function(mat, probs) {
  vals <- as.vector(abs(mat))
  quantile(vals, probs = probs, names = TRUE)
}

## Helper: create top-X% matrices (values < cutoff -> NA)
make_top_mats <- function(mat, cuts_named) {
  m_abs <- abs(mat)
  out <- lapply(names(cuts_named), function(nm) {
    thr <- cuts_named[[nm]]
    m <- m_abs
    m[m < thr] <- NA_real_
    m
  })
  names(out) <- names(cuts_named)
  out
}

## Apply to your three matrices
cuts_t0    <- get_thresholds(cor_t0_only, probs); 
# cuts_t0
#    99.95%   99.995% 
# 0.3940824 0.4536308 
cuts_t1    <- get_thresholds(cor_t1_only, probs);
# cuts_t1
#    99.95%   99.995% 
# 0.4116068 0.4730038 
cuts_cross <- get_thresholds(cor_cross,   probs)
# cuts_cross
#    99.95%   99.995% 
# 0.2882226 0.3352357 

top_t0    <- make_top_mats(cor_t0_only, cuts_t0)     # list of matrices: $`99.95th`, $`99.995th`
top_t0$99.95%   top_t0$99.995%
dim(top_t0$`99.95%`) # 21046 21356
dim(top_t0$`99.995%`) #21046 21356

top_t1    <- make_top_mats(cor_t1_only,  cuts_t1)    # same structure
top_cross <- make_top_mats(cor_cross,    cuts_cross) # same structure

## If you also want the raw vectors >= cutoff (for violin/boxplots later):
get_top_vectors <- function(mat, cuts_named) {
  v <- as.vector(abs(mat))
  lapply(cuts_named, function(thr) v[v >= thr])
}

topvec_t0    <- get_top_vectors(cor_t0_only, cuts_t0)       # $`99.95th`, $`99.995th`
topvec_t1    <- get_top_vectors(cor_t1_only, cuts_t1)
topvec_cross <- get_top_vectors(cor_cross, cuts_cross)

topvec_t0$99.95%   topvec_t0$99.995%  

# How many correlations survived each cutoff
length(topvec_t0[["99.95%"]]) # 224870
length(topvec_t0[["99.995%"]]) # 22473

head(topvec_t0[["99.95%"]])

summary(topvec_t0[["99.95%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3941  0.4021  0.4132  0.4206  0.4311  0.8534

summary(topvec_t1[["99.95%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.4116  0.4198  0.4311  0.4388  0.4496  0.8947 

summary(topvec_cross[["99.95%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.2882  0.2944  0.3029  0.3094  0.3170  0.8529 

summary(topvec_t0[["99.995%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.4537  0.4605  0.4700  0.4785  0.4866  0.8534 

summary(topvec_t1[["99.995%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.4730  0.4800  0.4898  0.4986  0.5061  0.8947 

summary(topvec_cross[["99.995%"]])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3352  0.3408  0.3488  0.3599  0.3627  0.8529

#save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_gene_pairs_all_genes.RData")
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_gene_pairs_all_genes.RData")

##########################
######## PLOTTING ########
##########################
## extract top entries with gene names
get_top_entries <- function(mat, cuts_named) {
  m_abs <- abs(mat)
  out <- lapply(cuts_named, function(thr) {
    # indices where correlation ≥ threshold
    keep_idx <- which(m_abs >= thr, arr.ind = TRUE)
    data.frame(
      blood_gene = rownames(mat)[keep_idx[, 1]],
      brain_gene = colnames(mat)[keep_idx[, 2]],
      correlation = mat[keep_idx]
    )
  })
  names(out) <- names(cuts_named)
  out
}

topentries_t0    <- get_top_entries(cor_t0_only, cuts_t0)
topentries_t1    <- get_top_entries(cor_t1_only, cuts_t1)
topentries_cross <- get_top_entries(cor_cross,   cuts_cross)

topentries_t0$99.95%   topentries_t0$99.995% 

##just looking at 99.95th percentile for now 
t0_top <- topentries_t0[["99.95%"]]
dim(t0_top)
#224870      3
head(t0_top)
          blood_gene         brain_gene correlation
1 ENSG00000149418.11 ENSG00000000003.14   0.4186721
2 ENSG00000154447.15 ENSG00000000003.14   0.4078119
3  ENSG00000240065.8 ENSG00000000003.14  -0.3980624
4  ENSG00000282420.1 ENSG00000000003.14  -0.4005615

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t0_top$blood_gene == t0_top$brain_gene, na.rm = TRUE); same_gene_rows
#322

## 2. How many unique blood genes
n_unique_blood <- length(unique(t0_top$blood_gene)); n_unique_blood
#20085

## 3. How many unique brain genes
n_unique_brain <- length(unique(t0_top$brain_gene)); n_unique_brain
#19747

t1_top <- topentries_t1[["99.95%"]]
dim(t1_top)
#224891      3
head(t1_top)
          blood_gene         brain_gene correlation
1 ENSG00000135842.17 ENSG00000000003.14   0.4455757
2 ENSG00000143226.13 ENSG00000000003.14   0.4257991
3 ENSG00000168228.15 ENSG00000000003.14   0.4138899
4 ENSG00000196663.16 ENSG00000000003.14   0.4163581

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t1_top$blood_gene == t1_top$brain_gene, na.rm = TRUE); same_gene_rows
# 310

## 2. How many unique blood genes
n_unique_blood <- length(unique(t1_top$blood_gene)); n_unique_blood
# 19447

## 3. How many unique brain genes
n_unique_brain <- length(unique(t1_top$brain_gene)); n_unique_brain
# 18910

#Which gene pairs are in T0 but not T1 and vice versa?
# Add a unique pair ID to each
t0_top$pair_id <- paste(t0_top$blood_gene, t0_top$brain_gene, sep = "_")
t1_top$pair_id <- paste(t1_top$blood_gene, t1_top$brain_gene, sep = "_")

# In T0 but not in T1
t0_only <- t0_top[!t0_top$pair_id %in% t1_top$pair_id, ]
dim(t0_only) #224186      4

# In T1 but not in T0
t1_only <- t1_top[!t1_top$pair_id %in% t0_top$pair_id, ]
dim(t1_only) #224207      4

overlap <- t0_top[t0_top$pair_id %in% t1_top$pair_id, ]
dim(overlap) #684   4

# Overlap at the gene level, not pair level
# Sometimes pairs don’t overlap but the genes do.
# Example: gene X might appear in T0 pairs and T1 pairs, but not with the same partner.

blood_overlap <- intersect(unique(t0_only$blood_gene), unique(t1_only$blood_gene))
brain_overlap <- intersect(unique(t0_only$brain_gene), unique(t1_only$brain_gene))

length(blood_overlap) #18499
length(brain_overlap) #17408

### HYPERGEOMETRIC PROBLEM 
N <- 21046L * 21356L
k <- floor(0.0005 * N)   # 224,729
K <- k

# P(X >= 684) where X ~ Hypergeometric(N, K, k)
pval <- phyper(q = 684 - 1, m = K, n = N - K, k = k, lower.tail = FALSE)

# Also report expectation & SD
mu <- k * (K / N)
sd <- sqrt(k * (K / N) * (1 - K / N) * ((N - k) / (N - 1)))

list(pval = pval, mean = mu, sd = sd)
# $pval
# [1] 2.184709e-291

# $mean
# [1] 112.3644

# $sd
# [1] 10.59491

###### OR 
# Parameters
n_total_pairs <- 21046 * 21356  # 449,374,576 total gene pairs
top_percent <- 0.0005  # 0.05%
n_top <- round(n_total_pairs * top_percent)  # number in top 0.05%

# Observed overlap
n_overlap_observed <- 684

# Simulation parameters
n_simulations <- 1000000

# Run simulations
set.seed(123)  # for reproducibility
null_overlaps <- replicate(n_simulations, {
  # Randomly sample n_top pairs from all possible pairs (set 1)
  sample1 <- sample(1:n_total_pairs, n_top, replace = FALSE)
  
  # Randomly sample n_top pairs from all possible pairs (set 2)
  sample2 <- sample(1:n_total_pairs, n_top, replace = FALSE)
  
  # Count overlap
  length(intersect(sample1, sample2))
})

# Analyze results
mean_overlap <- mean(null_overlaps); mean_overlap # 112.266
sd_overlap <- sd(null_overlaps); sd_overlap #10.43025
p_value <- sum(null_overlaps >= n_overlap_observed) / n_simulations; p_value #0 or very small


#######HYPERGEOMETRIC TEST | exact, no sim needed
# This is the exact analytical solution
# Total pairs, top pairs in each set, observed overlap
p_value_exact <- 1 - phyper(
  q = n_overlap_observed - 1,  # observed overlap - 1
  m = n_top,                    # pairs in T0 top set
  n = n_total_pairs - n_top,    # pairs NOT in T0 top set
  k = n_top                     # pairs selected for T1 top set
)
p_value_exact #0

# Expected overlap under null
expected_overlap <- n_top * (n_top / n_total_pairs)
expected_overlap #112.36 


#######################################################enrichment ~ 
# Enrichment of same-gene pairs within a single observed set
# Compare enrichment of same-gene pairs between two sets (e.g., T0 vs T1)
compare_same_gene_sets <- function(set_A, set_B, label_A = "T0", label_B = "T1", add_haldane = FALSE) {
  req_cols <- c("blood_gene", "brain_gene")
  stopifnot(all(req_cols %in% names(set_A)),
            all(req_cols %in% names(set_B)))
  
  # Count same vs different within each set
  a_same <- sum(set_A$blood_gene == set_A$brain_gene, na.rm = TRUE)
  a_diff <- nrow(set_A) - a_same
  
  b_same <- sum(set_B$blood_gene == set_B$brain_gene, na.rm = TRUE)
  b_diff <- nrow(set_B) - b_same
  
  # 2x2 table
  mat <- matrix(c(a_same, a_diff,
                  b_same, b_diff),
                nrow = 2, byrow = TRUE,
                dimnames = list(Set = c(label_A, label_B),
                                PairType = c("SameGene","DiffGene")))
  
  mat_for_test <- mat
  if (add_haldane) mat_for_test <- mat_for_test + 0.5
  
  ft <- fisher.test(mat_for_test)
  
  list(
    table = mat,
    counts = list(
      A = list(label = label_A, same = a_same, diff = a_diff, total = a_same + a_diff),
      B = list(label = label_B, same = b_same, diff = b_diff, total = b_same + b_diff)
    ),
    fisher = ft,
    odds_ratio = unname(ft$estimate),
    conf_int = ft$conf.int,
    p_value = ft$p.value
  )
}

res_t0_vs_t1 <- compare_same_gene_sets(t0_top, t1_top, label_A = "T0", label_B = "T1")
res_t0_vs_t1
# odds ratio 
#  1.038873 
# $p_value
# [1] 0.633173

res_t1_vs_t0 <- compare_same_gene_sets(t1_top, t0_top, label_A = "T1", label_B = "T0")
res_t1_vs_t0
# odds ratio 
#  0.9625816 
# $p_value
# [1] 0.633173





##########################
######## PLOTTING ########
##########################

## ---- 0) settings ----
paired_tests <- FALSE
p_adjust     <- "BH"
x_levels     <- c("T0_only", "T1_only", "Cross_timepoint")

## ---- 1) long format from top 99.95% vectors ----
library(dplyr)
library(tidyr)

cor_long <- list(
  T0_only        = topvec_t0[["99.9995%"]],
  T1_only        = topvec_t1[["99.9995%"]],
  Cross_timepoint = topvec_cross[["99.9995%"]]
) %>%
  enframe(name = "Group", value = "Correlation") %>%
  unnest(Correlation) %>%
  mutate(Group = factor(Group, levels = x_levels))

# cor_long <- list(
#   T0_only        = topvec_t0[["99.95%"]],
#   T1_only        = topvec_t1[["99.95%"]],
#   Cross_timepoint = topvec_cross[["99.95%"]]
# ) %>%
#   enframe(name = "Group", value = "Correlation") %>%
#   unnest(Correlation) %>%
#   mutate(Group = factor(Group, levels = x_levels))

## ---- 2) define comparisons ----
comparisons <- list(
  c("T0_only", "Cross_timepoint"),
  c("T1_only", "Cross_timepoint"),
  c("T0_only", "T1_only")
)

## ---- 3) Wilcoxon p-values ----
library(rstatix)

wilx_df <- cor_long %>%
  ungroup() %>%
  pairwise_wilcox_test(
    Correlation ~ Group,
    p.adjust.method = p_adjust,
    paired = paired_tests,
    comparisons = comparisons
  ) %>%
  select(group1, group2, p) %>%
  mutate(p_label = rstatix::p_format(p, accuracy = 0.001))

## ---- 4) Cohen’s d ----
library(purrr)
library(effsize)

compute_d <- function(g1, g2) {
  x <- cor_long %>% filter(Group == g1) %>% pull(Correlation)
  y <- cor_long %>% filter(Group == g2) %>% pull(Correlation)
  effsize::cohen.d(x, y, hedges.correction = FALSE, pooled = FALSE)$estimate |> as.numeric()
}

d_df <- wilx_df %>%
  mutate(
    d = map2_dbl(group1, group2, compute_d),
    d_label = sprintf("d = %.2f", d)
  )

## ---- 5) annotation data frame ----
yr  <- range(cor_long$Correlation, na.rm = TRUE)
pad <- 0.12 * diff(yr)

annot_df <- d_df %>%
  mutate(
    label      = paste0(d_label, ", p ", p_label),
    y.position = seq(yr[2] + pad, by = pad, length.out = n())
  )

lvls <- levels(cor_long$Group)
xpos <- setNames(seq_along(lvls), lvls)

annot_plot <- annot_df %>%
  mutate(
    xmin = xpos[group1],
    xmax = xpos[group2],
    xmid = (xmin + xmax) / 2,
    y    = y.position,
    ycap = y + pad * 0.18
  )

## ---- 6) ggplot (same as your scaffold) ----
library(ggplot2)

ylim_top <- max(yr[2], max(annot_plot$ycap, na.rm = TRUE)) + pad * 0.2

p <- ggplot(cor_long, aes(x = Group, y = Correlation, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian(ylim = c(yr[1], ylim_top), clip = "off") +
  theme_minimal(base_size = 12) +
  labs(title = "Top 99.95% Correlations", x = NULL, y = "Correlation") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 30, 10, 10)
  ) +
  geom_segment(data = annot_plot,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_plot,
               aes(x = xmin, xend = xmin, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_plot,
               aes(x = xmax, xend = xmax, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_label(data = annot_plot,
             aes(x = xmid, y = ycap, label = label),
             inherit.aes = FALSE, size = 3.5,
             label.size = 0, fill = "white", alpha = 0.92)
ggsave("/hpc/users/hoangd02/www/plots/lbp/timepoint_cor_99.9995_with_wilcoxontest.pdf", p, width = 10, height = 8)

ggsave("/hpc/users/hoangd02/www/plots/lbp/timepoint_cor_99.95_with_wilcoxontest.pdf", p, width = 10, height = 8)

###################    FORM5 NULL ~~ this has absolute cor_mat ~~ see if i did this for form3!
################### CODE FROM LBP_blood_brain_QC_20250603.rmd
#base_null <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form0.txt", data.table=FALSE)
form5_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form5_no_residID.txt", data.table=FALSE)
form5_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/null_summaries_blood_form5_with_ID.txt", data.table=FALSE)
##with residID do not exist 

#base_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form0.txt", data.table=FALSE)
form5_no_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_no_residID.txt", data.table=FALSE)
# 8       p99.95     0.23842920             0.133 ~224k
# 9      p99.995     0.27732301             0.098 ~22k
# 10    p99.9995     0.44724055             0.000 ~2k
# 11   p99.99995     0.74321385             0.000 ~224

form5_with_residID_summary <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_with_ID.txt", data.table=FALSE)
# 8       p99.95     0.22645702             0.458
# 9      p99.995     0.26274126             0.450
# 10    p99.9995     0.30767422             0.097
# 11   p99.99995     0.41535861             0.000

##in ~/qc2/
null_summaries_blood_form0.txt
null_summaries_blood_form3_no_residID.txt
null_summaries_blood_form3_with_residID.txt
null_summaries_blood_form5_no_residID.txt
null_summaries_blood_form5_with_ID.txt

comparison_of_null_and_real_data_summaries_blood_form0.txt
comparison_of_null_and_real_data_summaries_blood_form3_no_residID.txt
comparison_of_null_and_real_data_summaries_blood_form3_with_residID.txt
comparison_of_null_and_real_data_summaries_blood_form5_no_residID.txt
comparison_of_null_and_real_data_summaries_blood_form5_with_ID.txt










# Summary stats
cor_comparison <- data.frame(
  Combined = as.vector(cor_combined),
  T0_only = as.vector(cor_t0_only),
  T1_only = as.vector(cor_t1_only),
  Cross_timepoint = as.vector(cor_cross)
)

stats <- cor_comparison %>% 
  summarise_all(list(
    mean = mean, 
    median = median, 
    q99 = ~quantile(., 0.99),
    q99.5 = ~quantile(., 0.995),
    q99.95 = ~quantile(., 0.9995)
  ))
stats

####this code takes too long 
library(broom)
library(purrr)
# Function to perform pairwise comparisons
pairwise_tests <- function(data) {
  comparisons <- combn(names(data), 2, simplify = FALSE)
  
  results <- map_dfr(comparisons, function(pair) {
    group1 <- data[[pair[1]]]
    group2 <- data[[pair[2]]]
    
    # Wilcoxon test (non-parametric)
    wilcox_result <- wilcox.test(group1, group2)
    
    # Effect size (Cohen's d)
    cohens_d <- (mean(group1) - mean(group2)) / 
                sqrt((var(group1) + var(group2)) / 2)
    
    data.frame(
      comparison = paste(pair[1], "vs", pair[2]),
      wilcox_p = wilcox_result$p.value,
      cohens_d = cohens_d,
      mean_diff = mean(group1) - mean(group2),
      median_diff = median(group1) - median(group2)
    )
  })
  
  return(results)
}

# Run pairwise tests
test_results <- pairwise_tests(cor_comparison)
print(test_results)



##########i need to check whether this is equivalent to the no_swap_abs 

##no residID NOT ABSOLUTE VALUE - this is 73 blood samples and 146 brain samples (2 timepoints)
no_swap_abs  <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_0pct_nonabs_no_residID.rds")
dim(no_swap_abs) #dim(no_swap_abs)
no_swap_abs[1:5,1:5]

all_swap_abs <- readRDS("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap/spearman_swap_100pct_nonabs_no_residID.rds")


## ## ## ## ## ##
##how consistent the blood–brain correlations are across the two timepoints.
# Flatten upper triangles
vec_t0 <- cor_t0[upper.tri(cor_t0)]
vec_t1 <- cor_t1[upper.tri(cor_t1)]

# Compare similarity across timepoints
cor(vec_t0, vec_t1, method = "spearman")  # overall stability ~~ 0.002713632

## ## ## ## ## ##
genewise_stability <- sapply(1:nrow(cor_t0), function(i) {
  cor(cor_t0[i, ], cor_t1[i, ], method = "spearman")
})

# Assuming you already have genewise_stability as a numeric vector
df_stability <- data.frame(stability = genewise_stability)

plot <- ggplot(df_stability, aes(x = stability)) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "white",
    alpha = 0.7
  ) +
  labs(
    title = "Gene-wise stability across timepoints",
    x = "Spearman correlation (t0 vs t1 profiles)",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 14)

ggsave("/hpc/users/hoangd02/www/plots/lbp/Gene-wise_stability_across_timepoints.pdf", plot, width = 10, height = 6)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/Gene-wise_stability_across_timepoints.pdf

## ## ## ## ## ##
threshold_t0 <- quantile(abs(vec_t0), 0.9995)
threshold_t1 <- quantile(abs(vec_t1), 0.9995)

top_t0 <- which(abs(cor_t0) >= threshold_t0, arr.ind = TRUE)
top_t1 <- which(abs(cor_t1) >= threshold_t1, arr.ind = TRUE)

# Overlap between top sets
length(intersect(paste(top_t0[,1], top_t0[,2]),
                 paste(top_t1[,1], top_t1[,2]))) #673

save.image


################## SAME GENES ##############################
blood_genes <- rownames(cor_t0)   # 21,046 blood genes
brain_genes <- colnames(cor_t0)   # 21,356 brain genes
common_genes <- intersect(blood_genes, brain_genes)
length(common_genes)  # e.g., 17,533 overlapping genes
same_gene_cor <- sapply(common_genes, function(g) cor_t0[g, g])
summary(same_gene_cor)

# wrap into dataframe
df_same <- data.frame(correlation = same_gene_cor)

plot <- ggplot(df_same, aes(x = correlation)) +
  geom_histogram(
    bins = 50,              # same as breaks = 50
    fill = "steelblue",
    color = "white",
    alpha = 0.7
  ) +
  labs(
    title = "Same-gene blood brain correlations",
    x = "Spearman correlation",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)
ggsave("/hpc/users/hoangd02/www/plots/lbp/testtesttest.pdf", plot, width = 10, height = 6)

############################################
############# SAMPLE-LEVEL ################# NEW OCT 2 WITH CROSS TIME POINT START HERE
############################################
# brain full vs blood form5 no_residID
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
blood <- blood_form5_no_residID   # genes x blood_samples
brain <- brain_full_no_residID  # genes x brain_samples

##run the temporal component of dictionary
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)

two_pairs <- filter(dictionary, number_of_pair_brain == "2_pairs") #146, correct
dim(two_pairs) #146   8

two_pairs$abs_concordance <- NULL
two_pairs$concordance <- NULL

# --- Helper: correlate one blood–brain column pair on common genes ---
get_corr <- function(blood_id, brain_id, blood_mat, brain_mat) {
  common <- intersect(rownames(blood_mat), rownames(brain_mat))
  if (length(common) < 3) return(NA_real_)
  x <- blood_mat[common, blood_id, drop = TRUE]
  y <- brain_mat[common, brain_id, drop = TRUE]
  suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

# --- 1) Keep only rows that exist in matrices ---
two_pairs_clean <- two_pairs %>%
  filter(SAMPLE_ISMMS_blood %in% colnames(blood),
         SAMPLE_ISMMS_brain %in% colnames(brain)) %>%
  filter(timepoint %in% c(0, 1))

# --- 2) Same-timepoint correlations (T0_only, T1_only) ---
within_same <- two_pairs_clean %>%
  rowwise() %>%
  mutate(
    corr = get_corr(SAMPLE_ISMMS_blood, SAMPLE_ISMMS_brain, blood, brain),
    type = paste0("T", timepoint, "_only")
  ) %>%
  ungroup() %>%
  select(IID_ISMMS, type, corr) %>%
  distinct()

# --- 3) Cross-timepoint correlations (within person; avoid duplicates) ---
within_cross <- two_pairs_clean %>%
  inner_join(two_pairs_clean, by = "IID_ISMMS", suffix = c("_b", "_br")) %>%
  filter(timepoint_b != timepoint_br) %>%
  rowwise() %>%
  mutate(corr = get_corr(SAMPLE_ISMMS_blood_b, SAMPLE_ISMMS_brain_br, blood, brain)) %>%
  ungroup() %>%
  distinct(IID_ISMMS, corr) %>%
  mutate(type = "Cross_timepoint")

# --- 4) Long table (no pivot_wider) ---
cor_long <- bind_rows(within_same, within_cross)

table(cor_long$type)
# Cross_timepoint         T0_only         T1_only 
#             146              73              73 

######################################## NOV 6 ########################################
## just looking at signs 
library(dplyr)
library(tidyr)

summary_tbl <- cor_long %>%
  filter(type %in% c("T0_only", "T1_only")) %>%
  pivot_wider(names_from = type, values_from = corr) %>%
  mutate(
    T0_sign = if_else(T0_only > 0, "positive", "negative"),
    T1_sign = if_else(T1_only > 0, "positive", "negative")
  )

# (1) Counts of positives/negatives at each timepoint
t0_counts <- summary_tbl %>% count(T0_sign, name = "n_T0")
t1_counts <- summary_tbl %>% count(T1_sign, name = "n_T1")

# (2) Cross-tab of T0 vs T1 sign
cross_counts <- summary_tbl %>% count(T0_sign, T1_sign)

list(
  counts_T0 = t0_counts,
  counts_T1 = t1_counts,
  cross_table = cross_counts
)

$counts_T0
# A tibble: 2 × 2
  T0_sign   n_T0
  <chr>    <int>
1 negative    24
2 positive    49

$counts_T1
# A tibble: 2 × 2
  T1_sign   n_T1
  <chr>    <int>
1 negative    29
2 positive    44

$cross_table
# A tibble: 4 × 3
  T0_sign  T1_sign      n
  <chr>    <chr>    <int>
1 negative negative    10 / 73 = 13.7%
2 negative positive    14 / 73 = 19.2%
3 positive negative    19 / 73 = 26.0%
4 positive positive    30 / 73 = 41.1%

# Given how many people are positive at T0 and at T1 overall, 
# is there more concordance (same sign) than expected if signs were independent?

#Build a 2×2 contingency table
tab <- matrix(
  c(30, 19, 14, 10),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "T0_sign" = c("positive", "negative"),
    "T1_sign" = c("positive", "negative")
  )
)
tab
#           T1_sign
# T0_sign    positive negative
#   positive       30       19
#   negative       14       10

#proportion of patients whose sign stayed the same between T0 and T1
concordant <- sum(tab[1,1], tab[2,2])   # same sign
discordant <- sum(tab[1,2], tab[2,1])   # opposite sign
concordance_rate <- concordant / (concordant + discordant)
concordance_rate #0.5479452

chisq.test(tab)
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  tab
# X-squared = 2.6894e-31, df = 1, p-value = 1

fisher.test(tab) #being positive or negative at T0 does not necessarily tells you that you're positive or negative at T1
# 	Fisher's Exact Test for Count Data
# data:  tab
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.366536 3.391540
# sample estimates:
# odds ratio 
#   1.125951 

# The odds of being positive at T1 are about the same regardless of T0 sign — 
# completely consistent with chance.

vcd::assocstats(tab)  # for Cramér’s V and other association measures
#                       X^2 df P(> X^2)
# Likelihood Ratio 0.056091  1  0.81279
# Pearson          0.056237  1  0.81255

# Phi-Coefficient   : 0.028 
# Contingency Coeff.: 0.028 
# Cramer's V        : 0.028 

#About half the subjects retained the same sign across timepoints 
#very close to what you'd expect by random chance.

# There is no statistically significant association between the sign of correlations at T0 and T1.
# The proportion of individuals who stayed positive/negative (≈55%) is what you'd expect by chance if T0 and T1 signs were independent.

##In other words, knowing whether someone was “positive” or “negative” at T0 gives you essentially no information about what they’ll be at T1 
## at least at the level of sign direction.

###making a new cor_diff df to see which individuals have elevated vs deflated corr
### POSITIVE diff = higher T1 than T0 | T1 > T0
### NEGATIVE diff = lower T1 than T0

cor_diff <- cor_long %>%
  mutate(corr = as.numeric(corr)) %>%
  pivot_wider(
    names_from = type, values_from = corr,
    values_fn = ~ mean(as.numeric(.x), na.rm = TRUE)
  ) %>%
  mutate(diff = T1_only - T0_only)

head(cor_diff)
# A tibble: 6 × 5
  IID_ISMMS   T1_only T0_only Cross_timepoint    diff
  <chr>         <dbl>   <dbl>           <dbl>   <dbl>
1 PT-0018    0.0453   -0.0228        -0.00543  0.0680
2 PT-0021   -0.0352    0.0280         0.119   -0.0632
3 PT-0022   -0.0523    0.0916        -0.0494  -0.144 

summary(cor_diff$diff)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.42584 -0.13594 -0.03493 -0.01760  0.09984  0.56698

plot <- ggplot(cor_diff, aes(x = reorder(IID_ISMMS, diff), y = diff)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(diff, 2)), vjust = -0.1, size = 3) +
  labs(x = "Patient ID", y = "Difference in correlation", 
       title = "Difference in correlation per patient") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("/hpc/users/hoangd02/www/plots/lbp/test.pdf", plot, width = 18, height = 6)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/test.pdf

cor_diff <- as.data.frame(cor_diff)

###starting with the lowest T0
head(cor_diff[order(cor_diff$T0_only), ],10)
   IID_ISMMS     T1_only     T0_only Cross_timepoint       diff
19   PT-0051 -0.03014119 -0.35173931      0.26822319 0.32159812 #regress towards the mean ?
33   PT-0083  0.04278514 -0.25019047      0.07216529 0.29297561 #regress towards the mean ?
13   PT-0041 -0.02384490 -0.23447276      0.18130576 0.21062786 #regress towards the mean ?
70   PT-0183  0.35142438 -0.21555778     -0.09429534 0.56698216 ## strange
53   PT-0132 -0.09374664 -0.15142993      0.20065156 0.05768329 #regress towards the mean ?
72   PT-0193  0.02053329 -0.13151993      0.11126661 0.15205322 #regress towards the mean ?
7    PT-0030 -0.01148856 -0.11132791     -0.07074587 0.09983936 #regress towards the mean ?
12   PT-0040  0.15898845 -0.10786180      0.03539060 0.26685026
67   PT-0178  0.07510298 -0.09608550      0.06935169 0.17118847
65   PT-0169  0.22653374 -0.07840918      0.07080566 0.30494292

cor_diff[cor_diff$IID_ISMMS == "PT-0183", ] #BIG jump
#    IID_ISMMS   T1_only    T0_only Cross_timepoint      diff
# 70   PT-0183 0.3514244 -0.2155578     -0.09429534 0.5669822

dictionary[dictionary$IID_ISMMS == "PT-0183", ]
#     IID_ISMMS SAMPLE_ISMMS_blood SAMPLE_ISMMS_brain number_of_pair_brain
# 218   PT-0183   LBPSEMA4BLOOD508   LBPSEMA4BRAIN749              2_pairs
# 219   PT-0183   LBPSEMA4BLOOD027   LBPSEMA4BRAIN468              2_pairs
#     mymet_tissue_brain mymet_tissue_blood surgeryDate timepoint
# 218            R_Brain            R_Blood  2018-06-18         1
# 219            L_Brain            L_Blood  2018-05-21         0 

brain_metadata[brain_metadata$IID_ISMMS %in% c("PT-0183","PT-0051","PT-0132"), ] #STAR_Number_ for PT-0183 is higher
brain_metadata[brain_metadata$IID_ISMMS == "PT-0183", ]$mymet_rin_brain
#6.8 6.4

blood_metadata[blood_metadata$IID_ISMMS == "PT-0183", ]$mymet_rin_blood
#8.7 6.3

summary(brain_metadata$mymet_rin_brain)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  5.50    6.70    7.10    7.11    7.60    9.10

summary(blood_metadata$mymet_rin_blood)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  4.10    7.20    7.90    7.62    8.40   10.00


cor_diff[cor_diff$IID_ISMMS == "PT-0158", ]
#    IID_ISMMS    T1_only   T0_only Cross_timepoint       diff
# 61   PT-0158 -0.1206062 0.3052321     -0.01466998 -0.4258383


#########correlating T1_only, T0_only, and diff to all metadata variables
# join metrics to numeric metadata across subjects (BRAIN METADATA has both blood and brain)
joined <- cor_diff %>%
  select(IID_ISMMS, T1_only, T0_only, diff) %>%
  inner_join(brain_metadata %>% select(IID_ISMMS, where(is.numeric)),
             by = "IID_ISMMS")

# join metrics to numeric metadata across subjects BLOOD ONLY
# joined <- cor_diff %>%
#   select(IID_ISMMS, T1_only, T0_only, diff) %>%
#   inner_join(blood_only_metadata %>% select(IID_ISMMS, where(is.numeric)),
#              by = "IID_ISMMS")

# correlate each numeric metadata column with T0_only, T1_only, and diff
features <- setdiff(names(joined), c("IID_ISMMS","T1_only","T0_only","diff"))

cors_mat <- cor(
  joined[, features],
  joined[, c("T0_only","T1_only","diff")],
  use = "pairwise.complete.obs"  # handles NAs
)

# tidy df with feature names
cors_df <- data.frame(feature = rownames(cors_mat), cors_mat, row.names = NULL)
cors_df  # columns are feature, T0_only, T1_only, diff

###sort T0 from highest to lowest
head(cors_df[order(cors_df$diff, decreasing = TRUE), ],4) ##BRAIN
#                                  feature     T0_only    T1_only      diff
# 138                           GABA_brain -0.19204561 0.14635978 0.2358734
# 137                            GLU_brain -0.22371681 0.05102085 0.2019587
# 23  STAR_Uniquely_mapped_reads_pct_brain -0.08396409 0.13119338 0.1440828
# 115            mymet_rna_conc_ngul_brain -0.12840468 0.07401698 0.1431786

head(cors_df[order(cors_df$diff, decreasing = TRUE), ],4) ##BLOOD - LOOKS pretty good 
#                                             feature    T0_only      T1_only
# 130                       Monocyte_total_wilk_blood -0.1836395 -0.052867168
# 131                              Mono_DC_wilk_blood -0.1752020 -0.047767101
# 104            RNASeqMetrics_PCT_CODING_BASES_blood -0.1182940  0.005188722
# 24  STAR_pct_of_reads_mapped_to_multiple_loci_blood -0.0450354  0.095558568
#           diff
# 130 0.10803560
# 131 0.10469969
# 104 0.09350961
# 24  0.09263348

cor(brain_metadata[,c("GABA_brain","GLU_brain","ODC_brain")])
           GABA_brain  GLU_brain  ODC_brain
GABA_brain   1.000000  0.4299940 -0.4496800
GLU_brain    0.429994  1.0000000 -0.9882269
ODC_brain   -0.449680 -0.9882269  1.0000000


head(cors_df[order(cors_df$T0_only, decreasing = TRUE), ],10)
                                             feature   T0_only     T1_only
140                                        ODC_brain 0.2251452 -0.04495593
141                                         MG_brain 0.1786326  0.04097134
155        STAR_Number_of_reads_unmapped_other_blood 0.1632002 -0.09273624
252           RNASeqMetrics_MEDIAN_5PRIME_BIAS_blood 0.1477133 -0.06685436
274                         Monocyte_total_scp_blood 0.1429759  0.11191058
275                                Mono_DC_scp_blood 0.1420885  0.11530525
113 RNASeqMetrics_MEDIAN_5PRIME_TO_3PRIME_BIAS_brain 0.1354546 -0.04224841
277                              Monocytes_fig_blood 0.1282389  0.06052962
167           STAR_pct_of_reads_unmapped_other_blood 0.1214658 -0.09005433
80       InsertSizeMetrics_WIDTH_OF_10_PERCENT_brain 0.1208332 -0.14518005
           diff
140 -0.19935341
141 -0.11146144
155 -0.18116212
252 -0.15356041
274 -0.04099505
275 -0.03824844
113 -0.12920189
277 -0.06105340
167 -0.14765287
80  -0.18076572

head(cors_df[order(cors_df$T1_only, decreasing = TRUE), ],10)
                                               feature      T0_only   T1_only
104               RNASeqMetrics_PCT_CODING_BASES_brain  0.032529570 0.2504958
120                     GcBiasMetrics_GC_DROPOUT_brain  0.093548802 0.2254707
16                  STAR_Number_of_splices_AT_AC_brain -0.005194400 0.2130674
17       STAR_Number_of_splices_Annotated__sjdb__brain  0.029463519 0.2066961
19                  STAR_Number_of_splices_GT_AG_brain  0.029352266 0.2058820
21                  STAR_Number_of_splices_Total_brain  0.029094916 0.2057797
18                  STAR_Number_of_splices_GC_AG_brain  0.016496963 0.2019923
94                    RNASeqMetrics_CODING_BASES_brain  0.026944145 0.2005608
99  RNASeqMetrics_NUM_R2_TRANSCRIPT_STRAND_READS_brain  0.007834061 0.1882375
2                         FEATURECOUNTS_Assigned_brain  0.006041853 0.1774490
          diff
104 0.12781852
120 0.06596354
16  0.13381982
17  0.10346677
19  0.10305560
21  0.10318982
18  0.11050330
94  0.10165185
99  0.10873684
2   0.10353065

#######################################################################################
#this is correct, there are 2 cross-timepoint
sub <- subset(cor_long, type != "Cross_timepoint")

t0_sample_wise <- filter(sub, timepoint == "0")
summary(t0_sample_wise$corr)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.35174 -0.02818  0.06250  0.06171  0.16032  0.35356 
t1_sample_wise <- filter(sub, timepoint == "1")
summary(t1_sample_wise$corr)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.22020 -0.02384  0.02628  0.04412  0.08659  0.39867 

cor_long <- cor_long %>%
  mutate(type = factor(type, levels = c("T0_only","T1_only","Cross_timepoint")))

# ---- Paired Wilcoxon (align by IID_ISMMS; no pivot_wider) ----
pair_T0_T1 <- filter(cor_long, type == "T0_only") %>%
  inner_join(filter(cor_long, type == "T1_only"), by = "IID_ISMMS",
             suffix = c("_T0", "_T1"))
pair_T0_X  <- filter(cor_long, type == "T0_only") %>%
  inner_join(filter(cor_long, type == "Cross_timepoint"), by = "IID_ISMMS",
             suffix = c("_T0", "_X"))
pair_T1_X  <- filter(cor_long, type == "T1_only") %>%
  inner_join(filter(cor_long, type == "Cross_timepoint"), by = "IID_ISMMS",
             suffix = c("_T1", "_X"))

wilx_df <- tibble(
  group1 = c("T0_only","T0_only","T1_only"),
  group2 = c("T1_only","Cross_timepoint","Cross_timepoint"),
  p = c(
    wilcox.test(pair_T0_T1$corr_T0, pair_T0_T1$corr_T1, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(pair_T0_X$corr_T0,  pair_T0_X$corr_X,  paired = FALSE, exact = FALSE)$p.value,
    wilcox.test(pair_T1_X$corr_T1,  pair_T1_X$corr_X,  paired = , exact = FALSE)$p.value
  )
) %>%
  mutate(p_label = case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("= %.3f", p)
  ))

# ---- Cohen's d like your old code (unpaired, pooled = FALSE) ----
compute_d <- function(g1, g2) {
  x <- cor_long %>% filter(type == g1) %>% pull(corr)
  y <- cor_long %>% filter(type == g2) %>% pull(corr)
  effsize::cohen.d(x, y, hedges.correction = FALSE, pooled = FALSE)$estimate |> as.numeric()
}

library(dplyr)
library(ggplot2)
library(effsize)
library(purrr)
library(tibble)

d_df <- wilx_df %>%
  mutate(d = map2_dbl(group1, group2, compute_d),
         d_label = sprintf("d = %.2f", d))

# ---- Annotation rows (brackets + labels) ----
yr  <- range(cor_long$corr, na.rm = TRUE)
pad <- 0.12 * diff(yr)

lvls <- levels(cor_long$type)
xpos <- setNames(seq_along(lvls), lvls)

annot_df <- d_df %>%
  mutate(
    label      = paste0(d_label, ", p ", p_label),
    y.position = seq(yr[2] + pad, by = pad, length.out = n()),
    xmin       = xpos[group1],
    xmax       = xpos[group2],
    xmid       = (xmin + xmax) / 2,
    y          = y.position,
    ycap       = y + pad * 0.25
  )

ylim_top <- max(yr[2], max(annot_df$ycap, na.rm = TRUE)) + pad * 0.2

# ---- Plot (violin + box); colors map to 'type' like your old code ----
p <- ggplot(cor_long, aes(x = type, y = corr, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.05))  # 15% bottom, 5% top
  ) +
  coord_cartesian(ylim = c(yr[1], ylim_top), clip = "off") +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Spearman correlation") +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 10)
  )

# Brackets + labels (no warnings: use small data + inherit.aes = FALSE)
p <- p +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmin, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = annot_df,
               aes(x = xmax, xend = xmax, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_label(data = annot_df,
             aes(x = xmid, y = ycap, label = label),
             inherit.aes = FALSE, size = 3.5,
             label.size = 0, fill = "white", alpha = 0.92)

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_cross_timepoint.pdf", p, width = 10, height = 6)

### FOR WCPG 
##Violin and boxplot
p <- ggplot(cor_long, aes(x = type, y = corr, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  coord_cartesian(ylim = c(yr[1], ylim_top), clip = "off") +
  theme_minimal() +
  labs(x = NULL, y = "Spearman correlation") +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size = 16),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 30, 10, 10)
  )

# Brackets + labels
p <- p +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmin, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = annot_df,
               aes(x = xmax, xend = xmax, y = y, yend = y - pad*0.06),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_label(data = annot_df,
             aes(x = xmid, y = ycap, label = label),
             inherit.aes = FALSE,
             size = 6,                # larger label text (poster-ready)
             label.size = 0,
             fill = NA)
             #fill = "white", alpha = 0.92)

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_cross_timepoint_wcpg.pdf", p, width = 8, height = 5)

#### Boxplot and matching points between T0 and T1
library(dplyr)
library(tidyr)
library(ggplot2)

# Keep only T0/T1, drop NA, and require both timepoints per subject
df_pair <- cor_long %>%
  filter(type %in% c("T0_only", "T1_only")) %>%            # exclude cross-timepoint or others
  filter(!is.na(corr)) %>%
  group_by(IID_ISMMS) %>%
  filter(n_distinct(type) == 2) %>%                        # keep only IDs with both T0 & T1
  ungroup() %>%
  mutate(type = factor(type, levels = c("T0_only", "T1_only"),
                       labels = c("T0", "T1")))

# Add annotation bar coordinates
annot_df <- data.frame(
  xmin = 1,   # T0
  xmax = 2,   # T1
  y = 0.45,   # height of the bar (adjust depending on your data range)
  label = "d = 0.15, p = 0.32"
)

p <- ggplot(df_pair, aes(x = type, y = corr)) +
  # Boxplots (one per timepoint)
  geom_boxplot(aes(group = type),
               width = 0.35, outlier.shape = NA,
               fill = "grey90", color = "black") +
  # Connecting lines between T0 and T1 per subject
  geom_line(aes(group = IID_ISMMS),
            alpha = 0.25, linewidth = 0.6, color = "gray50") +
  # Points with transparency (easier to see overlap)
  geom_point(aes(group = IID_ISMMS),
             size = 2.8, alpha = 0.5, color = "black") +
  # Add bar and label above
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmin, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmax, xend = xmax, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_text(data = annot_df,
            aes(x = (xmin + xmax)/2, y = y + 0.02, label = label),
            inherit.aes = FALSE, size = 5.5) +
  labs(x = NULL, y = "Spearman correlation") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_t0_t1_matching_points_wcpg.pdf", p, width = 8, height = 5)

###combining T0 and T1 --> combined/within_person group to contrast with crosstimepoint
library(effsize)  # for Cohen's d

# 1. Combine T0 and T1 into one group
df_combined <- cor_long %>%
  filter(type %in% c("T0_only", "T1_only", "Cross_timepoint")) %>%
  mutate(group = ifelse(type %in% c("T0_only", "T1_only"),
                        "Combined",
                        "Cross_timepoint"))

# 2. Keep only individuals with both groups
df_pair <- df_combined %>%
  group_by(IID_ISMMS) %>%
  filter(all(c("Combined", "Cross_timepoint") %in% group)) %>%
  ungroup()

# 3. Paired Wilcoxon test
wilcox_test <- wilcox.test(
  corr ~ group,
  data = df_pair,
  paired = TRUE
)
#         Wilcoxon signed rank test with continuity correction
# data:  corr by group
# V = 5275, p-value = 0.8604

# 4. Compute Cohen’s d (paired)
d_result <- cohen.d(
  df_pair$corr,
  df_pair$group,
  paired = TRUE
)

# Format values for annotation
d_value <- round(d_result$estimate, 2)
p_value <- signif(wilcox_test$p.value, 3)

# 5. Plot with annotation bar
p <- ggplot(df_pair, aes(x = group, y = corr)) +
  geom_boxplot(aes(group = group),
               width = 0.35, fill = "grey90", color = "black", outlier.shape = NA) +
  geom_line(aes(group = IID_ISMMS),
            alpha = 0.25, linewidth = 0.6, color = "gray50") +
  geom_point(size = 2.8, alpha = 0.6, color = "black") +
  labs(x = NULL, y = "Spearman correlation") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# Annotation bar and label
annot_df <- data.frame(
  xmin = 1,
  xmax = 2,
  y = max(df_pair$corr, na.rm = TRUE) + 0.05,
  label = paste0("d = ", d_value, ", p = ", p_value)
)

p <- p +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmin, xend = xmin, y = y, yend = y - 0.02),
               linewidth = 0.8) +
  geom_segment(data = annot_df,
               aes(x = xmax, xend = xmax, y = y, yend = y - 0.02),
               linewidth = 0.8) +
  geom_text(data = annot_df,
            aes(x = (xmin + xmax)/2, y = y + 0.02, label = label),
            size = 5.5)

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_level_combined_cross_time_matching_points_wcpg.pdf", p, width = 8, height = 5)

##############################################################################
########################## NOV 5 2025 ########################################
#### investigating the null a bit more 
/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/

comparison_of_null_and_real_data_summaries_blood_form0.txt
comparison_of_null_and_real_data_summaries_blood_form3_no_residID.txt
comparison_of_null_and_real_data_summaries_blood_form3_with_residID.txt

comparison_of_null_and_real_data_summaries_blood_form5_no_residID.txt
comparison_of_null_and_real_data_summaries_blood_form5_no_residID_full_brain.txt
comparison_of_null_and_real_data_summaries_blood_form5_with_ID.txt
comparison_of_null_and_real_data_summaries_blood_form5_with_ID_full_brain.txt

form3_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_no_residID.txt", data.table=FALSE)
form3_no_residID
#            stat observed_value empirical_p_value
# 11        p99.5     0.20201327             0.063
# 12       p99.95     0.24839444             0.046
# 13      p99.995     0.29069225             0.018
# 14     p99.9995     0.43836985             0.000

form3_with_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form3_with_residID.txt", data.table=FALSE)
form3_with_residID
#            stat observed_value empirical_p_value
# 13      p99.995     0.26650653             0.273
# 14     p99.9995     0.31362213             0.060
# 15    p99.99995     0.41299296             0.000

######### THIS IS WITH BASE BRAIN #########
form5_no_residID_base <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_no_residID.txt", data.table=FALSE)
form5_no_residID_base
#            stat observed_value empirical_p_value
# 9      p99.995     0.27732301             0.098
# 10    p99.9995     0.44724055             0.000

form5_with_residID_base <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_with_ID.txt", data.table=FALSE)
form5_with_residID_base
#            stat observed_value empirical_p_value
# 10    p99.9995     0.30767422             0.097
# 11   p99.99995     0.41535861             0.000

######### THIS IS WITH FULL BRAIN #########
form5_no_residID_full <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_no_residID_full_brain.txt", data.table=FALSE)
form5_no_residID_full
#            stat observed_value empirical_p_value
# 9      p99.995     0.27210186             0.143
# 10    p99.9995     0.31611339             0.020 #this still translates to ~2,240 gene pairs 
# 11   p99.99995     0.59640267             0.000

form5_with_residID_full <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_with_ID.txt", data.table=FALSE)
form5_with_residID_full
#            stat observed_value empirical_p_value
# 9      p99.995     0.26274126             0.450
# 10    p99.9995     0.30767422             0.097
# 11   p99.99995     0.41535861             0.000 #this translates to ~224 gene pairs 

##redo the above with 99.9995th percentile for form5_no_residID_full
#and 99.99995th percentile for form5_with_residID_full

###

#code derived from LBP_blood_brain_QC_20250603.rmd
blood_form3_brain_baseline_spearman_cor_matrix.txt
blood_form3_brain_full_spearman_cor_matrix.txt

blood_form3_no_indivdualID_brain_full_spearman_cor_matrix.txt
blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt
####
blood_form5_brain_baseline_spearman_cor_matrix.txt
blood_form5_brain_full_spearman_cor_matrix.txt

blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt
blood_form5_no_indivdualID_brain_baseline_spearman_cor_matrix.txt

blood_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
blood_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
dim(blood_by_brain_cor) #21046 21356
row.names(blood_by_brain_cor) <- blood_by_brain_cor$V1
blood_by_brain_cor$V1 <- NULL

##Extracting the top 0.05% ~ 99.95th percentile
# Step 1. Compute the 99.95th percentile threshold
threshold <- quantile(as.numeric(as.matrix(blood_by_brain_cor)), probs = 0.9999995, na.rm = TRUE)
threshold 
#0.2184869 for form5_full_with_residID 99.95th
#0.2231395 for form5_full_no_residID 99.95th 

# Step 1. Compute the 99.995th percentile threshold for no_residID
threshold <- quantile(as.numeric(as.matrix(blood_by_brain_cor)), probs = 0.999995, na.rm = TRUE)
threshold 
#0.3069401 for form5_full_no_residID 99.9995th 

# Step 1. Compute the 99.9995th percentile threshold for with_residID
threshold <- quantile(as.numeric(as.matrix(blood_by_brain_cor)), probs = 0.9999995, na.rm = TRUE)
threshold 
#0.3333631 for form5_full_with_residID 99.99995th

# Step 2. Find the highly correlated pairs
mat <- as.matrix(blood_by_brain_cor)
blood_genes <- rownames(mat)
brain_genes <- colnames(mat)

# Get indices of correlations above threshold
high_idx <- which(mat >= threshold, arr.ind = TRUE)

# Step 3. Create a table of those pairs
top_pairs <- data.frame(
  blood_gene = blood_genes[high_idx[, 1]],
  brain_gene = brain_genes[high_idx[, 2]],
  correlation = mat[high_idx]
)

################### 99.95th PERCENTILE START #####################
############# form5_full_with_residID #############
####################################################
# Step 4. Sort descending by correlation - ALL GENES
top_pairs <- top_pairs[order(-top_pairs$correlation), ]
dim(top_pairs) #224730      3
head(top_pairs)
#                blood_gene         brain_gene correlation
# 15907   ENSG00000168209.5 ENSG00000090104.12   0.4450948
# 138690 ENSG00000111788.10  ENSG00000214826.5   0.4417699
# 103295  ENSG00000168209.5 ENSG00000178878.12   0.4352866
# 160508  ENSG00000238083.7  ENSG00000238083.7   0.4333312
# 101171  ENSG00000238083.7 ENSG00000176681.14   0.4213095

head(sort(table(top_pairs$blood_gene), decreasing = TRUE),10)
# ENSG00000135245.10 (lipid storage)  ENSG00000110887.8  ENSG00000279024.1  ENSG00000261251.1 
#                541                            464                448                331 
# ENSG00000082126.18  ENSG00000240859.2 ENSG00000102977.15  ENSG00000263089.1 
#                330                287                281                247 
#  ENSG00000231992.1 ENSG00000161395.14 
#                232                231 

#SAME GENE ONLY 
same_gene_pairs <- top_pairs[top_pairs$blood_gene == top_pairs$brain_gene, ]
dim(same_gene_pairs) #296   3
head(same_gene_pairs,10)
#                blood_gene         brain_gene correlation
# 160508  ENSG00000238083.7  ENSG00000238083.7   0.4333312
# 101165 ENSG00000176681.14 ENSG00000176681.14   0.4152370
# 222396  ENSG00000285534.1  ENSG00000285534.1   0.4093984
# 167476  ENSG00000247498.9  ENSG00000247498.9   0.3979235
# 142112  ENSG00000223496.3  ENSG00000223496.3   0.3965834

############# form5_full_no_residID #############
####################################################
top_pairs <- top_pairs[order(-top_pairs$correlation), ]
dim(top_pairs) #224738      3
head(top_pairs)
#                blood_gene         brain_gene correlation
# 153188 ENSG00000226259.10 ENSG00000226259.10   0.8780499
# 153969  ENSG00000226752.9  ENSG00000226752.9   0.8543236
# 65926  ENSG00000145736.14 ENSG00000145736.14   0.8382796

head(sort(table(top_pairs$blood_gene), decreasing = TRUE),10)
 ENSG00000228903.7 ENSG00000160439.16  ENSG00000105767.3  ENSG00000110887.8 
               824                601                519                335 
ENSG00000182010.11 ENSG00000157554.19 ENSG00000119777.20  ENSG00000279672.1 
               322                314                274                270 
ENSG00000079393.20 ENSG00000135845.10 
               252                249

#SAME GENE ONLY 
same_gene_pairs <- top_pairs[top_pairs$blood_gene == top_pairs$brain_gene, ]
dim(same_gene_pairs) #784   3
head(same_gene_pairs,10)
               blood_gene         brain_gene correlation
153188 ENSG00000226259.10 ENSG00000226259.10   0.8780499
153969  ENSG00000226752.9  ENSG00000226752.9   0.8543236
65926  ENSG00000145736.14 ENSG00000145736.14   0.8382796
168204  ENSG00000241945.8  ENSG00000241945.8   0.8358639
164887  ENSG00000237541.3  ENSG00000237541.3   0.8304804
206855  ENSG00000274602.5  ENSG00000274602.5   0.8253097
193680  ENSG00000265218.1  ENSG00000265218.1   0.8156110
208658  ENSG00000275895.7  ENSG00000275895.7   0.7910261
86402  ENSG00000164308.16 ENSG00000164308.16   0.7896608
163337  ENSG00000235833.1  ENSG00000235833.1   0.7893700

###### LOOK AT EXCEL FOR top_correlated_same_genes
################### 99.95th PERCENTILE ENDS #####################

################### 99.9995th PERCENTILE START ##################### 
#top 0.0005%
############# form5_full_no_residID #############
#################################################
top_pairs <- top_pairs[order(-top_pairs$correlation), ]
dim(top_pairs) #2248    3
head(sort(table(top_pairs$blood_gene), decreasing = TRUE),10)
#  ENSG00000228903.7 ENSG00000112576.12 ENSG00000160439.16  ENSG00000168209.5 
#                 37                 14                 12                 12 
#  ENSG00000203663.4  ENSG00000228049.7  ENSG00000267645.5 ENSG00000084070.12 
#                 12                 11                 11                 10 
# ENSG00000105808.17 ENSG00000168255.20 
#                 10                 10

#SAME GENE ONLY 
same_gene_pairs <- top_pairs[top_pairs$blood_gene == top_pairs$brain_gene, ]
dim(same_gene_pairs) #441   3 why does higher percentile yield more same genes hmmmmm

################### 99.99995th PERCENTILE START ##################### 
#top 0.00005%
############# form5_full_with_residID #############
###################################################
top_pairs <- top_pairs[order(-top_pairs$correlation), ]
dim(top_pairs) #2248    3
head(sort(table(top_pairs$blood_gene), decreasing = TRUE),10)
#  ENSG00000119138.4  ENSG00000166523.7  ENSG00000168209.5 ENSG00000096060.14 
#                  6                  6                  6                  5 
# ENSG00000109466.14 ENSG00000115590.14 ENSG00000112576.12  ENSG00000120129.6 
#                  4                  4                  3                  3 
# ENSG00000122025.14 ENSG00000123836.15 
#                  3                  3
#SAME GENE ONLY 
same_gene_pairs <- top_pairs[top_pairs$blood_gene == top_pairs$brain_gene, ]
dim(same_gene_pairs) #22  3

############ PERFORMING GO ANALYSIS FOR top_pairs | https://www.youtube.com/watch?v=JPwdqdo_tRg | https://geneontology.org/docs/go-enrichment-analysis/
#steps (if use online): 
#1. Paste gene names
#2. Select GO aspect (molecular function, biological process, cellular component)
#3. Select species
#4. Upload a reference (aka background) list at a later step 

############# form5_full_no_residID #############
top_pairs <- top_pairs[order(-top_pairs$correlation), ]
dim(top_pairs) #2248    3

library(clusterProfiler)
library(BiocManager)
library(org.Hs.eg.db) #homo sapiens
library(AnnotationDbi)

genes_to_test <- top_pairs$blood_gene
genes_to_test <- sub("\\..*", "", genes_to_test)

#The reference universe = all human genes that have a GO annotation in org.Hs.eg.db.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL", pvalueCutoff  = 0.10, pAdjustMethod = "BH")
GO_results_df <- as.data.frame(GO_results); GO_results_df
##nothing is enriched hmmmmm
##loosen the threshold, 1 significant result
           ONTOLOGY         ID          Description GeneRatio   BgRatio
GO:0045171       CC GO:0045171 intercellular bridge    13/868 102/22283
                 pvalue   p.adjust     qvalue
GO:0045171 0.0001653021 0.09356101 0.09356101
                                                                                                                                                                                                                    geneID
GO:0045171 ENSG00000134184/ENSG00000213366/ENSG00000134202/ENSG00000168765/ENSG00000173226/ENSG00000024422/ENSG00000100055/ENSG00000177483/ENSG00000243156/ENSG00000134201/ENSG00000040199/ENSG00000135596/ENSG00000062650
           Count
GO:0045171    13

barplot(GO_results, showCategory = 20)
dotplot(GO_results, showCategory = 20)
emapplot(GO_results)


# 2) See how many map to GO terms
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = genes_to_test,
                             keytype = "ENSEMBL",
                             columns = c("ENSEMBL","SYMBOL","GO"),
                             pvalueCutoff  = 0.10, 
                             pAdjustMethod = "BH")
n_input <- length(genes_to_test)
n_have_go <- length(unique(map$ENSEMBL[!is.na(map$GO)]))
cat("Input genes:", n_input, " | with GO terms:", n_have_go, "\n")
#Input genes: 2248  | with GO terms: 906 

# 3) Keep only those with GO annotation for enrichment input
genes_with_go <- unique(map$ENSEMBL[!is.na(map$GO)])

go_all <- enrichGO(
  gene          = genes_with_go,
  universe      = universe_all,           # or your ALL_GENES_TESTED
  OrgDb         = org.Hs.eg.db,           # <-- object, not string
  keyType       = "ENSEMBL",
  ont           = "ALL",                  # try "BP" first if still empty
  pvalueCutoff  = 0.10,                   # start a bit looser
  qvalueCutoff  = 0.20,
  pAdjustMethod = "BH",
  readable      = TRUE                    # converts to SYMBOL for readability
)
as.data.frame(go_all)[, c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")]
