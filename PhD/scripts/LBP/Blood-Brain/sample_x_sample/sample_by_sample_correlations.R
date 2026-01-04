##sample-level correlations
library(data.table)
library(dplyr)
library(ggplot2)
library(effsize)
library(ggsignif)

# Load data
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
table(brain_metadata$number_of_pair_brain) ##
head(dictionary)

# blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
# rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL
# blood_form <- blood_form5

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL
blood_form <- blood_form5_no_residID

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full)<- brain_full$V1
brain_full$V1 <- NULL
dim(brain_full) #21356   225

######################## only getting T0 for those individuals with repeated measures
dictionary$timepoint[dictionary$number_of_pair_brain == "1_pair"] <- 0
table(dictionary$timepoint)
#   0   1 
# 152  73
#this is correct

# Get common genes
blood_genes <- rownames(blood_form)
brain_genes <- rownames(brain_full)
common_genes <- intersect(blood_genes, brain_genes);
length(common_genes) #17533

# Subset both datasets to common genes
blood_matched <- blood_form[common_genes, ]
brain_matched <- brain_full[common_genes, ]

# Subset dictionary for timepoint == 0
dict_t0 <- subset(dictionary, timepoint == 0)

# Match sample IDs for blood and brain
blood_samples_t0 <- dict_t0$SAMPLE_ISMMS_blood
brain_samples_t0 <- dict_t0$SAMPLE_ISMMS_brain

# Subset matrices
blood_t0 <- blood_matched[, colnames(blood_matched) %in% blood_samples_t0]
brain_t0 <- brain_matched[, colnames(brain_matched) %in% brain_samples_t0]

dim(blood_t0)
#17533   152

dim(brain_t0)
#17533   152

# Compute Spearman correlation between blood samples and brain samples
sample_cor_matrix <- cor(blood_t0, brain_t0, method = "spearman")
#Rows = blood samples
#Columns = brain samples
#Each cell = correlation between one blood sample and one brain sample across all 17,533 genes

dim(sample_cor_matrix) #152 152
#do this again for all samples 

summary(as.vector(sample_cor_matrix))
#with residID
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.6594448 -0.0867099  0.0011819  0.0009162  0.0896992  0.6140214 

#no residID
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.6222952 -0.0830192  0.0006048  0.0011929  0.0854818  0.5922921

####### violin/boxplot of within vs between

# Extract within-person correlations (diagonal)
within_person <- diag(sample_cor_matrix)
summary(within_person)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.35174 -0.03666  0.06433  0.06275  0.16211  0.36916 

# Extract between-person correlations (all off-diagonal entries)
between_person <- sample_cor_matrix
diag(between_person) <- NA
between_person <- as.vector(between_person)
between_person <- between_person[!is.na(between_person)]

# Combine into one dataframe
df_cor <- data.frame(
  correlation = c(within_person, between_person),
  type = c(rep("Within-person", length(within_person)),
           rep("Between-person", length(between_person)))
)

# Make group order explicit so Cohen's d is "Within - Between"
df_cor$type <- factor(df_cor$type, levels = c("Between-person", "Within-person"))

# Tests
t_res <- t.test(correlation ~ type, data = df_cor); t_res     
# t = -3.347, df = 152.8, p-value = 0.001028          
w_res <- wilcox.test(correlation ~ type, data = df_cor, exact = FALSE); w_res  
# W = 1478010, p-value = 0.001155
# Effect size: Cohen's d

d_res <- effsize::cohen.d(correlation ~ type, data = df_cor, hedges.correction = FALSE)

# Combine label
p_lab <- sprintf("p = %.2g, d = %.2f", w_res$p.value, unname(d_res$estimate))

y_top <- max(df_cor$correlation, na.rm = TRUE)

p <- ggplot(df_cor, aes(x = type, y = correlation, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
  labs(title = "Blood Brain Correlations",
       x = "", y = "Spearman correlation") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Within-person","Between-person")),
              annotations = p_lab,
              y_position = y_top + 0.05, tip_length = 0.01, textsize = 4) +
  expand_limits(y = y_top + 0.10)

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_sample_within_between_people_violin_boxplot_no_residID.pdf", p, width = 8, height = 6)

### WCPG PLOT 
p <- ggplot(df_cor, aes(x = type, y = correlation, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.size = 0.8, alpha = 0.8, color = "black") +
  labs(
    x = NULL,
    y = "Spearman correlation"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("Within-person", "Between-person")),
    annotations = p_lab,
    y_position = y_top + 0.05,
    tip_length = 0.01,
    textsize = 6,
    vjust = 0.2
  ) +
  expand_limits(y = y_top + 0.10)

ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_sample_within_between_people_violin_boxplot_no_residID_wcpg.pdf", p, width = 8, height = 5)


#ggsave("/hpc/users/hoangd02/www/plots/lbp/sample_sample_within_between_people_violin_boxplot.pdf", p, width = 8, height = 6)

# Heatmap of the full 152 x 152 matrix
library(ComplexHeatmap)
library(circlize)
library(grid)

# 1) Define the mapping function with fixed limits
col_fun <- circlize::colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"))

# 2) (Optional but helpful) choose nice legend ticks
legend_breaks <- c(-0.7, -0.4, 0, 0.4, 0.7)

pdf("/hpc/users/hoangd02/www/plots/sample_sample_correlation_heatmap.pdf",
    width = 10, height = 8)

Heatmap(
  sample_cor_matrix,
  name = "rho",
  col = col_fun,                      # <- mapping with fixed limits
  heatmap_legend_param = list(        # <- force legend ticks
    at = legend_breaks,
    labels = legend_breaks
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Blood samples",
  column_title = "Brain samples",
  cell_fun = function(j, i, x, y, w, h, col) {
    if (i == j) grid.rect(x, y, w, h, gp = gpar(col = "black", lwd = 2, fill = NA))
  }
)
dev.off()

############################################### exploring within-person cor
head(df_cor)
#   correlation          type
# 1  0.09975948 Within-person
# 2  0.17753841 Within-person
# 3 -0.07189302 Within-person

within_person <- filter(df_cor, type == "Within-person")
summary(within_person$correlation)
##with resid ID
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.36709 -0.05088  0.03477  0.04106  0.14008  0.38457

##no resid ID
summary(within_person$correlation)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.35174 -0.03666  0.06433  0.06275  0.16211  0.36916 

# Put within-person correlations into a data frame
within_person$correlation <- as.numeric(within_person$correlation)

# Define bins
breaks <- c(-Inf, -0.3, -0.2, -0.1, 0,
             0.1, 0.2, 0.3, Inf)

labels <- c("< -0.3", "-0.3 to -0.2", "-0.2 to -0.1", "-0.1 to 0",
            "0 to <0.1", "0.1 to <0.2", "0.2 to <0.3", "0.3 and up")

# Assign bins
within_person$bin <- cut(within_person$correlation,
                         breaks = breaks, labels = labels, right = FALSE)

# Count per bin
bin_counts <- within_person %>%
  count(bin)

# Plot
plot <- ggplot(bin_counts, aes(x = bin, y = n, fill = bin)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "Within person Blood Brain Correlations",
       x = "Correlation bin", y = "Count") +
  theme(legend.position = "none")

ggsave("/hpc/users/hoangd02/www/plots/lbp/within_person_blood_brain_cor.pdf", plot, width = 8, height = 6)

ggsave("/hpc/users/hoangd02/www/plots/lbp/within_person_blood_brain_cor_no_residID.pdf", plot, width = 8, height = 6)

df<- read.csv("/sc/arion/projects/mscic1/results/Ariela/Ariela_data/Cognition/cognitiveDatExploration.csv")

dim(df)
# 9119   19
table(df$Tests)
table(df$tests_abbr) #this seems to be the broader umbrella of df$Tests
#there are some weird labeling for BAI and Beck Anxiety Inventory in df$Tests/df$tests_abbr

##lets just start with the smaller subset for now, can always go back and change --- 
# Get counts of each Tests value
test_counts <- table(df$Tests)

# Keep only those with count > 50
tests_over100 <- names(test_counts[test_counts > 100])
tests_over100
# [1] "Beck Anxiety Inventory"             "Becks Depression Inventory"        
# [3] "Boston Naming Test"                 "F-A-S phenomic verbal fluency test"

sub <- df[df$Tests %in% tests_over100, ]
table(sub$Tests)
            # Beck Anxiety Inventory         Becks Depression Inventory 
            #                    150                                114 
            #     Boston Naming Test F-A-S phenomic verbal fluency test 
            #                    220                                211

table(sub$tests_abbr)
# BAI BDI BNT FAS 
# 150 114 220 211 

table(sub$PT.ID)
PT-0036
PT-0045

test <- sub[sub$PT.ID == "PT-0036", ]; test #BNT
test <- sub[sub$PT.ID == "PT-0045", ]; test

library(dplyr)

# Prepare BNT rows from `sub`
bnt_df <- sub %>%
  filter(`tests_abbr` == "BNT") %>%
  mutate(
    `Raw.Score`  = suppressWarnings(as.numeric(`Raw.Score`)),
    Percentile   = suppressWarnings(as.numeric(gsub("[^0-9.]", "", Percentile))) # strip '%' etc.
  ) %>%
  filter(!is.na(`PT.ID`)) %>%
  group_by(`PT.ID`) %>%
  slice_max(order_by = `Raw.Score`, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(`PT.ID`,
         BNT_Raw.Score    = `Raw.Score`,
         BNT_Percentile   = Percentile)

# Join into dictionary
dictionary <- dictionary %>%
  left_join(bnt_df, by = c("IID_ISMMS" = "PT.ID"))

sum(is.na(dictionary$BNT_Raw.Score))
# [1] 76
sum(is.na(dictionary$BNT_Percentile))
# [1] 74


############################################### OLD ANALYSES
###############################################
################# ANALYSIS 1 ##################
###############################################
## plot correlations the diagonal of sample_cor_matrix (sample A-A) against 
#correlations of the non-paired sample in sample_cor_matrix (sample A- sample B, etc.)

# Pull paired correlations
paired_corrs <- mapply(function(blood, brain) {
  if (blood %in% rownames(sample_cor_matrix) &&
      brain %in% colnames(sample_cor_matrix)) {
    return(sample_cor_matrix[blood, brain])
  } else {
    return(NA_real_)
  }
}, dictionary$SAMPLE_ISMMS_blood, dictionary$SAMPLE_ISMMS_brain)

paired_df <- data.frame(
  type = "paired",
  corr = paired_corrs,
  IID_ISMMS = dictionary$IID_ISMMS
)
paired_df <- paired_df[!is.na(paired_df$corr), ]
dim(paired_df) #225   3

summary(paired_df$corr) ##base brain with blood form3_noID 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.267913 -0.092218 -0.025336 -0.008315  0.070044  0.306947 

summary(paired_df$corr) ##full brain with blood form3_noID 
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.36810 -0.02362  0.04656  0.05493  0.13863  0.46936 

# All correlations as a long dataframe
library(reshape2)
all_corrs <- melt(sample_cor_matrix, varnames = c("blood", "brain"), value.name = "corr")
dim(all_corrs) #50625     4
#225^2 = 50625

# Flag paired vs non-paired
all_corrs$type <- "non_paired"
for (i in seq_len(nrow(dictionary))) {
  b <- dictionary$SAMPLE_ISMMS_blood[i]
  br <- dictionary$SAMPLE_ISMMS_brain[i]
  all_corrs$type[all_corrs$blood == b & all_corrs$brain == br] <- "paired"
}

nonpaired_df <- filter(all_corrs, type == "non_paired")
dim(nonpaired_df) #50400     4
##225^2 -225 paired = 50400
summary(nonpaired_df$corr) ##base brain with blood form3_noID 
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.41389 -0.09920 -0.03983 -0.01752  0.06270  0.42958 

summary(nonpaired_df$corr) ##full brain with blood form3_noID 
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.635963 -0.077566  0.001123  0.001124  0.080332  0.597634

plot <- ggplot(all_corrs, aes(x = type, y = corr, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Blood Sample x Brain Sample Gene Expression Correlations Across Genes",
    x = "Type",
    y = "Spearman correlation"
  )
ggsave("/hpc/users/hoangd02/www/plots/sample_sample_correlation_box_violin_plot_full_brain.pdf", plot, width = 10, height = 6)

ggsave("/hpc/users/hoangd02/www/plots/sample_sample_correlation_box_violin_plot.pdf", plot, width = 10, height = 6)

###############################################
################# SPECIAL ANALYSIS ############ within cross-time point
###############################################
library(dplyr)
library(reshape2)
library(ggplot2)

## 0) Long form of all correlations (every blood × brain)
all_corrs <- melt(sample_cor_matrix, varnames = c("blood", "brain"), value.name = "corr")

## 1) Clean dictionary to just the pairs that exist in the matrix
dict_clean <- dictionary %>%
  distinct(IID_ISMMS, SAMPLE_ISMMS_blood, SAMPLE_ISMMS_brain, timepoint, number_of_pair_brain) %>%
  filter(SAMPLE_ISMMS_blood %in% rownames(sample_cor_matrix),
         SAMPLE_ISMMS_brain %in% colnames(sample_cor_matrix))

## 2) WITHIN: same timepoint (the exact A–A pairs from the dictionary)
within_same <- dict_clean %>%
  transmute(
    blood = SAMPLE_ISMMS_blood,
    brain = SAMPLE_ISMMS_brain,
    type  = "within_same_timepoint",
    IID_ISMMS
  )

## 3) WITHIN: cross timepoint (only for people with 2 pairs)
## Build all cross combinations within the same person: blood from row A × brain from row B, where timepoint differs
within_cross <- dict_clean %>%
  inner_join(dict_clean, by = "IID_ISMMS", suffix = c("_b", "_br")) %>%
  filter(timepoint_b != timepoint_br) %>%
  transmute(
    blood = SAMPLE_ISMMS_blood_b,
    brain = SAMPLE_ISMMS_brain_br,
    type  = "within_cross_timepoint",
    IID_ISMMS
  ) %>%
  distinct()  # avoid duplicates (each cross pairing can appear twice)

## 4) Combine labels and apply to the full table
labels_df <- bind_rows(within_same, within_cross) %>%
  select(blood, brain, type)

## Add a type column to all_corrs (default = between_people), then overwrite when matched
comparison_data <- all_corrs %>%
  left_join(labels_df, by = c("blood", "brain")) %>%
  mutate(type = ifelse(is.na(type), "between_people", type))

table(comparison_data$type)
        # between_people within_cross_timepoint  within_same_timepoint 
        #          50254                    146                    225
# Expect:
# - within_same_timepoint ≈ nrow(dict_clean)
# - within_cross_timepoint ≈ 2 * (# of IIDs with exactly two rows/timepoints)
# - between_people = total - (within_same + within_cross)

## 6) Plot the three distributions
comparison_data$type <- factor(
  comparison_data$type,
  levels = c("within_same_timepoint", "within_cross_timepoint", "between_people"),
  labels = c("Within (same timepoint)", "Within (cross timepoint)", "Between people")
)

p1<- ggplot(comparison_data, aes(x = type, y = corr, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sample x Sample Blood Brain Correlations Across Genes",
    x = NULL,
    y = "Spearman correlation"
  ) +
  guides(fill = "none")

ggsave("/hpc/users/hoangd02/www/plots/test10.pdf", p1, width = 10, height = 6)

## 7) Significance Testing
comparison_data %>%
  group_by(type) %>%
  summarise(n = n(), mean = mean(corr), median = median(corr),
            q1 = quantile(corr, 0.25), q3 = quantile(corr, 0.75))

# Wilcoxon tests (non-parametric) for differences
wilcox.test(corr ~ type, data = subset(comparison_data, type %in% c("Within (same timepoint)", "Between people")))
# 	Wilcoxon rank sum test with continuity correction

# data:  corr by type
# W = 6951775, p-value = 2.642e-09
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(corr ~ type, data = subset(comparison_data, type %in% c("Within (cross timepoint)", "Between people")))
# 	Wilcoxon rank sum test with continuity correction

# data:  corr by type
# W = 4669774, p-value = 1.174e-08
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(corr ~ type, data = subset(comparison_data, type %in% c("Within (same timepoint)", "Within (cross timepoint)")))
# 	Wilcoxon rank sum test with continuity correction

# data:  corr by type
# W = 15861, p-value = 0.5766
# alternative hypothesis: true location shift is not equal to 0


###############################################
################# ANALYSIS 2 ##################
###############################################
#comparing:
#Within-person correlations: Correlations between samples from the same individual (including cases where one person has 2 blood-brain pairs)
#Between-person correlations: Correlations between samples from different individuals
# Extract only TRUE paired correlations (same row in dictionary)
true_paired_corrs <- mapply(function(blood, brain) {
  if (blood %in% rownames(sample_cor_matrix) &&
      brain %in% colnames(sample_cor_matrix)) {
    return(sample_cor_matrix[blood, brain])
  } else {
    return(NA_real_)
  }
}, dictionary$SAMPLE_ISMMS_blood, dictionary$SAMPLE_ISMMS_brain)

true_paired_df <- data.frame(
  blood = dictionary$SAMPLE_ISMMS_blood,
  brain = dictionary$SAMPLE_ISMMS_brain,
  person_id = dictionary$IID_ISMMS,
  corr = true_paired_corrs,
  type = "within_person"
)
true_paired_df <- true_paired_df[!is.na(true_paired_df$corr), ]

# Mark these true pairs in all_corrs
all_corrs$type <- "between_people"
for(i in 1:nrow(true_paired_df)) {
  all_corrs$type[all_corrs$blood == true_paired_df$blood[i] & 
                 all_corrs$brain == true_paired_df$brain[i]] <- "within_person"
}

# Filter for comparison
comparison_data <- all_corrs[all_corrs$type %in% c("within_person", "between_people"), ]

# Density plot
p1 <- ggplot(comparison_data, aes(x = corr, fill = type)) +
  geom_density(alpha = 0.7) +
  labs(title = "Distribution of Blood-Brain Correlations",
       subtitle = "Within Person (True Pairs) vs Between People",
       x = "Correlation",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("within_person" = "lightblue", "between_people" = "red"),
                    labels = c("Between People", "Within Person")) +
  theme(legend.title = element_blank())
ggsave("/hpc/users/hoangd02/www/plots/test3.pdf", p1, width = 10, height = 6)
#ggsave("/hpc/users/hoangd02/www/plots/test.pdf", p1, width = 10, height = 6)

# Side-by-side comparison
p2 <- ggplot(comparison_data, aes(x = corr, fill = type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~type, scales = "free_y") +
  labs(title = "Distribution of Blood-Brain Correlations",
       x = "Correlation",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("within_person" = "lightblue", "between_people" = "red")) +
  theme(legend.position = "none")
ggsave("/hpc/users/hoangd02/www/plots/test2.pdf", p2, width = 10, height = 6)

# Statistical comparison
t.test(corr ~ type, data = comparison_data)
# 	Welch Two Sample t-test
# data:  corr by type
# t = -1.1527, df = 226.03, p-value = 0.2503
# alternative hypothesis: true difference in means between group between_people and group within_person is not equal to 0
# 95 percent confidence interval:
#  -0.02493653  0.00652999
# sample estimates:
# mean in group between_people  mean in group within_person 
#                 -0.017518041                 -0.008314771 

###############################################
################# ANALYSIS 3 ##################
###############################################
# Get individuals with 2 paired brain samples
two_pair_inds <- unique(brain_metadata$IID_ISMMS[brain_metadata$number_of_pair_brain == "2_pairs"])
n_two_pair <- length(two_pair_inds);  n_two_pair #73

# Extract correlations for matched pairs (same individual)
matched_correlations <- c()
sample_ids_blood <- colnames(blood_form)
sample_ids_brain <- colnames(v_brain$E)

# Find matched pairs based on individual ID
for (id in two_pair_inds) {
  # Get blood sample for this individual
  blood_cols <- which(grepl(id, sample_ids_blood))
  # Get brain samples for this individual  
  brain_cols <- which(brain_metadata$IID_ISMMS == id)
  
  if (length(blood_cols) == 1 & length(brain_cols) == 2) { ##this might be wrong ~~~
    # Get correlations between this person's blood and both brain samples
    cors <- sample_cor_matrix[blood_cols, brain_cols]
    matched_correlations <- c(matched_correlations, cors)
  }
}

# Extract off-diagonal correlations (unmatched pairs)
# Sample a subset to avoid memory issues
set.seed(123)
n_sample <- min(10000, length(sample_cor_matrix) - nrow(sample_cor_matrix))
off_diag_indices <- sample(which(upper.tri(sample_cor_matrix, diag = FALSE)), n_sample)
unmatched_correlations <- sample_cor_matrix[off_diag_indices]

# Create comparison dataframe for plotting
correlation_data <- data.frame(
  correlation = c(matched_correlations, unmatched_correlations),
  type = c(rep("Matched (same individual)", length(matched_correlations)),
           rep("Unmatched (different individuals)", length(unmatched_correlations)))
)

# Create plots
# 1. Density plot
density_plot <- ggplot(correlation_data, aes(x = correlation, fill = type, alpha = 0.7)) +
  geom_density() +
  labs(
    title = "Distribution of Sample x Sample Correlations",
    subtitle = "Blood vs Brain Expression Profiles",
    x = "Spearman Correlation",
    y = "Density",
    fill = "Comparison Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_alpha_identity()

# 2. Box plot
box_plot <- ggplot(correlation_data, aes(x = type, y = correlation, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(alpha = 0.3, width = 0.2, size = 0.5) +
  labs(
    title = "Sample x Sample Correlations: Matched vs Unmatched",
    subtitle = "Blood vs Brain Expression Profiles",
    x = "",
    y = "Spearman Correlation"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots
ggsave("/hpc/users/hoangd02/www/plots/sample_correlation_density.pdf", density_plot, width = 10, height = 6)
ggsave("/hpc/users/hoangd02/www/plots/sample_correlation_boxplot.pdf", box_plot, width = 8, height = 6)

# Save results
save(sample_cor_matrix, matched_correlations, unmatched_correlations, 
     file = "sample_correlation_results.RData")

