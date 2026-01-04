# https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html

library(data.table)
library(dplyr)
library(ggplot2)
library(variancePartition)
library(tidyr)
library(stringr)
library(tibble)

# Load data
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

head(v_brain$E)
head(v_blood$E)

##all forms
form3 <- ~ (1|IID_ISMMS) + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood 
form5 <- ~ (1|IID_ISMMS) + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood + RNASeqMetrics_PCT_INTRONIC_BASES_blood

form3 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood 
form5 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood

###form3 with no IDs
form3 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood 
varPart3 <- fitExtractVarPartModel(v_blood$E, form3, blood_metadata)
vp3 <- sortCols(varPart3)

pdf("/hpc/users/hoangd02/www/plots/lbp_blood_form3_noresidID_vpa.pdf", width = 15, height = 8)
plotVarPart(vp3)
dev.off()
https://hoangd02.u.hpc.mssm.edu/plots/lbp_blood_form3_noresidID_vpa.pdf

###form5 with no IDs
form5 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood
varPart5 <- fitExtractVarPartModel(v_blood$E, form5, blood_metadata)
vp5 <- sortCols(varPart5)

pdf("/hpc/users/hoangd02/www/plots/lbp_blood_form5_noresidID_vpa.pdf", width = 20, height = 8)
plotVarPart(vp5)
dev.off()

###form3 with IDs
form3 <- ~ (1|IID_ISMMS) + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood 
varPart3 <- fitExtractVarPartModel(v_blood$E, form3, blood_metadata)
vp3 <- sortCols(varPart3)

pdf("/hpc/users/hoangd02/www/plots/lbp_blood_form3_vpa.pdf", width = 15, height = 8)
plotVarPart(vp3)
dev.off()

####form5 with IDs 
form5 <- ~ (1|IID_ISMMS) + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood + RNASeqMetrics_PCT_INTRONIC_BASES_blood
varPart5 <- fitExtractVarPartModel(v_blood$E, form5, blood_metadata)
vp5 <- sortCols(varPart5)

pdf("/hpc/users/hoangd02/www/plots/lbp_blood_form5_vpa.pdf", width = 15, height = 8)
plotVarPart(vp5)
dev.off()

brain_full_form <- ~ (1|mymet_depletionbatch_brain) + mymet_rin_brain + (1|mymet_bank_brain) + (1|mymet_sex_brain) + RNASeqMetrics_INTRONIC_BASES_brain + ODC_brain + (1|mymet_tissue_brain)
#brain_full_form <- ~ (1|IID_ISMMS) + (1|mymet_depletionbatch_brain) + mymet_rin_brain + (1|mymet_bank_brain) + (1|mymet_sex_brain) + RNASeqMetrics_INTRONIC_BASES_brain + ODC_brain + (1|mymet_tissue_brain)
varPart <- fitExtractVarPartModel(v_brain$E, brain_full_form, brain_metadata)
vp <- sortCols(varPart)

#pdf("/hpc/users/hoangd02/www/plots/lbp_brain_form_full_vpa.pdf", width = 25, height = 8)
pdf("/hpc/users/hoangd02/www/plots/lbp_brain_form_full_no_residID_vpa.pdf", width = 25, height = 8)
plotVarPart(vp)
dev.off()


# Find all metadata variables ending in "_brain"
brain_vars <- grep("_brain$", colnames(brain_metadata), value = TRUE)

# Build formula string
form_str <- paste("~ (1|IID_ISMMS) +", paste(brain_vars, collapse = " + "))

# Convert to formula
form_all_brain <- as.formula(form_str)
form_all_brain

#filter out variables with zero variance 
# 1) Identify all *_brain columns
brain_vars <- grep("_brain$", colnames(brain_metadata), value = TRUE)

# 2) Standardize types (chars -> factors; leave numerics alone)
is_char <- vapply(brain_metadata[brain_vars], is.character, logical(1))
brain_metadata[brain_vars[is_char]] <- lapply(brain_metadata[brain_vars[is_char]], factor)

# 3) Find zero-variance variables (ignoring NA)
zero_var <- vapply(brain_metadata[brain_vars], function(x) {
  x_no_na <- x[!is.na(x)]
  if (length(x_no_na) == 0L) return(TRUE)        # all NA -> treat as zero variance
  if (is.numeric(x_no_na)) return(stats::var(x_no_na) == 0)
  if (is.factor(x_no_na) || is.logical(x_no_na)) return(nlevels(factor(x_no_na)) < 2)
  # fallback for other types
  length(unique(x_no_na)) < 2
}, logical(1))

if (any(zero_var)) {
  message("Dropping zero-variance variables: ",
          paste(names(zero_var)[zero_var], collapse = ", "))
}

brain_vars_clean <- setdiff(brain_vars, names(zero_var)[zero_var])

form_all_brain_clean <- as.formula(
  paste("~", paste(brain_vars_clean, collapse = " + "))
)

C <- canCorPairs(form_all_brain_clean, brain_metadata)
pdf("/hpc/users/hoangd02/www/plots/corr_matrix.pdf", width = 20, height = 20)
plotCorrMatrix(C)
dev.off()


# formula_test=formula(paste("~ (1|ceradsc) + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))
# formula_test0=formula(paste("~ 0 + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

# design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
# design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

# #Applying Voom Transformation 
# isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
# dge <- DGEList(counts=genedata[isexpr,]) 
# dge <- calcNormFactors(dge)

# #pdf("/sc/arion/projects/mscic1/results/jolie/my_voom_plot.pdf")
# v <- voom(dge, design, plot=TRUE)
# #dev.off()

# vp.ROSMAP <- fitExtractVarPartModel(v$E, formula_test, info_alltmp)
# percent_status_mean <- mean(vp.ROSMAP$ceradsc) * 100 #0.2282867%

## SEPT 4 2025 -- MAYBE PUT THIS IN A SEPARATE SCRIPT ! 
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

#form1 <- ~ (1|IID_ISMMS)
#form2 <- ~ (1|IID_ISMMS) + mymet_rin_blood
form <- ~1 
fit_blood <- dream(v_blood$E, form, blood_metadata)
resid_expr_blood <- residuals(fit_blood)

pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance

# form <- 
# fit_blood <- dream(v_blood$E, form, blood_metadata)
# resid_expr_blood <- residuals(fit_blood)

# pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
# pca_summary <- summary(pca_resid_blood)$importance

#####base blood vs base brain 
cor_matrix <- abs(cor(t(resid_expr_blood), t(v_brain$E), method = "spearman")) 
dim(cor_matrix) #21046 21356
cor_matrix_mean <- mean(as.matrix(cor_matrix)); cor_matrix_mean #0.05346138
cor_matrix_median <- median(as.matrix(cor_matrix)); cor_matrix_median #0.04457122
cor_matrix_90 <- quantile(as.matrix(cor_matrix), probs = 0.90); cor_matrix_90 #0.1112432
cor_matrix_95 <- quantile(as.matrix(cor_matrix), probs = 0.95); cor_matrix_95 #0.1336336
cor_matrix_99.95 <- quantile(as.matrix(cor_matrix), probs = 0.9995); cor_matrix_99.95 #0.2418879
cor_matrix_99.995 <- quantile(as.matrix(cor_matrix), probs = 0.99995); cor_matrix_99.995 #0.2867879
cor_matrix_99.9995 <- quantile(as.matrix(cor_matrix), probs = 0.999995); cor_matrix_99.9995 #0.4190052
cor_matrix_99.99995 <- quantile(as.matrix(cor_matrix), probs = 0.9999995); cor_matrix_99.99995 #0.7425941
cor_matrix_99.999995 <- quantile(as.matrix(cor_matrix), probs = 0.99999995); cor_matrix_99.999995 #0.7684807 
cor_matrix_max <- max(as.matrix(cor_matrix)); cor_matrix_max # 0.8395122

library(dplyr)
library(ggplot2)
library(tibble)

## --- Inputs you already have ---
# resid_expr_blood: genes x samples (G x N)
# v_brain$E       : genes x samples (G_brain x N)  -- N samples should match
# pca_resid_blood <- prcomp(t(resid_expr_blood), center=TRUE, scale.=TRUE)
# pca_summary <- summary(pca_resid_blood)$importance

## 1) Helper to residualize a matrix on given covariates (across samples)
# Y: N x G matrix (samples x genes)
# Z: N x p design matrix (include intercept)
residualize_matrix <- function(Y, Z) {
  M <- diag(nrow(Z)) - Z %*% solve(t(Z) %*% Z) %*% t(Z)
  M %*% Y
}

## 2) Labels for x-axis (Baseline, PC1, PC1+PC2, ..., PC1+...+PC10)
n_pcs <- 10
cum_labels <- c("Baseline",
                paste0("PC1", sapply(2:n_pcs, function(i) paste0("+PC", 2:i))))
cum_labels[2] <- "PC1"  # fix the PC1 label

## 3) Baseline correlation distribution (no PCs regressed out)
cm_vec_baseline <- as.numeric(as.matrix(abs(cor(t(resid_expr_blood),
                                               t(v_brain$E),
                                               method = "spearman"))))

perf_lines <- tibble(
  metric = c("Mean", "Median", "90th", "95th", "99.95th", "99.995th", "99.9995th", "Max"),
  value  = c(
    mean(cm_vec_baseline),
    median(cm_vec_baseline),
    quantile(cm_vec_baseline, 0.90),
    quantile(cm_vec_baseline, 0.95),
    quantile(cm_vec_baseline, 0.9995),
    quantile(cm_vec_baseline, 0.99995),
    quantile(cm_vec_baseline, 0.999995),
    max(cm_vec_baseline)
  )
)

## 4) For k = 0..10, regress out top k PCs from blood and recompute performance
# Scores are N x PCs (samples x PCs) since you did prcomp on t(resid_expr_blood)
scores <- pca_resid_blood$x  # N x min(N, G)

get_mean_abs_corr <- function(expr_blood_mat) {
  cm <- abs(cor(t(expr_blood_mat), t(v_brain$E), method = "spearman"))
  mean(as.numeric(cm))
}

N <- nrow(scores)
Y_blood <- t(resid_expr_blood)  # N x G

mean_series <- numeric(n_pcs + 1)

# k = 0 (baseline, no PCs removed)
mean_series[1] <- get_mean_abs_corr(resid_expr_blood)

# k = 1..10
for (k in 1:n_pcs) {
  Z <- cbind(Intercept = 1, scores[, 1:k, drop = FALSE])  # N x (k+1)
  Y_res <- residualize_matrix(Y_blood, Z)                 # N x G residuals
  expr_blood_k <- t(Y_res)                                # back to G x N
  mean_series[k + 1] <- get_mean_abs_corr(expr_blood_k)
}

pc_perf <- tibble(
  x = 0:n_pcs,
  label = cum_labels,
  mean_abs_spearman = mean_series
)

## 5) Plot: baseline threshold lines + trajectory of mean abs corr after removing PCs
plot <- ggplot() +
  # Horizontal reference lines from baseline distribution
  #geom_hline(data = perf_lines, aes(yintercept = value, color = metric), linewidth = 0.7) +
  # Trajectory of performance as you regress out more PCs
  geom_line(data = pc_perf, aes(x, mean_abs_spearman), linewidth = 1) +
  geom_point(data = pc_perf, aes(x, mean_abs_spearman), size = 2) +
  scale_x_continuous(
    name = "PCs regressed out from blood",
    breaks = 0:n_pcs,
    labels = cum_labels,
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  ylab("Mean |Spearman r| (blood vs brain genes)") +
  guides(color = guide_legend(title = "Baseline thresholds")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 12, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/test.pdf", plot, width = 10, height = 6)
https://hoangd02.u.hpc.mssm.edu/plots/lbp/test.pdf

############not just mean
# ---- 1) Stats helper: mean, median, tail quantiles, max ----
get_corr_stats <- function(expr_blood_mat) {
  cm <- abs(cor(t(expr_blood_mat), t(v_brain$E), method = "spearman"))
  v  <- as.numeric(cm)
  tibble(
    metric = c("Mean","Median","P90","P95","P99.95","P99.995","P99.9995","P99.99995","P99.999995","Max"),
    value  = c(
      mean(v, na.rm = TRUE),
      median(v, na.rm = TRUE),
      quantile(v, 0.90,      names = FALSE, na.rm = TRUE),
      quantile(v, 0.95,      names = FALSE, na.rm = TRUE),
      quantile(v, 0.9995,    names = FALSE, na.rm = TRUE),
      quantile(v, 0.99995,   names = FALSE, na.rm = TRUE),
      quantile(v, 0.999995,  names = FALSE, na.rm = TRUE),
      quantile(v, 0.9999995,  names = FALSE, na.rm = TRUE),
      quantile(v, 0.99999995,  names = FALSE, na.rm = TRUE),
      max(v, na.rm = TRUE)
    )
  )
}

# ---- 2) Labels: Baseline, PC1, PC1+PC2, ..., PC1+...+PCK ----
scores <- pca_resid_blood$x
n_pcs  <- min(15, ncol(scores))
labels_1toK <- sapply(1:n_pcs, function(k) paste(paste0("PC", 1:k), collapse = "+"))
cum_labels  <- c("Baseline", labels_1toK)

# ---- 3) Residualize utility (samples x genes input) ----
residualize_matrix <- function(Y, Z) {
  M <- diag(nrow(Z)) - Z %*% solve(t(Z) %*% Z) %*% t(Z)
  M %*% Y
}

# ---- 4) Compute stats for k = 0..n_pcs ----
N       <- nrow(scores)
Y_blood <- t(resid_expr_blood)  # N x G
stats_list <- vector("list", n_pcs + 1)

# k = 0 (Baseline)
stats_list[[1]] <- get_corr_stats(resid_expr_blood) %>%
  mutate(k = 0, label = "Baseline")

# k = 1..n_pcs
for (k in 1:n_pcs) {
  Z <- cbind(Intercept = 1, scores[, 1:k, drop = FALSE])  # N x (k+1)
  Y_res <- residualize_matrix(Y_blood, Z)                 # N x G
  expr_blood_k <- t(Y_res)                                # G x N
  stats_list[[k + 1]] <- get_corr_stats(expr_blood_k) %>%
    mutate(k = k, label = labels_1toK[k])
}

pc_perf <- bind_rows(stats_list)
write.csv(pc_perf,"/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/15PCs_by_performance_metrics_try2")
# ---- 5) Plot: one line per metric across k ----
# Make labels
lab_df <- pc_perf %>%
  mutate(lab = sprintf("%.3f", value))

# A small vertical offset so text sits just above the point (scale-aware)
y_off <- diff(range(pc_perf$value, na.rm = TRUE)) * 0.02

plot <- ggplot(pc_perf, aes(x = k, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(
    data = lab_df,
    aes(label = lab),
    nudge_y = y_off,          # put text directly above the point
    size = 3,
    show.legend = FALSE,
    check_overlap = TRUE      # drop a label if two collide
  ) +
  scale_x_continuous(
    name   = "PCs regressed out from blood",
    breaks = 0:n_pcs,
    labels = cum_labels
  ) +
  ylab("Absolute Spearman r summary") +
  guides(color = guide_legend(title = "Metric")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 12, hjust = 1),
        panel.grid.minor = element_blank())

ggsave("/hpc/users/hoangd02/www/plots/lbp/test3.pdf", plot, width = 10, height = 6)
https://hoangd02.u.hpc.mssm.edu/plots/lbp/test3.pdf
#ggsave("/hpc/users/hoangd02/www/plots/lbp/test2.pdf", plot, width = 10, height = 6)


#######################################
####base blood vs full brain no residID 
#######################################
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

cor_matrix <- abs(cor(t(resid_expr_blood), t(brain_full_no_residID), method = "spearman"))
dim(cor_matrix) #21046 21356

## 1) Build x-axis labels: "PC1", "PC1+PC2", ..., up to all PCs
# assumes `pca_summary` is like your printed table with a row "Cumulative Proportion"
cumprop <- as.numeric(pca_summary["Cumulative Proportion", ])
n_pcs   <- length(cumprop);n_pcs #225

cum_labels <- paste0("PC1", ifelse(2:n_pcs > 1, paste0("+PC", 2:n_pcs), ""))
# For readability, we’ll show every 5th label; you can change this.
breaks_show <- seq(1, n_pcs, by = max(1, floor(n_pcs/10)))
labels_show <- cum_labels[breaks_show]

## 2) Compute performance thresholds from your correlation matrix
cm_vec <- as.numeric(as.matrix(cor_matrix))

perf <- tibble(
  metric = c("Mean", "Median", "90th", "95th",
             "99.95th", "99.995th", "99.9995th", "99.99995th", "Max"),
  value  = c(
    mean(cm_vec),
    median(cm_vec),
    quantile(cm_vec, 0.90),
    quantile(cm_vec, 0.95),
    quantile(cm_vec, 0.9995),
    quantile(cm_vec, 0.99995),
    quantile(cm_vec, 0.999995),
    quantile(cm_vec, 0.9999995),
    max(cm_vec)
  )
)

## 3) Plot: horizontal lines for each performance metric across the PC axis
# (Optional) we can also annotate cumulative variance as a dashed line on a secondary axis,
# but to keep the y-axis as "performance", we only draw the performance lines here.

plot <- ggplot() +
  geom_hline(data = perf, aes(yintercept = value, color = metric), linewidth = 0.7) +
  # draw a vertical guide at “all PCs”
  geom_vline(xintercept = n_pcs, linetype = 3) +
  # Create an invisible rug of x positions so the axis renders properly
  geom_point(data = tibble(x = 1:n_pcs), aes(x, y = 0), alpha = 0) +
  scale_x_continuous(
    name   = "Cumulative PCs used",
    limits = c(1, n_pcs),
    breaks = breaks_show,
    labels = labels_show,
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  ylab("Performance metric (e.g., correlation)") +
  guides(color = guide_legend(title = "Thresholds")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/test2.pdf", plot, width = 10, height = 6)

##################################################################
##################################################################
##################################################################
form <- ~1 
fit_blood <- dream(v_blood$E, form, blood_metadata)
resid_expr_blood <- residuals(fit_blood)

pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance

pca = pca_resid_blood
resCor=canCorAllAgainstAll_Original(blood_only_metadata,as.data.frame(pca$x[,1:10]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor[,1:10]))
resCor=resCor[ordered_resCor,]

##head(resCor,20)- this is not all outputs
#                                                   PC1         PC2        PC3
# fastqDir_blood                              1.0000000 1.000000000 1.00000000
# SAMPLE_ISMMS_blood                          1.0000000 1.000000000 1.00000000
# LIMS_SEMA4_blood                            1.0000000 1.000000000 1.00000000
# ISM_SEMA4_blood                             1.0000000 1.000000000 1.00000000
# RSM_SEMA4_blood                             1.0000000 1.000000000 1.00000000
# IID_ISMMS                                   0.8551131 0.816056483 0.93242957
# mymet_rin_blood                             0.6311088 0.132517585 0.02419515
# RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood      0.5733436 0.283890762 0.10529995
# RNASeqMetrics_PCT_INTRONIC_BASES_blood      0.4996477 0.231864781 0.32771384
# Lymphocyte_total_wilk_blood                 0.4702496 0.102204789 0.81210824
# Lymphocyte_total_lm22_blood                 0.4689294 0.009721621 0.76570525
# RNASeqMetrics_PCT_MRNA_BASES_blood          0.4529274 0.458032960 0.19962742
# RNASeqMetrics_PCT_USABLE_BASES_blood        0.4501503 0.481321678 0.17978669
# STAR_Insertion_average_length_blood         0.4389597 0.778312492 0.06483234
# InsertSizeMetrics_WIDTH_OF_95_PERCENT_blood 0.4321672 0.077197561 0.16211442
# InsertSizeMetrics_WIDTH_OF_90_PERCENT_blood 0.4313449 0.024740045 0.16467761
# RNASeqMetrics_PCT_UTR_BASES_blood           0.4144708 0.460702308 0.23957936
# mymet_depletionbatch_blood                  0.3938933 0.334788244 0.33702357
# Lymphocyte_total_scp_blood                  0.3787107 0.106081395 0.77859559
# STAR_pct_of_reads_unmapped_other_blood      0.3721347 0.767935321 0.06172074
#                                                     PC4        PC5        PC6
# fastqDir_blood                              1.000000000 1.00000000 1.00000000
# SAMPLE_ISMMS_blood                          1.000000000 1.00000000 1.00000000
# LIMS_SEMA4_blood                            1.000000000 1.00000000 1.00000000
# ISM_SEMA4_blood                             1.000000000 1.00000000 1.00000000
# RSM_SEMA4_blood                             1.000000000 1.00000000 1.00000000
# IID_ISMMS                                   0.842650909 0.81326661 0.85661196
# mymet_rin_blood                             0.052874562 0.06153293 0.34348337
# RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood      0.143321679 0.17248512 0.18410073
# RNASeqMetrics_PCT_INTRONIC_BASES_blood      0.492406967 0.13402893 0.12354593
# Lymphocyte_total_wilk_blood                 0.101208464 0.06892320 0.08822096
# Lymphocyte_total_lm22_blood                 0.133648196 0.09901763 0.01287462
# RNASeqMetrics_PCT_MRNA_BASES_blood          0.326306292 0.18411636 0.04801426
# RNASeqMetrics_PCT_USABLE_BASES_blood        0.306824402 0.23563379 0.06715112
# STAR_Insertion_average_length_blood         0.010000392 0.10108029 0.11264761
# InsertSizeMetrics_WIDTH_OF_95_PERCENT_blood 0.050245318 0.01713785 0.10214457
# InsertSizeMetrics_WIDTH_OF_90_PERCENT_blood 0.022927311 0.02209864 0.21726882
# RNASeqMetrics_PCT_UTR_BASES_blood           0.320355867 0.08614417 0.05528395
# mymet_depletionbatch_blood                  0.431699774 0.73273513 0.46283011
# Lymphocyte_total_scp_blood                  0.033804123 0.12168765 0.11313826
# STAR_pct_of_reads_unmapped_other_blood      0.008342084 0.02468096 0.09741028

### CHECKING AT FORM 5 TO SEE WHAT ELSE IS HIGHLY CORRELATED 
form5 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood
fit_blood <- dream(v_blood$E, form5, blood_metadata)
resid_expr_blood <- residuals(fit_blood)

pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance

pca = pca_resid_blood
resCor=canCorAllAgainstAll_Original(blood_only_metadata,as.data.frame(pca$x[,1:10]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor[,1:10]))
resCor=resCor[ordered_resCor,]

info_all2 <- blood_only_metadata
summ=summary(pca)
SampleByVariable = t(resid_expr_blood) ###change this 
clonename<-rownames(SampleByVariable)

plot_pca_by_metadata(pca, summ, info_all2, clonename,
                     file_prefix = "blood", form_label = "form5_no_residID")

https://hoangd02.u.hpc.mssm.edu/plots/lbp/pca_plots_blood_form5_no_residID_2025-09-04.pdf

##################### ##################### ##################### THE RIGHT ONE
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
# Step 1: Start with intercept-only model
form0 <- ~ 1
fit_blood0 <- dream(v_blood$E, form0, blood_metadata)
resid_expr_blood0 <- residuals(fit_blood0)

# Step 2: PCA on residuals
pca_resid_blood <- prcomp(t(resid_expr_blood0), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance
pca_summary[1:3,1:15]
#                             PC1      PC2      PC3      PC4      PC5      PC6
# Standard deviation     67.58862 58.77984 43.69163 34.21867 27.94690 23.33321
# Proportion of Variance  0.21706  0.16417  0.09070  0.05564  0.03711  0.02587
# Cumulative Proportion   0.21706  0.38123  0.47193  0.52757  0.56468  0.59055
#                             PC7      PC8      PC9     PC10     PC11     PC12
# Standard deviation     21.00671 19.20182 18.76471 16.63870 15.24516 13.91267
# Proportion of Variance  0.02097  0.01752  0.01673  0.01315  0.01104  0.00920
# Cumulative Proportion   0.61151  0.62903  0.64576  0.65892  0.66996  0.67916
#                            PC13     PC14     PC15
# Standard deviation     12.98003 12.82771 11.98671
# Proportion of Variance  0.00801  0.00782  0.00683
# Cumulative Proportion   0.68716  0.69498  0.70181

# Extract sample-level PC scores
pc_scores <- as.data.frame(pca_resid_blood$x)
blood_metadata$PC1 <- pc_scores$PC1
blood_metadata$PC2 <- pc_scores$PC2
# ... add more as needed

# Step 3: Iteratively regress out PCs
form1 <- ~ PC1
fit_blood1 <- dream(v_blood$E, form1, blood_metadata)
resid_expr_blood1 <- residuals(fit_blood1)

form2 <- ~ PC1 + PC2
fit_blood2 <- dream(v_blood$E, form2, blood_metadata)
resid_expr_blood2 <- residuals(fit_blood2)

#streamline this process for 15 PCs: 
# assumes: library(variancePartition); library(limma)
# inputs:
#   v_blood$E        : genes x samples expression matrix
#   blood_metadata   : data.frame with one row per sample (rownames or a column matching colnames(v_blood$E))

# 1) Baseline residuals (~1)
form0 <- ~ 1
fit0  <- dream(v_blood$E, form0, blood_metadata)
resid0 <- residuals(fit0)  # genes x samples

# 2) PCA on baseline residuals
pca_resid <- prcomp(t(resid0), center = TRUE, scale. = TRUE)  # samples x genes -> PC scores per sample

# Keep top 15 PC scores and align to expression sample order
k_max <- 15
pc_scores <- as.data.frame(pca_resid$x[, 1:k_max, drop = FALSE])
# Ensure sample alignment:
pc_scores <- pc_scores[colnames(v_blood$E), , drop = FALSE]

# Rename to PC1..PC15
colnames(pc_scores) <- paste0("PC", seq_len(ncol(pc_scores)))

# Build augmented metadata without mutating your original
md_aug <- cbind(blood_metadata, pc_scores)

# 3) Iteratively regress out PCs 0..15
residuals_by_k <- vector("list", k_max + 1)
names(residuals_by_k) <- paste0("PC", 0:k_max)

# k = 0: Baseline already computed
residuals_by_k[[1]] <- resid0

# k = 1..15
for (k in 1:k_max) {
  form_k <- reformulate(termlabels = paste0("PC", 1:k))  # ~ PC1 + ... + PCk
  fit_k  <- dream(v_blood$E, form_k, md_aug)
  residuals_by_k[[k + 1]] <- residuals(fit_k)
}

# 4) (Optional) variance explained by PCs
pve <- (pca_resid$sdev^2) / sum(pca_resid$sdev^2)  # length = # of PCs available
pve_top15 <- pve[1:k_max]

# --- Outputs ---
# residuals_by_k[["PC0"]]  : baseline residuals (~1)
# residuals_by_k[["PC1"]]  : residuals after regressing PC1
# ...
# residuals_by_k[["PC15"]] : residuals after regressing PC1..PC15
# pve_top15                 : proportion variance explained for the top 15 PCs (for each one)

cor_matrix <- abs(cor(t(resid_expr_blood), t(v_brain$E), method = "spearman"))
dim(cor_matrix) #21046 21356

## computing correlations 
# Define which summary stats to compute
quantiles <- c(0.90, 0.95, 0.9995, 0.99995, 0.999995, 0.9999995)

# Function to compute correlation summary for one residual matrix
summarize_cor <- function(resid_expr, brain_expr, quantiles) {
  cor_mat <- abs(cor(t(resid_expr), t(brain_expr), method = "spearman"))
  flat <- as.vector(cor_mat)  # flatten to 1D
  stats <- c(
    mean     = mean(flat, na.rm = TRUE),
    median   = median(flat, na.rm = TRUE),
    quantile(flat, probs = quantiles, na.rm = TRUE),
    max      = max(flat, na.rm = TRUE)
  )
  return(stats)
}

# Apply to all residualized matrices (PC0..PC15)
cor_summary_list <- lapply(residuals_by_k, summarize_cor, brain_expr = v_brain$E, quantiles = quantiles)

# Convert list to a data.frame for easy plotting/comparison
cor_summary_df <- do.call(rbind, cor_summary_list)
cor_summary_df <- as.data.frame(cor_summary_df)
cor_summary_df$Residualization <- names(residuals_by_k)

# number of PCs you regressed out
n_pcs <- length(residuals_by_k) - 1

# Build cumulative labels
cumulative_labels <- c("Baseline",
                sapply(1:n_pcs, function(k) paste0("PC1 + … + PC", k)))

# Update in cor_summary_df
cor_summary_df$Residualization <- cumulative_labels
cor_summary_df$Residualization[2] <- "PC1"
cor_summary_df$Residualization[3] <- "PC1+PC2"

write.csv(cor_summary_df, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/15PCs_performance_metrics_of_449M")

# Pivot to long format (all numeric metric columns to one column)
plot_df <- cor_summary_df %>%
  pivot_longer(
    cols = c(mean, median, `90%`, `95%`, `99.95%`, `99.995%`, `99.9995%`, `99.99995%`, max),
    names_to = "metric",
    values_to = "value"
  )

# Factor Residualization to preserve order along x-axis
plot_df$Residualization <- factor(plot_df$Residualization,
                                  levels = cor_summary_df$Residualization)

library(scales)  # for number formatting

p1 <- ggplot(plot_df, aes(x = Residualization, y = value, color = metric, group = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = round(value, 3)),    # show values with 3 decimals
            vjust = -0.8, size = 3, show.legend = FALSE) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = "PCs regressed out from blood",
       y = "Absolute Spearman r summary",
       color = "Metric")

ggsave("/hpc/users/hoangd02/www/plots/lbp/15PCs_performance_metrics_449M.pdf", p1, width = 10, height = 6)
https://hoangd02.u.hpc.mssm.edu/plots/lbp/15PCs_performance_metrics_449M.pdf

##################### ALTOGETHER ##################### ##################### 
##################### ##################### ##################### 
### CHECKING AT FORM 5 TO SEE WHAT ELSE IS HIGHLY CORRELATED 
form5 <- ~ RNASeqMetrics_MEDIAN_3PRIME_BIAS_blood + mymet_rin_blood + STAR_Insertion_average_length_blood + Lymphocyte_total_wilk_blood + RNASeqMetrics_PCT_CODING_BASES_blood
fit_blood <- dream(v_blood$E, form5, blood_metadata)
resid_expr_blood <- residuals(fit_blood)

pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance

pca = pca_resid_blood
resCor=canCorAllAgainstAll_Original(blood_only_metadata,as.data.frame(pca$x[,1:15]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor[,1:15]))
resCor=resCor[ordered_resCor,]
resCor[1:10,1:13]

### UPDATING OLD PLOT
# ---- 1) Helper to read a correlation matrix file and compute stats ----
quantiles_named <- c(`90%`=0.90, `95%`=0.95, `99.95%`=0.9995, `99.995%`=0.99995,
                     `99.9995%`=0.999995, `99.99995%`=0.9999995)

compute_stats_from_file <- function(path) {
  x <- fread(path, data.table = FALSE)
  row.names(x) <- x$V1
  x$V1 <- NULL
  x <- abs(as.matrix(x))
  v <- as.numeric(x)

  c(
    Mean   = mean(v, na.rm = TRUE),
    Median = median(v, na.rm = TRUE),
    stats::quantile(v, probs = quantiles_named, names = TRUE, na.rm = TRUE)
  )
}

# ---- 2) Files (form0..form5) ----
files <- c(
  form0 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_brain_spearman_cor_matrix.txt",
  form1 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form1_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",
  form2 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form2_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",
  form3 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",
  form4 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form4_no_indivdualID_brain_baseline_spearman_cor_matrix.txt",
  form5 = "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_baseline_spearman_cor_matrix.txt"
)

# ---- 3) Compute stats for each form ----
stats_list <- lapply(files, compute_stats_from_file)
cor_stats  <- as.data.frame(do.call(rbind, stats_list))
cor_stats$Form <- rownames(cor_stats)

# Keep column order: Mean, Median, quantiles...
cor_stats <- cor_stats[, c("Form", "Mean", "Median", names(quantiles_named))]

# ---- 4) Long format for plotting ----
cor_stats_long <- cor_stats %>%
  mutate(Form = factor(Form, levels = names(files))) %>%
  pivot_longer(-Form, names_to = "Statistic", values_to = "Value")

# ---- 5) Plot (lines, points, and numeric labels) ----
p <- ggplot(cor_stats_long, aes(x = Form, y = Value, group = Statistic, color = Statistic)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.7, size = 3, show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    title = "Correlation Summary Statistics by Form (No Resid Individual ID, Base Brain)",
    x = "Form",
    y = "Absolute Spearman Correlation",
    color = "Statistic"
  )

ggsave(
  filename = "/hpc/users/hoangd02/www/plots/lbp/blood-brain_correlation_summary_by_form_20250904.pdf",
  plot = p,width = 10, height = 6)

########################## ARCHIVED ########################## 
form <- ~1 
fit_blood <- dream(v_blood$E, form, blood_metadata)
resid_expr_blood <- residuals(fit_blood)

pca_resid_blood <- prcomp(t(resid_expr_blood), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_resid_blood)$importance

pca = pca_resid_blood
resCor=canCorAllAgainstAll_Original(blood_only_metadata,as.data.frame(pca$x[,1:15]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor[,1:10]))
resCor=resCor[ordered_resCor,]
resCor[1:10,1:13]

blood_only_metadata = blood_only_metadata[,rownames(resCor)] #columns of blood_metadata2 is in the same order as rownames(resCor), which is descending corr order
#head(resCor,100)
info_all2 <- blood_only_metadata
summ=summary(pca)
SampleByVariable = t(resid_expr_blood) ###change this 
clonename<-rownames(SampleByVariable)

plot_pca_by_metadata(pca, summ, info_all2, clonename,
                     file_prefix = "blood", form_label = "form1")

https://hoangd02.u.hpc.mssm.edu/plots/lbp/pca_plots_blood_form0_2025-05-22.pdf
https://hoangd02.u.hpc.mssm.edu/plots/lbp/pca_plots_blood_form1_2025-05-22residualizedID.pdf
https://hoangd02.u.hpc.mssm.edu/plots/lbp/pca_plots_brain_form1_2025-05-22residualizedID.pdf


##################### ##################### ##################### THE RIGHT ONE
library(data.table)
library(dplyr)
library(ggplot2)
library(variancePartition)
library(tidyr)
library(stringr)
library(tibble)
library(limma)
library(glmnet)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
head(blood_only_metadata)
v_blood$E

brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

# ===========
# PARAMETERS
# ===========
ID_COL <- "IID_ISMMS"     # set to the subject grouping column in blood_only_metadata; set NULL if not needed
K_OUTER <- 5                 # outer CV folds
L_INNER <- 3                 # inner CV folds (per outer-train)
K_PC    <- 15                # number of PCs to consider when ranking covariates
ALPHAS  <- c(0, 0.5, 1)      # elastic-net alpha grid (ridge, elastic, lasso)
SEED    <- 123

# Candidate lists (EDIT to your data)
always_allow <- c(
  "mymet_rin_blood", "STAR_Insertion_average_length_blood ", "Lymphocyte_total_wilk_blood",
  "RNASeqMetrics_PCT_CODING_BASES_blood", "RNASeqMetrics_PCT_INTRONIC_BASES_blood"
)
never_adjust <- c(
  "mymet_sex_blood","mymet_age_blood"
)

# ==================================================
# 0) Helpers: folds, shuffle, metrics, ID filtering
# ==================================================
set.seed(SEED)

make_folds <- function(ids, K = 5, seed = 1L) {
  set.seed(seed)
  n <- length(ids)
  if (is.null(ids)) {
    fold_id <- sample(rep(1:K, length.out = n))
    split(seq_len(n), fold_id)
  } else {
    uids <- unique(ids)
    fold_u <- sample(rep(1:K, length.out = length(uids)))
    lapply(1:K, function(k) which(ids %in% uids[fold_u == k]))
  }
}

pair_shuffle <- function(idx) sample(idx, length(idx))

fisher_z   <- function(r) 0.5 * log((1 + r) / (1 - r))
inv_fisher <- function(z) (exp(2*z)-1)/(exp(2*z)+1)

primary_metric <- function(r_vec) {
  r_vec <- pmax(pmin(r_vec,  0.999999), -0.999999)
  median(fisher_z(r_vec), na.rm = TRUE)
}

# Secondary monitors on held-out
sec_metrics <- function(pred, truth) {
  # pred, truth: matrices [genes x samples] on held-out indices
  G <- nrow(truth)
  r_pear <- r_spear <- r2 <- nrmse <- mae <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    y <- truth[g, ]; yhat <- pred[g, ]
    if (sd(y) > 0 && sd(yhat) > 0) {
      r_pear[g] <- suppressWarnings(cor(yhat, y, method = "pearson"))
      r_spear[g] <- suppressWarnings(cor(yhat, y, method = "spearman"))
    }
    # R2
    sse <- sum((yhat - y)^2)
    sst <- sum((y - mean(y))^2)
    r2[g] <- if (sst > 0) 1 - sse/sst else NA_real_
    # Errors
    rmse <- sqrt(mean((yhat - y)^2))
    sd_y <- sd(y)
    nrmse[g] <- if (sd_y > 0) rmse/sd_y else NA_real_
    mae[g]   <- mean(abs(yhat - y))
  }
  list(
    median_fisher_z_pearson = primary_metric(r_pear),
    median_spearman         = median(r_spear, na.rm = TRUE),
    median_R2               = median(r2, na.rm = TRUE),
    median_NRMSE            = median(nrmse, na.rm = TRUE),
    median_MAE              = median(mae, na.rm = TRUE),
    r_per_gene              = r_pear
  )
}

# ID-like detector to avoid R^2 ~ 1 artifacts
flag_id_like <- function(x, nm, n, max_levels = 30) {
  nuniq <- dplyr::n_distinct(x)
  is_const <- nuniq <= 1
  unique_per_sample <- nuniq > 0.8 * n
  too_many_levels <- is.factor(x) && nlevels(x) > max_levels
  looks_like_path <- is.character(x) && any(grepl("[/\\\\]", x))
  name_matches <- grepl("(id|sample|barcode|lims|sema4|ismms|rsm|fastq|dir|path|uuid)",
                        nm, ignore.case = TRUE)
  is_const || unique_per_sample || too_many_levels || looks_like_path || name_matches
}

clean_meta_for_covsel <- function(meta, extra_drop = character()) {
  n <- nrow(meta)
  keep <- vapply(names(meta),
                 function(nm) !flag_id_like(meta[[nm]], nm, n),
                 logical(1))
  keep[names(meta) %in% extra_drop] <- FALSE
  meta[, keep, drop = FALSE]
}
 
# ===============================================
# 1) Rank metadata by PC association (TRAIN only)
# ===============================================
rank_covariates_by_pc <- function(Y_blood_train, meta_train, Kpc = 15,
                                  always_allow = character(),
                                  never_adjust = character()) {
  # Baseline residuals (~1) so metadata don't steer PCs
  fit0 <- variancePartition::dream(Y_blood_train, ~ 1, meta_train)
  R0   <- residuals(fit0)  # genes x samples

  PCA  <- prcomp(t(R0), center = TRUE, scale. = TRUE)
  Kpc  <- min(Kpc, ncol(PCA$x))
  scores <- as.data.frame(PCA$x[, 1:Kpc, drop = FALSE])
  w <- (PCA$sdev[1:Kpc]^2); w <- w / sum(w)

  # Candidate variables
  drop_names <- c("SampleID", ID_COL, never_adjust)
  cand_vars  <- setdiff(colnames(meta_train), drop_names)
  cand_vars  <- union(intersect(colnames(meta_train), always_allow), cand_vars)

  # per-PC R^2: numeric -> cor^2; factor -> one-way ANOVA R^2
  var_explained <- function(pc, z) {
    if (is.numeric(z)) {
      r <- suppressWarnings(cor(pc, z, use = "pair"))
      if (!is.finite(r)) r <- 0
      r^2
    } else {
      summary(lm(pc ~ z))$r.squared
    }
  }

  score_one_var <- function(z) {
    vk <- sapply(seq_len(Kpc), function(k) var_explained(scores[[k]], z))
    sum(w * vk, na.rm = TRUE)
  }

  tibble(
    var   = cand_vars,
    score = sapply(cand_vars, function(v) score_one_var(meta_train[[v]]))
  ) |>
    arrange(desc(score))
}

# ====================================
# 2) Build a small recipe library
# ====================================
build_recipes <- function(ranked_vars, meta_train, sizes = c(2, 4, 6), extras = list()) {
  # prune near-duplicate numeric variables
  prune_collinear <- function(df, vars, thr = 0.95) {
    keep <- character()
    for (v in vars) {
      if (!length(keep)) { keep <- v; next }
      if (is.numeric(df[[v]])) {
        keep_num <- keep[sapply(keep, \(k) is.numeric(df[[k]]))]
        if (!length(keep_num)) { keep <- c(keep, v); next }
        r <- suppressWarnings(cor(df[keep_num], df[[v]], use = "pair"))
        if (all(abs(r) < thr, na.rm = TRUE)) keep <- c(keep, v)
      } else {
        keep <- c(keep, v)
      }
    }
    keep
  }

  base_vars <- prune_collinear(meta_train, ranked_vars$var)
  recs <- list(R0 = ~ 1)
  for (s in sizes) {
    idx <- seq_len(min(s, length(base_vars)))
    if (length(idx) > 0) recs[[paste0("R", s)]] <- reformulate(base_vars[idx])
  }
  for (nm in names(extras)) recs[[nm]] <- extras[[nm]]
  recs
}

# ===================================================
# 3) Residualize with a recipe (fit on TRAIN only)
# ===================================================
residualize_with_recipe <- function(Y_blood, meta, recipe, train_idx, apply_idx) {
  fit <- variancePartition::dream(Y_blood[, train_idx, drop = FALSE],
                                  recipe, meta[train_idx, , drop = FALSE])
  beta <- fit$coefficients                        # (p x genes)
  X    <- model.matrix(recipe, meta[apply_idx, , drop = FALSE])  # N_apply x p
  Yhat <- t(X %*% beta)                           # genes x N_apply
  Y_blood[, apply_idx, drop = FALSE] - Yhat       # residuals
}

# =======================================================
# 4) Predictor: multi-output glmnet (mgaussian), tuned
# =======================================================
fit_predict_eval <- function(Y_blood_train_res, Y_brain_train,
                             Y_blood_val_res,   Y_brain_val,
                             alphas = c(0, 0.5, 1), seed = 1L) {
  set.seed(seed)
  X_tr <- t(Y_blood_train_res)   # samples x blood_genes
  Y_tr <- t(Y_brain_train)       # samples x brain_genes
  X_va <- t(Y_blood_val_res)
  # choose alpha by inner CV using cv.glmnet (tunes lambda)
  best <- list(score = -Inf, alpha = NA_real_, cvfit = NULL)
  for (a in alphas) {
    cvfit <- cv.glmnet(X_tr, Y_tr, family = "mgaussian", alpha = a, standardize = TRUE)
    # evaluate by predicting on validation
    pred <- predict(cvfit$glmnet.fit, newx = X_va, s = cvfit$lambda.min)[[1]]  # list -> array
    # Convert to [genes x samples]
    pred_mat <- t(pred)
    # Held-out metrics (primary + monitors)
    m <- sec_metrics(pred_mat, Y_brain_val)
    if (m$median_fisher_z_pearson > best$score) best <- list(score = m$median_fisher_z_pearson,
                                                             alpha = a, cvfit = cvfit, pred = pred_mat, metrics = m)
  }
  best
}

# ============================================================
# 5) Inner-CV: choose winning recipe (train-only, no leakage)
# ============================================================
select_recipe_innerCV <- function(Y_blood, meta_blood, Y_brain,
                                  train_idx, L = 3,
                                  ranked_vars_train,
                                  recipes,
                                  alphas = ALPHAS,
                                  seed = 1L) {
  set.seed(seed)
  # make inner folds on TRAIN indices (grouped by subject if available)
  if (!is.null(ID_COL) && ID_COL %in% colnames(meta_blood)) {
    ids_tr <- meta_blood[[ID_COL]][train_idx]
  } else ids_tr <- NULL
  inner_folds <- make_folds(ids_tr, K = L, seed = seed)

  recipe_scores <- tibble(recipe = names(recipes),
                          mean_primary = NA_real_, mean_delta = NA_real_)

  for (ri in seq_along(recipes)) {
    rec <- recipes[[ri]]
    primary_scores <- numeric(length(inner_folds))
    deltas <- numeric(length(inner_folds))

    for (f in seq_along(inner_folds)) {
      val_local_idx  <- train_idx[inner_folds[[f]]]
      trn_local_idx  <- setdiff(train_idx, val_local_idx)

      # residualize blood using rec (fit on inner-train; apply to both)
      Y_blood_tr_res <- residualize_with_recipe(Y_blood, meta_blood, rec, trn_local_idx, trn_local_idx)
      Y_blood_va_res <- residualize_with_recipe(Y_blood, meta_blood, rec, trn_local_idx, val_local_idx)

      # fit predictor + evaluate on inner-val (choose alpha inside)
      best <- fit_predict_eval(Y_blood_tr_res, Y_brain[, trn_local_idx, drop = FALSE],
                               Y_blood_va_res,   Y_brain[, val_local_idx, drop = FALSE],
                               alphas = alphas, seed = seed + f + ri)

      # primary metric on true pairing
      primary_scores[f] <- best$metrics$median_fisher_z_pearson

      # pair-shuffle guardrail on inner-val
      shuf <- sample(seq_along(val_local_idx))
      pred_shuf_metrics <- sec_metrics(best$pred, Y_brain[, val_local_idx[shuf], drop = FALSE])
      deltas[f] <- best$metrics$median_fisher_z_pearson - pred_shuf_metrics$median_fisher_z_pearson
    }

    recipe_scores$mean_primary[ri] <- mean(primary_scores, na.rm = TRUE)
    recipe_scores$mean_delta[ri]   <- mean(deltas, na.rm = TRUE)
  }

  # pick winner: highest mean_primary with mean_delta > 0; tie-break by fewer terms (parsimony)
  recipe_scores <- recipe_scores |>
    mutate(n_terms = sapply(recipes, function(r) ncol(model.matrix(r, meta_blood[train_idx, , drop = FALSE])) - 1)) |>
    arrange(desc(mean_primary), desc(mean_delta), n_terms)

  winner_name <- recipe_scores$recipe[1]
  list(winner = winner_name, table = recipe_scores)
}

# =====================================================
# 6) Outer CV loop: evaluate chosen recipe per fold
# =====================================================
run_outer_cv <- function(Y_blood, meta_blood, Y_brain,
                         K_outer = K_OUTER, L_inner = L_INNER,
                         id_col = ID_COL, seed = SEED) {

  if (!is.null(id_col) && id_col %in% colnames(meta_blood)) {
    ids_vec <- meta_blood[[id_col]]
  } else ids_vec <- NULL

  outer_folds <- make_folds(ids_vec, K = K_outer, seed = seed)

  results <- list()
  winner_names <- character(length(outer_folds))
  inner_tables <- vector("list", length(outer_folds))

  for (k in seq_along(outer_folds)) {
    test_idx  <- outer_folds[[k]]
    train_idx <- setdiff(seq_len(ncol(Y_blood)), test_idx)

    # Clean metadata (drop ID-like)
    meta_train_clean <- clean_meta_for_covsel(meta_blood[train_idx, , drop = FALSE])

    # Rank vars by PC association on TRAIN only
    ranked <- rank_covariates_by_pc(
      Y_blood_train = Y_blood[, train_idx, drop = FALSE],
      meta_train    = meta_train_clean,
      Kpc           = K_PC,
      always_allow  = always_allow,
      never_adjust  = never_adjust
    )

    # Build recipes on TRAIN only
    recipes <- build_recipes(
      ranked_vars = ranked,
      meta_train  = meta_train_clean,
      sizes       = c(2, 4, 6),
      extras      = list(
        R_batch_RIN = ~ batch + RIN,
        R_QC        = ~ batch + RIN + IntragenicRate + IntergenicRate + rRNARate
      )
    )

    # Inner-CV to select recipe
    sel <- select_recipe_innerCV(
      Y_blood = Y_blood, meta_blood = meta_blood, Y_brain = Y_brain,
      train_idx = train_idx, L = L_inner,
      ranked_vars_train = ranked, recipes = recipes,
      alphas = ALPHAS, seed = seed + k
    )
    winner <- sel$winner
    winner_names[k] <- winner
    inner_tables[[k]] <- sel$table

    # Residualize with winning recipe: fit on full TRAIN, apply to TRAIN & TEST
    Y_blood_train_res <- residualize_with_recipe(Y_blood, meta_blood, recipes[[winner]], train_idx, train_idx)
    Y_blood_test_res  <- residualize_with_recipe(Y_blood, meta_blood, recipes[[winner]], train_idx, test_idx)

    # Fit predictor on TRAIN (select alpha via CV), predict TEST
    best <- fit_predict_eval(Y_blood_train_res, Y_brain[, train_idx, drop = FALSE],
                             Y_blood_test_res,  Y_brain[, test_idx,  drop = FALSE],
                             alphas = ALPHAS, seed = seed + 100 + k)

    # Test metrics (true)
    m_true <- best$metrics

    # Test metrics (pair-shuffle)
    shuf <- sample(seq_along(test_idx))
    m_shuf <- sec_metrics(best$pred, Y_brain[, test_idx[shuf], drop = FALSE])

    results[[k]] <- list(
      fold = k,
      winner_recipe = winner,
      test_metrics_true = m_true,
      test_metrics_shuf = m_shuf,
      delta_primary = m_true$median_fisher_z_pearson - m_shuf$median_fisher_z_pearson
    )
    message(sprintf("Outer fold %d: winner=%s | primary(z) %.3f | Δ %.3f",
                    k, winner, m_true$median_fisher_z_pearson, results[[k]]$delta_primary))
  }

  list(results = results, winners = winner_names, inner_tables = inner_tables)
}

# =========================================================
# 7) Freeze global recipe & final CV with fixed recipe
# =========================================================
freeze_global_recipe <- function(winner_names) {
  tab <- sort(table(winner_names), decreasing = TRUE)
  as.character(names(tab)[1])
}

final_cv_with_fixed_recipe <- function(Y_blood, meta_blood, Y_brain,
                                       fixed_recipe, K_outer = K_OUTER,
                                       id_col = ID_COL, seed = SEED) {
  if (!is.null(id_col) && id_col %in% colnames(meta_blood)) ids_vec <- meta_blood[[id_col]] else ids_vec <- NULL
  folds <- make_folds(ids_vec, K = K_outer, seed = seed + 999)

  metrics <- vector("list", length(folds))
  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(ncol(Y_blood)), test_idx)

    # Fit residualization on TRAIN; apply to TRAIN & TEST
    Y_blood_train_res <- residualize_with_recipe(Y_blood, meta_blood, fixed_recipe, train_idx, train_idx)
    Y_blood_test_res  <- residualize_with_recipe(Y_blood, meta_blood, fixed_recipe, train_idx, test_idx)

    # Train predictor (alpha via CV) & predict TEST
    best <- fit_predict_eval(Y_blood_train_res, Y_brain[, train_idx, drop = FALSE],
                             Y_blood_test_res,  Y_brain[, test_idx,  drop = FALSE],
                             alphas = ALPHAS, seed = seed + 2000 + k)

    # True + shuffled
    m_true <- best$metrics
    shuf <- sample(seq_along(test_idx))
    m_shuf <- sec_metrics(best$pred, Y_brain[, test_idx[shuf], drop = FALSE])

    metrics[[k]] <- c(
      fold = k,
      primary_true = m_true$median_fisher_z_pearson,
      primary_shuf = m_shuf$median_fisher_z_pearson,
      delta = m_true$median_fisher_z_pearson - m_shuf$median_fisher_z_pearson,
      median_spearman = m_true$median_spearman,
      median_R2       = m_true$median_R2,
      median_NRMSE    = m_true$median_NRMSE,
      median_MAE      = m_true$median_MAE
    )
  }
  as.data.frame(do.call(rbind, metrics))
}

# ===========================
# 8) RUN THE WHOLE PIPELINE
# ===========================
# Assumes your objects are already aligned by sample order:
# v_blood$E, blood_only_metadata, brain_full_no_residID, brain_metadata

# Outer CV (nested selection of recipe per fold)
outer <- run_outer_cv(
  Y_blood   = v_blood$E,
  meta_blood= blood_only_metadata,
  Y_brain   = brain_full_no_residID,
  K_outer   = K_OUTER,
  L_inner   = L_INNER,
  id_col    = ID_COL,
  seed      = SEED
)

# Decide on a global recipe
global_name <- freeze_global_recipe(outer$winners)
message(sprintf("Global recipe (frozen): %s", global_name))

# Reconstruct that recipe from an inner-table (or rebuild from a ranked list on full data if needed)
# Easiest: rebuild recipes on ALL data using ranker, then pick by name:
meta_clean_all <- clean_meta_for_covsel(blood_only_metadata)
ranked_all <- rank_covariates_by_pc(
  Y_blood_train = v_blood$E,
  meta_train    = meta_clean_all,
  Kpc           = K_PC,
  always_allow  = always_allow,
  never_adjust  = never_adjust
)
recipes_all <- build_recipes(
  ranked_vars = ranked_all,
  meta_train  = meta_clean_all,
  sizes       = c(2, 4, 6),
  extras      = list(
    R_batch_RIN = ~ batch + RIN,
    R_QC        = ~ batch + RIN + IntragenicRate + IntergenicRate + rRNARate
  )
)
fixed_recipe <- recipes_all[[global_name]]

# Final CV with frozen global recipe (report these)
final_cv <- final_cv_with_fixed_recipe(
  Y_blood    = v_blood$E,
  meta_blood = blood_only_metadata,
  Y_brain    = brain_full_no_residID,
  fixed_recipe = fixed_recipe,
  K_outer    = K_OUTER,
  id_col     = ID_COL,
  seed       = SEED
)

print(final_cv)
# Summarize primary metric in r-space if you want:
median_r_z <- median(as.numeric(final_cv$primary_true), na.rm = TRUE)
cat("Final (median Fisher-z Pearson, test folds):", median_r_z, "\n",
    "Back-transformed r:", inv_fisher(median_r_z), "\n")