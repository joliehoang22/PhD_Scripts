library(ggpubr)
library(ggrastr)
library(effsize)
library(dplyr)
library(purrr)
library(tibble)
library(reshape2)
library(rstatix)

# Load all results
out_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap"

spearman_cross_timepoint_same_genes_20250923.rds
spearman_cross_timepoint_same_genes_abs_20250923.rds

spearman_same_timepoint_same_genes_combined_20250923.rds
spearman_same_timepoint_same_genes_combined_abs_20250923.rds

spearman_same_timepoint_same_genes_t0_only_20250923.rds
spearman_same_timepoint_same_genes_t0_only_ab_20250923.rds

spearman_same_timepoint_same_genes_t1_only_20250923.rds
spearman_same_timepoint_same_genes_t1_only_abs_20250923.rds

#cor_combined <- readRDS(file.path(out_dir, "spearman_same_timepoint_combined_20250923.rds"))
cor_t0_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_same_genes_t0_only_20250923.rds"))
cor_t1_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_same_genes_t1_only_20250923.rds"))
cor_cross <- readRDS(file.path(out_dir, "spearman_cross_timepoint_same_genes_20250923.rds"))

#dim for all is 21046 21356 ~~ it's just across different numbers of samples

# Basic dimensions and summaries
cat("Dimensions:\n")
cat("T0 only:", dim(cor_t0_only), "\n") 
cat("T1 only:", dim(cor_t1_only), "\n")
cat("Cross timepoint:", dim(cor_cross), "\n")

# Distribution summaries
summary(as.vector(cor_t0_only))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.798686 -0.079446 -0.001142 -0.001042  0.077255  0.853449 
summary(as.vector(cor_t1_only))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.827502 -0.085709 -0.001790 -0.001964  0.081914  0.894669 
summary(as.vector(cor_cross))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.778516 -0.054631  0.001398  0.001488  0.057565  0.852934

############ NON ABSOLUTE ############
######################################
cor_comparison <- data.frame(
  T0_only = as.vector(cor_t0_only),
  T1_only = as.vector(cor_t1_only),
  Cross_timepoint = as.vector(cor_cross)
)
cor_sample <- melt(cor_comparison)

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
## THIS IS FOR ALL SAME GENE PAIRS!!! not just 1 million 
##this is in one of the work_fri screen 


######################
###################### PLOTTING 
######################
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

paired_tests <- FALSE        # set TRUE if rows are matched by subject
p_adjust     <- "BH"       # use "BH" for FDR if you like
x_levels     <- c("T0_only","T1_only","Cross_timepoint")

## ---- 1) long format ----
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

## ---- 4) Cohen's d computed in the SAME order as wilcoxon ----
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
annot_df
#   group1  group2                 p p_label        d d_label   label   y.position
#   <chr>   <chr>              <dbl> <chr>      <dbl> <chr>     <chr>        <dbl>
# 1 T0_only Cross_timepoint 0        <0.001  -0.0307  d = -0.03 d = -0…       1.08
# 2 T1_only Cross_timepoint 0        <0.001  -0.0419  d = -0.04 d = -0…       1.27
# 3 T0_only T1_only         1.03e-46 <0.001   0.00756 d = 0.01  d = 0.…       1.46

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

ggsave("/hpc/users/hoangd02/www/plots/lbp/timepoint_SAME_GENES_cor_1million_per_group_with_wilcoxontest.pdf", p, width = 10, height = 8)
