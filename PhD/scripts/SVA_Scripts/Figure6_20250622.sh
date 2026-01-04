###### ###### ###### This figure includes 1-10 SV heatmap and Wilcoxon/Fisher DE enrichment 

###################################################################################
###### 1-10SV OF BP AND SZ from BP_and_SZ_1_SV_at_a_time_correlations.sh ##########
###################################################################################
https://hoangd02.u.hpc.mssm.edu/plots/bp_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values.pdf
https://hoangd02.u.hpc.mssm.edu/plots/sz_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values.pdf

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
library(tidyr)
library(ggrastr)
library(ggrepel)
library(patchwork)  
library(readr)

####### AD from MSBB code for publication 11112024.sh ######
############################################################
# Load the nsv=0 datasets
msbb_no_sva <- read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_no_sva_11.18.txt") #
rosmap_no_sva <- fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_11.18.txt", data.table = FALSE) #19476     7 |this code was used to make heatmaps - redo the heatmaps with this
#rosmap_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt", data.table = FALSE) #19476     7 | this code was used in DEA_NB.sh - FINAL VERSION | this has logFC_j so the correlation wouldn't work
rosmap_be <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_1.16.txt", data.table = FALSE) #26618 x 7 | FINAL VERSION

##check to see if this is identical -- small neglible difference, proceed with 11.18.txt one 
# rosmap_no_sva_ordered <- rosmap_no_sva[order(rosmap_no_sva$X),]
# rosmap_no_sva2_ordered <- rosmap_no_sva2[order(rosmap_no_sva2$X),]
# colnames(rosmap_no_sva2_ordered)[colnames(rosmap_no_sva2_ordered) == "logFC_j"] <- "logFC"
# rosmap_no_sva_ordered$X <- sub("\\..*", "", rosmap_no_sva_ordered$X)

# identical(rosmap_no_sva2_ordered, rosmap_no_sva_ordered) #FALSE
# identical(colnames(rosmap_no_sva_ordered), colnames(rosmap_no_sva2_ordered)) #TRUE
# identical(rownames(rosmap_no_sva_ordered), rownames(rosmap_no_sva2_ordered)) #FALSE

# all.equal(rosmap_no_sva_ordered, rosmap_no_sva2_ordered)
# # [1] "Attributes: < Component “row.names”: Mean relative difference: 0.0001010336 >"
# # [2] "Component “logFC”: Mean relative difference: 4.175885e-08"                    
# # [3] "Component “t”: Mean relative difference: 4.588178e-08"                        
# # [4] "Component “P.Value”: Mean relative difference: 7.316379e-08"                  
# # [5] "Component “adj.P.Val”: Mean relative difference: 8.195576e-08" 

# #the differences are extremely small, almost certainly due to floating point precision issues in R.

# Convert nsv=0 datasets to data frames for consistency
msbb_no_sva <- as.data.frame(msbb_no_sva)
rosmap_no_sva <- as.data.frame(rosmap_no_sva)

# Define paths for sva files
msbb_path <- "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/"
rosmap_path <- "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/"

# Create file names with full paths for MSBB and ROSMAP
msbb_files <- paste0(msbb_path, "msbb_sva_be_", 1:27, "_SVs.csv")
rosmap_files <- paste0(rosmap_path, "rosmap_sva_be_", 1:27, "_SVs.csv")

# Load sva datasets
msbb_list <- lapply(msbb_files, read.csv)
rosmap_list <- lapply(rosmap_files, read.csv)

# Convert sva datasets to data frames for consistency
msbb_list <- lapply(msbb_list, as.data.frame)
rosmap_list <- lapply(rosmap_list, as.data.frame)

# Add nsv=0 datasets to the beginning of each list
msbb_list <- c(list(msbb_no_sva), msbb_list)
rosmap_list <- c(list(rosmap_no_sva), rosmap_list)

# Confirm the lists have 28 elements each (nsv=0 to nsv=27)
cat("MSBB list length:", length(msbb_list), "\n")
cat("ROSMAP list length:", length(rosmap_list), "\n")

##genes are not aligned, check: 
# head(msbb_list[[2]])
# head(msbb_list[[1]])
# head(msbb_list[[22]])

# head(rosmap_list[[28]]) 
# head(rosmap_list[[2]])

tmp=lapply(rosmap_list,function(x){x$adj.P.Val<=0.05})
sapply(tmp,table)

##align the genes - MSBB
# Step 1: Identify the common set of genes across all datasets
common_genes <- Reduce(intersect, lapply(msbb_list, function(df) df$X))

# Step 2: Subset each dataset to include only the common genes
msbb_list <- lapply(msbb_list, function(df) df[df$X %in% common_genes, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
msbb_list <- lapply(msbb_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check <- all(sapply(msbb_list, function(df) identical(df$X, msbb_list[[1]]$X)))

if (alignment_check) {
  cat("Genes are successfully aligned across all datasets.\n")
} else {
  cat("Gene alignment failed. Check for issues in the datasets.\n")
}

##align the genes - ROSMAP
# Step 1: Identify the common set of genes across all datasets in rosmap_list
common_genes_rosmap <- Reduce(intersect, lapply(rosmap_list, function(df) df$X))

# Step 2: Subset each dataset in rosmap_list to include only the common genes
rosmap_list <- lapply(rosmap_list, function(df) df[df$X %in% common_genes_rosmap, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
rosmap_list <- lapply(rosmap_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check_rosmap <- all(sapply(rosmap_list, function(df) identical(df$X, rosmap_list[[1]]$X)))

if (alignment_check_rosmap) {
  cat("Genes are successfully aligned across all ROSMAP datasets.\n")
} else {
  cat("Gene alignment failed for ROSMAP datasets. Check for issues in the datasets.\n")
}

# Initialize matrices for Spearman correlation and overlap counts
cor_results_matrix <- matrix(NA, nrow = 28, ncol = 28)
overlap_counts_matrix <- matrix(NA, nrow = 28, ncol = 28)

# Loop through all combinations of MSBB and ROSMAP SVs (0 to 27)
for (i in 1:28) {
  for (j in 1:28) {
    # Merge the i-th MSBB dataset with the j-th ROSMAP dataset by 'X'
    sv_combo <- merge(
      msbb_list[[i]], 
      rosmap_list[[j]], 
      by = 'X', 
      all.x = TRUE, 
      all.y = TRUE, 
      suffixes = c('_msbb', '_rosmap')
    )
    
    # Calculate Spearman correlation if logFC columns exist
    if ("logFC_msbb" %in% names(sv_combo) && "logFC_rosmap" %in% names(sv_combo)) {
      valid_data <- !is.na(sv_combo$logFC_msbb) & !is.na(sv_combo$logFC_rosmap)
      if (sum(valid_data) > 0) {
        cor_test <- cor.test(
          sv_combo$logFC_msbb[valid_data], 
          sv_combo$logFC_rosmap[valid_data], 
          method = "spearman", 
          exact = FALSE
        )
        cor_results_matrix[i, j] <- cor_test$estimate
      }
    }
    
    # Count overlapping significant genes
    if ("adj.P.Val_msbb" %in% names(sv_combo) && "adj.P.Val_rosmap" %in% names(sv_combo)) {
      significant_genes <- sv_combo$adj.P.Val_msbb <= 0.05 & sv_combo$adj.P.Val_rosmap <= 0.05
      overlap_counts_matrix[i, j] <- sum(significant_genes, na.rm = TRUE)
    }
  }
}

# Convert matrices into long-form data frames
correlation_df <- as.data.frame(as.table(cor_results_matrix))
names(correlation_df) <- c("MSBB_SV", "ROSMAP_SV", "Correlation")
# head(correlation_df)

overlap_df <- as.data.frame(as.table(overlap_counts_matrix))
names(overlap_df) <- c("MSBB_SV", "ROSMAP_SV", "OverlapCount")

# Map SV=0 to SV=27 for the new dimension
correlation_df$MSBB_SV <- as.numeric(correlation_df$MSBB_SV) - 1
correlation_df$ROSMAP_SV <- as.numeric(correlation_df$ROSMAP_SV) - 1
overlap_df$MSBB_SV <- as.numeric(overlap_df$MSBB_SV) - 1
overlap_df$ROSMAP_SV <- as.numeric(overlap_df$ROSMAP_SV) - 1

# Merge correlation and overlap data frames
heatmap_df <- merge(correlation_df, overlap_df, by = c("MSBB_SV", "ROSMAP_SV"), all.x = TRUE)

# Ensure all combinations from SV=0 to SV=27 are present
all_svs <- expand.grid(MSBB_SV = 0:27, ROSMAP_SV = 0:27)
heatmap_df <- merge(all_svs, heatmap_df, by = c("MSBB_SV", "ROSMAP_SV"), all.x = TRUE)

# Calculate the total number of significant genes in MSBB
total_significant_msbb <- sum(msbb_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_msbb #682 DEGs

# Calculate the total number of significant genes in ROSMAP
total_significant_rosmap <- sum(rosmap_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_rosmap #1527 DEGs

# Function to calculate total significant genes for a given list of datasets
calculate_significant_genes <- function(data_list) {
  sapply(data_list, function(data) sum(data$adj.P.Val <= 0.05, na.rm = TRUE))
}

# Calculate the total number of significant genes for MSBB and ROSMAP for each SV
msbb_significant_genes <- calculate_significant_genes(msbb_list)
rosmap_significant_genes <- calculate_significant_genes(rosmap_list)

# Create custom axis labels for MSBB and ROSMAP
msbb_labels <- paste0(0:27, " SV (", msbb_significant_genes, " DEGs)")
rosmap_labels <- paste0(0:27, " SV (", rosmap_significant_genes, " DEGs)")

# Add labels to the heatmap data frame for plotting
heatmap_df$MSBB_Label <- factor(
  heatmap_df$MSBB_SV, 
  levels = 0:27, 
  labels = msbb_labels
)
heatmap_df$ROSMAP_Label <- factor(
  heatmap_df$ROSMAP_SV, 
  levels = 0:27, 
  labels = rosmap_labels
)

###############################################################################
# Subset data for SVs = 0 to 10
heatmap_df_subset1 <- subset(heatmap_df, MSBB_SV <= 10 & ROSMAP_SV <= 10)

# ps1 <- ggplot(heatmap_df_subset1, aes(x = MSBB_Label, y = ROSMAP_Label, fill = Correlation)) +
#   geom_tile(color = "white", size = 0.1) +
#   geom_text(
#     aes(label = ifelse(
#       is.na(Correlation), "", 
#       sprintf("%.2f\n(%s)", Correlation, OverlapCount)
#     )), 
#     color = "black", size = 3, lineheight = 0.9  # reduced from 4.5
#   ) +
#   scale_fill_gradient2(mid = "#FFFFC5", high = "#e0413f", midpoint = 0.15) +
#   #scale_fill_gradient2(mid = "#f9c646", high = "#e0413f", midpoint = 0.2) +
#   #scale_fill_gradient2(low = "#a0f600", high = "#d633a6") + #this is equivalent to most pink and low = white
#   coord_fixed() +
#   theme_minimal(base_size = 10) +  # reduced from 14
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 8),
#     axis.title.x = element_text(size = 10, face = "bold"),
#     axis.title.y = element_text(size = 10, face = "bold"),
#     legend.title = element_text(size = 9),
#     legend.text = element_text(size = 8),
#     legend.position = "bottom",
#     plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
#   ) +
#   labs(
#     title = "Spearman Correlation Heatmap for AD (SVs = 0 to 10)",
#     x = "MSBB",
#     y = "ROSMAP",
#     fill = bquote("Correlation "*rho)
#     #fill = bquote("Spearman's\ncorrelation"*rho)
#   )

# ps1 <- ggplot(heatmap_df_subset1, aes(x = MSBB_Label, y = ROSMAP_Label, fill = Correlation)) +
#   geom_tile(color = "white", size = 0.1) +
#   geom_text(
#     aes(label = ifelse(is.na(Correlation), "", sprintf("%.2f\n(%s)", Correlation, OverlapCount))),
#     color = "black", size = 4, lineheight = 0.9
#   ) +
#   scale_fill_gradient2(mid = "#FFFFC5", high = "#e0413f", midpoint = 0.15) +
#   coord_fixed() +
#   theme_minimal(base_size = 16) +
#   theme(
#     legend.position = "top", legend.margin = margin(t = 6, b = 6),
#     plot.margin = margin(5, 5, 5, 5),
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 18, face = "bold"),
#     axis.title.y = element_text(size = 18, face = "bold"),
#     legend.title = element_text(size = 16, face = "bold"),
#     legend.text = element_text(size = 14),
#     #plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
#   ) +
#   labs(
#     #title = "Spearman Correlation Heatmap for AD (SVs = 0 to 10)",
#     x = "MSBB",
#     y = "ROSMAP",
#     fill = bquote("Correlation "*rho)
#   )

##use this for legend
ps1 <- ggplot(heatmap_df_subset1, aes(x = MSBB_Label, y = ROSMAP_Label, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  geom_text(
    aes(label = ifelse(is.na(Correlation), "", sprintf("%.2f\n(%s)", Correlation, OverlapCount))),
    color = "black", size = 4, lineheight = 0.9) +
  scale_fill_gradient2(
    mid = "#FFFFC5", high = "#e0413f", midpoint = 0.15,
    guide = guide_colorbar(
      barwidth = 20, barheight = 0.8,
      title.position = "top", title.hjust = 0.5)) +
  coord_fixed() +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none", #change to top
    legend.box.margin = margin(t = 2, b = 2),       
    legend.margin = margin(t = 0, b = 0),      
    legend.text = element_text(size = 14),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)) +
  labs(
    x = "MSBB",
    y = "ROSMAP",
    fill = bquote("Correlation "*rho))

ggsave("/hpc/users/hoangd02/www/plots/rosmap_msbb_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_20250622.pdf", plot = ps1, width = 8, height = 8)
# ggsave("/hpc/users/hoangd02/www/plots/rosmap_msbb_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_20250519.pdf", plot = ps1, width = 8, height = 8)

###################################################################################
########### Fisher/Wilcoxon Enrichment Code in AD from DEA_NB.sh ##################
###################################################################################
##filtering out leek
library(rstatix)
library(dplyr)
library(MKdescr)

allResults <- allResults %>% filter(SVA != "leek")

# Calculate -log10(p-value), handle infinite values
allResults$log10pval <- -log10(allResults$pval)
allResults$log10pval[!is.finite(allResults$log10pval)] <- -log10(.Machine$double.xmin)

# MSBB pairwise test
stat.test1 <- allResults[allResults$data == "MSBB",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test1$data <- "MSBB"
attr(stat.test1, "args")$data <- allResults
stat.test1 <- stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)

# ROSMAP pairwise test
stat.test2 <- allResults[allResults$data == "ROSMAP",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test2$data <- "ROSMAP"
attr(stat.test2, "args")$data <- allResults
stat.test2 <- stat.test2 %>% add_xy_position(x = "data", dodge = 0.8)

# Combine stats
stat.test <- rbind(stat.test1, stat.test2)
#stat.test$data <- NULL
stat.test$p.adj.signif[stat.test$p.adj < 0.05] <- "*"
stat.test$p.adj.signif[stat.test$p.adj < 0.005] <- "**"
stat.test$p.adj.signif[stat.test$p.adj < 0.0005] <- "***"

# Create beeswarm plot, original 
# ggplot(allResults, aes(x = data, y = log10pval, color = SVA)) +
#   geom_beeswarm(dodge.width = 0.8, size = 1.5) +
#   #stat_pvalue_manual(stat.test, step.increase = 0.0, label = "p.adj") +
#   theme_bw() +
#   labs(y = expression(-log[10](p-value)),
#         x = "Datasets",
#         color = "Methods")

# Step 1: Prep data
allResults <- allResults %>%
  filter(SVA %in% c("be", "none")) %>%
  mutate(
    SVA = factor(SVA, levels = c("none", "be")),  # 1. Switch order
    facet_label = case_when(
      data == "ROSMAP" ~ "ROSMAP", ## double check this number
      data == "MSBB" ~ "MSBB" ## double check this number
    )
  )
allResults$pval <- ifelse(allResults$pval == 0, 1e-300, allResults$pval)

# Step 2: Create custom y-axis breaks
# breaks_rosmap <- 10^-(seq(0, 300, by = 20))  # 1, 1e-100, 1e-200, 1e-300
# breaks_msbb   <- 10^-(seq(0, 300, by = 20))   # 1, 1e-10, 1e-20, 1e-30

# # Plot - ryan approved (original)
# p2 <- ggplot(allResults, aes(x = SVA, y = pval, group = DElist, color = DElist)) +
#   geom_point(size = 4, shape = 16) +       
#   geom_line(alpha = 0.25) +
#   scale_y_neglog10( #significance scale 
#     breaks = c(breaks_rosmap, breaks_msbb),  
#     labels = scales::label_log()
#   ) +
#   facet_wrap(~facet_label, scales = "free_y") +  
#   scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
#   theme_bw() +
#   theme(
#     strip.text = element_text(size = 14, face = "bold"),  
#     axis.text.x = element_text(size = 14),  
#     axis.text.y = element_text(size = 14),                         
#     axis.title.x = element_text(size = 16, face = "bold"),         
#     axis.title.y = element_text(size = 16, face = "bold"),         
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 14)
#   ) +
#   labs(
#     y = "Fisher's Exact Test p-value",
#     x = "SVA Methods",
#     color = "DEG List"
#   )

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/AD_DEA_fisher_wilcoxon_beeswarm2.pdf", plot = p2, width = 14, height = 10)

allResults$pval <- ifelse(allResults$pval == 0, 1e-300, allResults$pval)

# Step 1: Define breaks
breaks_rosmap <- 10^-(c(seq(0, 100, by = 20), 200, 300))
breaks_msbb   <- 10^-(seq(0, 100, by = 20))

# Step 2: Split data
rosmap_data <- subset(allResults, grepl("ROSMAP", facet_label))
msbb_data   <- subset(allResults, grepl("MSBB", facet_label))

# Step 3: ROSMAP plot
p_rosmap <- ggplot(rosmap_data, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25, position = position_nudge(x = 0)) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5, position = position_nudge(x = 0)) +
  scale_y_neglog10(breaks = breaks_rosmap, labels = scales::label_log()) +
  facet_wrap(~facet_label, scales = "free_y") +
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

# Step 4: MSBB plot
p_msbb <- ggplot(msbb_data, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25, position = position_nudge(x = 0)) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5, position = position_nudge(x = 0)) +
  scale_y_neglog10(breaks = breaks_msbb, labels = scales::label_log()) +
  facet_wrap(~facet_label, scales = "free_y") +
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

# Step 5: Combine plots
# p2 <- p_rosmap + p_msbb

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/AD_DEA_fisher_wilcoxon_beeswarm_20250622.pdf", plot = p2, width = 14, height = 10)

combo1 <- ps1 + p_rosmap + p_msbb + plot_layout(widths = c(3, 2, 2))

# combo <- ps1 + p_rosmap + p_msbb + plot_layout(widths = c(1, 1.2, 1.2))
# combo <- ps1 / (p_rosmap + p_msbb) +
#   plot_layout(heights = c(1.3, 1.7), guides = "collect")
# combo <- ps1 / (p_rosmap + p_msbb) + plot_layout(heights = c(1.3, 1.7))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_fewSVs_AD_DEA_fisher_wilcoxon_beeswarm_20251029.pdf", plot = combo1, width = 21, height = 8)

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_fewSVs_AD_DEA_fisher_wilcoxon_beeswarm_20250623.pdf", plot = combo1, width = 21, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/sva/sim_results_oct2024/heatmap_fewSVs_AD_DEA_fisher_wilcoxon_beeswarm_20250622.pdf

########################################################
##################### FULL HEATMAP #####################
########################################################
ad_full <- ggplot(heatmap_df, aes(x = MSBB_Label, y = ROSMAP_Label, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  geom_text(
    aes(label = ifelse(
      is.na(Correlation), "", 
      sprintf("%.2f\n(%s)", Correlation, OverlapCount)
    )), 
    color = "black", size = 4, lineheight = 0.9  
  ) +
  scale_fill_gradient2(mid = "#FFFFC5", high = "#e0413f", midpoint = 0.15) +
  #scale_fill_gradient2(low = "white", high = "#A3392F") +
  coord_fixed() +
  theme_minimal(base_size = 16) +  # Increases overall theme text size
  theme(
    #plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
  ) +
  labs(
    #title = "Spearman Correlation Heatmap: MSBB vs ROSMAP",
    x = "MSBB",
    y = "ROSMAP",
    fill = bquote("Correlation "*rho))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_fewSVs_AD_full_20250623.pdf", plot = ad_full, width = 17, height = 17)

#################### BIPOLAR ###############################
############################################################
output_dir <- "/sc/arion/projects/mscic1/results/jolie/CMC/BP/"

hbcc_files <- paste0(output_dir, "bp_hbcc_sva_be_", 0:29, "_SVs.csv") #limiting factor with 29 SVs
mssm_files <- paste0(output_dir, "bp_mssm_sva_be_", 0:29, "_SVs.csv")

# Load sva datasets
hbcc_list <- lapply(hbcc_files, read.csv)
mssm_list <- lapply(mssm_files, read.csv)

# Convert sva datasets to data frames for consistency
hbcc_list <- lapply(hbcc_list, as.data.frame)
mssm_list <- lapply(mssm_list, as.data.frame)

# Confirm the lists have 30 elements each (nsv=0 to nsv=29)
cat("hbcc list length:", length(hbcc_list), "\n")
cat("mssm list length:", length(mssm_list), "\n")

##genes are not aligned, check: 
head(hbcc_list[[2]])
head(mssm_list[[29]])

tmp=lapply(mssm_list,function(x){x$adj.P.Val<=0.05})
sapply(tmp,table)

##align the genes - hbcc
# Step 1: Identify the common set of genes across all datasets
common_genes <- Reduce(intersect, lapply(hbcc_list, function(df) df$X))

# Step 2: Subset each dataset to include only the common genes
hbcc_list <- lapply(hbcc_list, function(df) df[df$X %in% common_genes, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
hbcc_list <- lapply(hbcc_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check <- all(sapply(hbcc_list, function(df) identical(df$X, hbcc_list[[1]]$X)))

if (alignment_check) {
  cat("Genes are successfully aligned across all datasets.\n")
} else {
  cat("Gene alignment failed. Check for issues in the datasets.\n")
}

##align the genes - mssm
# Step 1: Identify the common set of genes across all datasets in mssm_list
common_genes_mssm <- Reduce(intersect, lapply(mssm_list, function(df) df$X))

# Step 2: Subset each dataset in mssm_list to include only the common genes
mssm_list <- lapply(mssm_list, function(df) df[df$X %in% common_genes_mssm, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
mssm_list <- lapply(mssm_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check_mssm <- all(sapply(mssm_list, function(df) identical(df$X, mssm_list[[1]]$X)))

if (alignment_check_mssm) {
  cat("Genes are successfully aligned across all mssm datasets.\n")
} else {
  cat("Gene alignment failed for mssm datasets. Check for issues in the datasets.\n")
}

# Initialize a matrix to store Spearman correlation results (30x30)
cor_results_matrix <- matrix(NA, nrow = 30, ncol = 30)

# Loop through all combinations of hbcc and mssm SVs (0 to 29)
for (i in 1:30) {
  for (j in 1:30) {
    # Merge the i-th hbcc dataset with the j-th mssm dataset by 'X'
    sv_combo <- merge(
      hbcc_list[[i]], 
      mssm_list[[j]], 
      by = 'X', 
      all.x = TRUE, 
      all.y = TRUE, 
      suffixes = c('_hbcc', '_mssm')
    )
    
    # Check if logFC columns exist and have non-NA values
    if ("logFC_hbcc" %in% names(sv_combo) && "logFC_mssm" %in% names(sv_combo)) {
      valid_data <- !is.na(sv_combo$logFC_hbcc) & !is.na(sv_combo$logFC_mssm)
      if (sum(valid_data) > 0) {
        cor_test <- cor.test(
          sv_combo$logFC_hbcc[valid_data], 
          sv_combo$logFC_mssm[valid_data], 
          method = "spearman", 
          exact = FALSE
        )
        cor_results_matrix[i, j] <- cor_test$estimate
      }
    }
  }
}

# Convert the matrix into a long-form data frame
correlation_df <- as.data.frame(as.table(cor_results_matrix))
names(correlation_df) <- c("hbcc_SV", "mssm_SV", "Correlation")

# Map SV=0 to SV=29 for the new dimension
correlation_df$hbcc_SV <- as.numeric(correlation_df$hbcc_SV) - 1
correlation_df$mssm_SV <- as.numeric(correlation_df$mssm_SV) - 1

# Ensure all combinations from SV=0 to SV=29 are present
all_svs <- expand.grid(hbcc_SV = 0:29, mssm_SV = 0:29)
correlation_df <- merge(all_svs, correlation_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

##ADDING THE COUNT OF OVERLAPPING GENES ON THE TILES INSTEAD OF THE SPEARMAN RHO
# Initialize matrices for Spearman correlation and overlap counts
cor_results_matrix <- matrix(NA, nrow = 30, ncol = 30)
overlap_counts_matrix <- matrix(NA, nrow = 30, ncol = 30)

# Loop through all combinations of hbcc and mssm SVs (0 to 29)
for (i in 1:30) {
  for (j in 1:30) {
    # Merge the i-th hbcc dataset with the j-th mssm dataset by 'X'
    sv_combo <- merge(
      hbcc_list[[i]], 
      mssm_list[[j]], 
      by = 'X', 
      all.x = TRUE, 
      all.y = TRUE, 
      suffixes = c('_hbcc', '_mssm')
    )
    
    # Calculate Spearman correlation if logFC columns exist
    if ("logFC_hbcc" %in% names(sv_combo) && "logFC_mssm" %in% names(sv_combo)) {
      valid_data <- !is.na(sv_combo$logFC_hbcc) & !is.na(sv_combo$logFC_mssm)
      if (sum(valid_data) > 0) {
        cor_test <- cor.test(
          sv_combo$logFC_hbcc[valid_data], 
          sv_combo$logFC_mssm[valid_data], 
          method = "spearman", 
          exact = FALSE
        )
        cor_results_matrix[i, j] <- cor_test$estimate
      }
    }
    
    # Count overlapping significant genes
    if ("adj.P.Val_hbcc" %in% names(sv_combo) && "adj.P.Val_mssm" %in% names(sv_combo)) {
      significant_genes <- sv_combo$adj.P.Val_hbcc <= 0.05 & sv_combo$adj.P.Val_mssm <= 0.05
      overlap_counts_matrix[i, j] <- sum(significant_genes, na.rm = TRUE)
    }
  }
}

# Convert matrices into long-form data frames
correlation_df <- as.data.frame(as.table(cor_results_matrix))
names(correlation_df) <- c("hbcc_SV", "mssm_SV", "Correlation")

overlap_df <- as.data.frame(as.table(overlap_counts_matrix))
names(overlap_df) <- c("hbcc_SV", "mssm_SV", "OverlapCount")

# Map SV=0 to SV=29 for the new dimension
correlation_df$hbcc_SV <- as.numeric(correlation_df$hbcc_SV) - 1
correlation_df$mssm_SV <- as.numeric(correlation_df$mssm_SV) - 1
overlap_df$hbcc_SV <- as.numeric(overlap_df$hbcc_SV) - 1
overlap_df$mssm_SV <- as.numeric(overlap_df$mssm_SV) - 1

# Merge correlation and overlap data frames
heatmap_df <- merge(correlation_df, overlap_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

# Ensure all combinations from SV=0 to SV=29 are present
all_svs <- expand.grid(hbcc_SV = 0:29, mssm_SV = 0:29)
heatmap_df <- merge(all_svs, heatmap_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

# Calculate the total number of significant genes in hbcc
#hbcc_no_sva <- read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/BP/bp_hbcc_sva_be_0_SVs.csv")
hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp_11302024.txt",data.table=FALSE)
total_significant_hbcc <- sum(hbcc_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_hbcc #1303 DEGs

# Calculate the total number of significant genes in mssm
#mssm_no_sva <- read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/BP/bp_mssm_sva_be_0_SVs.csv")
mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_no_sva_bp_11302024.txt",data.table=FALSE)
total_significant_mssm <- sum(mssm_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_mssm #33 DEGs

# Function to calculate total significant genes for a given list of datasets
calculate_significant_genes <- function(data_list) {
  sapply(data_list, function(data) sum(data$adj.P.Val <= 0.05, na.rm = TRUE))
}

# Calculate the total number of significant genes for hbcc and mssm for each SV
hbcc_significant_genes <- calculate_significant_genes(hbcc_list)
mssm_significant_genes <- calculate_significant_genes(mssm_list)

# Create custom axis labels for hbcc and mssm
hbcc_labels <- paste0(0:29, " SV (", hbcc_significant_genes, " DEGs)")
mssm_labels <- paste0(0:29, " SV (", mssm_significant_genes, " DEGs)")

# Add labels to the heatmap data frame for plotting
heatmap_df$hbcc_Label <- factor(
  heatmap_df$hbcc_SV, 
  levels = 0:29, 
  labels = hbcc_labels
)
heatmap_df$mssm_Label <- factor(
  heatmap_df$mssm_SV, 
  levels = 0:29, 
  labels = mssm_labels
)

bp_full <- ggplot(heatmap_df, aes(x = hbcc_Label, y = mssm_Label, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(
    aes(label = ifelse(
      is.na(Correlation), "", 
      sprintf("%.2f\n(%s)", Correlation, OverlapCount)
    )), 
    color = "black", size = 4, lineheight = 0.9
  ) +
  scale_fill_gradient2(low = "white", mid = "#FFFFC5", high = "#e0413f", midpoint = 0.1,
                      breaks = c(0.1, 0.2, 0.3),) +
  #scale_fill_gradient2(low = "white", high = "#79AC78") +
  # scale_x_continuous(breaks = 0:29, expand = c(0, 0)) +
  # scale_y_continuous(breaks = 0:29, expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal(base_size = 16) +
  theme(
    #plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
    #panel.spacing = unit(0, "lines")
  ) +
  labs(
    #title = "Spearman Correlation Heatmap: CMC vs CMC-HBCC in Bipolar Disorder",
    x = "BP: CMC-HBCC",
    y = "BP: CMC",
    fill = bquote("Correlation "*rho))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_fewSVs_BP_full_20250623_test.pdf", plot = bp_full, width = 17, height = 17)

# ggsave("/hpc/users/hoangd02/www/plots/bp_mssm_hbcc_corr_0_to_29_heatmap_with_overlapping_genes_and_corr_values.pdf", plot=bp_full , width = 15, height = 15)

heatmap_df_subset1 <- subset(heatmap_df, hbcc_SV <= 10 & mssm_SV <= 10)

# ps1 <- ggplot(heatmap_df_subset1, aes(x = hbcc_Label, y = mssm_Label, fill = Correlation)) +
#   geom_tile(color = "white", size = 0.1) +
#   geom_text(
#     aes(label = ifelse(
#       is.na(Correlation), "", 
#       sprintf("%.2f\n(%s)", Correlation, OverlapCount)
#     )), 
#     color = "black", size = 3, lineheight = 0.9  # reduced from 4.5
#   ) +
#   scale_fill_gradient2(low = "#a0f600", high = "#d633a6") + #this is equivalent to most pink and low = white
#   #scale_fill_gradient2(mid = "#d633a6", high = "#a0f600") +
#   #scale_fill_gradient2(mid = "#a0f600", high = "#d633a6") +
#   coord_fixed() +
#   theme_minimal(base_size = 10) +  # reduced from 14
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 8),
#     axis.title.x = element_text(size = 10, face = "bold"),
#     axis.title.y = element_text(size = 10, face = "bold"),
#     legend.title = element_text(size = 9),
#     legend.text = element_text(size = 8),
#     plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
#   ) +
#   labs(
#     title = "Spearman Correlation Heatmap for BP (SVs = 0 to 10)",
#     x = "BP: CMC-HBCC Surrogate Variables",
#     y = "BP: CMC Surrogate Variables"
#   )

# ggsave("/hpc/users/hoangd02/www/plots/bp_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_test_20250519.pdf", plot = ps1, width = 8, height = 8)

# ggsave("/hpc/users/hoangd02/www/plots/bp_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values.pdf", plot = ps1, width = 8, height = 8)

#################### SCHIZOPHRENIA #########################
############################################################
output_dir <- "/sc/arion/projects/mscic1/results/jolie/CMC/"

hbcc_files <- paste0(output_dir, "sz_hbcc_sva_be_", 0:31, "_SVs.csv") #limiting factor with 31 SVs
mssm_files <- paste0(output_dir, "sz_mssm_sva_be_", 0:31, "_SVs.csv")

# Load sva datasets
hbcc_list <- lapply(hbcc_files, read.csv)
mssm_list <- lapply(mssm_files, read.csv)

# Convert sva datasets to data frames for consistency
hbcc_list <- lapply(hbcc_list, as.data.frame)
mssm_list <- lapply(mssm_list, as.data.frame)

# Confirm the lists have 32 elements each (nsv=0 to nsv=31)
cat("hbcc list length:", length(hbcc_list), "\n")
cat("mssm list length:", length(mssm_list), "\n")

##genes are not aligned, check: 
head(hbcc_list[[2]])
head(mssm_list[[29]])

tmp=lapply(mssm_list,function(x){x$adj.P.Val<=0.05})
sapply(tmp,table)

##align the genes - hbcc
# Step 1: Identify the common set of genes across all datasets
common_genes <- Reduce(intersect, lapply(hbcc_list, function(df) df$X))

# Step 2: Subset each dataset to include only the common genes
hbcc_list <- lapply(hbcc_list, function(df) df[df$X %in% common_genes, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
hbcc_list <- lapply(hbcc_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check <- all(sapply(hbcc_list, function(df) identical(df$X, hbcc_list[[1]]$X)))

if (alignment_check) {
  cat("Genes are successfully aligned across all datasets.\n")
} else {
  cat("Gene alignment failed. Check for issues in the datasets.\n")
}

##align the genes - mssm
# Step 1: Identify the common set of genes across all datasets in mssm_list
common_genes_mssm <- Reduce(intersect, lapply(mssm_list, function(df) df$X))

# Step 2: Subset each dataset in mssm_list to include only the common genes
mssm_list <- lapply(mssm_list, function(df) df[df$X %in% common_genes_mssm, ])

# Step 3: Sort each dataset by the 'X' column to ensure consistent order
mssm_list <- lapply(mssm_list, function(df) df[order(df$X), ])

# Step 4: Verify alignment
alignment_check_mssm <- all(sapply(mssm_list, function(df) identical(df$X, mssm_list[[1]]$X)))

if (alignment_check_mssm) {
  cat("Genes are successfully aligned across all mssm datasets.\n")
} else {
  cat("Gene alignment failed for mssm datasets. Check for issues in the datasets.\n")
}

# Initialize a matrix to store Spearman correlation results 
cor_results_matrix <- matrix(NA, nrow = 32, ncol = 32)

# Loop through all combinations of hbcc and mssm SVs (0 to 31)
for (i in 1:32) {
  for (j in 1:32) {
    # Merge the i-th hbcc dataset with the j-th mssm dataset by 'X'
    sv_combo <- merge(
      hbcc_list[[i]], 
      mssm_list[[j]], 
      by = 'X', 
      all.x = TRUE, 
      all.y = TRUE, 
      suffixes = c('_hbcc', '_mssm')
    )
    
    # Check if logFC columns exist and have non-NA values
    if ("logFC_hbcc" %in% names(sv_combo) && "logFC_mssm" %in% names(sv_combo)) {
      valid_data <- !is.na(sv_combo$logFC_hbcc) & !is.na(sv_combo$logFC_mssm)
      if (sum(valid_data) > 0) {
        cor_test <- cor.test(
          sv_combo$logFC_hbcc[valid_data], 
          sv_combo$logFC_mssm[valid_data], 
          method = "spearman", 
          exact = FALSE
        )
        cor_results_matrix[i, j] <- cor_test$estimate
      }
    }
  }
}

# Convert the matrix into a long-form data frame
correlation_df <- as.data.frame(as.table(cor_results_matrix))
names(correlation_df) <- c("hbcc_SV", "mssm_SV", "Correlation")

# Map SV=0 to SV=31 for the new dimension
correlation_df$hbcc_SV <- as.numeric(correlation_df$hbcc_SV) - 1
correlation_df$mssm_SV <- as.numeric(correlation_df$mssm_SV) - 1

# Ensure all combinations from SV=0 to SV=31 are present
all_svs <- expand.grid(hbcc_SV = 0:31, mssm_SV = 0:31)
correlation_df <- merge(all_svs, correlation_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

##ADDING THE COUNT OF OVERLAPPING GENES ON THE TILES INSTEAD OF THE SPEARMAN RHO
# Initialize matrices for Spearman correlation and overlap counts
cor_results_matrix <- matrix(NA, nrow = 32, ncol = 32)
overlap_counts_matrix <- matrix(NA, nrow = 32, ncol = 32)

# Loop through all combinations of hbcc and mssm SVs (0 to 31)
for (i in 1:32) {
  for (j in 1:32) {
    # Merge the i-th hbcc dataset with the j-th mssm dataset by 'X'
    sv_combo <- merge(
      hbcc_list[[i]], 
      mssm_list[[j]], 
      by = 'X', 
      all.x = TRUE, 
      all.y = TRUE, 
      suffixes = c('_hbcc', '_mssm')
    )
    
    # Calculate Spearman correlation if logFC columns exist
    if ("logFC_hbcc" %in% names(sv_combo) && "logFC_mssm" %in% names(sv_combo)) {
      valid_data <- !is.na(sv_combo$logFC_hbcc) & !is.na(sv_combo$logFC_mssm)
      if (sum(valid_data) > 0) {
        cor_test <- cor.test(
          sv_combo$logFC_hbcc[valid_data], 
          sv_combo$logFC_mssm[valid_data], 
          method = "spearman", 
          exact = FALSE
        )
        cor_results_matrix[i, j] <- cor_test$estimate
      }
    }
    
    # Count overlapping significant genes
    if ("adj.P.Val_hbcc" %in% names(sv_combo) && "adj.P.Val_mssm" %in% names(sv_combo)) {
      significant_genes <- sv_combo$adj.P.Val_hbcc <= 0.05 & sv_combo$adj.P.Val_mssm <= 0.05
      overlap_counts_matrix[i, j] <- sum(significant_genes, na.rm = TRUE)
    }
  }
}

# Convert matrices into long-form data frames
correlation_df <- as.data.frame(as.table(cor_results_matrix))
names(correlation_df) <- c("hbcc_SV", "mssm_SV", "Correlation")

overlap_df <- as.data.frame(as.table(overlap_counts_matrix))
names(overlap_df) <- c("hbcc_SV", "mssm_SV", "OverlapCount")

# Map SV=0 to SV=31 for the new dimension
correlation_df$hbcc_SV <- as.numeric(correlation_df$hbcc_SV) - 1
correlation_df$mssm_SV <- as.numeric(correlation_df$mssm_SV) - 1
overlap_df$hbcc_SV <- as.numeric(overlap_df$hbcc_SV) - 1
overlap_df$mssm_SV <- as.numeric(overlap_df$mssm_SV) - 1

# Merge correlation and overlap data frames
heatmap_df <- merge(correlation_df, overlap_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

# Ensure all combinations from SV=0 to SV=31 are present
all_svs <- expand.grid(hbcc_SV = 0:31, mssm_SV = 0:31)
heatmap_df <- merge(all_svs, heatmap_df, by = c("hbcc_SV", "mssm_SV"), all.x = TRUE)

# Calculate the total number of significant genes in hbcc
sz_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_no_sva.txt",data.table=FALSE)#
total_significant_hbcc <- sum(sz_hbcc_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_hbcc # 2065 DEGs

# Calculate the total number of significant genes in mssm
sz_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_no_sva.txt",data.table=FALSE)
total_significant_mssm <- sum(sz_mssm_no_sva$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_mssm # 233 DEGs

# Function to calculate total significant genes for a given list of datasets
calculate_significant_genes <- function(data_list) {
  sapply(data_list, function(data) sum(data$adj.P.Val <= 0.05, na.rm = TRUE))
}

# Calculate the total number of significant genes for hbcc and mssm for each SV
hbcc_significant_genes <- calculate_significant_genes(hbcc_list)
mssm_significant_genes <- calculate_significant_genes(mssm_list)

# Create custom axis labels for hbcc and mssm
hbcc_labels <- paste0(0:31, " SV (", hbcc_significant_genes, " DEGs)")
mssm_labels <- paste0(0:31, " SV (", mssm_significant_genes, " DEGs)")

# Add labels to the heatmap data frame for plotting
heatmap_df$hbcc_Label <- factor(
  heatmap_df$hbcc_SV, 
  levels = 0:31, 
  labels = hbcc_labels
)
heatmap_df$mssm_Label <- factor(
  heatmap_df$mssm_SV, 
  levels = 0:31, 
  labels = mssm_labels
)

#############################################################################
## new code with correlation values and DEGs
# Assume OverlapCount is already in correlation_df
# If not, you need to add a column named OverlapCount before plotting
sz_full <- ggplot(heatmap_df, aes(x = hbcc_Label, y = mssm_Label, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(
    aes(label = ifelse(
      is.na(Correlation), "", 
      sprintf("%.2f\n(%s)", Correlation, OverlapCount)
    )), 
    color = "black", size = 4, lineheight = 0.9
  ) +
  scale_fill_gradient2(mid = "#FFFFC5", high = "#e0413f",midpoint=0.1,
                        breaks = c(0.1, 0.2, 0.3)) +
  # scale_x_continuous(breaks = 0:29, expand = c(0, 0)) +
  # scale_y_continuous(breaks = 0:29, expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal(base_size = 16) +
  theme(
    #plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
    #panel.spacing = unit(0, "lines")
  ) +
  labs(
    #title = "Spearman Correlation Heatmap: CMC vs CMC-HBCC in Schizophrenia",
    x = "SZ: CMC-HBCC",
    y = "SZ: CMC",
    fill = bquote("Correlation "*rho))

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_fewSVs_SZ_full_20250623.pdf", plot = sz_full, width = 17, height = 17)

### FINAL
pdf("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/heatmap_ad_bp_sz.pdf", width = 17, height = 17)

print(wrap_plots(ad_full, ncol = 1))
print(wrap_plots(bp_full, ncol = 1))
print(wrap_plots(sz_full, ncol = 1))

dev.off()

# ggsave("/hpc/users/hoangd02/www/plots/sz_mssm_hbcc_corr_0_to_31_heatmap_with_overlapping_genes_and_corr_values.pdf", plot=plot, width = 15, height = 15)

heatmap_df_subset1 <- subset(heatmap_df, hbcc_SV <= 10 & mssm_SV <= 10)

# ps1 <- ggplot(heatmap_df_subset1, aes(x = hbcc_Label, y = mssm_Label, fill = Correlation)) +
#   geom_tile(color = "white", size = 0.1) +
#   geom_text(
#     aes(label = ifelse(
#       is.na(Correlation), "", 
#       sprintf("%.2f\n(%s)", Correlation, OverlapCount)
#     )), 
#     color = "black", size = 3, lineheight = 0.9  # reduced from 4.5
#   ) +
#   scale_fill_gradient2(low = "#a0f600", high = "#d633a6") + #this is equivalent to most pink and low = white
#   #scale_fill_gradient2(low = "white", high = "#734F96") +
#   coord_fixed() +
#   theme_minimal(base_size = 10) +  # reduced from 14
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 8),
#     axis.title.x = element_text(size = 10, face = "bold"),
#     axis.title.y = element_text(size = 10, face = "bold"),
#     legend.title = element_text(size = 9),
#     legend.text = element_text(size = 8),
#     plot.margin = margin(t = 10, r = 10, b = 10, l = 20)
#   ) +
#   labs(
#     title = "Spearman Correlation Heatmap for SZ (SVs = 0 to 10)",
#     x = "SZ: CMC-HBCC Surrogate Variables",
#     y = "SZ: CMC Surrogate Variables"
#   )

# ggsave("/hpc/users/hoangd02/www/plots/sz_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_20250519.pdf", plot = ps1, width = 8, height = 8)

# ggsave("/hpc/users/hoangd02/www/plots/sz_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values.pdf", plot = ps1, width = 8, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/rosmap_msbb_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_20250519.pdf

https://hoangd02.u.hpc.mssm.edu/plots/sz_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_20250519.pdf

https://hoangd02.u.hpc.mssm.edu/plots/bp_mssm_hbcc_corr_0_to_10_heatmap_with_overlapping_genes_and_cor_values_test_20250519.pdf








































