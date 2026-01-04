###### ###### ###### This figure real data logFC 
###################################################################################
###### from SVA_all_code.sh ##########
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
library(tidyr)
library(ggrastr)
library(ggrepel)
library(patchwork)  

####################### AD #################################
##################### ROSMAP ###############################
rosmap_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt", data.table = FALSE) #19476     7
rosmap_be <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_1.16.txt", data.table = FALSE) #19476     7
rosmap_leek <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_leek_1.16.txt", data.table = FALSE) #19476     7

names(rosmap_no_sva)[names(rosmap_no_sva) == "logFC_j"] <- "logFC"

####################### MSBB ############################## 
msbb_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_no_sva_1.11.txt", data.table = FALSE) #26618 x 7
msbb_be <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_be_1.11.txt", data.table = FALSE) #26618 x 7
msbb_leek <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_leek_1.11.txt", data.table = FALSE) #26618 x 7

names(msbb_no_sva)[names(msbb_no_sva) == "logFC_j"] <- "logFC"

############################################################
ad_no_sva_combo <-merge(msbb_no_sva,rosmap_no_sva,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_msbb','_rosmap'))
#ad_cor <-cor.test(ad_no_sva_combo$logFC_j_msbb, ad_no_sva_combo$logFC_j_rosmap, method = "spearman") #0.4603241 
ad_cor <-cor.test(ad_no_sva_combo$logFC_msbb, ad_no_sva_combo$logFC_rosmap, method = "spearman",alternative="greater"); ad_cor #0.4603241 

ad_sva_be_combo <- merge(msbb_be,rosmap_be,by='X', all.x=TRUE, all.y=TRUE,suffixes=c('_msbb','_rosmap'))
ad_cor_sva <- cor.test(ad_sva_be_combo$logFC_msbb, ad_sva_be_combo$logFC_rosmap, method = "spearman",alternative="greater"); ad_cor_sva #0.1157661 

# Data
ad_cor_values <- data.frame(
  Method = c("No SVA", "SVA 'BE'"),
  Spearman = c(ad_cor$estimate, ad_cor_sva$estimate),
  Disease = "AD"
)

# Factor to control shape and order
ad_cor_values$Method <- factor(ad_cor_values$Method, levels = c("No SVA", "SVA 'BE'"))

# Plot
ad <- ggplot(ad_cor_values, aes(x = Spearman, y = Disease, shape = Method)) +
  geom_line(aes(group = Disease), linewidth = 1, color = "gray50") +  # connecting line
  geom_point(size = 4) +
  scale_shape_manual(values = c(16, 17)) +  # circle for No SVA, triangle for SVA 'BE'
  xlim(0, 0.5) +
  theme_minimal(base_size = 14) +
  labs(x = "Spearman Correlation", y = NULL, shape = "Adjustment") +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )
ad

ggsave(filename = "/hpc/users/hoangd02/www/plots/real_data_test.pdf",
  plot = ad, width = 8, height = 8) 

####################### SZ #################################
####################### CMC ###############################
sz_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_no_sva.txt",data.table=FALSE) #19086     7
sz_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_sva_be.txt",data.table=FALSE)
sz_mssm_no_sva$X <- sub("\\..*", "", sz_mssm_no_sva$X)
sz_mssm_sva_be$X <- sub("\\..*", "", sz_mssm_sva_be$X)

##################### CMC-HBCC ###############################
sz_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_no_sva.txt",data.table=FALSE) #19086     7
sz_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_sva_be.txt",data.table=FALSE)
sz_hbcc_no_sva$X <- sub("\\..*", "", sz_hbcc_no_sva$X)
sz_hbcc_sva_be$X <- sub("\\..*", "", sz_hbcc_sva_be$X)

####################### BP #################################
####################### CMC ###############################
bp_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp_11302024.txt",
				 data.table=FALSE) #19017     7
bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp.txt",data.table=FALSE)
bp_hbcc_no_sva$X <- sub("\\..*", "", bp_hbcc_no_sva$X)
bp_hbcc_sva_be$X <- sub("\\..*", "", bp_hbcc_sva_be$X)

##################### CMC-HBCC ###############################s
bp_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_no_sva_bp_11302024.txt",
			     data.table=FALSE) #19017     7
bp_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_sva_be_bp.txt",data.table=FALSE)
bp_mssm_no_sva$X <- sub("\\..*", "", bp_mssm_no_sva$X)
bp_mssm_sva_be$X <- sub("\\..*", "", bp_mssm_sva_be$X)


####################### ASD #################################
####################### UCLA ###############################
ucla<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_without_sva_09232024.txt",
	  data.table=FALSE)  #20850 x 7
ucla_sva<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_with_sva_09232024.txt",data.table=FALSE) 

####################### YALE ###############################
##from ASD code for publication 09242024.sh
yale<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_without_sva_09232024.txt",
	  data.table=FALSE) #18255 x 7
yale_sva<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_with_sva_09232024.txt",data.table=FALSE)


# Function to merge and calculate Spearman correlation
get_spearman_logFC <- function(df1, df2, suffix1, suffix2, sva = FALSE) {
  combo <- merge(df1, df2, by = "X", all = TRUE, suffixes = c(paste0("_", suffix1), paste0("_", suffix2)))
  
  # Pick correct column names
  col1 <- paste0("logFC_", suffix1)
  col2 <- paste0("logFC_", suffix2)
  
  # If columns are named like logFC_j_msbb instead, adjust here
  if (!col1 %in% names(combo)) col1 <- grep(paste0("logFC.*", suffix1), names(combo), value = TRUE)[1]
  if (!col2 %in% names(combo)) col2 <- grep(paste0("logFC.*", suffix2), names(combo), value = TRUE)[1]
  
  cor_result <- cor.test(combo[[col1]], combo[[col2]], method = "spearman", alternative = "greater")
  return(c(cor = cor_result$estimate, p.value = cor_result$p.value))
}

ad_no_sva_result <- get_spearman_logFC(msbb_no_sva,rosmap_no_sva, "msbb", "rosmap")
ad_sva_result    <- get_spearman_logFC(msbb_be,rosmap_be, "msbb", "rosmap")

bp_no_sva_result <- get_spearman_logFC(bp_mssm_no_sva, bp_hbcc_no_sva, "mssm", "hbcc")
bp_sva_result    <- get_spearman_logFC(bp_mssm_sva_be, bp_hbcc_sva_be, "mssm", "hbcc")

sz_no_sva_result <- get_spearman_logFC(sz_mssm_no_sva, sz_hbcc_no_sva, "mssm", "hbcc")
sz_sva_result    <- get_spearman_logFC(sz_mssm_sva_be, sz_hbcc_sva_be, "mssm", "hbcc")

asd_no_sva_result <- get_spearman_logFC(ucla, yale, "ucla", "yale")
asd_sva_result    <- get_spearman_logFC(ucla_sva, yale_sva, "ucla", "yale")

cor_values <- data.frame(
  Method = rep(c("No SVA", "SVA 'BE'"), times = 4),
  Spearman = c(
    ad_no_sva_result["cor.rho"], ad_sva_result["cor.rho"],
    asd_no_sva_result["cor.rho"], asd_sva_result["cor.rho"],
    bp_no_sva_result["cor.rho"], bp_sva_result["cor.rho"],
    sz_no_sva_result["cor.rho"], sz_sva_result["cor.rho"]
  ),
  Disease = rep(c("AD","ASD", "BP", "SZ"), each = 2)
)

cor_values$Method <- factor(cor_values$Method, levels = c("No SVA", "SVA 'BE'"))
cor_values$Disease <- factor(cor_values$Disease, levels = c("SZ", "BP", "ASD","AD"))

## LINE PLOT
all <- ggplot(cor_values, aes(x = Spearman, y = Disease, color = Method)) +
  geom_line(aes(group = Disease), linewidth = 1, color = "gray50") +
  geom_point(size = 4, shape = 16) +
  scale_x_continuous(
    limits = c(-0.1, 0.5),
    breaks = seq(-0.1, 0.5, by = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Spearman Correlation of the log"[2]*"FC between Datasets"),
    y = NULL, color = "SVA Approaches"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()  # ← removes lines *between* 0.1 ticks
  )

## BAR PLOT IDEA
bar_plot_all <- ggplot(cor_values, aes(x = Method, y = Spearman, fill = Method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Disease, nrow = 2) +
  scale_y_continuous(
    limits = c(-0.1, 0.5),
    breaks = seq(-0.1, 0.5, by = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    y = expression("Spearman Correlation of the log"[2]*"FC between Datasets"),
    x = "SVA Approach"
  ) +
  theme(panel.grid.minor.y = element_blank())

############################ ###### PD ################# ################# ################# 
############################# THIS IS FROM SVA_all_code.sh #############################

no_sva_spearman_cor_matrix <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/no_sva_spearman_cor_matrix_final_final.txt", header = TRUE, sep = "\t")
sva_spearman_cor_matrix <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/sva_spearman_cor_matrix_final_final.txt", header = TRUE, sep = "\t")

rownames(no_sva_spearman_cor_matrix) <- no_sva_spearman_cor_matrix$X
rownames(sva_spearman_cor_matrix) <- sva_spearman_cor_matrix$X

no_sva_spearman_cor_matrix$X <- NULL
sva_spearman_cor_matrix$X <- NULL

rownames(no_sva_spearman_cor_matrix) <- gsub("_none$", "", rownames(no_sva_spearman_cor_matrix))
rownames(sva_spearman_cor_matrix) <- gsub("_sva$", "", rownames(sva_spearman_cor_matrix))

#make sure they have the same col name
colnames(no_sva_spearman_cor_matrix) <- sub("_none$", "", colnames(no_sva_spearman_cor_matrix))
colnames(sva_spearman_cor_matrix) <- sub("_sva$", "", colnames(sva_spearman_cor_matrix))

# Convert correlation matrix to long format for ggplot
sva_cor_matrix_long <- melt(as.matrix(sva_spearman_cor_matrix))
no_sva_cor_matrix_long <- melt(as.matrix(no_sva_spearman_cor_matrix))

##repeat the same thing for no_sva_spearman_cor_matrix AND sva_spearman_cor_matrix
names(sva_cor_matrix_long) <- c("Dataset1", "Dataset2", "value")
names(no_sva_cor_matrix_long) <- c("Dataset1", "Dataset2", "value")

order<-c("GSE20292", "GSE8397", "GSE20164", "GSE20163","GSE24378", "GSE7621", "GSE20141", "GSE49036")
sva_cor_matrix_long$Dataset1 <- factor(sva_cor_matrix_long$Dataset1, levels = order)
sva_cor_matrix_long$Dataset2 <- factor(sva_cor_matrix_long$Dataset2, levels = order)

no_sva_cor_matrix_long$Dataset1 <- factor(no_sva_cor_matrix_long$Dataset1, levels = order)
no_sva_cor_matrix_long$Dataset2 <- factor(no_sva_cor_matrix_long$Dataset2, levels = order)

# Generate heatmap
# heatmap_plot <- ggplot(no_sva_cor_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = round(value, 3)), vjust = 1) +
#   scale_fill_gradient2(low = "blue", high = "red", mid= "white",
#                        limit = c(-1, 1), space = "Lab", 
#                        name = "Spearman\nCorrelation") +
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
#         axis.text.y = element_text(size = 12), 
#         plot.title = element_text(hjust = 0.5, size = 20),
#         legend.title = element_text(size = 14),  
#         legend.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 16),  
#         axis.title.y = element_text(size = 16)) + 
#   coord_fixed() +
#   labs(title = "Dataset Replicability without SVA") + 
#   xlab("Datasets") + ylab("Datasets"); heatmap_plot 

##ggsave("/hpc/users/hoangd02/www/plots/GEO/no_sva_spearman_cor_matrix_heatmap_final.png", plot = heatmap_plot, width = 10, height = 10)

##MAKING THE DIFFERENCE MATRIX
rownames(no_sva_spearman_cor_matrix) <- no_sva_spearman_cor_matrix$X
rownames(sva_spearman_cor_matrix) <- sva_spearman_cor_matrix$X

no_sva_spearman_cor_matrix$X <- NULL
sva_spearman_cor_matrix$X <- NULL

rownames(no_sva_spearman_cor_matrix) <- gsub("_none$", "", rownames(no_sva_spearman_cor_matrix))
rownames(sva_spearman_cor_matrix) <- gsub("_sva$", "", rownames(sva_spearman_cor_matrix))

#make sure they have the same col name
colnames(no_sva_spearman_cor_matrix) <- sub("_none$", "", colnames(no_sva_spearman_cor_matrix))
colnames(sva_spearman_cor_matrix) <- sub("_sva$", "", colnames(sva_spearman_cor_matrix))

common_columns <- intersect(colnames(sva_spearman_cor_matrix), colnames(no_sva_spearman_cor_matrix))
sva_spearman_cor_matrix <- sva_spearman_cor_matrix[, common_columns]
no_sva_spearman_cor_matrix <- no_sva_spearman_cor_matrix[, common_columns]

sva_matrix <- as.matrix(sva_spearman_cor_matrix)
no_sva_matrix <- as.matrix(no_sva_spearman_cor_matrix)

# Create the difference matrix
library(seriation)
difference_matrix <- no_sva_matrix - sva_matrix
rowdist <-dist(difference_matrix)
coldist <-dist(t(difference_matrix))
roworder<-seriate(rowdist)
colorder<-seriate(coldist)

#PLOT
difference_matrix_long <- melt(difference_matrix)

# Rename columns to match the expected format
names(difference_matrix_long) <- c("Dataset1", "Dataset2", "value")

# Ensure the levels of Dataset1 and Dataset2 match the order of dataset names
dataset_order <- rownames(difference_matrix)
difference_matrix_long$Dataset1 <- factor(difference_matrix_long$Dataset1, levels = rownames(difference_matrix)[unlist(roworder)])
difference_matrix_long$Dataset2 <- factor(difference_matrix_long$Dataset2, levels = colnames(difference_matrix)[unlist(colorder)])

difference_matrix_long$Dataset1 <- factor(
  difference_matrix_long$Dataset1,
  levels = c("4", "1", "7", "6", "3", "2", "5", "8"),
  labels = c("GSE20292", "GSE8397", "GSE20164", "GSE20163", 
             "GSE24378", "GSE7621", "GSE20141", "GSE49036"))

cor_values$Disease <- factor(cor_values$Disease, levels = c("AD", "ASD", "BP", "SZ"))
cor_values$Method <- factor(cor_values$Method, levels = c("No SVA", "SVA 'BE'"))

## BAR PLOT IDEA
# Plot
bar_plot_all <- ggplot(cor_values, aes(x = Method, y = Spearman, fill = Method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(
      label = sprintf("%.2f", Spearman),
      vjust = ifelse(Spearman >= 0, -0.4, 1.2)
    ),
    position = position_dodge(width = 0.9),
    size = 4
  ) +
  facet_wrap(~ Disease, nrow = 1) +
  scale_y_continuous(
    limits = c(-0.15, 0.55),  # lower limit adjusted for negative labels
    breaks = seq(-0.1, 0.5, by = 0.1)
  ) +
  scale_fill_manual(
    name = "SVA Approaches",
    values = c("No SVA" = "navy", "SVA 'BE'" = "lightblue")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    y = expression("Correlation of the log"[2]*"FC between Datasets"),
    x = "SVA Approach"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )
#bar_plot_all

# Generate the heatmap
heatmap_plot <- ggplot(difference_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = 1) +
  scale_fill_gradient2(
    high = "red", mid = "white", low = "blue",
    limit = c(-0.2, 0.2), space = "Lab",
    name = expression("Correlation " * rho)
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 25, vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  ) +  
  coord_fixed() +
  labs(x = "PD Datasets", y = "PD Datasets")
  #heatmap_plot

#ggsave("/hpc/users/hoangd02/www/plots/GEO/difference_no_sva_spearman_cor_matrix_heatmap_final.pdf", plot = heatmap_plot, width = 10, height = 10)

combo <- bar_plot_all + heatmap_plot + plot_layout(widths = c(2, 2))  # bar_plot_all is 1/3 the width
ggsave("/hpc/users/hoangd02/www/plots/figure5_deg_reproducibility.pdf", plot = combo, width = 12, height = 6)

####################### ARCHIVED ###############################
############### CROSS-DISEASE COMPARISONS ########################
# List of datasets (no SVA and SVA)
all_datasets <- list(
  `AD:ROSMAP_no_sva` = rosmap_no_sva,
  `AD:MSBB_no_sva` = msbb_no_sva,
  `BP:MSSM_no_sva` = bp_mssm_no_sva,
  `BP:HBCC_no_sva` = bp_hbcc_no_sva,
  `SZ:MSSM_no_sva` = sz_mssm_no_sva,
  `SZ:HBCC_no_sva` = sz_hbcc_no_sva,

  `AD:ROSMAP_sva` = rosmap_be,
  `AD:MSBB_sva` = msbb_be,
  `BP:MSSM_sva` = bp_mssm_sva_be,
  `BP:HBCC_sva` = bp_hbcc_sva_be,
  `SZ:MSSM_sva` = sz_mssm_sva_be,
  `SZ:HBCC_sva` = sz_hbcc_sva_be
)

# Base names without SVA suffix
dataset_base_names <- names(all_datasets)[grepl("_no_sva$", names(all_datasets))]
dataset_base_names <- gsub("_no_sva$", "", dataset_base_names)

# All unique pairs (including cross-disease)
dataset_pairs <- combn(dataset_base_names, 2, simplify = FALSE)

get_cross_spearman <- function(base1, base2) {
  df1_no <- all_datasets[[paste0(base1, "_no_sva")]]
  df2_no <- all_datasets[[paste0(base2, "_no_sva")]]
  df1_sva <- all_datasets[[paste0(base1, "_sva")]]
  df2_sva <- all_datasets[[paste0(base2, "_sva")]]

  # Subset to X + logFC and merge
  df1_no <- df1_no[, c("X", "logFC")]
  df2_no <- df2_no[, c("X", "logFC")]
  df1_sva <- df1_sva[, c("X", "logFC")]
  df2_sva <- df2_sva[, c("X", "logFC")]
  names(df1_no)[2] <- names(df1_sva)[2] <- "logFC_1"
  names(df2_no)[2] <- names(df2_sva)[2] <- "logFC_2"

  merge_no <- merge(df1_no, df2_no, by = "X")
  merge_sva <- merge(df1_sva, df2_sva, by = "X")

  spearman_no <- if (nrow(merge_no) > 0) cor(merge_no$logFC_1, merge_no$logFC_2, method = "spearman") else NA
  spearman_sva <- if (nrow(merge_sva) > 0) cor(merge_sva$logFC_1, merge_sva$logFC_2, method = "spearman") else NA

  data.frame(
    Dataset_Pair = rep(paste(base1, "vs", base2), 2),
    Method = c("No SVA", "SVA 'BE'"),
    Spearman = c(spearman_no, spearman_sva),
    Disease = "Cross"
  )
}


# Apply to all pairs
cross_cor_df <- do.call(rbind, lapply(dataset_pairs, function(p) get_cross_spearman(p[1], p[2])))

all_cor_df <- rbind(cross_cor_df, cor_values)  # cor_values = your previous AD/BP/SZ summary rows

a <- ggplot(cross_cor_df, aes(x = Spearman, y = Dataset_Pair, color = Method)) +
  geom_line(aes(group = Dataset_Pair), linewidth = 0.7, alpha = 0.6) +
  geom_point(size = 3, shape = 16) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Spearman Correlation of the log"[2]*"(FC) between Datasets"),
    y = NULL,
    color = "Adjustment"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave(filename = "/hpc/users/hoangd02/www/plots/real_data_test3.pdf",
  plot = a, width = 8, height = 8) 

##ordering
get_pair_type <- function(pair) {
  if (grepl("^AD:.* vs AD:.*", pair)) return("A_AD-AD")
  if (grepl("^BP:.* vs BP:.*", pair)) return("B_BP-BP")
  if (grepl("^SZ:.* vs SZ:.*", pair)) return("C_SZ-SZ")
  if (grepl("^AD:.* vs BP:.*|^BP:.* vs AD:.*", pair)) return("D_AD-BP")
  if (grepl("^AD:.* vs SZ:.*|^SZ:.* vs AD:.*", pair)) return("E_AD-SZ")
  if (grepl("^BP:.* vs SZ:.*|^SZ:.* vs BP:.*", pair)) return("F_BP-SZ")
  return("G_Other")
}

# Assign and sort
pair_order_df <- data.frame(Dataset_Pair = unique(cross_cor_df$Dataset_Pair))
pair_order_df$group <- sapply(pair_order_df$Dataset_Pair, get_pair_type)

# Apply the new factor order
pair_order_df <- pair_order_df %>%
  arrange(group, Dataset_Pair)

cross_cor_df$Dataset_Pair <- factor(
  cross_cor_df$Dataset_Pair,
  levels = rev(pair_order_df$Dataset_Pair)
)

# Plot again
a <- ggplot(cross_cor_df, aes(x = Spearman, y = Dataset_Pair, color = Method)) +
  geom_line(aes(group = Dataset_Pair), linewidth = 0.7, alpha = 0.6) +
  geom_point(size = 3, shape = 16) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Spearman Correlation of the log"[2]*"(FC) between Datasets"),
    y = NULL,
    color = "Adjustment"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave(filename = "/hpc/users/hoangd02/www/plots/real_data_test3.pdf",
       plot = a, width = 8, height = 8)


####################### PD #################################
##PD datasets: 28 unique pairs (8 choose 2 = 8x7/2 = 28)
# You're comparing dataset 1 vs 2, 1 vs 3, ..., 1 vs 8
# Then 2 vs 3, 2 vs 4, ..., 2 vs 8
# And so on...

GSE8397_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #22283     7
GSE8397_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE7621_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #54318     7
GSE7621_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE24378_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #61359     7
GSE24378_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20292_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #22283     7
GSE20292_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t") ##no X, just probeID
GSE20292_sva$X <- sub("_.*", "", GSE20292_sva$probeID)
GSE20292_sva$probeID <- NULL
GSE20292_sva <- GSE20292_sva[, c("X", setdiff(names(GSE20292_sva), "X"))]

GSE20141_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #54675     7
GSE20141_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")
  
GSE20163_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #22283     7
GSE20163_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")
  
GSE20164_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #22283     7
GSE20164_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")
  
GSE49036_none = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t") #54675     7
GSE49036_sva = read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t") 

pd_datasets <- list(
  GSE8397_none = GSE8397_none,
  GSE7621_none = GSE7621_none,
  GSE24378_none = GSE24378_none,
  GSE20292_none = GSE20292_none,
  GSE20141_none = GSE20141_none,
  GSE20163_none = GSE20163_none,
  GSE20164_none = GSE20164_none,
  GSE49036_none = GSE49036_none,
  
  GSE8397_sva = GSE8397_sva,
  GSE7621_sva = GSE7621_sva,
  GSE24378_sva = GSE24378_sva,
  GSE20292_sva = GSE20292_sva,
  GSE20141_sva = GSE20141_sva,
  GSE20163_sva = GSE20163_sva,
  GSE20164_sva = GSE20164_sva,
  GSE49036_sva = GSE49036_sva
)

clean_pair_name <- function(name) {
  gsub("_none|_sva", "", name)
}

pd_no_sva_combo <- merge(GSE8397_none,GSE49036_none,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_GSE8397','_GSE49036'))
cor_result <- cor.test(pd_no_sva_combo$logFC_GSE8397, pd_no_sva_combo$logFC_GSE49036, method = "spearman", alternative="greater")
cor_result
#-0.01245489

#pd_no_sva_combo will contain all unique X values from both datasets, and rows that don’t exist in both datasets will have NA values in one of the logFC columns
## If either logFC_GSE8397 or logFC_GSE49036 has NA values (which is likely after the full join), cor.test will automatically exclude those pairs.
# sum(is.na(pd_no_sva_combo$logFC_GSE8397) | is.na(pd_no_sva_combo$logFC_GSE49036))
# ad_sva_be_combo<-merge(msbb_sva_be,rosmap_sva_be,by='X', all.x=TRUE, all.y=TRUE,suffixes=c('_msbb','_rosmap'))
# #cor.test(sva_be_combo$logFC_msbb, sva_be_combo$logFC_rosmap, method = "spearman")  #r=0.1157661, p-value < 2.2e-16

# Get all unique PD dataset base names (e.g., GSE8397)
dataset_names <- unique(clean_pair_name(names(pd_datasets)))
pair_grid <- combn(dataset_names, 2, simplify = FALSE)

# For each unique pair, compute Spearman for _none and _sva
get_pairwise_results <- function(dataset1, dataset2) {
  df1_none <- pd_datasets[[paste0(dataset1, "_none")]]
  df2_none <- pd_datasets[[paste0(dataset2, "_none")]]
  df1_sva  <- pd_datasets[[paste0(dataset1, "_sva")]]
  df2_sva  <- pd_datasets[[paste0(dataset2, "_sva")]]

  # Subset logFC columns
  df1_none <- df1_none[, c("X", "logFC")]
  df2_none <- df2_none[, c("X", "logFC")]
  df1_sva  <- df1_sva[, c("X", "logFC")]
  df2_sva  <- df2_sva[, c("X", "logFC")]

  names(df1_none)[2] <- "logFC_1"
  names(df2_none)[2] <- "logFC_2"
  names(df1_sva)[2] <- "logFC_1"
  names(df2_sva)[2] <- "logFC_2"

  merged_none <- merge(df1_none, df2_none, by = "X")
  merged_sva  <- merge(df1_sva, df2_sva, by = "X")

  cor_none <- if (nrow(merged_none) > 0) cor(merged_none$logFC_1, merged_none$logFC_2, method = "spearman") else NA
  cor_sva  <- if (nrow(merged_sva) > 0)  cor(merged_sva$logFC_1,  merged_sva$logFC_2,  method = "spearman") else NA

  data.frame(
    Dataset_Pair = paste(dataset1, "vs", dataset2),
    Method = c("No SVA", "SVA 'BE'"),
    Spearman = c(cor_none, cor_sva),
    Disease = "PD"
  )
}

# Compute all pairwise results
pd_pairs_df <- do.call(rbind, lapply(pair_grid, function(p) get_pairwise_results(p[1], p[2])))

# Assume cor_values already contains AD, BP, SZ summary correlations
cor_values$Dataset_Pair <- paste0(cor_values$Disease, "_summary")
cor_values_fmt <- cor_values[, c("Dataset_Pair", "Method", "Spearman", "Disease")]

cor_all <- rbind(cor_values_fmt, pd_pairs_df)

a <- ggplot(cor_all, aes(x = Spearman, y = Dataset_Pair, color = Method)) +
  geom_line(aes(group = Dataset_Pair), linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 3, shape = 16) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Spearman Correlation of the log"[2]*"(FC) between Datasets"),
    y = NULL,
    color = "Adjustment"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave("/hpc/users/hoangd02/www/plots/real_data_test1.pdf", plot = a, width = 8, height = 8)

####### reordering the lines
pair_diff <- cor_all %>%
  group_by(Dataset_Pair) %>%
  summarize(
    abs_diff = abs(diff(Spearman)),  # Assumes exactly 2 rows per pair
    .groups = "drop"
  )

cor_all <- left_join(cor_all, pair_diff, by = "Dataset_Pair")

cor_all$Dataset_Pair <- factor(cor_all$Dataset_Pair,
                               levels = cor_all %>%
                                 distinct(Dataset_Pair, abs_diff) %>%
                                 arrange(abs_diff) %>%
                                 pull(Dataset_Pair))

a <- ggplot(cor_all, aes(x = Spearman, y = Dataset_Pair, color = Method)) +
  geom_line(aes(group = Dataset_Pair), linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 3, shape = 16) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Spearman Correlation of the log"[2]*"(FC) between Datasets"),
    y = NULL,
    color = "Adjustment"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank()
  ) + xlim(-0.1, 0.5)

ggsave("/hpc/users/hoangd02/www/plots/real_data_test2.pdf", plot = a, width = 8, height = 8)














