###### ###### ###### This figure shows cosine angles and ROSMAP shuffling  
###################################################################################
############################## from cosine_angle_calculation.sh ###################
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

rosmap_angles<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/sva_angles_contrast_ROSMAP2.csv", row.names = 1) #48 SVs
# rosmap_angles_matrix <- as.matrix(rosmap_angles)
# rosmap_angles_matrix

msbb_angles <-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/sva_design_angles_MSBB_only_plaquemean.csv", row.names = 1) #27 SVs

bp_hbcc_angles <-read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_BP_HBCC_only_dx2.csv", row.names = 1) #29 SVs

bp_mssm_angles <-read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_BP_MSSM_only_dx.csv", row.names = 1) #36 SVs

sz_hbcc_angles <-read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_SZ_HBCC_only_dx.csv", row.names = 1) #31 SVs

sz_mssm_angles <-read.csv("/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_SZ_MSSM_only_dx.csv", row.names = 1) #50 SVs

##summary stats for angles (depicted in boxplot)
column_means <- colMeans(sz_mssm_angles)
summary_stats <- summary(column_means); summary_stats

# #rosmap_angles
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 79.75   86.31   88.27   87.60   89.18   89.94 

# #msbb_angles
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 81.49   84.94   87.13   86.64   88.80   89.75

# #bp_hbcc_angles 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 74.16   83.76   87.46   85.85   88.69   89.99 

# #bp_mssm_angles
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 76.78   85.70   86.63   86.53   88.92   89.83 

# #sz_hbcc_angles
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 76.82   86.09   87.51   86.51   88.91   89.91 

# #sz_mssm_angles
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 80.21   86.75   88.14   87.78   89.28   89.94 


# Combine the datasets into a single dataframe
dataframes <- list(
  rosmap_angles = rosmap_angles,
  msbb_angles = msbb_angles,
  bp_hbcc_angles = bp_hbcc_angles,
  bp_mssm_angles = bp_mssm_angles,
  sz_hbcc_angles = sz_hbcc_angles,
  sz_mssm_angles = sz_mssm_angles
)

# Add a dataset identifier and reshape to long format
combined_data <- bind_rows(lapply(names(dataframes), function(name) {
  df <- dataframes[[name]]
  df$dataset <- name  # Add a dataset column
  df  # Return the modified dataframe
}), .id = "id") %>%
  pivot_longer(
    cols = starts_with("V"),  # All columns that start with "V"
    names_to = "Variable",
    values_to = "Value"
  )

dataset_order <- c("rosmap_angles", "msbb_angles",  "bp_mssm_angles", "bp_hbcc_angles",
  "sz_mssm_angles", "sz_hbcc_angles")

dataset_labels <- c("AD: ROSMAP", "AD: MSBB", "BP: CMC", "BP: CMC-HBCC", "SZ: CMC", "SZ: CMC-HBCC")

custom_colors <- c(
  "rosmap_angles" = "#A3392F",    
  "msbb_angles" = "#E37C72",   
  "bp_mssm_angles" = "#BBD8A3",   
  "bp_hbcc_angles" = "#79AC78",
  "sz_mssm_angles" = "#C4A5D6",   
  "sz_hbcc_angles" = "#734F96"     
)

p1<-ggplot(combined_data, aes(x = dataset, y = Value, fill = dataset)) +
  geom_boxplot() +
  scale_x_discrete(
    limits = dataset_order, 
    labels = dataset_labels  
  ) +
  scale_fill_manual(values = custom_colors) +  
  labs(
    x = "Datasets",
    y = "Angles between SVs and Phenotype of Interest"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),  
    legend.position = "none"  
  )

ggsave("/hpc/users/hoangd02/www/plots/cosine_angles.pdf", plot = p1, width = 9, height = 7) 

#ggsave("/hpc/users/hoangd02/www/plots/cosine_angles_combo_plot2.pdf", plot = p1, width = 8, height = 5) 

###################################################################################
############################## from ROSMAP_random_shuffling.R #####################
## Load in ROSMAP data, 1 SV at at time from Figure6_20250622.sh
rosmap_sva_be_0_SVs <- fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_11.18.txt", data.table = FALSE) #19476     7

rosmap_path <- "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/"
rosmap_files <- paste0(rosmap_path, "rosmap_sva_be_", 1:48, "_SVs.csv")
rosmap_list <- lapply(rosmap_files, read.csv)
rosmap_list <- lapply(rosmap_list, as.data.frame)
cat("ROSMAP list length:", length(rosmap_list), "\n")

rosmap_list <- c(list(rosmap_sva_be_0_SVs), rosmap_list)

head(rosmap_list[[1]])

deg_counts <- sapply(rosmap_list, function(res) {
  sum(res$adj.P.Val < 0.05, na.rm = TRUE)
})

deg_summary <- data.frame(
  n_sv = 0:48,
  DEGs = deg_counts
)

shuffled_res <- read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/random_shuffling_combined_results_with_seed.csv")
dim(shuffled_res) #5100    5

# p2 <- ggplot(shuffled_res, aes(x = factor(n_sv), y = DEGs, fill = Method)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("SVA" = "lightblue", "No SVA" = "navy")) +
#   scale_x_discrete(
#     breaks = as.character(c(0, 10, 20, 30, 40, 50)),
#     labels = c("0", "10", "20", "30", "40", "50")
#   ) +
#   labs(
#     x = "Number of Surrogate Variables",
#     y = "Number of DEGs",
#     fill = "Method"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 0.5),
#     legend.position = "none"
#   )

deg_summary$Source <- "Unshuffled Data"
shuffled_res$Source <- "Shuffled Data"

p3 <- ggplot(shuffled_res, aes(x = factor(n_sv), y = DEGs)) +
  geom_boxplot(aes(fill = Source)) +
  geom_point(
    data = deg_summary,
    aes(x = factor(n_sv), y = DEGs, color = Source, shape = Source),
    size = 2.5
  ) +
  scale_fill_manual(
    values = c("Shuffled Data" = "gray80"),
    name = "Data Type"
  ) +
  scale_color_manual(
    values = c("Unshuffled Data" = "red"),
    name = "Data Type"
  ) +
  scale_shape_manual(
    values = c("Unshuffled Data" = 17),
    name = "Data Type"
  ) +
  scale_x_discrete(
    breaks = as.character(c(0, 10, 20, 30, 40, 50)),
    labels = c("0", "10", "20", "30", "40", "50")
  ) +
  labs(
    x = "Number of Surrogate Variables",
    y = "Number of DEGs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

#ggsave("/hpc/users/hoangd02/www/plots/rosmap_shuffled_boxplot.pdf", plot = p3, width = 9, height = 7) 

combo <- p1 + p3

ggsave("/hpc/users/hoangd02/www/plots/cosine_angles_and_rosmap_shuffled_boxplot.pdf", plot = combo, width = 14, height = 6) 

# a <- ggplot(shuffled_res, aes(x = Method, y = DEGs, fill = Method)) +
#   geom_boxplot() +
#   facet_wrap(~n_sv, scales = "free_y") +
#   theme_minimal() +
#   labs(title = "Distribution of DEGs with and without SVA across varying n_sv",
#        x = "Method",
#        y = "Number of DEGs",
#        fill = "Method")

# ggsave("/hpc/users/hoangd02/www/plots/test.pdf", plot = a, width = 9, height = 7) 



















