##Ryan's original code to calculate the cosine angle: 
design_dotprod <- t(design) %*% design
design_norm <- sqrt(colSums(design^2))
design_norm_prod_mat <- design_norm %*% t(design_norm)
design_cos_value_mat <- design_dotprod / design_norm_prod_mat
design_cos_value_mat[] <- pmin(1,pmax(-1, design_cos_value_mat))
design_angle_mat_rads <- acos(design_cos_value_mat)
design_angle_mat_degs <- design_angle_mat_rads / pi * 180


##converting Ryan's original code into a function: 
#note: (1) change second "design" to design_sva in step 1 to compute design_dotprod; (2) there will be two norms, one is norm of the original design and the second is the norm of the sva design; (3) change design_norm to sva_design_norm in step 3; step 4 to step 7: should be the same; final output: the rows should list the covariates in the original design and the columns should be the different SVs


compute_design_angles <- function(design, sva_design, supplement_angles_over_90 = TRUE) {
  # Step 1: Calculate the dot product of the original design and the SVA design
  design_dotprod <- t(design) %*% sva_design
  
  # Step 2: Compute the norm of each column in both the original and SVA designs
  design_norm <- sqrt(colSums(design^2))
  sva_design_norm <- sqrt(colSums(sva_design^2))
  
  # Step 3: Compute the product of norms for each pair of original and SVA columns
  design_norm_prod_mat <- outer(design_norm, sva_design_norm)

  #The outer function in R computes the outer product of two vectors, returning a matrix where each element is the product of the corresponding elements from the two input vectors.
  
  # Step 4: Calculate the cosine similarity matrix
  design_cos_value_mat <- design_dotprod / design_norm_prod_mat
  
  # Step 5: Clamp the values between -1 and 1 to avoid numerical errors
  design_cos_value_mat[] <- pmin(1, pmax(-1, design_cos_value_mat))
  
  # Step 6: Compute the angle matrix in radians
  design_angle_mat_rads <- acos(design_cos_value_mat)
  
  # Step 7: Convert the angles from radians to degrees
  design_angle_mat_degs <- design_angle_mat_rads / pi * 180

  if (supplement_angles_over_90) {
  	design_angle_mat_degs[] <- pmin(design_angle_mat_degs, 180 - design_angle_mat_degs)
  }
  
  # Assign row and column names for clarity
  rownames(design_angle_mat_degs) <- colnames(design)
  colnames(design_angle_mat_degs) <- colnames(sva_design)
  
  # Return the result
  return(design_angle_mat_degs)
}


############################################ STEPS ############################################
#1. prepare original design and SVA design matrices - TRY ON MSBB!!!!
#run sva, extract SVs (svobj$sv) and bind them in your original design matrix
library(sva)

status="PlaqueMean"
formula_test=formula(paste("~ ", status," + PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))
formula_test0=formula(paste("~ PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))

design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

#Applying Voom Transformation 
isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata) 
dge <- DGEList(counts=genedata[isexpr,]) 
dge <- calcNormFactors(dge) 
dge <- dge[, colnames(dge) %in% rownames(design)] 
 
v <- voom(dge, design, plot=TRUE)

n.sv = num.sv(v$E,design,method="be") #27 SVs for "be" 0 SV for "leek"

if (n.sv > 0) {
  svobj <- sva(v$E, design, design0, n.sv = n.sv)
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  design_sv <- cbind(design, temp)
} else {
  design_sv <- design
}

#2. call the compute_design_angles function above
if (n.sv > 0) {
  sva_angles <- compute_design_angles(design, svobj$sv)
} else {
  message("No surrogate variables estimated; skipping angle computation.")
  sva_angles <- NULL
}
 
#3. output results
if (!is.null(sva_angles)) {
  print(sva_angles)  
  write.csv(sva_angles, "sva_design_angles_MSBB.csv", row.names = TRUE) 
}

#subset only for "PlaqueMean"
plaquemean_angles_matrix <- sva_angles["PlaqueMean", , drop = FALSE]

# Check the result
write.csv(plaquemean_angles_matrix, "sva_design_angles_MSBB_only_plaquemean.csv", row.names = TRUE) 

################################for ROSMAP##########################################
##WITH A CONTRAST 
#1. prepare original design and SVA design matrices - TRY ON MSBB!!!!
#run sva, extract SVs (svobj$sv) and bind them in your original design matrix
library(sva)

# Filter lowly expressed genes and calculate normalization factors
##for regular DEA, use this: 
formula_test=formula(paste("~ 0 + ceradsc + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

formula_test0=formula(paste("~ 0 + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
dge <- DGEList(counts=genedata[isexpr,]) 
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

# Fit the model with the original design
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = design)
contrast_design <- design %*% contrast.matrix

# Add a name to the contrast vector for clarity
colnames(contrast_design) <- "ceradsc1_ceradsc4"

n.sv = num.sv(v$E,design,method="be") #48 SVs for "be" 0 SV for "leek"
if (n.sv > 0) {
  sva_angles_contrast <- compute_design_angles(contrast_design, svobj$sv)
} else {
  message("No surrogate variables estimated; skipping angle computation.")
  sva_angles_contrast <- NULL
}

# View or save the results
if (!is.null(sva_angles_contrast)) {
  print(sva_angles_contrast)  # Print the angles matrix
  write.csv(sva_angles_contrast, "sva_angles_contrast_ROSMAP2.csv", row.names = TRUE)  # Save to a file
}

################################SZ##########################################
#see SZ code for publication.R

MSSM_Penn_Pitt_formula_test <- formula(paste("~ Dx + Reported_Gender + RIN +","scale(IntronicRate) + scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +", "Institution * (ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3)", sep = ""))
MSSM_Penn_Pitt_formula_test0 <- formula(paste("~ Reported_Gender + RIN +","scale(IntronicRate) + scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +", "Institution * (ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3)", sep = ""))

design.MSSM=model.matrix(formula(paste0("~",gsub("1 |","",as.character(MSSM_Penn_Pitt_formula_test)[2],fixed=TRUE))),METADATA.MSSM)
design0.MSSM=model.matrix(formula(paste0("~",gsub("1 |","",as.character(MSSM_Penn_Pitt_formula_test0)[2],fixed=TRUE))), METADATA.MSSM)

NIMH_HBCC_formula_test <- formula(paste("~ Dx + Reported_Gender + RIN +",
    "scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +",
    "ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3",
    sep = ""))
NIMH_HBCC_formula_test0 <- formula(paste("~ Reported_Gender + RIN +",
    "scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +",
    "ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3",
    sep = ""))

design.HBCC=model.matrix(formula(paste0("~",gsub("1 |","",as.character(NIMH_HBCC_formula_test)[2],fixed=TRUE))),METADATA.HBCC)
design0.HBCC=model.matrix(formula(paste0("~",gsub("1 |","",as.character(NIMH_HBCC_formula_test0 )[2],fixed=TRUE))), METADATA.HBCC)

v.MSSM <- voom(dge.MSSM, design.MSSM, plot=FALSE) 
v.HBCC <- voom(dge.HBCC, design.HBCC, plot=FALSE) 

##RUN THIS ONE AT A TIME FOR v.MSSM AND v.HBCC 
n.sv = num.sv(v.MSSM$E, design.MSSM ,method="be") #50 SVs for "be" 0 SV for "leek"
n.sv = num.sv(v.HBCC$E, design.HBCC ,method="be") #31 SVs for "be" 0 SV for "leek"

#change out v and design 
if (n.sv > 0) {
  svobj <- sva(v.HBCC$E, design.HBCC, design0.HBCC, n.sv = n.sv)
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  design_sv <- cbind(design.HBCC, temp)
} else {
  design_sv <- design
}

if (n.sv > 0) {
  sva_angles <- compute_design_angles(design.HBCC, svobj$sv)
} else {
  message("No surrogate variables estimated; skipping angle computation.")
  sva_angles <- NULL
}
 
#3. output results
if (!is.null(sva_angles)) {
  print(sva_angles)  
  write.csv(sva_angles, "/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_SZ_HBCC.csv", row.names = TRUE) 
}

#subset only for "DxSCZ"
DxSCZ_angles_matrix <- sva_angles["DxSCZ", , drop = FALSE]; DxSCZ_angles_matrix

# Check the result
write.csv(DxSCZ_angles_matrix, "/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_SZ_HBCC_only_dx.csv", row.names = TRUE) 


################################BP##########################################
##it's a bit weird when im not producing the same results 

#see BP code for publication.R
MSSM_Penn_Pitt_formula_test <- formula(paste("~ Dx + Reported_Gender + RIN +","scale(IntronicRate) + scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +", "Institution * (ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3)", sep = ""))
MSSM_Penn_Pitt_formula_test0 <- formula(paste("~ Reported_Gender + RIN +","scale(IntronicRate) + scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +", "Institution * (ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3)", sep = ""))

design.MSSM=model.matrix(formula(paste0("~",gsub("1 |","",as.character(MSSM_Penn_Pitt_formula_test)[2],fixed=TRUE))),METADATA.MSSM) #dim: 336 x 22
design0.MSSM=model.matrix(formula(paste0("~",gsub("1 |","",as.character(MSSM_Penn_Pitt_formula_test0)[2],fixed=TRUE))), METADATA.MSSM)

NIMH_HBCC_formula_test <- formula(paste("~ Dx + Reported_Gender + RIN +",
    "scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +",
    "ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3",
    sep = ""))
NIMH_HBCC_formula_test0 <- formula(paste("~ Reported_Gender + RIN +",
    "scale(IntragenicRate) + scale(IntergenicRate) + scale(rRNARate) +",
    "ageOfDeath + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3",
    sep = ""))

design.HBCC=model.matrix(formula(paste0("~",gsub("1 |","",as.character(NIMH_HBCC_formula_test)[2],fixed=TRUE))),METADATA.HBCC)
design0.HBCC=model.matrix(formula(paste0("~",gsub("1 |","",as.character(NIMH_HBCC_formula_test0 )[2],fixed=TRUE))), METADATA.HBCC)
v.MSSM <- voom(dge.MSSM, design.MSSM, plot=FALSE) 
v.HBCC <- voom(dge.HBCC, design.HBCC, plot=FALSE) 

##RUN THIS ONE AT A TIME FOR v.MSSM AND v.HBCC 
n.sv = num.sv(v.MSSM$E, design.MSSM ,method="be") #36 SVs for "be" 0 SV for "leek"
n.sv = num.sv(v.HBCC$E, design.HBCC ,method="be") #29 SVs for "be" 0 SV for "leek"

#change out v and design 
if (n.sv > 0) {
  svobj <- sva(v.HBCC$E, design.HBCC, design0.HBCC, n.sv = n.sv)
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  design_sv <- cbind(design.HBCC, temp)
} else {
  design_sv <- design.HBCC
}

if (n.sv > 0) {
  sva_angles <- compute_design_angles(design.HBCC, svobj$sv)
} else {
  message("No surrogate variables estimated; skipping angle computation.")
  sva_angles <- NULL
}
 
#3. output results
if (!is.null(sva_angles)) {
  print(sva_angles)  
  write.csv(sva_angles, "/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_BP_HBCC2.csv", row.names = TRUE) 
}

#subset only for "DxSCZ"
DxBP_angles_matrix <- sva_angles["DxBP", , drop = FALSE]; DxBP_angles_matrix

# Check the result
write.csv(DxBP_angles_matrix, "/sc/arion/projects/mscic1/results/jolie/CMC/sva_design_angles_BP_HBCC_only_dx2.csv", row.names = TRUE) 



#sva_angles_contrast_ROSMAP.csv has the original non-supplementary angles

# #Prepare the data for plotting
# sv_indices <- paste0("SV", 1:48)
# cosine_angles <- as.numeric(sva_angles_contrast)
# plot_data <- data.frame(SV = sv_indices, Angle = cosine_angles)

# #Bar plot
# library(ggplot2)
# plot<-ggplot(plot_data, aes(x = SV, y = Angle)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   theme_minimal() +
#   labs(title = "Cosine Angles Between Primary Variable of Interest and Surrogate Variables",
#        x = "Surrogate Variables (SVs)",
#        y = "Cosine Angle (degrees)") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/cosine_angles_SVs_primary_var_ROSMAP.pdf", plot = plot, width = 8, height = 5)

# # Dot plot with line
# plot<-ggplot(plot_data, aes(x = as.numeric(gsub("SV", "", SV)), y = Angle)) +
#   geom_point(size = 3, color = "blue") +          # Add dots
#   geom_line(color = "gray", size = 1) +          # Connect dots with a line
#   scale_x_continuous(breaks = 1:48) +            # Ensure all SVs are shown on the x-axis
#   theme_minimal() +
#   labs(title = "Cosine Angles Between Primary Variable of Interest and Surrogate Variables",
#        x = "Surrogate Variables (SVs)",
#        y = "Cosine Angle (degrees)") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/cosine_angles_SVs_primary_var_ROSMAP2.pdf", plot = plot, width = 8, height = 5)


###reading results in
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

#rosmap_angles
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
79.75   86.31   88.27   87.60   89.18   89.94 

#msbb_angles
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
81.49   84.94   87.13   86.64   88.80   89.75

#bp_hbcc_angles 
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
74.16   83.76   87.46   85.85   88.69   89.99 

#bp_mssm_angles
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
76.78   85.70   86.63   86.53   88.92   89.83 

#sz_hbcc_angles
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
76.82   86.09   87.51   86.51   88.91   89.91 

#sz_mssm_angles
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
80.21   86.75   88.14   87.78   89.28   89.94 


library(dplyr)
library(tidyr)
library(ggplot2)

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

dataset_order <- c("rosmap_angles", "msbb_angles", "bp_hbcc_angles",
  "bp_mssm_angles", "sz_hbcc_angles", "sz_mssm_angles")

dataset_labels <- c("AD: ROSMAP", "AD: MSBB", "BP: HBCC", "BP: MPP", "SZ: HBCC", "SZ: MPP")

# Create the boxplot
p1<-ggplot(combined_data, aes(x = dataset, y = Value, fill = dataset)) +
  geom_boxplot() +
  scale_x_discrete(
    limits = dataset_order,  # Order of datasets
    labels = dataset_labels  # Custom labels for datasets
  ) +
  labs(
    x = "Datasets",
    y = "Angles"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Remove legend for clarity
  )

ggsave("/hpc/users/hoangd02/www/plots/cosine_angles_combo_plot2.pdf", plot = p1, width = 8, height = 5) 
























