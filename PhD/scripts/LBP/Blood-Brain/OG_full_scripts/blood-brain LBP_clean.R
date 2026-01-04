##Blood - living brain transcriptomics directory
/sc/arion/projects/mscic1/results/jolie/LBP/LBP_and_Blood_CompiledData 

##what are these files??
lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMI_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_18JAN2022.RDS
lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMICellTypes_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_11FEB2022.RDS
lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS
lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_forQC_530LIVINGSamples_23FEB2022.RDS
lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_18FEB2022.RDS

##Blood - living brain transcriptomics 
/sc/arion/projects/mscic1/results/jolie/LBP/metadata_covar_lbp_blood.RDS

##Blood - living brain transcriptomics 
/sc/arion/projects/mscic1/results/jolie/LBP/living_postmortem_files.RDS

## blood-brain folder 
#screen section cor_play

########### from 3. raw_blood_GE.sh ####START HERE 2-13-2025
rm(list=ls())
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)

metadata <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")          
dim(metadata) ##[1] 530 170   
metadata$SAMPLE_ISMMS
metadata$IID_ISMMS ##finding the paired samples aka overlapping blood and brain samples with the same IID_ISMMS
 
##filter out brain and blood 
table(metadata$mymet_brain)
# blood brain 
#   243   287 

table(metadata$mymet_tissue)
# L_Blood L_Brain R_Blood R_Brain 
#     128     157     115     130 

metadata$IID_ISMMS ##list the IDs 
uniqueN(metadata[,.(IID_ISMMS)]) ####[1] 172

brain_metadata <- metadata[metadata$mymet_brain == "brain", ]
dim(brain_metadata)
#287 170
uniqueN(brain_metadata[,.(IID_ISMMS)]) ####[1] 171 (number of people)

blood_metadata <- metadata[metadata$mymet_brain == "blood", ]
dim(blood_metadata)
#243 170
uniqueN(blood_metadata[,.(IID_ISMMS)]) ####[1] 155 (number of people)

common_samples <- merge(brain_metadata, blood_metadata, by = "IID_ISMMS")
head(common_samples)
uniqueN(common_samples[,.(IID_ISMMS)]) ####[1] 154 (number of people)

colnames(common_samples) <- gsub("\\.x$", "_brain", colnames(common_samples))
colnames(common_samples) <- gsub("\\.y$", "_blood", colnames(common_samples))

common_samples[,c("IID_ISMMS","SAMPLE_ISMMS_brain","BARCODE_ISMMS_brain","STAR_Average_mapped_length_brain","GLU_brain","GLU_blood","STAR_Average_mapped_length_blood","bbstatus_brain","bbstatus_blood","mymet_tissue_brain","mymet_tissue_blood")]

table(common_samples$IID_ISMMS) #values are either 1, 2, or 4 ##################################
one_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 1)  
one_sample <- as.data.frame(one_sample)  
dim(one_sample) #37 339
table(one_sample$IID_ISMMS) #37 people 
# PT-0019 PT-0025 PT-0036 PT-0037 PT-0039 PT-0045 PT-0053 PT-0058 PT-0062 PT-0064 
#       1       1       1       1       1       1       1       1       1       1 
# PT-0068 PT-0070 PT-0074 PT-0078 PT-0079 PT-0080 PT-0090 PT-0095 PT-0099 PT-0106 
#       1       1       1       1       1       1       1       1       1       1 
# PT-0117 PT-0130 PT-0131 PT-0140 PT-0144 PT-0148 PT-0149 PT-0150 PT-0155 PT-0171 
#       1       1       1       1       1       1       1       1       1       1 
# PT-0173 PT-0187 PT-0189 PT-0190 PT-0191 PT-0194 PT-0196 
#       1       1       1       1       1       1       1 

test <- filter(one_sample,IID_ISMMS == "PT-0019")
test3<- filter(common_samples, SAMPLE_ISMMS_blood == "LBPSEMA4BLOOD326") #LBPSEMA4BLOOD326

##check whether mymet_tissue_brain and mymet_tissue_blood start with the same letter (either "R" or "L") for the same IID_ISMMS in one_sample
one_sample$same_prefix <- substr(one_sample$mymet_tissue_brain, 1, 1) == 
                              substr(one_sample$mymet_tissue_blood, 1, 1)

# Count TRUE vs FALSE cases
table(one_sample$same_prefix)
# FALSE  TRUE 
#     1    36 

##there was one sample (PT-0117) that was mismatched (they have L_Brain and R_Blood so i changed it to L_Blood)
one_sample$mymet_tissue_blood[one_sample$IID_ISMMS == "PT-0117" & one_sample$mymet_tissue_blood == "R_Blood"] <- "L_Blood"
one_sample$mymet_timepoint_blood[one_sample$IID_ISMMS == "PT-0117" & one_sample$mymet_timepoint_blood == "right"] <- "left"

two_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 2)  
two_sample <- as.data.frame(two_sample)  
dim(two_sample) #74 339
table(two_sample$IID_ISMMS) #37 people 

######### INVESTIGATING TWO_SAMPLE CASES FURTHER ##################
### 8 INDIVIDUALS WITH 1 BRAIN SAMPLE AND 2 BLOOD SAMPLES
# Count distinct tissue values per ID
id_counts <- aggregate(mymet_tissue_brain ~ IID_ISMMS, data = two_sample, function(x) length(unique(x)))

# Get IDs where mymet_tissue_brain has only one unique value
valid_ids <- id_counts$IID_ISMMS[id_counts$mymet_tissue_brain == 1]

# Subset original data
one_brain_two_blood_samples <- two_sample[two_sample$IID_ISMMS %in% valid_ids, ]
table(one_brain_two_blood_samples$IID_ISMMS)
# PT-0032 PT-0046 PT-0054 PT-0071 PT-0100 PT-0116 PT-0124 PT-0177 #8 people
#       2       2       2       2       2       2       2       2

# Keep only matching L_Brain-L_Blood and R_Brain-R_Blood pairs
paired_one_brain_two_blood_samples <- one_brain_two_blood_samples %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))

# Check the new row count
nrow(paired_one_brain_two_blood_samples) #8

################################################################## FINAL: paired_one_brain_two_blood_samples

### 29 INDIVIDUALS WITH 2 BRAIN SAMPLES AND 1 BLOOD SAMPLE
table(two_sample$IID_ISMMS)

# Count distinct tissue values per ID
id_counts <- aggregate(mymet_tissue_brain ~ IID_ISMMS, data = two_sample, function(x) length(unique(x)))

# Get IDs where mymet_tissue_brain has only one unique value
valid_ids <- id_counts$IID_ISMMS[id_counts$mymet_tissue_brain == 2]

# Subset original data
two_brain_one_blood_sample <- two_sample[two_sample$IID_ISMMS %in% valid_ids, ]
table(two_brain_one_blood_sample$IID_ISMMS)
# PT-0020 PT-0056 PT-0084 PT-0088 PT-0091 PT-0094 PT-0102 PT-0109 PT-0111 PT-0115 
#       2       2       2       2       2       2       2       2       2       2 
# PT-0118 PT-0119 PT-0120 PT-0121 PT-0136 PT-0138 PT-0139 PT-0143 PT-0145 PT-0147 
#       2       2       2       2       2       2       2       2       2       2 
# PT-0154 PT-0157 PT-0161 PT-0163 PT-0172 PT-0174 PT-0175 PT-0179 PT-0184 
#       2       2       2       2       2       2       2       2       2

# Keep only matching L_Brain-L_Blood and R_Brain-R_Blood pairs
paired_two_brain_one_blood_sample <- two_brain_one_blood_sample %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))

# Check the new row count
nrow(paired_two_brain_one_blood_sample) #29

################################################################## OR FINAL: two_sample


################## INVESTIGATING FOUR_SAMPLE CASES FURTHER ####################################
four_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 4)  
four_sample <- as.data.frame(four_sample)  
dim(four_sample) #320 339
table(four_sample$IID_ISMMS) #80 people 

test<-filter(four_sample, IID_ISMMS == "PT-0024")
test[,c("IID_ISMMS","mymet_tissue_blood","mymet_tissue_brain","mymet_extractionkit_blood","mymet_extractionkit_brain","Monocyte_total_scp_blood","GABA_brain")]

# Keep only matching L_Brain-L_Blood and R_Brain-R_Blood pairs
paired_four_sample <- four_sample %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))

nrow(paired_four_sample) #160

one_sample$same_prefix <- NULL
two_sample$same_prefix <- NULL

# common_samples_with_dup <- rbind(one_sample,two_sample,paired_four_sample)
# nrow(one_sample) #37
# nrow(two_sample) #74
# nrow(paired_four_sample) #160
# dim(common_samples_with_dup) #271 339 (37 + 74 + 160 = 271)

common_samples_without_dup <- rbind(one_sample, paired_one_brain_two_blood_samples, paired_two_brain_one_blood_sample, paired_four_sample)
nrow(one_sample) #37
nrow(paired_one_brain_two_blood_samples) #8
nrow(paired_two_brain_one_blood_sample) #29
nrow(paired_four_sample) #160
dim(common_samples_without_dup) #234 339 (37 + 8 + 29 + 160 = 234)

##drop 1 sample with ID (IID_ISMMS == "PT-0117") because brain and blood sample with this ID might be collected from different time point
common_samples_without_dup <- common_samples_without_dup[common_samples_without_dup$IID_ISMMS != "PT-0117", ]
dim(common_samples_without_dup) #233 339

length(unique(common_samples_without_dup$IID_ISMMS)) #153 people

# table(common_samples_with_dup$IID_ISMMS) #values are either 2 or 1
# overlapping_IDs <- names(table(common_samples_with_dup$IID_ISMMS)) 
# length(overlapping_IDs)#154 IDs

brain_metadata <-common_samples_without_dup[common_samples_without_dup$mymet_brain_brain == "brain", ] 
blood_metadata <-common_samples_without_dup[common_samples_without_dup$mymet_brain_blood == "blood", ] 

identical(brain_metadata$SAMPLE_ISMMS_brain, blood_metadata$SAMPLE_ISMMS_brain) #TRUE
identical(brain_metadata$SAMPLE_ISMMS_blood, blood_metadata$SAMPLE_ISMMS_blood) #TRUE

##good sanity check here: 
dim(blood_metadata) #233 339
dim(brain_metadata) #233 339

length(unique(blood_metadata$IID_ISMMS)) #153
length(unique(brain_metadata$IID_ISMMS)) #153

# ## remove samples with RIN < 6 -- should have done this in the beginning! but thats okay!
# sum(blood_metadata$mymet_rin_blood < 6, na.rm = TRUE) #28 samples - with what PT ID?
# sum(brain_metadata$mymet_rin_brain < 6, na.rm = TRUE) #18 samples - with what PT ID?

# blood_metadata <- blood_metadata[blood_metadata$mymet_rin_blood >= 6, ]
# brain_metadata <- brain_metadata[brain_metadata$mymet_rin_brain >= 6, ]

# dim(blood_metadata) #205 339
# dim(brain_metadata) #215 339

# overlapping_IIDs <- intersect(blood_metadata$IID_ISMMS, brain_metadata$IID_ISMMS)
# length(overlapping_IIDs) #141 people

# # Subset blood_metadata and ensure correct order
# blood_metadata_matched <- blood_metadata[match(overlapping_IIDs, blood_metadata$IID_ISMMS), ]
# # Subset brain_metadata and ensure correct order
# brain_metadata_matched <- brain_metadata[match(overlapping_IIDs, brain_metadata$IID_ISMMS), ]

# dim(blood_metadata_matched) #141 339
# dim(brain_metadata_matched) #141 339

# identical(brain_metadata$IID_ISMMS, blood_metadata$IID_ISMMS)

# table(blood_metadata$mymet_tissue_blood)
# # L_Blood L_Brain R_Blood R_Brain 
# #     112       0      93       0 

# table(blood_metadata$mymet_tissue_brain)
# # L_Blood L_Brain R_Blood R_Brain 
# #       0     112       0      93

# length(unique(blood_metadata$IID_ISMMS)) #143
# length(unique(brain_metadata2$IID_ISMMS)) #141

# missing_in_brain <- setdiff(unique(blood_metadata$IID_ISMMS), unique(brain_metadata$IID_ISMMS))
# print(missing_in_brain)  # Shows which samples are missing from brain_metadata
# #"PT-0139" "PT-0163" #why dont these people have brain data?

# brain_metadata[brain_metadata$IID_ISMMS %in% c("PT-0198", "PT-0163"), ]

############################ EXPRESSION #################
raw_count <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_530LIVINGsamples_15FEB2022.RDS")                                                                             
dim(raw_count) ##[1] 58929   530  

blood_ge <- raw_count[, colnames(raw_count) %in% blood_metadata$SAMPLE_ISMMS_blood]
dim(blood_ge)
#58929   233

brain_ge <- raw_count[, colnames(raw_count) %in% brain_metadata$SAMPLE_ISMMS_brain]
dim(brain_ge)
#58929   233

## OR -- this is for ALL samples, not just paired samples deriving from blood_metadata
blood_ge <- raw_count[, grep("BLOOD", colnames(raw_count))]
#58929   243

write.table(blood_ge, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_raw_ALL_samples_not_paired_samples.txt", sep="\t", quote=FALSE, row.names=TRUE)

identical(blood_metadata$IID_ISMMS, brain_metadata$IID_ISMMS) #TRUE
identical(brain_metadata, blood_metadata) #TRUE
#blood_metadata[,c("IID_ISMMS","SAMPLE_ISMMS_brain","SAMPLE_ISMMS_blood")]
#brain_metadata[,c("IID_ISMMS","SAMPLE_ISMMS_brain","SAMPLE_ISMMS_blood")]

identical(colnames(brain_ge), brain_metadata$SAMPLE_ISMMS_brain) #FALSE
# Reorder columns of brain_ge to match the order of brain_metadata$SAMPLE_ISMMS_brain
brain_ge <- brain_ge[, match(brain_metadata$SAMPLE_ISMMS_brain, colnames(brain_ge))]
identical(colnames(brain_ge), brain_metadata$SAMPLE_ISMMS_brain) #TRUE
rownames(brain_metadata) <- brain_metadata$SAMPLE_ISMMS_blood

identical(colnames(blood_ge), blood_metadata$SAMPLE_ISMMS_blood) #FALSE
# Reorder columns of brain_ge to match the order of brain_metadata$SAMPLE_ISMMS_brain
blood_ge <- blood_ge[, match(blood_metadata$SAMPLE_ISMMS_blood, colnames(blood_ge))]
identical(colnames(blood_ge), blood_metadata$SAMPLE_ISMMS_blood) #TRUE
rownames(blood_metadata) <- blood_metadata$SAMPLE_ISMMS_blood
identical(colnames(blood_ge),rownames(blood_metadata)) #TRUE

##
save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline.RData")

## globin genes removal - manually mapping because biomaRt SUCKS
#https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000161544;r=17:76527356-76551175
globin_genes <- c("CYGB", "HBA1", "HBA2", 
                  "HBB", "HBD", "HBE1", 
                  "HBG1", "HBG2", "HBM", 
                  "HBQ1", "HBZ", "MB") #12 genes

globin_ensembl_ids <- c("ENSG00000161544.10", "ENSG00000206172.8", "ENSG00000188536.13", 
                         "ENSG00000244734.4", "ENSG00000223609.11", "ENSG00000213931.7", 
                         "ENSG00000213934.9", "ENSG00000196565.15", "ENSG00000206177.7", 
                         "ENSG00000086506.3", "ENSG00000130656.6", "ENSG00000198125.13")
#12 genes

# Grab rows where row names match the IDs in globin_ensembl_ids
matched_rows <- blood_ge[rownames(blood_ge) %in% globin_ensembl_ids, ]
matched_rows #seems like we have 9 globin genes in our data

# Get existing gene IDs in matched_rows
existing_genes <- rownames(matched_rows)

# Find missing genes
missing_genes <- setdiff(globin_ensembl_ids, existing_genes)
missing_genes
#we dont have "ENSG00000223609.11" "ENSG00000213934.9"  "ENSG00000130656.6"
#aka "HBD", "HBG1" or "HBZ"

dim(blood_ge)
#58929   233

# Remove globin genes from blood_ge
blood_ge <- blood_ge[!rownames(blood_ge) %in% globin_ensembl_ids, ]

dim(blood_ge)
#58920   233

####### BLOOD
# Convert raw counts to DGEList
dge_blood <- DGEList(counts = blood_ge)

# Convert to CPM
cpm_blood <- cpm(dge_blood)

# Filtering: Keep genes with CPM >= 1 in at least 10% of samples
min_samples <- ceiling(0.1 * ncol(blood_ge))  # 10% of samples
keep_genes <- rowSums(cpm_blood >= 1) >= min_samples
dge_blood <- dge_blood[keep_genes, , keep.lib.sizes = FALSE]  # Apply filtering

# Apply TMM normalization for composition bias
dge_blood <- calcNormFactors(dge_blood, method = "TMM")

# Define formula (no covariates, just an intercept)
formula <- ~ 1

# Run voomWithDreamWeights with corrected metadata
voom_blood <- voomWithDreamWeights(dge_blood, formula = formula, data = blood_metadata)

# Extract log2 CPM expression values
blood_log2_cpm <- voom_blood$E

dim(blood_log2_cpm)
#21046   233

####### BRAIN
# Convert raw counts to DGEList
dge_brain <- DGEList(counts = brain_ge)

# Convert to CPM
cpm_brain <- cpm(dge_brain)

# Filtering: Keep genes with CPM >= 1 in at least 10% of samples
min_samples <- ceiling(0.1 * ncol(brain_ge))  # 10% of samples
keep_genes <- rowSums(cpm_brain >= 1) >= min_samples
dge_brain <- dge_brain[keep_genes, , keep.lib.sizes = FALSE]  # Apply filtering

# Apply TMM normalization for composition bias
dge_brain <- calcNormFactors(dge_brain, method = "TMM")

# Define formula (no covariates, just an intercept)
formula <- ~ 1

# Run voomWithDreamWeights with corrected metadata
voom_brain <- voomWithDreamWeights(dge_brain, formula = formula, data = brain_metadata)

# Extract log2 CPM expression values
brain_log2_cpm <- voom_brain$E
dim(brain_log2_cpm)
#21360   233

#### FINDING COMMON GENES
## Brain: 21360 genes
## Blood: 21046 genes
## Common: 17533 genes

# Genes that are in Brain but not blood: 21360 - 17533 = 3827 genes 
# Genes that are in Blood but not brain: 21046 - 17533 = 3513 genes 
# Genes that are in Either (AKA Number of genes expressed in at least one (union)): 21360 + 3513 = 24873 genes
# Brain gene + Genes that are in Blood but not brain

## Common: 17533 genes
common_genes <- intersect(rownames(blood_log2_cpm), rownames(brain_log2_cpm))
num_overlapping <- length(common_genes)#17539 genes
num_overlapping

# Extract gene lists
brain_genes <- rownames(brain_log2_cpm)
blood_genes <- rownames(blood_log2_cpm)

# Percentage of overlapping genes in each dataset
percent_brain <- (num_overlapping / length(brain_genes)) * 100; percent_brain #82.08333 %
percent_blood <- (num_overlapping / length(blood_genes)) * 100; percent_blood #83.30799 %

# Get gene names from both datasets
blood_genes <- rownames(blood_log2_cpm)
brain_genes <- rownames(brain_log2_cpm)

# Find overlapping genes
common_genes <- intersect(rownames(blood_log2_cpm), rownames(brain_log2_cpm))

# Compute the percentage of overlapping genes
percent_overlap <- (length(common_genes) / length(unique(c(rownames(blood_log2_cpm), rownames(brain_log2_cpm))))) * 100
#denominator = 24873 genes
length(unique(c(rownames(blood_log2_cpm), rownames(brain_log2_cpm))))
# Genes that are in Either: 21360 + 3513 = 24873 genes

percent_overlap
#70.49767
#This is the percentage of genes expressed in both tissues relative to the total gene pool detected in either dataset


###### PCs on brain_log2_cpm and blood_log2_cpm
## see pca_blood-brain LBP_clean.sh, resulting in blood_technical_PCs and brain_technical_PCs
blood_technical_PCs 
brain_technical_PCs

#add the following biological covariates: mymet_age_brain, mymet_sex_brain, mymet_rin_brain, mymet_seqbatch_brain, and mymet_depletionbatch_brain

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_after_pca_before_residualizing_out_covariates.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_after_pca_before_residualizing_out_covariates.RData")

library(variancePartition)

# Define model for blood
# blood_formula <- as.formula("~ mymet_age_blood + mymet_sex_blood + mymet_depletionbatch_blood +
#                              PC1 + PC2 + PC3 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")

# # Define model for brain
# brain_formula <- as.formula("~ mymet_age_brain + mymet_sex_brain + mymet_depletionbatch_brain +
#                              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7")

# mymet_depletionbatch_brain
# #mymet_rin_brain + mymet_seqbatch_brain + mymet_depletionbatch_brain 

# # Run VPA for blood
# blood_varPart <- fitExtractVarPartModel(blood_log2_cpm, blood_formula, data = blood_metadata)

# # Run VPA for brain
# brain_varPart <- fitExtractVarPartModel(brain_log2_cpm, brain_formula, data = brain_metadata)

# # Generate variance partitioning plot
# blood_varPart_plot <- plotVarPart(blood_varPart)
# brain_varPart_plot <- plotVarPart(brain_varPart)

# # Save the plots
# ggsave(filename = "/hpc/users/hoangd02/www/plots/blood_varPart_plot2.png",
#        plot = blood_varPart_plot, 
#        width = 24, height = 6, dpi = 300)

# # Save the plot
# ggsave(filename = "/hpc/users/hoangd02/www/plots/brain_varPart_plot2.png",
#        plot = brain_varPart_plot, 
#        width = 20, height = 6, dpi = 300)

# # Convert categorical variables to numeric factors
# blood_metadata$mymet_seqbatch_blood <- as.numeric(as.factor(blood_metadata$mymet_seqbatch_blood))
# blood_metadata$mymet_depletionbatch_blood <- as.numeric(as.factor(blood_metadata$mymet_depletionbatch_blood))
# blood_metadata$mymet_sex_blood <- ifelse(blood_metadata$mymet_sex_blood == 1, 1, 0) #accidentally put 1 as male here whereas 0 should be male
# blood_metadata$mymet_sex_blood <- 1 - blood_metadata$mymet_sex_blood
# table(blood_metadata$mymet_sex_blood)
#  #  0   1 
#  # 93 140 
#  #0 is female and 1 is male, 140 M and 93 F

# # Recreate the numeric covariates dataset
# numeric_covariates <- blood_metadata[, c("mymet_age_blood","mymet_sex_blood", "mymet_rin_blood",
#                                          "mymet_seqbatch_blood","mymet_depletionbatch_blood",
#                                          "PC1", "PC2", "PC3", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

# # Compute Spearman correlation matrix
# cor_matrix<- cor(numeric_covariates, use = "pairwise.complete.obs", method = "spearman")

# # Visualize correlation heatmap
# library(corrplot)
# # Open PNG device
# png("/hpc/users/hoangd02/www/plots/blood_corrplot.png", width = 7, height = 7, units = "in", res = 300)

# # Generate the correlation plot
# corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7, diag = FALSE)

# # Close PNG device
# dev.off()

# ### cor for brain
# # Convert categorical variables to numeric factors for brain metadata
# brain_metadata$mymet_seqbatch_brain <- as.numeric(as.factor(brain_metadata$mymet_seqbatch_brain))
# brain_metadata$mymet_depletionbatch_brain <- as.numeric(as.factor(brain_metadata$mymet_depletionbatch_brain))
# brain_metadata$mymet_sex_brain <- ifelse(brain_metadata$mymet_sex_brain == "M", 1, 0)
# table(brain_metadata$mymet_sex_brain)
#  #  0   1 
#  # 93 140

# # Define numeric covariates including the PCs for brain
# numeric_covariates_brain <- brain_metadata[, c("mymet_age_brain", "mymet_sex_brain","mymet_rin_brain", "mymet_seqbatch_brain", "mymet_depletionbatch_brain","PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")]

# # Compute correlation matrix using Spearman (for ranked correlations)
# cor_matrix_brain<-cor(numeric_covariates_brain, use = "pairwise.complete.obs", method = "spearman")

# # Open PNG device to save the plot
# png("/hpc/users/hoangd02/www/plots/brain_corrplot.png", width = 7, height = 7, units = "in", res = 300)

# # Generate the correlation plot
# corrplot(cor_matrix_brain, method = "color", type = "upper", tl.cex = 0.7, diag = FALSE)

# # Close PNG device
# dev.off()


# # Load required libraries
# library(variancePartition)
# library(ggplot2)

# #Step 1: Define the Covariates to Regress Out
# blood_covariates_to_remove <- c("PC1", "PC2", "PC3", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
#                                 "mymet_age_blood", "mymet_sex_blood","mymet_depletionbatch_blood")

# brain_covariates_to_remove <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
#                                 "mymet_age_brain", "mymet_sex_brain","mymet_depletionbatch_blood") 

# ## STEP 2: Create the Design Matrix
# blood_design <- model.matrix(~ ., data = blood_metadata[, blood_covariates_to_remove])
# brain_design <- model.matrix(~ ., data = brain_metadata[, brain_covariates_to_remove])

# #Step 3: Residualize Gene Expression
# blood_residuals <- apply(blood_log2_cpm, 1, function(y) {
#   lm(y ~ blood_design)$residuals
# })
# blood_residuals <- t(blood_residuals)  # Transpose back to match original structure
# dim(blood_residuals)

# brain_residuals <- apply(brain_log2_cpm, 1, function(y) {
#   lm(y ~ brain_design)$residuals
# })
# brain_residuals <- t(brain_residuals)  # Transpose back to match original structure
# dim(brain_residuals)

################# CORRELATIONS 
############BEFORE RESIDUALIZING FOR COVARIATES
# Subset both matrices to only include common genes
blood_subset <- blood_log2_cpm[common_genes, ]
brain_subset <- brain_log2_cpm[common_genes, ]
#dim of each is 17533 x 233

######################## MAY 22 2025 START ###########################
# Subset genes
blood_sub <- blood_subset[1:100, ]
brain_sub <- blood_subset[1:100, ]

#test for subset | should i expect r=1 for the same gene in 2 different tissues?
cor_matrix_sub <- cor(t(blood_sub), t(brain_sub), method = "spearman")
dim(cor_matrix_sub)  # Should be 100 x 100
cor_matrix_sub[1:5, 1:5]  # View first few correlations

#test for all 
cor_matrix_all <- cor(t(blood_subset), t(brain_subset), method = "spearman")
dim(cor_matrix_all)  # Should be 100 x 100
cor_matrix_all[1:5, 1:5]  # View first few correlations

write.table(cor_matrix_all, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/cor_matrix_all_20250522.txt", sep="\t", quote=FALSE, row.names=TRUE)

library(data.table)
cor_matrix <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/cor_matrix_all_20250522.txt", data.table=FALSE)
dim(cor_matrix) #17533 x 17534
cor_matrix[1:5,1:5]

######################## MAY 22 2025 END ###########################
blood_subset_resid <- blood_residuals[common_genes, ]
brain_subset_resid <- brain_residuals[common_genes, ]
#dim of each is 17533 x 233

#this is before removing globin genes
# write.table(blood_subset, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_cpm_20250307.txt", sep="\t", quote=FALSE, row.names=TRUE)
# write.table(brain_subset, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/brain_expression_cpm_20250307.txt", sep="\t", quote=FALSE, row.names=TRUE)

#this is after removing globin genes and before residualization of covariates
write.table(blood_subset, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_cpm_nonresidualized_20250308.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(brain_subset, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/brain_expression_cpm_nonresidualized_20250308.txt", sep="\t", quote=FALSE, row.names=TRUE)

#this is after removing globin genes and after residualization of covariates
write.table(blood_subset_resid, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_cpm_residualized_20250309.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(brain_subset_resid, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/brain_expression_cpm_residualized_20250309.txt", sep="\t", quote=FALSE, row.names=TRUE)

blood_subset <- read.table("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_cpm_nonresidualized_20250308.txt", 
                         sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

brain_subset <- read.table("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/brain_expression_cpm_nonresidualized_20250308.txt", 
                         sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

#colnames(blood_subset_resid)<-colnames(brain_subset_resid) #double check if this is true

### PER ONE GENE 
# Compute Spearman correlation for each gene across samples
#AKA the correlation between expression values of a given gene across paired samples (blood vs. brain)
#this answers "Are genes consistently co-expressed in both blood and brain across individuals?"
#this helps identify genes whose expression is stable across tissues
#High correlation (r ≈ 1): Gene expression is similar in blood & brain across individuals.
#Low correlation (r ≈ 0): Gene expression varies independently between tissues.
gene_correlation_matrix <- mapply(function(x, y) cor(x, y, method = "spearman"), 
                                  as.data.frame(t(blood_subset)), #transpose so that rows = samples and columns = genes
                                  as.data.frame(t(brain_subset)))


#mapply() applies the function element-wise over corresponding columns in the two df, it takes two vectors (x and y) representing a gene's expression across samples and computes
#x = expression of a gene across samples in blood
#y = expression of the same gene across samples in brain
#for example:
x = c(5.2, 6.1, 4.9, 5.5, 5.8)  # GeneA expression in blood
y = c(5.1, 6.0, 4.8, 5.6, 5.7)  # GeneA expression in brain
Sample  GeneA (Blood) GeneA (Brain)
S1          5.2         5.1
S2          6.1         6.0
#a single Spearman correlation coefficient for each gene is generated

# Convert to matrix and assign row names
gene_correlation_matrix <- as.matrix(gene_correlation_matrix)
rownames(gene_correlation_matrix) <- common_genes

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/retracing_my_step.RData")


# View the first few correlations
head(gene_correlation_matrix)

summary(gene_correlation_matrix)
       # V1          
 # Min.   :-0.28555  
 # 1st Qu.:-0.02153  
 # Median : 0.02653  
 # Mean   : 0.03751  
 # 3rd Qu.: 0.08081  
 # Max.   : 0.82265 
##there are genes that are very similar in blood and brain 
##there are genes that are not 
# Define highly correlated genes (including both positive and negative correlations)
highly_correlated_genes <- rownames(gene_correlation_matrix)[abs(gene_correlation_matrix[, 1]) > 0.5]

# Define lowly correlated genes (close to zero)
lowly_correlated_genes <- rownames(gene_correlation_matrix)[abs(gene_correlation_matrix[, 1]) < 0.1]

# Print the number of highly correlated and weakly correlated genes
length(highly_correlated_genes)  # 131 Genes with strong positive or negative correlation
length(lowly_correlated_genes)   # 13654 Genes with weak correlation

##show some version of GO here

pdf("/hpc/users/hoangd02/www/plots/Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples_20250308.pdf")
hist(abs(gene_correlation_matrix[, 1]), breaks = 50, col = "blue",
     main = "Distribution of Absolute Gene Correlations",
     xlab = "Absolute Pearson Correlation")
dev.off()

#old versions
https://hoangd02.u.hpc.mssm.edu/plots/Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples.pdf


############AFTER RESIDUALIZING FOR COVARIATES
# gene_correlation_matrix <- mapply(function(x, y) cor(x, y, method = "spearman"), 
#                                   as.data.frame(t(blood_subset_resid)), #transpose so that rows = samples and columns = genes
#                                   as.data.frame(t(brain_subset_resid)))
# summary(gene_correlation_matrix)
# #pearson
# # Min.   :-0.38513  
# # 1st Qu.:-0.01984  
# # Median : 0.03498  
# # Mean   : 0.05183  
# # 3rd Qu.: 0.09425  
# # Max.   : 0.87528 

# #spearman
# #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # -0.38759 -0.02033  0.03356  0.04893  0.09288  0.84020

# length(gene_correlation_matrix)
# #17533

# # Define highly correlated genes (including both positive and negative correlations)
# highly_correlated_genes <- rownames(gene_correlation_matrix)[abs(gene_correlation_matrix[, 1]) > 0.5]
# lowly_correlated_genes <- rownames(gene_correlation_matrix)[abs(gene_correlation_matrix[, 1]) < 0.1]

# length(highly_correlated_genes)  # 141 Genes with strong positive or negative correlation
# length(lowly_correlated_genes)   # 14289 Genes with weak correlation

# pdf("/hpc/users/hoangd02/www/plots/Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples_residualized_20250308.pdf")
# #hist(abs(gene_correlation_matrix[, 1]), breaks = 50, col = "blue",
# hist(gene_correlation_matrix[, 1], breaks = 50, col = "blue",
#      main = "Distribution of Absolute Gene Correlations",
#      xlab = "Absolute Pearson Correlation")
# dev.off()

# pdf("/hpc/users/hoangd02/www/plots/Distribution_of_Gene_Correlations_with_233pairedsamples_residualized_spearman_with_direction_20250308.pdf")
# #hist(abs(as.numeric(gene_correlation_matrix)), breaks = 50, col = "blue",
# hist(as.numeric(gene_correlation_matrix), breaks = 50, col = "blue",
#      main = "Distribution of Gene Correlations",
#      xlab = "Spearman Correlation")
# dev.off()


# #old versions
# https://hoangd02.u.hpc.mssm.edu/plots/Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples.pdf
# Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples_residualized_20250308.pdf (Pearson)
# Distribution_of_Absolute_Gene_Correlations_with_233pairedsamples_residualized_spearman_20250308.pdf


########### SOME EQTL
#MAIN GTEX Paper: https://www.science.org/doi/10.1126/science.aaz1776
#downloaded GTEx_Analysis_v8_eQTL.tar from https://www.gtexportal.org/home/downloads/adult-gtex/qtl
https://genetics.opentargets.org/variant/7_105640741_G_A

#load data
eqtl_brain <- read.table("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/Brain_Frontal_Cortex_BA9.v8.egenes.txt", header = TRUE, sep = "\t")
eqtl_blood <- read.table("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/Whole_Blood.v8.egenes.txt", header = TRUE, sep = "\t")
dim(eqtl_blood)
#20315    33
dim(eqtl_brain)
#24676    33

##all the gene IDs are unique
length(unique(eqtl_brain$gene_id)) #24676
length(unique(eqtl_blood$gene_id)) #20315
##top SNPs for each gene and for all genes that are expressed- confirm!

##the raw datasets contain all genes with multiple SNPs regulating that gene

# Only genes present in BOTH will be retained
eqtl_merged <- merge(eqtl_brain, eqtl_blood, by = "gene_id", suffixes = c("_brain", "_blood"))
dim(eqtl_merged)  
#18280    65

##i deduce all shared genes but with multiple SNPs 

##same 

variant_pos_blood and variant_pos_brain
101217027 101558812

head(eqtl_merged)
             gene_id gene_name_brain gene_chr_brain gene_start_brain
1 ENSG00000000003.14          TSPAN6           chrX        100627109
2 ENSG00000000419.12            DPM1          chr20         50934867
3 ENSG00000000457.13           SCYL3           chr1        169849631
4 ENSG00000000460.16        C1orf112           chr1        169662007
5 ENSG00000000938.12             FGR           chr1         27612064
6 ENSG00000000971.15             CFH           chr1        196651878
  gene_end_brain strand_brain num_var_brain beta_shape1_brain beta_shape2_brain
1      100639991            -          4243           1.03647           561.422
2       50958555            -          8010           1.03573          1046.030
3      169894267            -          8208           1.04586           584.843
4      169854080            +          8572           1.06801           595.303
5       27635277            -          4887           1.08925           446.391
6      196747504            +          6058           1.02143           475.943
  true_df_brain pval_true_df_brain        variant_id_brain tss_distance_brain
1       127.406        0.000306625  chrX_101558812_G_T_b38             918821
2       121.118        0.000164714  chr20_50150782_G_A_b38            -807773
3       120.089        0.000158221  chr1_169796199_G_A_b38             -98068
4       118.679        0.001888670 chr1_170368611_T_TG_b38             706604
5       118.415        0.000990476   chr1_27340744_T_C_b38            -294533
6       124.216        0.000158904  chr1_196150180_A_G_b38            -501698

sub2<- eqtl_merged[,c("gene_id","variant_id_brain","variant_id_blood")]

sub_eqtl_merged ## IS THE CORRECT ONE THAT REPRESENT A GENE - SNP PAIRING

sub_eqtl_merged<-eqtl_merged[eqtl_merged$variant_id_brain == eqtl_merged$variant_id_blood,]
sub<- sub_eqtl_merged[,c("gene_id","variant_id_brain","variant_id_blood")]
dim(sub_eqtl_merged) #422 genes that share eQTL

#Add an Indicator for Whether a Gene Has an eQTL
# Convert gene_correlation_matrix to a proper data frame
gene_correlation_df <- data.frame(
  Gene = rownames(gene_correlation_matrix),  # Extract row names correctly
  Correlation = as.numeric(gene_correlation_matrix[,1])  # Convert the first column to numeric values
)

# Add eQTL presence indicator (1 = has eQTL, 0 = no eQTL)
#gene_correlation_df$has_eQTL <- ifelse(gene_correlation_df$Gene %in% eqtl_merged$gene_id, "Yes", "No")
gene_correlation_df$has_eQTL <- ifelse(gene_correlation_df$Gene %in% sub_eqtl_merged$gene_id, "Yes", "No")

# Check how many genes have eQTLs
table(gene_correlation_df$has_eQTL)
#    No   Yes 
# 17385   148 
#0.9915588 -- 99%

library(ggplot2)

has_eqtl_plot <-ggplot(gene_correlation_df, aes(x = Correlation, fill = has_eQTL)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("Yes" = "red", "No" = "blue")) +
  labs(title = "Distribution of Gene Correlations with eQTL Overlay",
       x = "Spearman Correlation",
       y = "Count",
       fill = "Has eQTL?") +
  theme_bw()
has_eqtl_plot

ggsave(filename = "/hpc/users/hoangd02/www/plots/gene_correlation_with_eqtl_plot.pdf",
       plot = has_eqtl_plot, 
       width = 7, height = 7, dpi = 300)

### Modify the Histogram with Matching Effect Sizes
# Compute Effect Size Agreement:
# If slope_brain ≈ slope_blood, the SNP has a similar genetic effect on both tissues.
# If slope_brain and slope_blood have opposite signs, the eQTL has discordant effects.

#Categorize Genes by Effect Similarity:
sub_eqtl_merged$Effect_Agreement <- ifelse(
  sign(sub_eqtl_merged$slope_brain) == sign(sub_eqtl_merged$slope_blood), "Same Direction", "Opposite Direction")

table(sub_eqtl_merged$Effect_Agreement)
# Opposite Direction     Same Direction 
#               15                407  

#Merge with Correlation Data
gene_correlation_df <- merge(
  gene_correlation_df, sub_eqtl_merged[, c("gene_id", "Effect_Agreement")],
  by.x = "Gene", by.y = "gene_id", all.x = TRUE
)
gene_correlation_df$Effect_Agreement[is.na(gene_correlation_df$Effect_Agreement)] <- "No eQTL"

table(gene_correlation_df$Effect_Agreement)
# No eQTL     Opposite Direction     Same Direction 
# 17385                  3                145

#Plot Histogram with Effect Agreement:
# plot<-ggplot(gene_correlation_df, aes(x = Correlation, fill = Effect_Agreement)) +
#   geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
#   scale_fill_manual(values = c("Same Direction" = "red", "Opposite Direction" = "blue", "No eQTL" = "gray")) +
#   labs(title = "Distribution of Gene Correlations with eQTL Effect Agreement",
#        x = "Spearman Correlation",
#        y = "Count",
#        fill = "eQTL Effect") +
#   theme_bw()

plot <- ggplot(gene_correlation_df, aes(x = Correlation, fill = Effect_Agreement)) +
  geom_histogram(binwidth = 0.05, color = "black") +  # Keep black borders
  facet_wrap(~Effect_Agreement, scales = "fixed") +  # Fix x and y axes
  scale_fill_manual(values = c("Same Direction" = "red", "Opposite Direction" = "blue", "No eQTL" = "gray")) +
  labs(title = "Gene Correlation Distribution by eQTL Effect Agreement",
       x = "Spearman Correlation",
       y = "Count",
       fill = "eQTL Effect Agreement") +
  theme_bw()
plot 

ggsave(filename = "/hpc/users/hoangd02/www/plots/gene_correlation_with_eqtl_direction_faceted.pdf",
       plot = plot, 
       width = 9, height = 9, dpi = 300)

plot <- ggplot(gene_correlation_df, aes(x = Effect_Agreement, y = Correlation, fill = Effect_Agreement)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Same Direction" = "red", "Opposite Direction" = "blue", "No eQTL" = "gray")) +
  labs(
    title = "Gene Correlation Distribution by eQTL Effect Direction",
    x = "eQTL Effect Agreement",
    y = "Spearman Correlation",
    fill = "Effect Agreement"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
#plot

table(gene_correlation_df$Effect_Agreement)
# No eQTL Opposite Direction     Same Direction 
# 17385                  3                145

oppo_direction <- filter(gene_correlation_df, gene_correlation_df$Effect_Agreement == "Opposite Direction")

ggsave(filename = "/hpc/users/hoangd02/www/plots/gene_correlation_with_eqtl_boxplot4.pdf",
       plot = plot, 
       width = 12, height = 7, dpi = 300)

##same directions and magnitude 
# Categorize SNPs into 4 Groups
# We'll classify each SNP into one of the following:

# Same Direction & Similar Magnitude → (slope_brain ≈ slope_blood and sign matches)
# Same Direction & Different Magnitude → (sign matches, but magnitudes differ)
# Opposite Direction & Similar Magnitude → (sign differs, but magnitudes are similar)
# Opposite Direction & Different Magnitude → (sign differs, and magnitudes differ)

# Define threshold for magnitude similarity (10% difference)
threshold <- 0.05

# Compute relative magnitude difference
#REMEMBER TO RUN THIS!!!!
#eqtl_merged$Magnitude_Diff <- abs(eqtl_merged$slope_brain - eqtl_merged$slope_blood) /
#                             ((abs(eqtl_merged$slope_brain) + abs(eqtl_merged$slope_blood)) / 2)
# This tells us how different the slopes are relative to their overall size.
# Small values (~0.0-0.1) mean the slopes are very similar in magnitude.
# Large values (>0.5-1.0) mean the slopes differ significantly.

#Categorize SNPs --- RUN LINES 808 AND 809
sub_eqtl_merged$Effect_Category <- ifelse(
  sign(sub_eqtl_merged$slope_brain) == sign(sub_eqtl_merged$slope_blood) & sub_eqtl_merged$Magnitude_Diff < threshold, "Same Direction & Similar Magnitude",
  ifelse(sign(sub_eqtl_merged$slope_brain) == sign(sub_eqtl_merged$slope_blood) & sub_eqtl_merged$Magnitude_Diff >= threshold, "Same Direction & Different Magnitude",
  ifelse(sign(sub_eqtl_merged$slope_brain) != sign(sub_eqtl_merged$slope_blood) & sub_eqtl_merged$Magnitude_Diff < threshold, "Opposite Direction & Similar Magnitude",
         "Opposite Direction & Different Magnitude")))

gene_correlation_df <- merge(
  gene_correlation_df, sub_eqtl_merged[, c("gene_id", "Effect_Category")],
  by.x = "Gene", by.y = "gene_id", all.x = TRUE
)

# Replace NA with "No eQTL"
gene_correlation_df$Effect_Category[is.na(gene_correlation_df$Effect_Category)] <- "No eQTL"


table(gene_correlation_df$Effect_Category)
#                                  No eQTL 
#                                    11872 
# Opposite Direction & Different Magnitude 
#                                     2491 
#     Same Direction & Different Magnitude 
#                                     2952 
#       Same Direction & Similar Magnitude 
#                                      218 

# Identify genes that should be "Opposite Direction & Similar Magnitude"
subset_check <- eqtl_merged[
  sign(eqtl_merged$slope_brain) != sign(eqtl_merged$slope_blood) &  # Opposite signs
  eqtl_merged$Magnitude_Diff < threshold,  # Similar magnitude
]

# Check how many genes meet this condition
nrow(subset_check) #0 

#there are truly no genes in the Opposite Direction & Similar Magnitude category

## distributions - different or not?
# Extract correlation values for each category
same_dir_sim_mag <- gene_correlation_df$Correlation[gene_correlation_df$Effect_Category == "Same Direction & Similar Magnitude"]
same_dir_diff_mag <- gene_correlation_df$Correlation[gene_correlation_df$Effect_Category == "Same Direction & Different Magnitude"]
opp_dir_diff_mag <- gene_correlation_df$Correlation[gene_correlation_df$Effect_Category == "Opposite Direction & Different Magnitude"]
no_eqtl <- gene_correlation_df$Correlation[gene_correlation_df$Effect_Category == "No eQTL"]


# Perform pairwise Kolmogorov-Smirnov tests
ks.test(same_dir_sim_mag, same_dir_diff_mag)   # Compare Same Dir & Similar Mag vs. Same Dir & Diff Mag
#   Asymptotic two-sample Kolmogorov-Smirnov test
# data:  same_dir_sim_mag and same_dir_diff_mag
# D = 0.16154, p-value = 5.005e-05
# alternative hypothesis: two-sided
ks.test(same_dir_sim_mag, opp_dir_diff_mag)    # Compare Same Dir & Similar Mag vs. Opp Dir & Diff Mag
# data:  same_dir_sim_mag and opp_dir_diff_mag
# D = 0.20076, p-value = 1.92e-07
# alternative hypothesis: two-sided
ks.test(same_dir_sim_mag, no_eqtl)             # Compare Same Dir & Similar Mag vs. No eQTL
# data:  same_dir_sim_mag and no_eqtl
# D = 0.22786, p-value = 4.44e-10
# alternative hypothesis: two-sided
# Warning message:
# In ks.test.default(same_dir_sim_mag, no_eqtl) :
#   p-value will be approximate in the presence of ties
ks.test(same_dir_diff_mag, opp_dir_diff_mag)   # Compare Same Dir & Diff Mag vs. Opp Dir & Diff Mag
#   Asymptotic two-sample Kolmogorov-Smirnov test
# data:  same_dir_diff_mag and opp_dir_diff_mag
# D = 0.061341, p-value = 7.684e-05
# alternative hypothesis: two-sided
ks.test(same_dir_diff_mag, no_eqtl)            # Compare Same Dir & Diff Mag vs. No eQTL
# D = 0.092461, p-value < 2.2e-16
ks.test(opp_dir_diff_mag, no_eqtl)             # Compare Opp Dir & Diff Mag vs. No eQTL
# D = 0.036442, p-value = 0.008433

#visualize
plot <- ggplot(gene_correlation_df, aes(x = Correlation, fill = Effect_Category)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  facet_wrap(~Effect_Category, scales = "fixed") +  # Separate facets for each category
  scale_fill_manual(values = c(
    "Same Direction & Similar Magnitude" = "red",
    "Same Direction & Different Magnitude" = "orange",
    "Opposite Direction & Similar Magnitude" = "blue",
    "Opposite Direction & Different Magnitude" = "purple",
    "No eQTL" = "gray"
  )) +
  labs(x = "Spearman Correlation",
       y = "Count",
       fill = "eQTL Effect Category") +
  theme_bw() + 
  theme(legend.position = "none",  # Remove the legend
  strip.text = element_text(size = 15),  # Increase facet text size
  axis.title = element_text(size = 18),  # Increase axis title size
  axis.text = element_text(size = 16))  # Increase axis text size


ggsave(filename = "/hpc/users/hoangd02/www/plots/gene_correlation_with_eqtl_direction_magnitude_with_no_legend.png",
       plot = plot, 
       width = 11, height = 11, dpi = 300)


#See if strongly correlated genes are functionally related
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Remove version numbers (e.g., ENSG00000005889.15 → ENSG00000005889)
highly_correlated_genes_clean <- sub("\\..*", "", highly_correlated_genes)
length(highly_correlated_genes_clean)
# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = highly_correlated_genes_clean,  # Ensembl IDs
                       column = "SYMBOL",  # Convert to gene symbols
                       keytype = "ENSEMBL",  # Input format is Ensembl
                       multiVals = "first")  # Use the first match

# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = highly_correlated_genes_clean,  
#                        column = "SYMBOL",  
#                        keytype = "ENSEMBL",  
#                        multiVals = "list")  # Returns a list of gene symbols

gene_symbols_df <- data.frame(
  Ensembl_ID = rep(names(gene_symbols), sapply(gene_symbols, length)),
  Gene_Symbol = unlist(gene_symbols)
)
dim(gene_symbols_df)
[1] 152   2

# Remove any missing values
gene_symbols_df2 <- na.omit(gene_symbols_df) #there are no NAs

# Check the first few mapped genes
head(gene_symbols) ##so few!

head(gene_symbols_df2)
# Run GO enrichment analysis with converted gene symbols
enrich_result <- enrichGO(gene = gene_symbols_df2$Gene_Symbol,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",  # Now using gene symbols
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH",
                          readable = TRUE)

# View the top enriched pathways
head(enrich_result)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

functional <- dotplot(enrich_result, showCategory = 5, title = "GO Biological Process Enrichment")
functional

ggsave(filename = "/hpc/users/hoangd02/www/plots/go_enrichment_for_corr_above_0.5.png",
       plot = functional, 
       width = 9, height = 6, dpi = 300)

#barplot(enrich_result, showCategory = 20, title = "GO Biological Process Enrichment")
# dotplot(enrich_result, showCategory = 5, title = "GO Biological Process Enrichment")

# emapplot(enrich_result, showCategory = 30)  # Adjust category number if needed
#cnetplot(enrich_result, showCategory = 10, circular = TRUE, colorEdge = TRUE)


############################# non-residualized genes results

                   ID
GO:0019882 GO:0019882
GO:0048002 GO:0048002
GO:0002399 GO:0002399
GO:0002503 GO:0002503
GO:0000028 GO:0000028
GO:0002396 GO:0002396
                                                          Description GeneRatio
GO:0019882                        antigen processing and presentation      6/45
GO:0048002     antigen processing and presentation of peptide antigen      5/45
GO:0002399                      MHC class II protein complex assembly      3/45
GO:0002503 peptide antigen assembly with MHC class II protein complex      3/45
GO:0000028                           ribosomal small subunit assembly      3/45
GO:0002396                               MHC protein complex assembly      3/45
             BgRatio       pvalue     p.adjust       qvalue
GO:0019882 108/18903 2.054759e-07 0.0001068474 8.781389e-05
GO:0048002  64/18903 4.175673e-07 0.0001085675 8.922755e-05
GO:0002399  16/18903 6.908516e-06 0.0008981070 7.381204e-04
GO:0002503  16/18903 6.908516e-06 0.0008981070 7.381204e-04
GO:0000028  19/18903 1.189459e-05 0.0010377962 8.529256e-04
GO:0002396  20/18903 1.397033e-05 0.0010377962 8.529256e-04
                                                 geneID Count
GO:0019882 KDM5D/ERAP2/HLA-DQB1/HLA-DRB5/HLA-G/HLA-DQA1     6
GO:0048002       ERAP2/HLA-DQB1/HLA-DRB5/HLA-G/HLA-DQA1     5
GO:0002399                   HLA-DQB1/HLA-DRB5/HLA-DQA1     3
GO:0002503                   HLA-DQB1/HLA-DRB5/HLA-DQA1     3
GO:0000028                      RRP7A/PWP2/LOC102724159     3
GO:0002396                   HLA-DQB1/HLA-DRB5/HLA-DQA1     3

##HLA subtypes!

############################# residualized genes results
                   ID
GO:0048002 GO:0048002
GO:0000028 GO:0000028
GO:0019882 GO:0019882
GO:0002478 GO:0002478
GO:0002399 GO:0002399
GO:0002503 GO:0002503
                                                                Description
GO:0048002           antigen processing and presentation of peptide antigen
GO:0000028                                 ribosomal small subunit assembly
GO:0019882                              antigen processing and presentation
GO:0002478 antigen processing and presentation of exogenous peptide antigen
GO:0002399                            MHC class II protein complex assembly
GO:0002503       peptide antigen assembly with MHC class II protein complex
           GeneRatio   BgRatio       pvalue     p.adjust       qvalue
GO:0048002      6/86  64/18903 4.509803e-07 0.0004189607 0.0003911661
GO:0000028      4/86  19/18903 1.469055e-06 0.0006823762 0.0006371061
GO:0019882      6/86 108/18903 9.812905e-06 0.0030387296 0.0028371347
GO:0002478      4/86  40/18903 3.220496e-05 0.0074796019 0.0069833913
GO:0002399      3/86  16/18903 4.878152e-05 0.0075530057 0.0070519253
GO:0002503      3/86  16/18903 4.878152e-05 0.0075530057 0.0070519253
                                                 geneID Count
GO:0048002 ERAP2/HLA-DQB1/HLA-DRB5/HLA-G/HLA-A/HLA-DQA1     6
GO:0000028                RRP7A/RPS28/PWP2/LOC102724159     4
GO:0019882 ERAP2/HLA-DQB1/HLA-DRB5/HLA-G/HLA-A/HLA-DQA1     6
GO:0002478             HLA-DQB1/HLA-DRB5/HLA-A/HLA-DQA1     4
GO:0002399                   HLA-DQB1/HLA-DRB5/HLA-DQA1     3
GO:0002503                   HLA-DQB1/HLA-DRB5/HLA-DQA1     3

#################### BRAIN GENIE ## FOR NON-QC'ed DATASETS ########################################
##to run BRAINGENIE, blood and brain GE matrices must have the SAME number of GENES but can have different samples

#https://github.com/hessJ/BrainGENIE?tab=readme-ov-file
#https://github.com/hessJ/BrainGENIE/wiki/How-to-run-BrainGENIE
# in bash
module load git
git ls-remote https://github.mountsinai.org/Beckmann-lab/hello.git

git clone https://github.com/hessJ/BrainGENIE
cd BrainGENIE

#download the trained models (GTEx v8) from Zenodo:
wget https://zenodo.org/record/6350240/files/normalized_expression_dat_gtexv8.tar.gz
tar -xvf normalized_expression_dat_gtexv8.tar.gz

#/hpc/users/hoangd02/BrainGENIE 

## in R
version #R version 4.2.0
setwd("/hpc/users/hoangd02/BrainGENIE")  
source("braingenie_methods.R")
library(plyr)

# --- Load blood transcriptome for new sample that you want to impute brain region transcriptome
#new_dat = data.frame(fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_cpm.txt"))
new_dat=blood_subset  
dim(new_dat)  #17537   233  
new_dat=blood_subset_resid
dim(new_dat)  #17533   233
sum(is.na(new_dat)) #0

new_dat <- as.matrix(new_dat)  # Convert data to matrix
mode(new_dat) <- "numeric"  # Ensure numeric format

# --- BrainGENIE functions
# load paired blood-brain transcriptome data for GTExv8
#load_expr_data(path_to_data = "~/BrainGENIE/normalized_expression_dat_gtexv8/Brain_Frontal_Cortex_BA9/")
#/hpc/users/hoangd02

blood_expr = blood_subset
brain_expr = brain_subset

colnames(new_dat) <- colnames(brain_expr)
identical(colnames(brain_expr),colnames(new_dat)) #TRUE

# re-train LR models using genes from new samples 
# Note: the parameter `tissue` is a character string appended to the output file for record-keeping purposes
retrain_gtex(gene_list = rownames(new_dat), output = "~/BrainGENIE/", 
             tissue = "PFC", ncomps = 40, 
             prop_for_test_set = 0.0, n_folds = 5)

#tissue = "Frontal_Cortex_BA9", ncomps = 40, #ncomps = 2,
##for gene_list use rownames not colnames LOL 

# load LR prediction performance from cross-validation
perf = load_cv_performance()
dim(perf)
#17533    12

##mean cross-val R2
summary(perf$Rsq)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01517 0.04934 0.06650 0.10009 0.47542 

mean_train_R2 <- mean(perf$train.cor^2, na.rm=TRUE)  # Mean Training R² across genes
print(paste("Mean Training R²:", mean_train_R2)) #0.492387626722538

##or
perf$train.cor^2 #Mean Training R2 
perf$Rsq #mean CV R2
perf$Cor #mean CV Pearson's r

# Mean Training R² across all genes
mean_train_R2 <- mean(perf$train.cor^2, na.rm=TRUE)

# Mean CV R² across all genes
mean_CV_R2 <- mean(perf$Rsq, na.rm=TRUE)

# Max CV R² across all genes
max_CV_R2 <- max(perf$Rsq, na.rm=TRUE)

print(paste("Mean Training R²:", mean_train_R2))
print(paste("Mean CV R²:", mean_CV_R2))
print(paste("Max CV R²:", max_CV_R2))
#LBP Blood on GTEx Brain -- WRONG
# [1] "Mean Training R²: 0.492387626722538"
# [1] "Mean CV R²: 0.066503804029847"
# [1] "Max CV R²: 0.475423036583265"

#LBP Blood on LBP Brain residualized
# [1] "Mean Training R²: 0.26144158433924"
# [1] "Mean CV R²: 0.0125430618430637"
# [1] "Max CV R²: 0.514924966039937"

#LBP Blood on LBP Brain nonresidualized
# [1] "Mean Training R²: 0.249969054551564"
# [1] "Mean CV R²: 0.0116273147230426"
# [1] "Max CV R²: 0.94731791906899"

# filter genes that were significantly predicted via cross-validation
perf = perf[(perf$Cor >= 0.1 & perf$fdr < 0.05), ]
dim(perf)
#2 PCs with non-residualized dataset: 3274   12 - not many genes that are significant.....
#5 PCs with residualized dataset: 267  12
#10 PCs with residualized dataset: 321  12
#20 PCs with residualized dataset: 783  12
#40 PCs with residualized dataset: 1335   12

dim(perf)
#40 PCs with non-residualized dataset: 1070   12

perf$
perf$gene         perf$SE           perf$FisherZ      perf$pval
perf$Cor          perf$SD           perf$sample.size  perf$fdr
perf$train.cor    perf$Rsq          perf$Zscore       perf$tissue

table(perf$sample.size)
#The total number of GTEx v8 samples for Frontal Cortex (BA9) is around 125 samples.
#only 25 samples were used for testing in each fold, each fold likely had ~25 samples for testing (125 / 5 = 25)

# Mean Training R² across all genes
mean_train_R2 <- mean(perf$train.cor^2, na.rm=TRUE)

# Mean CV R² across all genes
mean_CV_R2 <- mean(perf$Rsq, na.rm=TRUE)

# Mean Pearson's r cor
mean_pearson_r <- mean(perf$Cor, na.rm=TRUE)

print(paste("Mean Training R²:", mean_train_R2))
print(paste("Mean CV R²:", mean_CV_R2))
print(paste("Mean Pearson's r:", mean_pearson_r))

## nonresidualized trained on LBP blood ###use this for proposal 
# [1] "Mean Training R²: 0.351655694576075"
# [1] "Mean CV R²: 0.080360501414165"
# [1] "Mean Pearson's r: 0.258867415104499"

## residualized
# [1] "Mean Training R²: 0.32174272319786"
# [1] "Mean CV R²: 0.0609393033883516"
# [1] "Mean Pearson's r: 0.239787269018631"

# [1] "Mean Training R²: 0.32174272319786"
##This means that, on average, the model explains 53% of the variance in brain gene expression in the training set.

# [1] "Mean CV R²: 0.0609393033883516"
#This indicates that, on average, the model explains only 10.8% of the variance in brain gene expression in the test (validation) set.

# [1] "Mean Pearson's r: 0.239787269018631"
# This is the mean correlation between predicted and actual brain gene expression values across CV folds.
# A Pearson’s r of 0.318 suggests a weak to moderate positive correlation between predicted and actual values.


# ##mean cross-val R2
# summary(perf$Rsq)
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.03158 0.06371 0.09348 0.10804 0.13864 0.47542 

# mean_train_R2 <- mean(perf$train.cor^2, na.rm=TRUE)  # Mean Training R² across genes
# print(paste("Mean Training R²:", mean_train_R2)) #0.530431205472619

# mean_CV_r <- mean(perf$Cor, na.rm=TRUE)  # Mean across genes
# print(paste("Mean CV Pearson’s r:", mean_CV_r)) #0.210388033394538

###keeping only the predicted brain gene expression that correlates at least moderately with the actual brain expression (cor that survives multiple testing correction)

# run PCA on full GTEx data
blood.pca = fit_pca(gene_list = rownames(new_dat))

head(blood.pca$pca$x)
dim(blood.pca$pca$x) #233 233

length(blood.pca$genes) # 17533

# train LR models using full GTEx data
trained.models = fit_lr_weights_in_gtex(pca_model = blood.pca, gene_list = perf$gene, tissue = 'Frontal_Cortex_BA9', n_comps = 40)
#Training LR models for: 1335 genes

# obtain PCs in new samples based on PCA solution from GTEx
fit.pca.to.new.samples = predict_pca(dat = t(new_dat), pca_model = blood.pca)
#dat expectation: rows = samples and columns = genes so i have to transpose this cuz the original new_dat has the opposite

# run imputation step on the PCA solution in new samples based on LR model weights obtained from full GTEx
imputed = impute_gxp(pca_in_new_sample = fit.pca.to.new.samples, trained_model = trained.models, scale = TRUE)
dim(imputed) 
#233 1335

imputed_transposed <- t(imputed)
dim(imputed_transposed) #3277  233

write.table(imputed_transposed, "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/imputed_brain_GE.txt", sep="\t", quote=FALSE, row.names=TRUE)

dim(brain_subset_resid) #17533   233
common_genes = intersect(rownames(imputed_transposed), rownames(brain_subset_resid))
length(common_genes) #1335

# Subset both matrices to only include common genes
imputed_transposed_df <- imputed_transposed[common_genes, ]
dim(imputed_transposed_df)
#1335  233
brain_subset_resid <- brain_subset_resid[common_genes, ]
dim(brain_subset_resid)
#1335  233

identical(rownames(imputed_transposed_df),rownames(brain_subset_resid)) #TRUE

##correlation across samples for each gene
gene_correlation_matrix <- mapply(function(x, y) cor(x, y, method = "spearman"), 
                                  as.data.frame(t(imputed_transposed_df)), #transpose so that rows = samples and columns = genes
                                  as.data.frame(t(brain_subset_resid)))
# Convert to matrix and assign row names
gene_correlation_matrix <- as.matrix(gene_correlation_matrix)
summary(gene_correlation_matrix)
 #       V1        
 # Min.   :0.2318  
 # 1st Qu.:0.4558  
 # Median :0.4833  
 # Mean   :0.4862  
 # 3rd Qu.:0.5111  
 # Max.   :0.8098 

 gene_correlation_df <- data.frame(
  Gene = names(gene_correlation_matrix),  # Extract names instead of rownames
  Correlation = as.numeric(gene_correlation_matrix)  # Convert to numeric values
)

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/braingenie_run_20250309.RData")
##old, this has no eqtl
save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/braingenie_run_20250308.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/braingenie_run_20250308.RData")

hist(as.numeric(gene_correlation_matrix), breaks = 50, col = "blue",
     main = "Distribution of Gene Correlations",
     xlab = "Spearman Correlation")
dev.off()

######################################## NON-RESIDUALIZED
dim(brain_subset) #17537   233
common_genes = intersect(rownames(imputed_transposed), rownames(brain_subset))
length(common_genes) #883 genes

# Subset both matrices to only include common genes
imputed_transposed_df <- imputed_transposed[common_genes, ]
dim(imputed_transposed_df)
#883 233
brain_subset_df <- brain_subset[common_genes, ]
dim(brain_subset_df)
#883 233

identical(rownames(imputed_transposed_df),rownames(brain_subset_df)) #TRUE

##correlation across samples for each gene
gene_correlation_matrix <- mapply(function(x, y) cor(x, y, method = "pearson"), 
                                  as.data.frame(t(imputed_transposed_df)), #transpose so that rows = samples and columns = genes
                                  as.data.frame(t(brain_subset_df)))
# Convert to matrix and assign row names
gene_correlation_matrix <- as.matrix(gene_correlation_matrix)
summary(gene_correlation_matrix)
 #       V1          
 # Min.   :-0.20969  
 # 1st Qu.:-0.05441  
 # Median :-0.01907  
 # Mean   :-0.01621  
 # 3rd Qu.: 0.01873  
 # Max.   : 0.19763 

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/braingenie_run_20250307.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/braingenie_run_20250307.RData")


##correlation across genes for each sample
cor_results = apply(imputed_matched, 1, function(sample) cor(sample, real_brain_matched[rownames(imputed_matched) == rownames(real_brain_matched),], use = "pairwise.complete.obs"))

# Summary statistics
summary(cor_results)

library(ggplot2)

cor_df = data.frame(Sample = names(cor_results), Correlation = cor_results)

ggplot(cor_df, aes(x = Correlation)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Correlation of Imputed vs. Real Brain Expression",
       x = "Pearson Correlation",
       y = "Count")





gene_correlation_matrix <- mapply(function(x, y) cor(x, y, method = "pearson"), 
                                  as.data.frame(t(blood_subset)), #transpose so that rows = samples and columns = genes
                                  as.data.frame(t(brain_subset)))

#mapply() applies the function element-wise over corresponding columns in the two df, it takes two vectors (x and y) representing a gene's expression across samples and computes
#x = expression of a gene across samples in blood
#y = expression of the same gene across samples in brain
#for example:
x = c(5.2, 6.1, 4.9, 5.5, 5.8)  # GeneA expression in blood
y = c(5.1, 6.0, 4.8, 5.6, 5.7)  # GeneA expression in brain
Sample  GeneA (Blood) GeneA (Brain)
S1          5.2         5.1
S2          6.1         6.0
#cor(x, y, method = "pearson")

#a single Pearson correlation coefficient for each gene is generated



# View the first few correlations
head(gene_correlation_matrix)

mean(gene_correlation_matrix)
#0.02869604

summary(gene_correlation_matrix)







##i have to make sure that samples match!!!



### PER ONE INDIVIDUAL
#Compute Pearson Correlation for Each Sample Across Genes
#This means calculating the correlation between gene expression profiles of a sample across all genes
#this answer: "Do individuals have similar gene expression patterns in blood and brain?"
#identify whether global gene expression profiles are correlated between tissues

#High correlation (r ≈ 1): The overall expression profile in blood is similar to the brain for that individual.
#Low correlation (r ≈ 0): The individual's gene expression differs significantly between the tissues.

#Compute Pearson Correlation for Each Sample Across Genes
# Step 1: Find common genes and subset both datasets
common_genes <- intersect(rownames(blood_log2_cpm), rownames(brain_log2_cpm))
blood_subset <- blood_log2_cpm[common_genes, ]
brain_subset <- brain_log2_cpm[common_genes, ]

# Step 2: Ensure sample IDs match in both matrices
common_samples <- intersect(colnames(blood_subset), colnames(brain_subset))

SAMPLE_ISMMS_brain and SAMPLE_ISMMS_blood

##this is hard to do because they have different colnames!!!!

# Subset both matrices to have the same samples
blood_subset <- blood_subset[, common_samples]
brain_subset <- brain_subset[, common_samples]

# Step 3: Create an empty square matrix to store correlations
sample_correlation_matrix <- matrix(NA, nrow = length(common_samples), ncol = length(common_samples),
                                    dimnames = list(common_samples, common_samples))

# Step 4: Compute Pearson correlation for each sample across genes
for (sample in common_samples) {
    sample_correlation_matrix[sample, sample] <- cor(blood_subset[, sample], brain_subset[, sample], method = "pearson")
}

# View the first few rows
head(sample_correlation_matrix)


















































