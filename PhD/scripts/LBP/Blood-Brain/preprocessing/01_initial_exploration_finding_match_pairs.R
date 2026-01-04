##code retrieve from LBP_clean_blood_brain.ipynb

rm(list=ls())
library(data.table)
library(ggplot2)
library(readxl)
library(biomaRt)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(matrixStats)

#We also set the random seed for reproducibility:
set.seed(2025)

## Metadata, which has everything 
metadata <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")          
dim(metadata) #530 170 
metadata[1:5, 1:10] #there are 6 ID variables
#uniqueN(metadata[,.(IID_ISMMS)]) #there are 172 people total 
length(unique(metadata$IID_ISMMS)) #there are 172 people total 

## just checking whether blood and brain samples pairing is correct
sub <- metadata[,c("IID_ISMMS","SAMPLE_ISMMS","LIMS_SEMA4")] #"LBPSEMA4BLOOD395_0"  in "LIMS_SEMA4"
sub_blood <- sub[grepl("BLOOD", SAMPLE_ISMMS)] #243  3
sub_brain <- sub[grepl("BRAIN", SAMPLE_ISMMS)] #287  3

a <- merge(sub_blood, sub_brain, by="IID_ISMMS") #431   5

## based on previous deep dive in the sample space in blood-brain LBP_clean.sh, I've identified a mislabel; 
#this sample will later be dropped, see Sample Space Exploration
## there was one sample (PT-0117) that was mismatched (they have L_Brain and R_Blood so I'll changed it to L_Blood)
metadata[metadata$IID_ISMMS == "PT-0117", c("mymet_tissue", "mymet_brain", "mymet_timepoint")]
idx <- which(metadata$IID_ISMMS == "PT-0117" & metadata$mymet_tissue == "R_Blood")
metadata$mymet_tissue[idx] <- "L_Blood"
metadata$mymet_timepoint[idx] <- "left"
## Check again to see if the mislabeling is fixed
metadata[metadata$IID_ISMMS == "PT-0117", c("mymet_tissue", "mymet_brain", "mymet_timepoint")]

## How many blood and brain samples do I have?
table(metadata$mymet_brain)
# blood brain 
#   243   287 --> sum = 530
table(metadata$mymet_tissue)
# L_Blood L_Brain R_Blood R_Brain 
#     129     157     114     130 

brain_metadata <- metadata[metadata$mymet_brain == "brain", ]
dim(brain_metadata) #287 x 170 columns
#uniqueN(brain_metadata[,.(IID_ISMMS)]) #171 people with brain samples | this means that one person has only blood sample then 
length(unique(brain_metadata$IID_ISMMS)) #171 people with brain samples

blood_metadata <- metadata[metadata$mymet_brain == "blood", ]
dim(blood_metadata) #243 x 170
#unique(blood_metadata[,.(IID_ISMMS)]) #155 people with blood samples

blood_metadata <- as.data.frame(blood_metadata)
rownames(blood_metadata) <- blood_metadata$SAMPLE_ISMMS

brain_metadata <- as.data.frame(brain_metadata)
rownames(brain_metadata) <- brain_metadata$SAMPLE_ISMMS

## Common samples, which have individuals who have BOTH blood and brain samples
common_samples <- merge(brain_metadata, blood_metadata, by = "IID_ISMMS")
dim(common_samples) #431 339
#uniqueN(common_samples[,.(IID_ISMMS)])  #there are 154 people with both blood and brain samples
length(unique(common_samples$IID_ISMMS))

colnames(common_samples) <- gsub("\\.x$", "_brain", colnames(common_samples))
colnames(common_samples) <- gsub("\\.y$", "_blood", colnames(common_samples))

tbl <- table(common_samples$IID_ISMMS)
sum(tbl == 1)  # number of IDs with 1 sample = 37
sum(tbl == 2)  # number of IDs with 2 samples = 37
sum(tbl == 4)  # number of IDs with 4 samples = 73 x4 | 80

## Deep dive into sample space
############################################
############################################
one_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 1)  
one_sample <- as.data.frame(one_sample)  
dim(one_sample) #37 339 | 37 people with only 1 sample

one_sample[,c("IID_ISMMS","mymet_tissue_blood","mymet_tissue_brain")]
table(one_sample$mymet_tissue_blood)
table(one_sample$mymet_tissue_brain)
#L_Blood L_Brain R_Blood R_Brain 
#      0      26       0      11 

##check whether mymet_tissue_brain and mymet_tissue_blood start with the same letter (either "R" or "L") for the same IID_ISMMS in one_sample
one_sample$same_prefix <- substr(one_sample$mymet_tissue_brain, 1, 1) == substr(one_sample$mymet_tissue_blood, 1, 1)

# Count TRUE vs FALSE cases
table(one_sample$same_prefix) ## ALL TRUE! Good.

one_sample$IID_ISMMS #PT-0019, 25, 140, 196

############################################
############################################
two_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 2)  
two_sample <- as.data.frame(two_sample)  
dim(two_sample) #74 339
table(two_sample$IID_ISMMS) #37 people with 2 samples | PT-0020, 32, 163

### 8 INDIVIDUALS WITH 1 BRAIN SAMPLE AND 2 BLOOD SAMPLES
# Count distinct tissue values per ID
id_counts <- aggregate(mymet_tissue_brain ~ IID_ISMMS, data = two_sample, function(x) length(unique(x)))
# Get IDs where mymet_tissue_brain has only one unique value
valid_ids <- id_counts$IID_ISMMS[id_counts$mymet_tissue_brain == 1]

# Subset original data
one_brain_two_blood_samples <- two_sample[two_sample$IID_ISMMS %in% valid_ids, ]
# make sure the L_Brain-L_Blood and R_Brain-R_Blood pairing is accurate:
paired_one_brain_two_blood_samples <- one_brain_two_blood_samples %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))

table(paired_one_brain_two_blood_samples$IID_ISMMS) #8 people
# PT-0032 PT-0046 PT-0054 PT-0071 PT-0100 PT-0116 PT-0124 PT-0177 
#       2       2       2       2       2       2       2       2

##investigation in pairing
sub <- subset(common_samples, IID_ISMMS == "PT-0032")
sub[,c("IID_ISMMS","mymet_tissue_blood","mymet_tissue_brain","RNASeqMetrics_MEDIAN_5PRIME_BIAS_blood","RNASeqMetrics_MEDIAN_5PRIME_BIAS_brain")]
dictionary2 <- common_samples[,c("IID_ISMMS","SAMPLE_ISMMS_blood","SAMPLE_ISMMS_brain","mymet_tissue_brain","mymet_tissue_blood")] #can also do this with blood ID 
"number_of_pair_brain"

dictionary2_matched <- dictionary2[dictionary2$mymet_tissue_brain == gsub("_Blood", "_Brain", dictionary2$mymet_tissue_blood), ]

#   IID_ISMMS mymet_tissue_blood mymet_tissue_brain
#45   PT-0032            L_Blood            L_Brain
#46   PT-0032            R_Blood            L_Brain
#   RNASeqMetrics_MEDIAN_5PRIME_BIAS_blood
#45                               1.025763
#46                               0.939081
#   RNASeqMetrics_MEDIAN_5PRIME_BIAS_brain
#45                               0.965238
#46                               0.965238

table(paired_one_brain_two_blood_samples$mymet_tissue_brain)
#L_Blood L_Brain R_Blood R_Brain 
#      0       6       0       2

### 29 INDIVIDUALS WITH 2 BRAIN SAMPLES AND 1 BLOOD SAMPLE
# Count distinct tissue values per ID
id_counts <- aggregate(mymet_tissue_brain ~ IID_ISMMS, data = two_sample, function(x) length(unique(x)))

# Get IDs where mymet_tissue_brain has only one unique value
valid_ids <- id_counts$IID_ISMMS[id_counts$mymet_tissue_brain == 2]

# Subset original data
two_brain_one_blood_sample <- two_sample[two_sample$IID_ISMMS %in% valid_ids, ]
# make sure the L_Brain-L_Blood and R_Brain-R_Blood pairing is accurate:
paired_two_brain_one_blood_sample <- two_brain_one_blood_sample %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))

table(two_brain_one_blood_sample$IID_ISMMS) #29 people
# PT-0020 PT-0056 PT-0084 PT-0088 PT-0091 PT-0094 PT-0102 PT-0109 PT-0111 PT-0115 
#       2       2       2       2       2       2       2       2       2       2 
# PT-0118 PT-0119 PT-0120 PT-0121 PT-0136 PT-0138 PT-0139 PT-0143 PT-0145 PT-0147 
#       2       2       2       2       2       2       2       2       2       2 
# PT-0154 PT-0157 PT-0161 PT-0163 PT-0172 PT-0174 PT-0175 PT-0179 PT-0184 
#       2       2       2       2       2       2       2       2       2
sub <- subset(common_samples, IID_ISMMS == "PT-0020")
sub[,c("IID_ISMMS","mymet_tissue_blood","mymet_tissue_brain","RNASeqMetrics_MEDIAN_5PRIME_BIAS_blood","RNASeqMetrics_MEDIAN_5PRIME_BIAS_brain")]
#  IID_ISMMS mymet_tissue_blood mymet_tissue_brain
#6   PT-0020            L_Blood            L_Brain
#7   PT-0020            L_Blood            R_Brain
#  RNASeqMetrics_MEDIAN_5PRIME_BIAS_blood RNASeqMetrics_MEDIAN_5PRIME_BIAS_brain
#6                               0.949469                               0.712104
#7                               0.949469                               0.599933

table(two_brain_one_blood_sample$mymet_tissue_blood)
#L_Blood L_Brain R_Blood R_Brain 
#     28       0      30       0

############################################
############################################
four_sample <- common_samples %>%
  group_by(IID_ISMMS) %>%
  filter(n() == 4)  
four_sample <- as.data.frame(four_sample)  

four_sample[1:5,1:5]

dim(four_sample) #320 339
table(four_sample$IID_ISMMS) #80 people 

##investigation in pairing
test2 <- subset(metadata, IID_ISMMS == "PT-0198")
test2[1:5,1:5]

test<- subset(four_sample, IID_ISMMS == "PT-0198")
test[1:5,1:5]

test[,c("IID_ISMMS","SAMPLE_ISMMS.x","mymet_extractiondate.x","mymet_tissue.x","SAMPLE_ISMMS.y","mymet_extractiondate.y","mymet_tissue.y")]

#LBPSEMA4BRAIN709 and LBPSEMA4BLOOD078 left
#LBPSEMA4BRAIN288 and LBPSEMA4BLOOD379 right

LBPSEMA4BRAIN709   PT-0198   LBPSEMA4BLOOD078   LBPSEMA4BRAIN709
LBPSEMA4BRAIN288   PT-0198   LBPSEMA4BLOOD379   LBPSEMA4BRAIN288

# make sure the L_Brain-L_Blood and R_Brain-R_Blood pairing is accurate:
paired_four_sample <- four_sample %>%
  filter((mymet_tissue_brain == "L_Brain" & mymet_tissue_blood == "L_Blood") |
         (mymet_tissue_brain == "R_Brain" & mymet_tissue_blood == "R_Blood"))
nrow(paired_four_sample) #160

one_sample$same_prefix <- NULL
two_sample$same_prefix <- NULL

############################################
############################################
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

brain_metadata <-common_samples_without_dup[common_samples_without_dup$mymet_brain_brain == "brain", ] 
blood_metadata <-common_samples_without_dup[common_samples_without_dup$mymet_brain_blood == "blood", ] 

identical(brain_metadata$SAMPLE_ISMMS_brain, blood_metadata$SAMPLE_ISMMS_brain) #TRUE
identical(brain_metadata$SAMPLE_ISMMS_blood, blood_metadata$SAMPLE_ISMMS_blood) #TRUE

##good sanity check here: 
dim(blood_metadata) #233 339
dim(brain_metadata) #233 339
identical(brain_metadata, blood_metadata) #TRUE | this is because I the columns are "duplicated" with underscore
#_blood and _brain

length(unique(blood_metadata$IID_ISMMS)) #153
length(unique(brain_metadata$IID_ISMMS)) #153

########################################################
################### Expression data ####################
########################################################
raw_count <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_530LIVINGsamples_15FEB2022.RDS")                                                                                   
dim(raw_count) ##58929 x 530  

blood_ge <- raw_count[, colnames(raw_count) %in% blood_metadata$SAMPLE_ISMMS_blood]
dim(blood_ge) #58929 x 233

brain_ge <- raw_count[, colnames(raw_count) %in% brain_metadata$SAMPLE_ISMMS_brain]
dim(brain_ge) #58929 x 233

identical(colnames(blood_ge), blood_metadata$SAMPLE_ISMMS_blood) #FALSE
# Reorder columns of brain_ge to match the order of brain_metadata$SAMPLE_ISMMS_brain
blood_ge <- blood_ge[, match(blood_metadata$SAMPLE_ISMMS_blood, colnames(blood_ge))]
identical(colnames(blood_ge), blood_metadata$SAMPLE_ISMMS_blood) #TRUE
rownames(blood_metadata) <- blood_metadata$SAMPLE_ISMMS_blood
identical(colnames(blood_ge),rownames(blood_metadata)) #TRUE

identical(colnames(brain_ge), brain_metadata$SAMPLE_ISMMS_brain) #FALSE
# Reorder columns of brain_ge to match the order of brain_metadata$SAMPLE_ISMMS_brain
brain_ge <- brain_ge[, match(brain_metadata$SAMPLE_ISMMS_brain, colnames(brain_ge))]
identical(colnames(brain_ge), brain_metadata$SAMPLE_ISMMS_brain) #TRUE
rownames(brain_metadata) <- brain_metadata$SAMPLE_ISMMS_brain
identical(rownames(brain_metadata), colnames(brain_ge)) #TRUE

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_20250522.RData")


