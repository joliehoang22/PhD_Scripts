## code derived from LBP_blood_brain_QC_20250603.rmd and all_form3_all_brain_cor.R
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
library(purrr)
library(caret)
library(tidyr)

set.seed(2025)

## Do I have to re-QC everything if i subset it to only individuals with 2 pairs? I think so. 
# Then maybe I might have more than 73 people because previously some of them were excluded for QC purposes...

## Research question: Is blood day 1 more correlated to brain day 30 than blood day 30 to brain day 30?

# let's just calculate the correlation for 1 patient (or a few) without re-QC first then decide how to proceed.

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250811.RData")
dim(dictionary) #225 x 6
head(dictionary)
sub <- filter(dictionary, IID_ISMMS == "PT-0024"); sub
table(brain_metadata$number_of_pair_brain) 

# need to pull surgery day -- this is from lbp_allBatches_Blood_CellTypePropEstimates_CorrelationWithLabs_09FEB2022_JH.R in Lora_scripts
mblood <- readRDS(file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_SurgeryDates_ClinicalLabsDates_Mapping_182BloodSamples_08FEB2022.RDS")
dim(mblood) #182 x 5
colnames(mblood)
#[1] "iid"         "LBPSEMA4_ID" "surgeryDate" "labDate"     "dateDif"

filter(mblood, LBPSEMA4_ID == "LBPSEMA4BLOOD326")

filter(mblood, LBPSEMA4_ID == "LBPSEMA4BLOOD557")

dictionary <- dictionary %>%
  left_join(mblood %>% select(LBPSEMA4_ID, surgeryDate),
            by = c("SAMPLE_ISMMS_blood" = "LBPSEMA4_ID"))

sum(is.na(dictionary$surgeryDate)) #58 (58/225 = 26%)
#not all samples have surgeryDate here -- maybe this is already filtered for samples that have labs

##try a different metadata
dop <- as.data.table(read_excel("/sc/arion/projects/psychgen/lbp/data/emr/lbp_lel2021_clinical/lel2021_anesthesia_data_pull_11APR2022.xlsx"), na=c(""))
dop <- dop[,c("subject_id","mrn","brain_surgery_date")]
dim(dop) #551   3
filter(dop,subject_id=="PT-0025")

##this is where i can get the timepoint! 
myfile <- "/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_forSEMA4_awcFormatted.tsv"
dt <- fread(myfile, na=c("na","","NA"))[grep("brain", tissue, ignore.case=T)]
dt <- dt[,c("iid","LBPSEMA4_ID","phe","tissue","timepoint")]
#can i assume that left always have go first? 
#maybe try to find the sample ID first then if can't assume that left get the earlier date?

dictionary <- dictionary %>%
  left_join(dt %>% select(LBPSEMA4_ID, timepoint),
            by = c("SAMPLE_ISMMS_brain" = "LBPSEMA4_ID"))

##switching the time points so 0 == first timepoint and 1 == second timepoint
dictionary$timepoint <- 1 - dictionary$timepoint

# 1) For each subject, compute earliest and latest surgery dates
dictionary <- dictionary %>%
  mutate(
    surgeryDate = as.Date(surgeryDate),
    timepoint   = as.integer(timepoint)
  )

dop <- dop %>%
  mutate(brain_surgery_date = as.Date(brain_surgery_date))  # convert from POSIXct to Date

# 1) For each subject, compute earliest and latest surgery dates (robust to all-NA)
dop_two <- dop %>%
  group_by(subject_id) %>%
  summarise(
    early_date = if (all(is.na(brain_surgery_date))) as.Date(NA) else min(brain_surgery_date, na.rm = TRUE),
    late_date  = if (all(is.na(brain_surgery_date))) as.Date(NA) else max(brain_surgery_date, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Join and fill ONLY for 2_pairs:
#    timepoint == 0 -> early_date; timepoint == 1 -> late_date
dictionary <- dictionary %>%
  left_join(dop_two, by = c("IID_ISMMS" = "subject_id")) %>%
  mutate(
    surgeryDate = case_when(
      number_of_pair_brain == "2_pairs" & timepoint == 0 ~ early_date,
      number_of_pair_brain == "2_pairs" & timepoint == 1 ~ late_date,
      TRUE                                                ~ surgeryDate
    )
  ) %>%
  select(-early_date, -late_date)

#there are still 22 NAs
dictionary[is.na(dictionary$surgeryDate), ]

filter(dop, subject_id == "PT-0025")
filter(dop, subject_id == "PT-0148")
##there are a few 1 pair with 2 surgery dates
filter(dop, subject_id == "PT-0058")
filter(dop, subject_id == "PT-0070")

# Keep only subjects with exactly one distinct surgery date (and drop NA dates)
single_date <- dop %>%
  filter(!is.na(brain_surgery_date)) %>%
  group_by(subject_id) %>%
  filter(n_distinct(brain_surgery_date) == 1) %>%
  summarise(surgeryDate = first(brain_surgery_date), .groups = "drop")

# Join and fill only where dictionary$surgeryDate is NA
dictionary <- dictionary %>%
  left_join(single_date, by = c("IID_ISMMS" = "subject_id"), suffix = c("", ".from_dop")) %>%
  mutate(
    surgeryDate = if_else(is.na(surgeryDate), surgeryDate.from_dop, surgeryDate)
  ) %>%
  select(-surgeryDate.from_dop)

sum(is.na(dictionary$surgeryDate)) #16

dictionary[is.na(dictionary$surgeryDate), ]
filter(dop, subject_id == "PT-0162")

##filling in the final surgery dates -- matched by timepoint; if 1 = later date; 0 = earlier date
# Build per-subject earliest and latest surgery dates for subjects with 2+ dates
multi_dates <- dop %>%
  filter(!is.na(brain_surgery_date)) %>%
  group_by(subject_id) %>%
  summarise(
    n_dates = n_distinct(brain_surgery_date),
    early   = min(brain_surgery_date),
    late    = max(brain_surgery_date),
    .groups = "drop"
  ) %>%
  filter(n_dates >= 2)

# Join and fill only remaining NAs in dictionary:
# timepoint == 0 -> early; timepoint == 1 -> late
dictionary <- dictionary %>%
  left_join(multi_dates, by = c("IID_ISMMS" = "subject_id")) %>%
  mutate(
    surgeryDate = if_else(
      is.na(surgeryDate) & timepoint == 0 & !is.na(early), early,
      if_else(is.na(surgeryDate) & timepoint == 1 & !is.na(late),  late,  surgeryDate)
    )
  ) %>%
  select(-n_dates, -early, -late)

# Check how many still missing
sum(is.na(dictionary$surgeryDate))

filter(dictionary, IID_ISMMS == "PT-0162")
filter(dictionary, IID_ISMMS == "PT-0175")

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")

##################################################
## now try some correlations 
tail(dictionary)
#    IID_ISMMS SAMPLE_ISMMS_blood SAMPLE_ISMMS_brain number_of_pair_brain
#220   PT-0192   LBPSEMA4BLOOD572   LBPSEMA4BRAIN405              2_pairs
#221   PT-0192   LBPSEMA4BLOOD515   LBPSEMA4BRAIN146              2_pairs

#222   PT-0193   LBPSEMA4BLOOD328   LBPSEMA4BRAIN010              2_pairs
#223   PT-0193   LBPSEMA4BLOOD672   LBPSEMA4BRAIN172              2_pairs

#224   PT-0198   LBPSEMA4BLOOD078   LBPSEMA4BRAIN709              2_pairs
#225   PT-0198   LBPSEMA4BLOOD379   LBPSEMA4BRAIN288              2_pairs
#    mymet_tissue_brain mymet_tissue_blood surgeryDate timepoint
#220            L_Brain            L_Blood  2018-10-18         0
#221            R_Brain            R_Blood  2018-11-26         1

#222            L_Brain            L_Blood  2018-10-18         0 ##this and the person above have the same surgery dates?
#223            R_Brain            R_Blood  2018-11-26         1

#224            L_Brain            L_Blood  2019-02-11         0
#225            R_Brain            R_Blood  2019-03-11         1


### LETS TRY IID_ISMMS == PT-0192 FIRST !
#code derived from LBP_blood_brain_QC_20250603.rmd
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1
blood_form3$V1 <- NULL
blood_form3[1:5,1:5]
dim(blood_form3) #21046   225

brain_base <- v_brain$E

pt_0192_day1_blood <- blood_form3[ , "LBPSEMA4BLOOD572", drop = FALSE]; dim(pt_0192_day1_blood) #21046     1
pt_0192_day1_brain <- brain_base[ , "LBPSEMA4BRAIN405", drop = FALSE]; dim(pt_0192_day1_brain) #21356     1

pt_0192_day39_blood <- blood_form3[ , "LBPSEMA4BLOOD515", drop = FALSE]; dim(pt_0192_day39_blood) #21046     1
pt_0192_day39_brain <- brain_base[ , "LBPSEMA4BRAIN146", drop = FALSE]; dim(pt_0192_day39_brain) #21356     1

#                  comparison spearman_correlation
#1   Day1 blood vs Day1 brain           -0.1847637
#2  Day1 blood vs Day39 brain           -0.1545536
#3 Day39 blood vs Day39 brain            0.2906551

#############
############# change this to 0193
#############
pt_0192_day1_blood <- blood_form3[ , "LBPSEMA4BLOOD328", drop = FALSE]; dim(pt_0192_day1_blood) #21046     1
pt_0192_day1_brain <- brain_base[ , "LBPSEMA4BRAIN010", drop = FALSE]; dim(pt_0192_day1_brain) #21356     1

pt_0192_day39_blood <- blood_form3[ , "LBPSEMA4BLOOD672", drop = FALSE]; dim(pt_0192_day39_blood) #21046     1
pt_0192_day39_brain <- brain_base[ , "LBPSEMA4BRAIN172", drop = FALSE]; dim(pt_0192_day39_brain) #21356     1

#                  comparison spearman_correlation
#1   Day1 blood vs Day1 brain           0.22194969
#2  Day1 blood vs Day39 brain           0.26809810
#3 Day39 blood vs Day39 brain          -0.04880653

#############
############# change this to 0198
#############

#                  comparison spearman_correlation
#1   Day1 blood vs Day1 brain           0.06444818
#2  Day1 blood vs Day39 brain           0.05691988
#3 Day39 blood vs Day39 brain           0.26610617

pt_0192_day1_blood <- blood_form3[ , "LBPSEMA4BLOOD078", drop = FALSE]; dim(pt_0192_day1_blood) #21046     1
pt_0192_day1_brain <- brain_base[ , "LBPSEMA4BRAIN709", drop = FALSE]; dim(pt_0192_day1_brain) #21356     1

pt_0192_day39_blood <- blood_form3[ , "LBPSEMA4BLOOD379", drop = FALSE]; dim(pt_0192_day39_blood) #21046     1
pt_0192_day39_brain <- brain_base[ , "LBPSEMA4BRAIN288", drop = FALSE]; dim(pt_0192_day39_brain) #21356     1

# Match genes for Day 1 blood vs Day 1 brain
common_genes_day1 <- intersect(rownames(pt_0192_day1_blood), rownames(pt_0192_day1_brain))
cor_day1 <- cor(
  pt_0192_day1_blood[common_genes_day1, , drop = TRUE],
  pt_0192_day1_brain[common_genes_day1, , drop = TRUE],
  method = "spearman"
)

# Match genes for Day 1 blood vs Day 39 brain
common_genes_day39 <- intersect(rownames(pt_0192_day1_blood), rownames(pt_0192_day39_brain))
cor_day39 <- cor(
  pt_0192_day1_blood[common_genes_day39, , drop = TRUE],
  pt_0192_day39_brain[common_genes_day39, , drop = TRUE],
  method = "spearman"
)

# Match genes for Day 39 blood vs Day 39 brain
common_genes_day39_pair <- intersect(rownames(pt_0192_day39_blood), rownames(pt_0192_day39_brain))
cor_day39_pair <- cor(
  pt_0192_day39_blood[common_genes_day39_pair, , drop = TRUE],
  pt_0192_day39_brain[common_genes_day39_pair, , drop = TRUE],
  method = "spearman"
)

# Combine into one table
results <- data.frame(
  comparison = c(
    "Day1 blood vs Day1 brain",
    "Day1 blood vs Day39 brain",
    "Day39 blood vs Day39 brain"
  ),
  spearman_correlation = c(cor_day1, cor_day39, cor_day39_pair)
)
results


### BLOOD FORM 3 AND BASELINE BRAIN 
cor_matrix <- cor(t(blood_form3), t(v_brain$E), method = "spearman")
dim(cor_matrix)



##################################################
######### NO INDIVIDUAL | Brain baseline ######### 
##################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1; blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2483944
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2906922
