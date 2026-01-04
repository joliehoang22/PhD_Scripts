## code derived from LBP_blood_brain_QC_20250729.R

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
library(furrr)      
library(future)    
library(purrr)      
set.seed(2025)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250729.RData")
head(dictionary)

table(dictionary$IID_ISMMS)

no_go <- c("PT-0003","PT-0007","PT-0008","PT-0014","PT-0050","PT-0066",
            "PT-0067", "PT-0075","PT-0077","PT-0107","PT-0129","PT-0137",
            "PT-0151","PT-0168","PT-0176","PT-0216","PT-0259","PT-0417","PT-0428")

#good to go! no one has that ID 
