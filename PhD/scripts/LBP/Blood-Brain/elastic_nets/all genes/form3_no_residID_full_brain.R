#code derived from OF_blood-brain_elastic_net_form3_base_brain_with_residID.R
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
library(furrr)      # for future_map
library(future)     # for plan()
library(purrr)      # for map-style helpers

set.seed(2025)

elastic_full_no_residID <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_full_brain_predicted_from_blood.txt", data.table=FALSE)
summary(elastic_full_no_residID$r2_train)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.061   0.161   0.264   0.321   0.431   0.999   11733
summary(elastic_full_no_residID$r2_test) ##############highest
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.003   0.013   0.040   0.040   0.898   11733 

elastic_net_model_form3_objects_all_full_brain_genes.rds