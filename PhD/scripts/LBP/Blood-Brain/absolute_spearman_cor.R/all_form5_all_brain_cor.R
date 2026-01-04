## code derived from LBP_blood_brain_QC_20250603.rmd
# for empirical p values see OG_null_with_empirical_p.R

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
##################################################
######### NO INDIVIDUAL | Brain baseline ######### 
##################################################
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1; blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2384292
blood_form5_99.995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.99995); blood_form5_99.995 #0.277323 

##########################################################
######### Blood WITH INDIVIDUAL | Brain baseline ######### 
##########################################################
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1; blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.226457
blood_form5_99.995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.99995); blood_form5_99.995 #0.2627413

##################################################
######### NO INDIVIDUAL | Brain full ############# 
##################################################
blood_form5_by_brain_cor  <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1; blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2337137
blood_form5_99.995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.99995); blood_form5_99.995 #0.2721019

##########################################################
######### Blood WITH INDIVIDUAL | Brain full ############# 
##########################################################
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1; blood_form5_by_brain_cor$V1 <- NULL
blood_form5_by_brain_cor <- abs(blood_form5_by_brain_cor)
blood_form5_99.95 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.9995); blood_form5_99.95 #0.2304362
blood_form5_99.995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.99995); blood_form5_99.995 #0.2677602

##########################################################
##########################################################
form5_full_brain_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_full_brain_no_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form5_full_brain_no_residID$r2_test, 0.9995, na.rm = TRUE) #0.8347526

form5_full_brain_with_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_full_brain_with_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form5_full_brain_with_residID$r2_test, 0.9995, na.rm = TRUE) #0.480376

form5_base_brain_with_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_base_brain_with_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form5_base_brain_with_residID$r2_test, 0.9995, na.rm = TRUE) #0.7753061 

form5_base_brain_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form5_base_brain_no_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form5_base_brain_no_residID$r2_test, 0.9995, na.rm = TRUE) #0.9140312

