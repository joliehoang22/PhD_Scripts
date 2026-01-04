###SVA REPRODUCIBILITY
library(variancePartition)
library(preprocessCore)
library(BiocParallel)
library(assertthat)
library(data.table)
library(tidyverse)
library(seriation)
library(GEOquery)
library(reshape2)
library(stringr)
library(ggplot2)
library(foreach)
library(limma)
library(edgeR)
library(dplyr)
library(readr)
library(sva)

##GENERAL STRUCTURE: 
#1. LOAD AND FORMAT METADATA AND GENE EXPRESSION DATA
#2. VOOM AND LIMMA ANALYSES PREPARATION AND APPLICATION TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
#4. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "BE" AND "LEEK" METHODS
#5. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
#6. SUMMARIZE AND SAVE DEA WITH SVA RESULTS

##PAY ATTENTION TO CONTRAST: makeContrasts("disease"-"control", levels=design_sv)

##GO ANALYSES 
#MSBB code for Jolie 05082024.sh

#put separate codes --> lab github --> clone and publicize
#write functions for the same analyses needed (DEA, SVA)


##################################################
##########ALZHEIMER'S DISEASE - ROSMAP############
##################################################
#we also explored 2 subsets of ROSMAP code and do within dataset DEA reproducibility
#ROSMAP halves code 05272024.sh ???

##ROSMAP code for publication 05082024.sh 
#1 question: 
#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
results$X<-rownames(results)
names(results)[1] <- "logFC_j"

#######might have some missing results line here!!!
write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

###use rosmap_no_sva_1.11.txt NOT rosmap_no_sva_1.16.txt 
#apparently the logFC signs are flipped in these two, probably due to makeContrasts(ceradsc1-ceradsc4, levels=design_sv), but the correlation between rosmap_no_sva_1.11.txt and msbb_no_sva_1.11.txt is POSITIVE 0.46 and rosmap_no_sva_1.16.txt and msbb_no_sva_1.11.txt is NEGATIVE 0.46 

###use rosmap_sva_be_1.16.txt NOT rosmap_sva_be_1.11.txt, cor btw rosmap_sva_be_1.16.txt and msbb_sva_be_1.11.txt is POSITIVE 0.12 and rosmap_sva_be_1.11.txt and msbb_sva_be_1.11.txt is NEGATIVE 0.12 and 


##################################################
##########ALZHEIMER'S DISEASE - MSBB##############
##################################################

##MSBB code for publication 05082024.sh
#there are no makeContrast like ROSMAP because I used PlaqueMean - not setting contrast by disease-HC
status="PlaqueMean"
results<-topTable(fit2,coef="PlaqueMean",number=Inf) # Extract top results

#1 question
fit2 <- eBayes(fit) 
results<-topTable(fit2,coef="PlaqueMean",number=Inf) 
##might have some missing results line here!!!

##################################################
############SCHIZOPHRENIA - CMC-MPP###############
##################################################

##################################################
############SCHIZOPHRENIA - CMC-HBCC##############
##################################################

#SZ Code.R is the same as cmc_SZ_psychencode_sva.R, use SZ code for publication.R


##################################################
#########BIPOLAR DISEASE - CMC-HBCC###############
##################################################

##################################################
#########BIPOLAR DISEASE - CMC-MSSM###############
##################################################

#cmc_BP_psychencode_sva.R, use BP code for publication.R

##################################################
#####AUTISM SPECTRUM DISORDER - UCLA-ASD##########
##################################################

##################################################
#####AUTISM SPECTRUM DISORDER - YALE-ASD##########
##################################################

##ASD code for publication 09242024.sh - DONE

##################################################
########PARKINSON'S DISEASE - GSE8397#############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE8397#############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE7621#############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE24378############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE20292############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE20141############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE20163############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE20164############
##################################################

##################################################
########PARKINSON'S DISEASE - GSE49036############
##################################################

#GEO.R

##################################################
##################################################
#LOAD ALL DATASETS FOR CORRELATIONS AND PLOT 

###### ALZHEIMER'S DISEASE ######
#################################

#in SVA Part 2 DEA.sh 
library(ggplot2)
library(readr)
library(data.table)
library(patchwork)

msbb_no_sva <-read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_no_sva_1.11.txt")
rosmap_no_sva <-read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt") 

msbb_sva_be<-read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_be_1.11.txt") 

rosmap_sva_be <-read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_1.16.txt") 

ad_no_sva_combo<-merge(msbb_no_sva,rosmap_no_sva,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_msbb','_rosmap'))
#cor.test(no_sva_combo$logFC_j_msbb, no_sva_combo$logFC_j_rosmap, method = "spearman", exact=FALSE)
#r=0.4603241, p-value < 2.2e-16

ad_sva_be_combo<-merge(msbb_sva_be,rosmap_sva_be,by='X', all.x=TRUE, all.y=TRUE,suffixes=c('_msbb','_rosmap'))
#cor.test(sva_be_combo$logFC_msbb, sva_be_combo$logFC_rosmap, method = "spearman")  #r=0.1157661, p-value < 2.2e-16

########## SCHIZOPHRENIA ########
#################################
sz_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_no_sva.txt",data.table=FALSE)#
sz_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_no_sva.txt",data.table=FALSE)#

sz_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_sva_be.txt",data.table=FALSE)# 
sz_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_sva_be.txt",data.table=FALSE)#

sz_no_sva_combo<-merge(sz_mssm_no_sva, sz_hbcc_no_sva,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_mssm','_hbcc')) #19086 genes
#cor.test(no_sva_combo$logFC_mssm, no_sva_combo$logFC_hbcc) #na.omit = default; r=0.382888;  p-value < 2.2e-16

sz_sva_be_combo<-merge(sz_mssm_sva_be, sz_hbcc_sva_be,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_mssm','_hbcc')) #19086 genes
#cor.test(sva_be_combo$logFC_mssm, sva_be_combo$logFC_hbcc) #na.omit = default; r=0.1255918;  p-value < 2.2e-16

#further plotting code is in SZ code for publication.R

######## BIPOLAR DISEASE ########
#################################
bp_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp.txt",data.table=FALSE)#
bp_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_no_sva_bp.txt",data.table=FALSE)#

bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp.txt",data.table=FALSE)# 
bp_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_sva_be_bp.txt",data.table=FALSE)#

bp_no_sva_combo<-merge(bp_hbcc_no_sva, bp_mssm_no_sva,by='X', all.x=TRUE, all.y=TRUE,suffixes=c('_hbcc','_mssm'))
#rho=0.2960745, p-value < 2.2e-16

bp_sva_be_combo<-merge(bp_hbcc_sva_be, bp_mssm_sva_be,by='X', all.x=TRUE, all.y=TRUE,suffixes=c('_hbcc','_mssm'))
#cor.test(sva_combo$logFC_hbcc, sva_combo$logFC_mssm, method = "spearman",alternative="greater") #rho=0.08719992, p-value < 2.2e-16

# Define individual plots
ad_plot_no_sva <- ggplot(ad_no_sva_combo, aes(x = logFC_j_msbb, y = logFC_j_rosmap)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Alzheimer's Disease - without SVA",
    x = "MSBB logFC",
    y = "ROSMAP logFC",
    #fill = "Density level" 
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") + 
  geom_smooth(method = 'lm') + 
  xlim(-0.06, 0.04) + ylim(-0.6, 1) 

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/ad_plot_no_sva.pdf", plot = ad_plot_no_sva , width = 5, height = 5) 

ad_plot_sva_be <- ggplot(ad_sva_be_combo, aes(x = logFC_msbb, y = logFC_rosmap)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Alzheimer's Disease - with SVA",
    x = "MSBB logFC",
    y = "ROSMAP logFC"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") +
  geom_smooth(method = 'lm') +
  xlim(-0.06, 0.04) + ylim(-0.6, 1) 

bp_plot_no_sva<- ggplot(bp_no_sva_combo, aes(x = logFC_mssm, y = logFC_hbcc)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Bipolar Disorder - without SVA",
    x = "MPP logFC",
    y = "HBCC logFC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") + 
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Identity line
  xlim(-1, 1) + ylim(-1, 1.5)


bp_plot_sva_be<- ggplot(bp_sva_be_combo, aes(x = logFC_mssm, y = logFC_hbcc)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Bipolar Disorder - with SVA",
    x = "MPP logFC",
    y = "HBCC logFC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") + 
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Identity line
  xlim(-1, 1) + ylim(-1, 1.5)

sz_plot_no_sva<- ggplot(sz_no_sva_combo, aes(x = logFC_mssm, y = logFC_hbcc)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Schizophrenia - without SVA",
    x = "MPP logFC",
    y = "HBCC logFC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") + 
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Identity line
  xlim(-1, 1) + ylim(-1, 1.5)


sz_plot_sva_be<- ggplot(sz_sva_be_combo, aes(x = logFC_mssm, y = logFC_hbcc)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
  scale_fill_viridis_c(guide = "none") + 
  theme_bw() +
  labs(
    title = "Schizophrenia - with SVA",
    x = "MPP logFC",
    y = "HBCC logFC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none") + 
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Identity line
  xlim(-1, 1) + ylim(-1, 1.5)

# Combine plots
combo_disease <- 
  (ad_plot_no_sva + ad_plot_sva_be + 
   bp_plot_no_sva + bp_plot_sva_be + 
   sz_plot_no_sva + sz_plot_sva_be) +
  plot_layout(ncol = 2, nrow = 3) + 
  plot_annotation(tag_levels = 'A')

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_disease_sprearman_cor_plot2.pdf", plot = combo_disease , width = 10, height = 12) 

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/combo_disease_cor_plot2.png", plot = combo_disease , width = 10, height = 12) 


######## AUTISM SPECTRUM DISORDER ########
##########################################
yale<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_without_sva_09232024.txt",data.table=FALSE) #18255 x 7
yale_sva<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_with_sva_09232024.txt",data.table=FALSE) #18255 x 7

ucla<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_without_sva_09232024.txt",data.table=FALSE) #20850 x 7
ucla_sva<-fread("/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_with_sva_09232024.txt",data.table=FALSE) #20850 x 7

no_sva_combo<-merge(yale,ucla,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_yale','_ucla'))
x<-cor.test(no_sva_combo$logFC_yale, no_sva_combo$logFC_ucla, method = "spearman", exact=FALSE,alternative = "greater");x  #na.omit = default; r= -0.09316965;  p-value < 2.2e-16

sva_combo<-merge(yale_sva,ucla_sva,by='X', all.x=TRUE, all.y=TRUE, suffixes=c('_yale_sva','_ucla_sva'))
y<-cor.test(sva_combo$logFC_yale_sva, sva_combo$logFC_ucla_sva, method = "spearman", exact=FALSE, alternative = "greater")#na.omit = default; r= -0.04837939;  p-value = 1.012e-10

###### PARKINSON'S DISEASE ######
#################################
GSE8397_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE8397_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE7621_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE7621_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE24378_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE24378_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20292_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20292_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20141_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20141_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20163_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20163_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20164_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20164_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE49036_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE49036_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")


combo_no_sva<-merge(GSE8397_none,GSE7621_none, by='X', all.x=TRUE, all.y=TRUE) #22283x13
combo_sva<-merge(GSE8397_sva,GSE7621_sva, by='X', all.x=TRUE, all.y=TRUE) #22283x13

cor.test(combo_no_sva$logFC.x, combo_no_sva$logFC.y,method="spearman") #r=0.2091094
cor.test(combo_sva$logFC.x, combo_sva$logFC.y,method="spearman") #r=0.6351632 UNEXPECTED RESULTS!!!!

# Create lists for _none and _sva datasets
datasets <- list(
  GSE8397_none = GSE8397_none,
  GSE7621_none = GSE7621_none,
  GSE24378_none = GSE24378_none,
  GSE20292_none = GSE20292_none,
  GSE20141_none = GSE20141_none,
  GSE20163_none = GSE20163_none,
  GSE20164_none = GSE20164_none,
  GSE49036_none = GSE49036_none
)

datasets <- list(
  GSE8397_sva = GSE8397_sva,
  GSE7621_sva = GSE7621_sva,
  GSE24378_sva = GSE24378_sva,
  GSE20292_sva = GSE20292_sva,
  GSE20141_sva = GSE20141_sva,
  GSE20163_sva = GSE20163_sva,
  GSE20164_sva = GSE20164_sva,
  GSE49036_sva = GSE49036_sva
)

###single cor: 
combo_no_sva<-merge(GSE8397_none,GSE7621_none, by='X', all.x=TRUE, all.y=TRUE) #22283x13
combo_sva<-merge(GSE8397_sva,GSE7621_sva, by='X', all.x=TRUE, all.y=TRUE) #22283x13

cor.test(combo_no_sva$logFC.x, combo_no_sva$logFC.y,method="spearman") #r=0.009507782 
cor.test(combo_sva$logFC.x, combo_sva$logFC.y,method="spearman") #r=-0.009109724 UNEXPECTED RESULTS!!!!

calculate_spearman_correlations <- function(datasets) {
  dataset_names <- names(datasets)
  n <- length(datasets)
  cor_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(dataset_names, dataset_names))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      data_i <- datasets[[i]]
      data_j <- datasets[[j]]
      merged_data <- merge(data_i, data_j, by = 'X', all.x = TRUE, all.y = TRUE)
      cor_value <- cor(merged_data$logFC.x, merged_data$logFC.y, method = "spearman", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
    }
  }
  return(cor_matrix)
}

spearman_cor_matrix <- calculate_spearman_correlations(datasets)
diag(spearman_cor_matrix) <- 1
spearman_cor_matrix

# write.table(spearman_cor_matrix , file = "/sc/arion/projects/mscic1/results/jolie/GEO/no_sva_spearman_cor_matrix_final_final.txt", sep = "\t", row.names = TRUE, col.names = NA)
# write.table(spearman_cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/GEO/sva_spearman_cor_matrix_final_final.txt", sep = "\t", row.names = TRUE, col.names = NA)

##START HERE FOR PLOTTING PD
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
heatmap_plot <- ggplot(no_sva_cor_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid= "white",
                       limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16)) + 
  coord_fixed() +
  labs(title = "Dataset Replicability without SVA") + 
  xlab("Datasets") + ylab("Datasets"); heatmap_plot 


ggsave("/hpc/users/hoangd02/www/plots/GEO/no_sva_spearman_cor_matrix_heatmap_final.png", plot = heatmap_plot, width = 10, height = 10)

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

# Generate the heatmap
heatmap_plot <- ggplot(difference_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid= "white",  # purple to yellow
                       limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        axis.title.x = element_text(size = 16),  # Adjust x-axis title size
        axis.title.y = element_text(size = 16)) +  # Adjust y-axis title size
  coord_fixed() +
  labs(title = "Difference between LogFC Correlations Without SVA and With SVA", x = "Datasets", y = "Datasets"); heatmap_plot

ggsave("/hpc/users/hoangd02/www/plots/GEO/difference_no_sva_spearman_cor_matrix_heatmap_final.pdf", plot = heatmap_plot, width = 10, height = 10)

#ggsave("/hpc/users/hoangd02/www/plots/GEO/difference_no_sva_spearman_cor_matrix_heatmap_final.png", plot = heatmap_plot, width = 10, height = 10)


