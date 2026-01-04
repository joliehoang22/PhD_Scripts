module load R/4.0.3
R

rm(list=ls())
options(stringsAsFactors=F)
library(data.table)
library(openxlsx)
###############################################################################################################
labs <- read.xlsx("/sc/arion/projects/psychgen/lbp/data/emr/msd494/MSD494LBPMRN112018.xlsx", sheet = 11, colNames = T, rowNames = F)

mblood <- readRDS(file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_SurgeryDates_ClinicalLabsDates_Mapping_182BloodSamples_08FEB2022.RDS")

livlabs4 <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_CBC_ClinicalLabResults_09FEB2022.RDS")

labcast2 <- readRDS(file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_CBC_ClinicalLabResults_LongFormat_dcast_09FEB2022.RDS")

###############################################################################################################
#take only the lab dates from mblood within 3 days of the surgery date

mblood[, dateDif2:=tstrsplit(dateDif, split=" ", fixed=T, keep=1L)]
mblood$dateDif2 <- abs(as.numeric(mblood$dateDif2))

summary(mblood$dateDif2)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 0.00   10.00   13.50   34.25   18.00  972.00

mblood2 <- mblood[dateDif2<4, ]


#subset for only the labdates dates in mblood
labdates <- unique(mblood2$labDate)

livlabs5 <- livlabs4[newDate2 %in% labdates, ]
labcast3 <- labcast2[newDate2 %in% labdates, ]


unique(labcast3$ORDER_CODE)
#[1] 203-CBC & PLT & DIFF 202-CBC & PLT  

cbc <- labcast2[ORDER_CODE=="203-CBC & PLT & DIFF", ]

 apply(is.na(cbc2), 2, sum)

#keep only percentage values

cbc2 <- cbc[UNIT_OF_MEASURE=="%", ]

cbcx <- labcast2[UNIT_OF_MEASURE=="%", ]

#keep only CBC & PLT & DIFF order codes

cbc.order <- unique(labs[labs$ORDER_CODE=="203-CBC & PLT & DIFF", ]$RESULT_CODE) 

cbc3 <- cbc2[, colnames(cbc2) %in% c( "subject_id", "newDate2", "rInd", "ORDER_CODE", "UNIT_OF_MEASURE", intersect(cbc.order, colnames(cbc2))), with=FALSE ]

length(which(is.na(cbc3$`2638-PLATELET`)))

length(which(is.na(cbc3$`2638-PLATELET`)))
#[1] 159

apply(is.na(cbc3), 2, sum)

# subject_id         newDate2        rInd       ORDER_CODE  UNIT_OF_MEASURE 20243-LYMPHOCYTE   20244-MONOCYTE 20245-NEUTROPHIL 
#   0                0                0                0                0                8                8                8 
# 20246-EOSINOPHIL   20247-BASOPHIL    2638-PLATELET 
#                8                8              159

cbc3 <- cbc3[, 1:(ncol(cbc3)-1)] #remove the platelet column as it's all NAs

cbc3[rowSums(is.na(cbc3)) > 0,]

cbc3[rowSums(is.na(cbc3)) == 5,]

cbc4 <- cbc3[! rowSums(is.na(cbc3)) > 0,] #this is what we're left with as far as lab data, at the end of everything for only the CBC DIFF PLT


apply(is.na(cbcx), 2, sum) #cbcx has both of the CBC order codes 

cbcx2 <- cbcx[ , !colSums(is.na(cbcx)) > 151, with=FALSE] 
#cbcx2 <- cbcx[ , colSums(is.na(cbcx)) < 152]

apply(is.na(cbcx2), 2, sum) #cbcx has both of the CBC order codes 

###############################################################################################################
#Cibersort outputs

#1) Lake -- brain single cell
lakeout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_lake/CIBERSORTx_Results.txt")

dim(lakeout)
#[1] 533   9

colnames(lakeout)
# [1] "Mixture"     "GLU"         "GABA"        "AST"         "ODC"        
# [6] "MG"          "P-value"     "Correlation" "RMSE"    

#2) LM22 -- blood facs
lm22out <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_LM22/output2/CIBERSORTx_Results.txt")

dim(lm22out)
#[1] 243  26

colnames(lm22out)
#  [1] "Mixture"                      "B cells naive"               
#  [3] "B cells memory"               "Plasma cells"                
#  [5] "T cells CD8"                  "T cells CD4 naive"           
#  [7] "T cells CD4 memory resting"   "T cells CD4 memory activated"
#  [9] "T cells follicular helper"    "T cells gamma delta"         
# [11] "T cells regulatory (Tregs)"   "NK cells resting"            
# [13] "NK cells activated"           "Monocytes"                   
# [15] "Macrophages M0"               "Macrophages M1"              
# [17] "Macrophages M2"               "Dendritic cells resting"     
# [19] "Dendritic cells activated"    "Mast cells resting"          
# [21] "Mast cells activated"         "Eosinophils"                 
# [23] "Neutrophils"                  "P-value"                     
# [25] "Correlation"                  "RMSE"        

#From Noam:
# mer_LM22$Lymphocyte=mer_LM22$`B cells naive` + mer_LM22$`B cells memory` + mer_LM22$`Plasma cells` + mer_LM22$`T cells CD8` + mer_LM22$`T cells CD4 naive` + mer_LM22$`T cells CD4 memory resting` + mer_LM22$`T cells CD4 memory activated` + mer_LM22$`T cells follicular helper` + mer_LM22$`T cells gamma delta` + mer_LM22$`T cells regulatory (Tregs)`


#lymphocytes
lm22out$Lymphocyte_total=lm22out$`B cells naive` + lm22out$`B cells memory` + lm22out$`Plasma cells` + lm22out$`T cells CD8` + lm22out$`T cells CD4 naive` + lm22out$`T cells CD4 memory resting` + lm22out$`T cells CD4 memory activated` + lm22out$`T cells follicular helper` + lm22out$`T cells gamma delta` + lm22out$`T cells regulatory (Tregs)` + lm22out$`NK cells resting` + lm22out$`NK cells activated`

#Also monocytes is their own category not inclusive of macrophages and dendritic cells 

lm22out$Mono_Macro_DC <- lm22out$Monocytes + lm22out$`Macrophages M0` + lm22out$`Macrophages M1` + lm22out$`Macrophages M2` + lm22out$`Dendritic cells resting` + lm22out$`Dendritic cells activated` 


#3) Wilk -- single cell
wilkout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_wilk/output5/CIBERSORTx_Results.txt")

dim(wilkout)
#[1] 243  24

#PB = Peripheral Blood B Cell

colnames(wilkout)
#  [1] "Mixture"                   "RBC"                      
#  [3] "IgM PB"                    "IgG PB"                   
#  [5] "IgA PB"                    "CD14 Monocyte"            
#  [7] "CD8m T"                    "CD4m T"                   
#  [9] "CD4n T"                    "B"                        
# [11] "Platelet"                  "IFN-stim CD4 T"           
# [13] "NK"                        "Neutrophil"               
# [15] "CD16 Monocyte"             "Proliferative Lymphocytes"
# [17] "gd T"                      "pDC"                      
# [19] "Developing Neutrophil"     "SC & Eosinophil"          
# [21] "DC"                        "P-value"                  
# [23] "Correlation"               "RMSE"  

#From Noam:
# mer_wilk$Lymphocyte=mer_wilk$`CD8m T` + mer_wilk$`CD4m T` + mer_wilk$`B` + mer_wilk$`IFN-stim CD4 T` + mer_wilk$`Proliferative Lymphocytes` + mer_wilk$`gd T` + mer_wilk$`IgM PB` + mer_wilk$`IgG PB` + mer_wilk$`IgA PB`
# mer_wilk$Monocyte=mer_wilk$`CD14 Monocyte` + mer_wilk$`CD16 Monocyte`

#CD4nT +  NK missing from lymphocytes above

wilkout$Lymphocyte_total=wilkout$`CD8m T` + wilkout$`CD4m T` + wilkout$`B` + wilkout$`IFN-stim CD4 T` + wilkout$`Proliferative Lymphocytes` + wilkout$`gd T` + wilkout$`IgM PB` + wilkout$`IgG PB` + wilkout$`IgA PB` + wilkout$`NK` + wilkout$`CD4n T`

#gM PB + IgG PB + IgA PB + CD8mT + CD4nT + CD4nT + B + NK + gd T + IFN-stim CD4 T+ Proliferative Lymphocytes = lymphocytes

wilkout$Monocyte_total=wilkout$`CD14 Monocyte` + wilkout$`CD16 Monocyte`

wilkout$Mono_DC=wilkout$`CD14 Monocyte` + wilkout$`CD16 Monocyte` + wilkout$`pDC` + wilkout$DC


#4) SCP424
scpout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_SCP424/output4/CIBERSORTx_Results.txt")

dim(scpout)
#[1] 243  14

colnames(scpout)
#  [1] "Mixture"                     "CD4+_T_cell"                
#  [3] "Cytotoxic_T_cell"            "Natural_killer_cell"        
#  [5] "CD16+_monocyte"              "CD14+_monocyte"             
#  [7] "Megakaryocyte"               "B_cell"                     
#  [9] "Dendritic_cell"              "Plasmacytoid_dendritic_cell"
# [11] "Unassigned"                  "P-value"                    
# [13] "Correlation"                 "RMSE"      

#From Noam:
# mer_SCP424$Lymphocyte=mer_SCP424$`CD4+_T_cell` + mer_SCP424$`Cytotoxic_T_cell` + mer_SCP424$`B_cell`
# mer_SCP424$Monocyte=mer_SCP424$`CD16+_monocyte` + mer_SCP424$`CD14+_monocyte` 

scpout$Lymphocyte_total=scpout$`CD4+_T_cell` + scpout$`Cytotoxic_T_cell` + scpout$`B_cell` + scpout$`Natural_killer_cell`
scpout$Monocyte_total=scpout$`CD16+_monocyte` + scpout$`CD14+_monocyte`

scpout$Mono_DC=scpout$`CD16+_monocyte` + scpout$`CD14+_monocyte` + scpout$`Dendritic_cell` + scpout$`Plasmacytoid_dendritic_cell`


#5) Cibersort2
figout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_Fig2b-WholeBlood_RNAseq/output3/CIBERSORTx_Results.txt")

colnames(figout)
# [1] "Mixture"     "T cells CD8" "Monocytes"   "T cells CD4" "NKT cells"  
# [6] "B cells"     "NK cells"    "P-value"     "Correlation" "RMSE"  


dim(figout)
#[1] 243  10

#From Noam:
#mer_NSCLC_PBMC$Lymphocyte=mer_NSCLC_PBMC$`T cells CD8` + mer_NSCLC_PBMC$`T cells CD4` + mer_NSCLC_PBMC$`B cells`

figout$lymphocyte_total=figout$`T cells CD8` + figout$`T cells CD4` + figout$`B cells` + figout$`NKT cells` + figout$`NK cells`

###############################################################################################################
#now, need to merge cbc4 with all of the lymphocyte and monocyte cell type estimate columns by sample ID

#first, merge mblood2 with cbc4 via labDate
colnames(cbc4)[2] <- "labDate" 
colnames(cbc4)[1] <- "iid" 

cbc5 <- merge(cbc4, mblood2[, c(1,2,4)])

length(intersect(cbc4$labDate, mblood2$labDate))

cbc4[as.character(labDate) %in% intersect(as.character(cbc4$labDate), as.character(mblood2$labDate)),] #14 rows intersect based on the lab date


colnames(cbcx2)[2] <- "labDate" 
colnames(cbcx2)[1] <- "iid" 

cbcx2[as.character(labDate) %in% intersect(as.character(cbcx2$labDate), as.character(mblood2$labDate)),] #15 rows intersect based on the lab date

cbcx3 <- merge(cbcx2, mblood2[, c(1,2,4)]) #why is this only 11 now, i don't understand; does it have something to do with the date conversion

intersect(paste(cbcx2$iid, cbcx2$labDate), paste(mblood2$iid, mblood2$labDate))

# [1] "PT-0021 2014-06-12" "PT-0029 2014-08-04" "PT-0036 2014-10-02" "PT-0044 2014-11-10" "PT-0048 2015-02-24" "PT-0054 2015-03-28"
#  [7] "PT-0057 2015-04-28" "PT-0110 2016-08-25" "PT-0130 2017-03-10" "PT-0144 2017-06-29" "PT-0172 2018-03-16"

#ok yes so based on both, only 11 intersect -- this checks out unfortunately 

############################################################################################################
lm22out$Mixture [1]
lm22out$Lymphocyte_total [27]
lm22out$Mono_Macro_DC [28]
lm22out$Monocytes [14]

wilkout$Mixture [1]
wilkout$Lymphocyte_total [25]
wilkout$Monocyte_total [26]
wilkout$Mono_DC [27]

scpout$Mixture [1]
scpout$Lymphocyte_total [15]
scpout$Monocyte_total [16]
scpout$Mono_DC [17]

figout$Mixture [1]
figout$lymphocyte_total [11]
figout$Monocytes [3]

cellcounts <- merge(lm22out[, c(1, 27, 14, 28)], wilkout[, c(1, 25, 26, 27)], by="Mixture", suffixes = c("_lm22", "_wilk"))

cellcounts <- merge(cellcounts, scpout[, c(1, 15, 16, 17)], by="Mixture", suffixes = c("", "_scp"))

cellcounts3 <- merge(cellcounts, figout[, c(1, 11, 3)], by="Mixture", suffixes = c("", "_fig"))

colnames(cellcounts3)[3] <-"Monocytes_lm22"
colnames(cellcounts3)[4] <-"Mono_Macro_DC_lm22"
colnames(cellcounts3)[6] <-"Monocyte_total_wilk"
colnames(cellcounts3)[7] <-"Mono_DC_wilk"
colnames(cellcounts3)[8] <-"Lymphocyte_total_scp"
colnames(cellcounts3)[11] <-"Lymphocyte_total_fig"


countmrg <- merge(cbcx3, cellcounts3, by.x="LBPSEMA4_ID", by.y="Mixture", all.x=TRUE)

saveRDS(countmrg, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_Cibersortx_CellCountEstimates_CBCLabs_Merge_11samples_11FEB2022.RDS")


saveRDS(cellcounts, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_Cibersortx_CellCountEstimates_4References_243samples_11FEB2022.RDS")

#################################################################################################################
apply(is.na(countmrg), 2, sum) 

sapply(countmrg, class)

#for(i in colnames(countmrg[, 7:27])){if (class(countmrg[[i]]) == "character") countmrg[[i]] <- as.numeric(countmrg[[i]])}

for(i in colnames(countmrg)[7:27]){if (class(countmrg[[i]]) == "character") countmrg[[i]] <- as.numeric(countmrg[[i]])} #they need to be numeri to do the cor

mrgcor<- cor(countmrg[, 7:27], method="spearman")



#get rid of the LBPSEMA4BLOOD556 and do that corr separately

count2 <- countmrg[! LBPSEMA4_ID=="LBPSEMA4BLOOD556",]
apply(is.na(count2), 2, sum)

count3 <- count2[, !colSums(is.na(count2)) == 10, with=FALSE]
apply(is.na(count3), 2, sum)
sapply(count3, class)

countcor <- cor(count3[, 7:ncol(count3)], method="spearman")


cp /sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_wgcna/lbp_allBatches_RAPiD_11BloodSamples_CibersortxDeconvolution_CBCLabCounts_Cors_heatmap_11FEB2022.pdf /hpc/users/hoangd02/www/plots/lbp/

pdf("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_wgcna/lbp_allBatches_RAPiD_11BloodSamples_CibersortxDeconvolution_CBCLabCounts_Cors_heatmap_11FEB2022.pdf", height=15,width=15)
#plotCorrMatrix(xcat_moduleTraitPvalue2 ,margins = c(25, 25))
heatmap(countcor, margins = c(25, 25), cexRow = 0.5, cexCol = 0.5, col=heat.colors(12))
legend(x="right", legend="col")
dev.off()

scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_wgcna/lbp_allBatches_RAPiD_11BloodSamples_CibersortxDeconvolution_CBCLabCounts_Cors_heatmap_11FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_QC/


#############
#scp seems to be the best one based on the heatmap

find /sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/ -name '*Covariates*'

lbpcov <- readRDS(file = "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMI_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_18JAN2022.RDS")

fc <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")

colnames(cellcounts3)[1] <- "SAMPLE_ISMMS"
covmrg <- merge(lbpcov, cellcounts3, by="SAMPLE_ISMMS", all=TRUE)

#merge with brain estimates
colnames(lakeout)[1] <- "SAMPLE_ISMMS"
covmrg2 <- merge(covmrg, lakeout[, 1:6], by="SAMPLE_ISMMS", all=TRUE)

saveRDS(covmrg2, file = "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMICellTypes_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_11FEB2022.RDS")


#cbc4 <- cbc3[! rowSums(is.na(cbc3)) > 0,] 

# de.compare <- merge(hde, pm.de, by="gene", suffixes = c(".harvardpmide", ".livpmde"))
# idx <-  cbcx2[as.character(labDate) %in% intersect(as.character(cbcx2$labDate), as.character(mblood2$labDate)),]$iid #15 rows intersect based on the lab date

# idx[!idx %in% cbcx3$iid]
# #[1] PT-0039 PT-0111 PT-0134

# cbcx2[iid %in% idx[!idx %in% cbcx3$iid],]

# mblood2[iid %in% idx[!idx %in% cbcx3$iid], ]

# apply(is.na(cbc3), 2, sum)

# cbc5[, colSums(is.na(cbc5)) > 0]

# apply(is.na(cbc5), 2, sum)


# #cbc5[rowSums(is.na(cbc5)) > 0,]

# library(arsenal)
# #comparedf(cbc4, cbc5)
# summary(comparedf(cbc4, cbc5))


# cbc6 <- merge(cbc4, mblood2[, c(1,2, 4)], by = c("iid", "labDate"))

# cbc6[which(duplicated(cbc6[, c(1:10)])),]

