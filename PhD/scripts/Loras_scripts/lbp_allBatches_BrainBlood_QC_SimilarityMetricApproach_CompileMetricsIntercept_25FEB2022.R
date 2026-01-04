module load R/4.0.3 
R

rm(list=ls())
options(stringsAsFactors=F)
##################################################
# setup
  library(data.table)
  library(edgeR)
  library(limma)
  library(foreach)
  library(parallel)
  library(doMC)
  library(variancePartition)

  options(cores = detectCores())
  registerDoMC(16)
#############################################################################################################
# success <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/keepIntercepts/success", header=F)
covs <- fread("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_CovsForSelfNonSelfMetricComparisons_22FEB2022.csv")
cv <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")

#############################################################################################################
#Merge everything

#Test
cov1 <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/Covariate1.tsv")

cov2 <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/Covariate2.tsv")

dt <- merge(cov1, cov2, all.x=T, all.y=T)


setwd("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/keepIntercepts/")
files <- paste0("Covariate", seq(1, 145, by=1), ".tsv")

  count <- 0 
  for (i in files) {
      count <- count + 1
      print(count)
      path <- paste0("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/keepIntercepts/", i)
      if (!file.exists(path)) {
          cat ("No output for covariate", i,"\n")
      } else { 
          fc <- fread(path)
          if (count == 1) {
           dt <- fc
      } else {
        dt <- merge(dt, fc, all.x=T, all.y=T)
        }
      }
    } 


#ADD THE ONE WITH NO COVS TO THIS AND SAVE
##################################################
nw2 <- readRDS(file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_Vobj_Rho_22FEB2022.RDS")

  mean1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho)
  mean2 <-mean(nw2[class1=="not_samePersonSameSurgery"]$rho)
  mean3 <- mean(nw2[class2=="samePerson"]$rho)
  mean4 <- mean(nw2[class2=="not_samePerson"]$rho)

# value we want to maximize
  metric1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho) / mean(nw2[class1=="not_samePersonSameSurgery"]$rho) #[1] 1.005123
  metric2 <- mean(nw2[class2=="samePerson"]$rho) / mean(nw2[class2=="not_samePerson"]$rho) #[1] 1.000946


master <- data.table(Covariate="Vobject", "meanRho_SamePersonSameSurg"= mean1, "meanRho_NotSamePersonSameSurg"=mean2, "meanRho_SamePerson"=mean3, "meanRho_NotSamePerson"=mean4, "mean_samePersonSurg_notSameSame"=metric1, "mean_samePerson_notSamePerson"=metric2)

#"PC1cor"=NA, "PC2cor"=NA, "PC3cor"=NA, "PC4cor"=NA, "PC5cor"=NA,

master2 <- merge(master, dt, by=intersect(names(master), names(dt)), all.x=T, all.y=T)

master3 <- master2[c(128, 1:127, 129:146), ]


#ADD THE VOBJECT PC CORS WITH THE COVS AND MERGE
##################################################
vobjpc <- readRDS(file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Vobject_CovPCCorrs_25FEB2022.RDS")


vobjpc$Covariate <- rownames(vobjpc)
master4 <- merge(vobjpc, master3, by="Covariate", all.y=T)

dim(master4)
#[1] 146  12

colnames(master4)[2] <- "vobj-PC1_cor"
colnames(master4)[3] <- "vobj-PC2_cor"
colnames(master4)[4] <- "vobj-PC3_cor"
colnames(master4)[5] <- "vobj-PC4_cor"
colnames(master4)[6] <- "vobj-PC5_cor"

master5 <- as.data.table(master4)

#MAKE THE SUBTRACTION METRICS FROM VOBJ
##################################################
master
#    Covariate meanRho_SamePersonSameSurg meanRho_NotSamePersonSameSurg
# 1:   Vobject                  0.5611774                     0.5583173
#    meanRho_SamePerson meanRho_NotSamePerson mean_samePersonSurg_notSameSame
# 1:          0.5588487             0.5583208                        1.005123
#    mean_samePerson_notSamePerson
# 1:                      1.000946

colnames(master) #master here is just the vobject values
# [1] "Covariate"                       "meanRho_SamePersonSameSurg"     
# [3] "meanRho_NotSamePersonSameSurg"   "meanRho_SamePerson"             
# [5] "meanRho_NotSamePerson"           "mean_samePersonSurg_notSameSame"
# [7] "mean_samePerson_notSamePerson"  

master5[, vobj_diff_samePersonSurg_notSameSame:=master$mean_samePersonSurg_notSameSame-mean_samePersonSurg_notSameSame]
master5[, vobj_diff_mean_samePerson_notSamePerson:=master$mean_samePersonSurg_notSameSame-mean_samePerson_notSamePerson]


saveRDS(master5, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_AllCovariates_Rho_SimilarityMetrics_BasedOnResidualsWithIntercept_25FEB2022.RDS")


#############################################################################################################
apply(master5[,2:12], 2, range)
apply(master5[2:146,2:12], 2, range)
apply(master5[,7:12], 2, range)

apply(master5, 2, range)

