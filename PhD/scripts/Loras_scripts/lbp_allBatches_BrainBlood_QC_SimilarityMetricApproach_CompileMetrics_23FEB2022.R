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
success <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/success", header=F)
covs <- fread("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_CovsForSelfNonSelfMetricComparisons_22FEB2022.csv")
cv <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")

test <- data.table(paste0("cov", seq(1, 145, by=1), ".stdout"))

test[!V1 %in% success$V1]$V1
#[1] "cov6.stdout"   "cov126.stdout" "cov144.stdout" "cov145.stdout"

#The error for all of these 
# Warning message:
# Partial NA coefficients for 23202 probe(s) 
# Error in prcomp.default(SampleByVariable, scale = T) : 
#   cannot rescale a constant/zero column to unit variance
# Calls: prcomp -> prcomp.default
# Execution halted

covs[6]
#                                covs
# 1: FASTQC_Per_tile_sequence_quality

covs[126]
#                   covs
# 1: mymet_extractiondate

covs[144]
#          covs
# 1: rapidBatch

covs[145]
#          covs
# 1: s4newbatch


#Fix these factors by dropping the empty factor levels in order to get rid of the error during the jobs

table(cv$FASTQC_Per_tile_sequence_quality)
# FAIL PASS WARN 
#    0  529    1
cv$FASTQC_Per_tile_sequence_quality <- droplevels(cv$FASTQC_Per_tile_sequence_quality)

table(cv$mymet_extractiondate)
cv$mymet_extractiondate <- droplevels(cv$mymet_extractiondate)

table(cv$rapidBatch)
cv$rapidBatch <- droplevels(cv$rapidBatch)

table(cv$s4newbatch)
cv$s4newbatch <- droplevels(cv$s4newbatch)


saveRDS(cv, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_forQC_530LIVINGSamples_23FEB2022.RDS")


#############################################################################################################
#Merge everything

#Test
cov1 <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/Covariate1.tsv")

cov2 <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/Covariate2.tsv")

dt <- merge(cov1, cov2, all.x=T, all.y=T)


setwd("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/")
files <- paste0("Covariate", seq(1, 145, by=1), ".tsv")

  count <- 0 
  for (i in files) {
      count <- count + 1
      print(count)
      path <- paste0("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/", i)
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
nw2 <- readRDS(file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_Vobj_Rho_22FEB2022.RDS")

  mean1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho)
  mean2 <-mean(nw2[class1=="not_samePersonSameSurgery"]$rho)
  mean3 <- mean(nw2[class2=="samePerson"]$rho)
  mean4 <- mean(nw2[class2=="not_samePerson"]$rho)

# value we want to maximize
  metric1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho) / mean(nw2[class1=="not_samePersonSameSurgery"]$rho) #[1] 1.005123
  metric2 <- mean(nw2[class2=="samePerson"]$rho) / mean(nw2[class2=="not_samePerson"]$rho) #[1] 1.000946


master <- data.table(Covariate="Vobject", "PC1cor"=NA, "PC2cor"=NA, "PC3cor"=NA, "PC4cor"=NA, "PC5cor"=NA, "meanRho_SamePersonSameSurg"= mean1, "meanRho_NotSamePersonSameSurg"=mean2, "meanRho_SamePerson"=mean3, "meanRho_NotSamePerson"=mean4, "mean_samePersonSurg_notSameSame"=metric1, "mean_samePerson_notSamePerson"=metric2)


master2 <- merge(master, dt, by=intersect(names(master), names(dt)), all.x=T, all.y=T)

master3 <- master2[c(128, 1:127, 129:146), ]

saveRDS(master3, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_AllCovariates_Rho_SimilarityMetrics_BasedOnResiduals_23FEB2022.RDS")


#############################################################################################################
apply(master3[,2:12], 2, range)
apply(master3[2:146,2:12], 2, range)
apply(master3[,7:12], 2, range)

#      meanRho_SamePersonSameSurg meanRho_NotSamePersonSameSurg
# [1,]                 -0.7487446                    -0.7589946
# [2,]                  0.6076121                     0.5947678
#      meanRho_SamePerson meanRho_NotSamePerson mean_samePersonSurg_notSameSame
# [1,]         -0.7547386            -0.7589942                      -0.3788882
# [2,]          0.5977120             0.5947804                    2661.7974260
#      mean_samePerson_notSamePerson
# [1,]                   -1715.06129
# [2,]                      25.19967



colclasses <- lapply(cv, class)

classes <- do.call(rbind, colclasses)
classes2 <- as.data.table(classes, keep.rownames="Covariate")
colnames(classes2)[2] <- "Class"

classes3 <- classes2[Covariate %in% covs$covs, ]

factors <- classes3[Class=="factor", ]
numeric <- classes3[Class=="numeric", ]

factors2 <- as.data.table(factors$Covariate)
numeric2 <- as.data.table(numeric$Covariate)
colnames(factors2)[1] <- "covs"
colnames(numeric2)[1] <- "covs"


fwrite(factors2, quo=F, row=F, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_FactorCovsForSelfNonSelfMetricComparisons_25FEB2022.csv")
fwrite(numeric2, quo=F, row=F, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_NumericCovsForSelfNonSelfMetricComparisons_25FEB2022.csv")



