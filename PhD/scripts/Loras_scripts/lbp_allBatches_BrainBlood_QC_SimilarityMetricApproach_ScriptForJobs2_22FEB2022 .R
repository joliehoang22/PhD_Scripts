test <- fread("/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/Covariate1.tsv")



#+NAME: WIP_compare_correlation_matrices_from_liv_and_pm_residuals_PERM
#+BEGIN_SRC shell

# setup 
  module load R/4.0.3
  RSCRIPT=/sc/arion/work/liharl02/lbp_allBatches_BloodBrain_QC_similarityMetrics_22FEB2022.r
  SCRATCH=/sc/arion/projects/psychgen/lbp/scratch/lbp_allBatches_BloodBrain_QC/similarityMetric_QCapproach/

# test run
  cd ${SCRATCH} 
  i=1
  Rscript ${RSCRIPT} ${i}

  #run the r script command which is Rscript ${RSCRIPT} ${i}
  #if that runs and the jobs don't it's the jobs fault

# the failed ones -- define the numbers as an array
  cd ${SCRATCH}
  ARRAY=(6 126 144 145)
  echo ${ARRAY[*]}
  
  for i in "${ARRAY[@]}"
  do mybsub psychgen cov${i} 5000 1:00 premium 16 "Rscript ${RSCRIPT} ${i}"
  done

  #run the r script command which is Rscript ${RSCRIPT} ${i}
  #if that runs and the jobs don't it's the jobs fault

# send jobs
  cd ${SCRATCH}
  for i in {1..145} #or for i in i through nrow covs 
  do mybsub psychgen cov${i} 5000 1:00 premium 16 "Rscript ${RSCRIPT} ${i}" #maybe change from premium to private if not working fast enough; cov is being defined here as what the file will be called; the 5000 is memory=50GB; make 50000 if it fails; once it starts get killed after 1 hour if it doesn't finish; 16 is the cores/nodes/thing for parallelization; Rscript is the unix command for an r script; i after the name of the R script is the argument that's being taken by the R script 
  done

# check 
  cd ${SCRATCH}
  ls cov*stdout | sort | uniq  > sent
  grep Success cov*stdout | awk -F":" '{print $1}' | sort | uniq > success
  comm -23 sent success | sort | uniq > fail
  wc -l sent success fail
  ##10000 sent
  ##10000 success
  ##   0 fail
  grep Fail cov*stdout | awk -F":" '{print $1}' | sort | uniq > fail
 




# combine
  cd ${SCRATCH}
  head -1 perm1.tsv > combined.tsv
  cat perm*.tsv | grep iteration -v >> combined.tsv


  # clean up
  rm ${SCRATCH}/perm*

  #kill a job by job id
  bkill -J cov145

#############################################################################################################
# Approach
# Make a script that takes as input expression and covariate data and the current formula that has been selected
# Script would then, as a job for each covariate, calculate the relationship between the covariate and PC2 on the current set of residuals or vobject [use Noam's function and calculate all covs to all pcs and subset for the one], and the metric
# Then [for every covariate], it would calculate the residuals including that covariate added to the formula, redo the metric, and return the covariate's name, it's relationship to the PCs, and the metric
# Would have 154 of those for each round; 10 min for all of them to run
# Read the 154 outputs into R and put them together

# THE SCRIPT WOULD BE FOR JUST ONE COVARIATE, AND THEN IN THE JOB WE WOULD SPECIFY FOR I IN ALLCOVS

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
canCorAllAgainstAll_Original <- function(X, Y = X,minimum_intersect=0) {
    # Compute canonical correlation of all columns of X against all columns of Y,
    # similar to variancePartition::canCorPairs.
    library(stringr)
    X_formulas <- lapply(str_c("~", colnames(X)), as.formula)
    X_varList <- lapply(X_formulas, function(xf) model.matrix.lm(xf, X, na.action = "na.pass")[,-1, drop = FALSE])
    Y_formulas <- lapply(str_c("~", colnames(Y)), as.formula)
    Y_varList <- lapply(Y_formulas, function(yf) model.matrix.lm(yf, Y, na.action = "na.pass")[,-1, drop = FALSE])
    XY_cc <- matrix(nrow = ncol(X), ncol = ncol(Y), data = 0,
                    dimnames = list(colnames(X), colnames(Y)))
    for (ix in seq_along(X_varList)) {
        keep1 = apply(X_varList[[ix]], 1, function(x) !any(is.na(x)))
        for (iy in seq_along(Y_varList)) {
            keep2 = apply(Y_varList[[iy]], 1, function(x) !any(is.na(x)))
            keep = keep1 & keep2
            if(sum(keep)>minimum_intersect){
            fit <- cancor(X_varList[[ix]][keep, , drop = FALSE], Y_varList[[iy]][keep, , drop = FALSE])
            # Using root-mean-square to summarize, as discussed with Gabriel Hoffman
            XY_cc[ix,iy] <- sqrt(mean(fit$cor^2))
        }else{
            XY_cc[ix,iy] <- NA
        }
        }
    }
    return(XY_cc)
    }
#############################################################################################################
args <- commandArgs(trailingOnly=TRUE)
covNumber <- as.integer(args[[1]]) #145 covs total
outFile <- paste0("Covariate", covNumber, ".tsv")

#start at the seventh column of the covs table; 7:153, 170

# covs <- colnames(cv)
# covs <- covs[7:170]
# covs <- covs[c(164, 1:163)]
# covs <- covs[1:148]
# covs <- covs[c(1:124, 127:148)]
# covs <- covs[c(1:131, 133:146)]
# covs2<-data.table(covs)

# fwrite(covs2, "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_CovsForSelfNonSelfMetricComparisons_22FEB2022.csv")

# data
  fc <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_530LIVINGsamples_15FEB2022.RDS")
  cv <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")
  nw <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_18FEB2022.RDS")
  covs <- fread("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_CovsForSelfNonSelfMetricComparisons_22FEB2022.csv")

  covOfInterest <- covs[covNumber]$covs

  isexpr <- rowSums(cpm(fc)>=1) >= 0.1*ncol(fc)
  geneExpr <- DGEList( fc[isexpr,] )
  geneExpr <- calcNormFactors( geneExpr )
  
  geneExpr <- geneExpr[, match(cv$SAMPLE_ISMMS, colnames(geneExpr))]
  vobj <- voom(geneExpr)
  dim(geneExpr) #[1] 23202   530
  
  # identical(cv$SAMPLE_ISMMS, colnames(geneExpr))
  rownames(cv) <- cv$SAMPLE_ISMMS

  form <- as.formula(~ get(covOfInterest)) 
  design <- model.matrix(form, cv)
  expr.data <- vobj$E
  res.lmgroup <- lmFit(expr.data, design)
  res = residuals(res.lmgroup, vobj) #residuals

#SampleByVariable=t(cov(vobj$E))  
SampleByVariable=t(cov(res))  
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)

#The covariates-PCs correlation table
resCor=canCorAllAgainstAll_Original(cv,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]
all <- as.data.frame(resCor, keep.rownames=TRUE)
#all$cov <- rownames(all)

covCors <- all[rownames(all)==eval(covOfInterest), ] #need to pull out only the corr for the selected variable


# number of brain blood pairs (self)
  dt <- cv[,.(iid=IID_ISMMS, sid=SAMPLE_ISMMS, timepoint=mymet_timepoint, bbstatus)]
  ct <- merge( dt[bbstatus=="Brain",.(iid, brain=sid, timepoint)],
              dt[bbstatus=="Blood",.(iid, blood=sid, timepoint)], 
              by=c("iid","timepoint") )
  nrow(ct) #[1] 233 ... why only 233? I remember there were 244 pairs sequenced

# all brain blood pairs (self and nonself) -- 108345 pairs
  # ct[,class:="self"]
  # nw <- copy(ct)
  # for (i in 1:nrow(ct)) {
  #     print(i)
  #     curBrain <- ct[i]$brain
  #     curBlood <- ct[i]$blood
  #     addBrain <- ct[brain!=curBrain]
  #     addBlood <- copy(addBrain)
  #     addBrain[,brain:=curBrain] #makes every single brain be the current brain as all of those pairs are nonself
  #     addBlood[,blood:=curBlood]
  #     addBrain[,class:="nonself"]
  #     addBlood[,class:="nonself"]
  #     nw <- rbind(nw, addBrain, addBlood)  
  # }
  
  nw[,timepoint:=NULL]
  nw[,iid:=NULL]
  
# self vs nonself rho
  nw2 <- foreach(i = 1:nrow(nw), .combine = rbind )%dopar%{ #parallelizes;  instead of running one iteration of the loop at a time it will run 20 at a time or however many cores you set
      if (i %% 100 == 0 ) cat("\r",i," of ", nrow(nw),"\t\t")
      brn <- nw[i]$brain
      bld <- nw[i]$blood
      cls <- nw[i]$class
      vc1 <- res[, brn]
      vc2 <- res[, bld]
      rho <- cor(vc1, vc2, method="spearman")
      add <- data.table("brain"=brn, "blood"=bld, "class"=cls, "rho"=rho)
  }
  nw2 <- merge(merge(nw2, ct[,.(iid_brain=iid, brain)]), ct[,.(iid_blood=iid, blood)], by="blood" ) 
  nw2[class=="self",class1:="samePersonSameSurgery"]
  nw2[class=="nonself",class1:="not_samePersonSameSurgery"]
  nw2[iid_brain==iid_blood, class2:="samePerson"]
  nw2[iid_brain!=iid_blood, class2:="not_samePerson"]
  nw2[,class:=NULL]

  mean1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho)
  mean2 <-mean(nw2[class1=="not_samePersonSameSurgery"]$rho)
  mean3 <- mean(nw2[class2=="samePerson"]$rho)
  mean4 <- mean(nw2[class2=="not_samePerson"]$rho)

# value we want to maximize
  metric1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho) / mean(nw2[class1=="not_samePersonSameSurgery"]$rho) #[1] 1.005123
  metric2 <- mean(nw2[class2=="samePerson"]$rho) / mean(nw2[class2=="not_samePerson"]$rho) #[1] 1.000946


  # saveRDS(nw2, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_Vobj_Rho_22FEB2022.RDS")
  
  #return the covariate's name, it's relationship to the PCs, and the metric


  #create the master table that you rbind the rest of them to -- 

  # master <- data.table(Covariate=covs, "PC1cor"=123, "PC2cor"=123, "PC3cor"=123, "PC4cor"=123, "PC5cor"=123, "mean_samePersonSurg_notSameSame"=12345, "mean_samePerson_notSamePerson"=12345)

  # master[cov==eval(covOfInterest)]$PC1cor <- covCors$PC1
  # master[cov==eval(covOfInterest)]$PC2cor <- covCors$PC2
  # master[cov==eval(covOfInterest)]$PC3cor <- covCors$PC3
  # master[cov==eval(covOfInterest)]$PC4cor <- covCors$PC4
  # master[cov==eval(covOfInterest)]$PC5cor <- covCors$PC5
  # master[cov==eval(covOfInterest)]$mean_samePersonSurg_notSameSame <- mean1
  # master[cov==eval(covOfInterest)]$mean_samePerson_notSamePerson <- mean2


master <- data.table(Covariate=eval(covOfInterest), "PC1cor"=covCors$PC1, "PC2cor"=covCors$PC2, "PC3cor"=covCors$PC3, "PC4cor"=covCors$PC4, "PC5cor"=covCors$PC5, "meanRho_SamePersonSameSurg"= mean1, "meanRho_NotSamePersonSameSurg"=mean2, "meanRho_SamePerson"=mean3, "meanRho_NotSamePerson"=mean4, "mean_samePersonSurg_notSameSame"=metric1, "mean_samePerson_notSamePerson"=metric2)


fwrite(master, sep='\t', quo=F, row=F, file=outFile)


saveRDS(all, "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Vobject_CovPCCorrs_25FEB2022.RDS")
