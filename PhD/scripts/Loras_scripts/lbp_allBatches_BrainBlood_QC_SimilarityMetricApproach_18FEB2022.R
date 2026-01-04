# Approach
# Make a script that takes as input expression and covariate data and the current formula that has been selected
# Script would then, as a job for each covariate, calculate the relationship between the covariate and PC2 on the current set of residuals or vobject [use Noam's function and calculate all covs to all pcs and subset for the one], and the metric
# Then [for every covariate], it would calculate the residuals including that covariate added to the formula, redo the metric, and return the covariate's name, it's relationship to the PCs, and the metric
# Would have 154 of those for each round; 10 min for all of them to run
# Read the 154 outputs into R and put them together

# THE SCRIPT WOULD BE FOR JUST ONE COVARIATE, AND THEN IN THE JOB WE WOULD SPECIFY FOR I IN ALLCOVS

module load R/4.0.3 #very important that you load this R version, dream will not work with the regular version that auto loads on minerva
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
  setwd("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData")

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
covOfInterest <- #HOWEVER THIS GETS PULLED FROM UNIX 

# data
  fc <- readRDS("lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_530LIVINGsamples_15FEB2022.RDS")
  cv <- readRDS("lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")
  nw <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_18FEB2022.RDS")

  isexpr <- rowSums(cpm(fc)>=1) >= 0.1*ncol(fc)
  geneExpr <- DGEList( fc[isexpr,] )
  geneExpr <- calcNormFactors( geneExpr )
  identical(cv$SAMPLE_ISMMS, colnames(geneExpr))
  rownames(cv) <- cv$SAMPLE_ISMMS
  vobj <- voom(geneExpr)
  dim(geneExpr) #[1] 23202   530


cellcounts <- tail(colnames(cv), n=17) #CONFIRM WHAT IS CORRECT HERE
cellcounts <- cellcounts[1:16] #remove the last one, bbstatus

# remove <- c("SAMPLE_ISMMS", cellcounts)

# form=as.formula(paste("~",paste(colnames(cv)[!colnames(cv) %in% remove],collapse="+")))
# C = canCorPairs(form,cv[,!colnames(cv) %in% remove, with=F])

form <- ~ (1|IID_ISMMS) #this is the formula you will use for the residuals
vobj = voomWithDreamWeights( geneExpr, form, cv, BPPARAM = MulticoreParam(5))
identical(cv$SAMPLE_ISMMS, colnames(vobj$E))

library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 20) #this is to make the functions run -- you want this number and the number below in MultiCoreParam() to multiply to 100 
#Sys.setenv(OMP_NUM_THREADS = 10)

resfit = dream(vobj, form, cv, BPPARAM = MulticoreParam(5), computeResiduals = TRUE) #calculating the residuals
resfit = lmfit(vobj, computeResiduals = TRUE) #calculating the residuals
res <- residuals(resfit)


DATA <- cv[, !colnames(cv) %in% cellcounts, with=FALSE] #not sure if this is needed here
plotpath = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/" 
vobj <- vobjAllGenes

SampleByVariable=t(cov(vobj$E)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)

all <- as.data.frame(resCor, keep.rownames=TRUE)
#need to pull out only the corr for the selected variable


# C2<-as.data.frame(C[rownames(C), "bbstatus"])
# colnames(C2)[1] <- "bbstatus_corr"
# C2$covs <- rownames(C2)
# all <- merge(C2, all, by.x="covs", by.y=0)
# all <- all[order(-all$PC1),] #ordered by highest correlation with PC1
# head(all,40) 


# number of brain blood pairs (self)
  dt <- cv[,.(iid=IID_ISMMS, sid=SAMPLE_ISMMS, timepoint=mymet_timepoint, bbstatus)]
  ct <- merge( dt[bbstatus=="Brain",.(iid, brain=sid, timepoint)],
              dt[bbstatus=="Blood",.(iid, blood=sid, timepoint)], 
              by=c("iid","timepoint") )
  nrow(ct) #[1] 233 ... why only 233? I remember there were 244 pairs sequenced

# all brain blood pairs (self and nonself) -- 108345 pairs
  ct[,class:="self"]
  nw <- copy(ct)
  for (i in 1:nrow(ct)) {
      print(i)
      curBrain <- ct[i]$brain
      curBlood <- ct[i]$blood
      addBrain <- ct[brain!=curBrain]
      addBlood <- copy(addBrain)
      addBrain[,brain:=curBrain] #makes every single brain be the current brain as all of those pairs are nonself
      addBlood[,blood:=curBlood]
      addBrain[,class:="nonself"]
      addBlood[,class:="nonself"]
      nw <- rbind(nw, addBrain, addBlood)  
  }
  
  nw[,timepoint:=NULL]
  nw[,iid:=NULL]

  #saveRDS(nw, file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_18FEB2022.RDS")
  
# self vs nonself rho
  nw2 <- foreach(i = 1:nrow(nw), .combine = rbind )%dopar%{ #parallelizes;  instead of running one iteration of the loop at a time it will run 20 at a time or however many cores you set
      if (i %% 100 == 0 ) cat("\r",i," of ", nrow(nw),"\t\t")
      brn <- nw[i]$brain
      bld <- nw[i]$blood
      cls <- nw[i]$class
      vc1 <- vobj$E[, brn]
      vc2 <- vobj$E[, bld]
      rho <- cor(vc1, vc2, method="spearman")
      add <- data.table("brain"=brn, "blood"=bld, "class"=cls, "rho"=rho)
  }
  nw2 <- merge(merge(nw2, ct[,.(iid_brain=iid, brain)]), ct[,.(iid_blood=iid, blood)], by="blood" ) 
  nw2[class=="self",class1:="samePersonSameSurgery"]
  nw2[class=="nonself",class1:="not_samePersonSameSurgery"]
  nw2[iid_brain==iid_blood, class2:="samePerson"]
  nw2[iid_brain!=iid_blood, class2:="not_samePerson"]
  nw2[,class:=NULL]

# value we want to maximize
  mean(nw2[class1=="samePersonSameSurgery"]$rho) / mean(nw2[class1=="not_samePersonSameSurgery"]$rho) #[1] 1.005123
  mean(nw2[class2=="samePerson"]$rho) / mean(nw2[class2=="not_samePerson"]$rho) #[1] 1.000946


  saveRDS(nw2, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BLOODandBRAIN_Pairs_SelfNonSelf_Vobj_Rho_22FEB2022.RDS")
  #return the covariate's name, it's relationship to the PCs, and the metric

  mean1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho)
  mean2 <-mean(nw2[class1=="not_samePersonSameSurgery"]$rho)
  mean3 <- mean(nw2[class2=="samePerson"]$rho)
  mean4 <- mean(nw2[class2=="not_samePerson"]$rho)

# value we want to maximize
  metric1 <- mean(nw2[class1=="samePersonSameSurgery"]$rho) / mean(nw2[class1=="not_samePersonSameSurgery"]$rho) #[1] 1.005123
  metric2 <- mean(nw2[class2=="samePerson"]$rho) / mean(nw2[class2=="not_samePerson"]$rho) #[1] 1.000946


master <- data.table(Covariate=eval(covOfInterest), "PC1cor"=covCors$PC1, "PC2cor"=covCors$PC2, "PC3cor"=covCors$PC3, "PC4cor"=covCors$PC4, "PC5cor"=covCors$PC5, "meanRho_SamePersonSameSurg"= mean1, "meanRho_NotSamePersonSameSurg"=mean2, "meanRho_SamePerson"=mean3, "meanRho_NotSamePerson"=mean4, "mean_samePersonSurg_notSameSame"=metric1, "mean_samePerson_notSamePerson"=metric2)

