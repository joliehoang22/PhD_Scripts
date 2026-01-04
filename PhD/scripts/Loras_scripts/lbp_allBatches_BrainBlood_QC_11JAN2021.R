#1. REMOVE MISLABELED SAMPLES AS BEFORE -- ALEX MISLABELING ANALYSIS
#2. CELL TYPE PROPORTIONS -- Cibersort a la Noam 
#3. checking for outliers via pca, not other methods
#4. Add PMI -- done
#5. Add Bank -- done


module load R/4.0.3
R

rm(list=ls())
options(stringsAsFactors=F)
##################################################
#load libraries
##################################################
library(limma)
library(edgeR)
library(Glimma)
library(variancePartition)#, lib.loc="~/.Rlib")
library(ggplot2)
library(gridExtra)
library(grid)
library(doParallel)
registerDoParallel(20)
library(sp)
library(CellMix)
library(biomaRt)
library(gsubfn)
library(data.table)
library(sp)
library(Matrix)
###############################################################################################################

#Run the functions first

multiplot_same_legend <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  #get_legend_info
  mylegend<-g_legend(plots[[1]])

  numPlots = length(plots)

  for(i in 1:numPlots){
    plots[[i]] <- plots[[i]] + theme(legend.position="none")
  }
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
    if(is.null(mylegend)==F){
      grid.draw(mylegend)
    }    
  }
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg)>0){
    legend <- tmp$grobs[[leg]]
  }else{
    legend <- c()
  }
  return(legend)}


###############################################################################################################

# oldcov <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7MAY2021.RDS")

# oldfc <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_onlyBRAIN_532samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7MAY2021.RDS")

mv /sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS /sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/

mv /sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_BLOODandBRAIN_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS /sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/

lbpcov <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_BLOODandBRAIN_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")

fc <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")

###############################################################################################################

table(lbpcov$mymet_race)

table(lbpcov$mymet_ethnicity)

lbpcov[mymet_race=="White", mymet_race:="white"]

lbpcov$mymet_race<-droplevels(lbpcov$mymet_race)

# American Indian           asian    asian indian           black         chinese 
#               1              24               4               7               3 
#           other       pakistani         Unknown           white 
#              22               3             207             505


table(lbpcov$mymet_sex)

#  F   M 
# 300 476


fc$Geneid <- NULL
lbpcov <- lbpcov[match(colnames(fc), lbpcov$SAMPLE_ISMMS),]
identical(lbpcov$SAMPLE_ISMMS, colnames(fc))


# isexpr <- rowSums(cpm(fc)>=1) >= 0.1*ncol(fc)
# vobj <- voom(calcNormFactors(DGEList(counts=fc[isexpr,])))

# identical(lbpcov$SAMPLE_ISMMS, colnames(vobj$E))
#[1] TRUE
#lbpcov <- lbpcov[match(colnames(vobj$E), lbpcov$SAMPLE_ISMMS),]

dim(lbpcov)
#[1] 776 155

########################################################################################
#Add PMI and Bank calculated the last time around -- PMI for pm brain samples only

oldcov <- readRDS("/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/lbp_allBatches_RAPiDMYMET_QCMetrics_PlusCellTypeProportions_PlusPMI-and-Bank_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7JUN2021.RDS")

oldcov2 <- readRDS("/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/lbp_allBatches_RAPiD_MYMET_QCMetrics_PlusCellTypeProportions_PlusPMI-and-Bank_FINAL_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")

# [159] "Bank"                                                          
# [160] "cold_pmi_CORRECTED"

######################
#Bank
#######################

#All Blood is LIVING
#All Living brain is LIVING
#The PM banks are the same as before

allsamples <- data.frame(lbpcov$SAMPLE_ISMMS, lbpcov$mymet_tissue)
bank <- data.frame(oldcov$SAMPLE_ISMMS, oldcov$Bank)

allbank <- merge(allsamples, bank, by.x="lbpcov.SAMPLE_ISMMS", by.y="oldcov.SAMPLE_ISMMS", all.x=TRUE)

dim(allbank)
#[1] 776   2

colnames(allbank)[1] <- "SAMPLE_ISMMS"
colnames(allbank)[2] <- "mymet_tissue"
colnames(allbank)[3] <- "Bank"

oldcov[!oldcov$SAMPLE_ISMMS %in% allbank$SAMPLE_ISMMS]$SAMPLE_ISMMS
#[1] character(0)

allbank <- as.data.table(allbank)

allbank[mymet_tissue=="R_Blood", Bank:="LIVING"]
allbank[mymet_tissue=="L_Blood", Bank:="LIVING"]

allbank[is.na(Bank)]$SAMPLE_ISMMS
#[1] "LBPSEMA4BRAIN107"

lbpcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN107"]$mymet_living
# [1] 1
# Levels: 0 1

allbank[SAMPLE_ISMMS=="LBPSEMA4BRAIN107", Bank:="LIVING"]

length(which(is.na(allbank$Bank)))
#[1] 0

allbank2 <- allbank[, c(1,3)]


lbpcov2 <- merge(lbpcov, allbank2, by="SAMPLE_ISMMS")

identical(lbpcov2$Bank, allbank2$Bank)
#[1] TRUE


######################
#PMI
#######################

final.everything <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")

finalcov <- final.everything$covariates 

finalcov[mymet_postmortem==1 & cold_pmi_CORRECTED==0]$SAMPLE_ISMMS #dead brain with pmi of 0?
#[1] "LBPSEMA4BRAIN535"
#According to the PMI spreadsheet, this should be NA for PMI -- need to update final.everything; this is a Miami sample so not used in the PMI analyses phew

finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN535", cold_pmi_CORRECTED:=NA]

finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN535"]$IID_ISMMS
#"HBHJ_17_001_FC"

is.na(finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN535"]$cold_pmi_CORRECTED)
#[1] TRUE


final.everything$covariates[SAMPLE_ISMMS=="LBPSEMA4BRAIN535", cold_pmi_CORRECTED]
#[1] 0

final.everything$covariates[SAMPLE_ISMMS=="LBPSEMA4BRAIN535", cold_pmi_CORRECTED:=NA]


xx <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")$covariates
library(arsenal)
summary(comparedf(xx, final.everything$covariates))

# Table: Differences detected

# var.x                var.y                 ..row.names..  values.x   values.y    row.x   row.y
# -------------------  -------------------  --------------  ---------  ---------  ------  ------
# cold_pmi_CORRECTED   cold_pmi_CORRECTED              355  0          NA            355     355


#NEED TO INFORM ALEX OF THE MIAMI SAMPLE ISSUE AND RE-SAVE FINAL.EVERYTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#save the corrected final.everything, with everything identical to the old one except for this value 
saveRDS(final.everything, file="/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")

#THIS is where everything should be saved, NOT psychgen2
saveRDS(final.everything, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")
#################################################################################################################################

find /dir/to/search -name "lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS" -print

find /sc/arion/projects/psychgen/lbp/files/files_from_psychgen2 -name "lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS" -print

find /sc/arion/projects/psychgen/lbp/files/files_from_psychgen2 -wholename "*lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS*" -print


find /sc/arion/projects/psychgen/lbp/files/ -name "Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_batches_multiomics_SEP202032_09082020_2.csv" -print


find . -name '*.pl'

find . -name "*lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS*"

find /sc/arion/projects/psychgen/lbp/files -wholename "*lbp_allBatches_230-PostmortemBrainSamples_Harvard_Columbia_Miami_PMI_time-formatted-as-decimal_7JUN2021.RDS*" -print


#################################################################################################################################


#allsamples <- data.frame(lbpcov$SAMPLE_ISMMS, lbpcov$mymet_tissue)
pmi <- data.frame(finalcov$SAMPLE_ISMMS, finalcov$cold_pmi_CORRECTED)

allpmi <- merge(allsamples, pmi, by.x="lbpcov.SAMPLE_ISMMS", by.y="finalcov.SAMPLE_ISMMS", all.x=TRUE)

dim(allpmi)
#

colnames(allpmi)[1] <- "SAMPLE_ISMMS"
colnames(allpmi)[2] <- "mymet_tissue"
colnames(allpmi)[3] <- "cold_pmi_CORRECTED"

pmstatus <- data.frame(lbpcov2$SAMPLE_ISMMS, lbpcov2$mymet_postmortem, lbpcov2$Bank)

allpmi <- merge(allpmi, pmstatus, by.x="SAMPLE_ISMMS", by.y="lbpcov2.SAMPLE_ISMMS")
colnames(allpmi)[4] <- "mymet_postmortem"
colnames(allpmi)[5] <- "Bank"


allpmi[is.na(allpmi$cold_pmi_CORRECTED),]

nrow(allpmi[is.na(allpmi$cold_pmi_CORRECTED),])
#[1] 274

allpmi <- as.data.table(allpmi)

table(allpmi$mymet_postmortem)

    # living postmortem 
    #    530        246

# allpmi[mymet_tissue=="R_Blood", cold_pmi_CORRECTED:=0]
# allpmi[mymet_tissue=="L_Blood", cold_pmi_CORRECTED:=0]

allpmi[Bank=="LIVING", cold_pmi_CORRECTED:=0]


nrow(allpmi[is.na(allpmi$cold_pmi_CORRECTED),])
#[1] 19

table(allpmi[is.na(allpmi$cold_pmi_CORRECTED)]$Bank)

# COLUMBIA  HARVARD    MIAMI   LIVING 
#       18        0        1        0 

pminas2 <- allpmi[is.na(allpmi$cold_pmi_CORRECTED), ]

#     SAMPLE_ISMMS mymet_tissue cold_pmi_CORRECTED mymet_postmortem     Bank
#  1: LBPSEMA4BRAIN007      R_Brain                 NA       postmortem COLUMBIA
#  2: LBPSEMA4BRAIN043      R_Brain                 NA       postmortem COLUMBIA
#  3: LBPSEMA4BRAIN176      R_Brain                 NA       postmortem COLUMBIA
#  4: LBPSEMA4BRAIN178      R_Brain                 NA       postmortem COLUMBIA
#  5: LBPSEMA4BRAIN245      R_Brain                 NA       postmortem COLUMBIA
#  6: LBPSEMA4BRAIN258      R_Brain                 NA       postmortem COLUMBIA
#  7: LBPSEMA4BRAIN285      L_Brain                 NA       postmortem COLUMBIA
#  8: LBPSEMA4BRAIN400      R_Brain                 NA       postmortem COLUMBIA
#  9: LBPSEMA4BRAIN406      L_Brain                 NA       postmortem COLUMBIA
# 10: LBPSEMA4BRAIN407      R_Brain                 NA       postmortem COLUMBIA
# 11: LBPSEMA4BRAIN474      L_Brain                 NA       postmortem COLUMBIA
# 12: LBPSEMA4BRAIN490      L_Brain                 NA       postmortem COLUMBIA
# 13: LBPSEMA4BRAIN535      R_Brain                 NA       postmortem    MIAMI
# 14: LBPSEMA4BRAIN570      R_Brain                 NA       postmortem COLUMBIA
# 15: LBPSEMA4BRAIN598      L_Brain                 NA       postmortem COLUMBIA
# 16: LBPSEMA4BRAIN639      L_Brain                 NA       postmortem COLUMBIA
# 17: LBPSEMA4BRAIN711      L_Brain                 NA       postmortem COLUMBIA
# 18: LBPSEMA4BRAIN717      L_Brain                 NA       postmortem COLUMBIA
# 19: LBPSEMA4BRAIN779      L_Brain                 NA       postmortem COLUMBIA

#####################################################################################################
colna <- finalcov[Bank=="COLUMBIA" & is.na(cold_pmi_CORRECTED)]$IID_ISMMS 
length(colna)
#[1] 15 --> this checks out, this is the correct number of patient IDs that have NAs in Columbia per the Postmortem_Brain_Samples_Chart_LEL_8APR2021

colnaids <- finalcov[Bank=="COLUMBIA" & is.na(cold_pmi_CORRECTED)]$SAMPLE_ISMMS 

pminas[! pminas$SAMPLE_ISMMS %in% colnaids]

#        SAMPLE_ISMMS mymet_tissue cold_pmi_CORRECTED mymet_postmortem     Bank
# 1: LBPSEMA4BRAIN007      R_Brain                 NA       postmortem COLUMBIA
# 2: LBPSEMA4BRAIN245      R_Brain                 NA       postmortem COLUMBIA
# 3: LBPSEMA4BRAIN535      R_Brain                 NA       postmortem    MIAMI
# 4: LBPSEMA4BRAIN598      L_Brain                 NA       postmortem COLUMBIA


finalcov[Bank=="COLUMBIA" & is.na(cold_pmi_CORRECTED)]$SAMPLE_ISMMS

#  [1] "LBPSEMA4BRAIN043" "LBPSEMA4BRAIN176" "LBPSEMA4BRAIN178" "LBPSEMA4BRAIN258"
#  [5] "LBPSEMA4BRAIN285" "LBPSEMA4BRAIN400" "LBPSEMA4BRAIN406" "LBPSEMA4BRAIN407"
#  [9] "LBPSEMA4BRAIN474" "LBPSEMA4BRAIN490" "LBPSEMA4BRAIN570" "LBPSEMA4BRAIN639"
# [13] "LBPSEMA4BRAIN711" "LBPSEMA4BRAIN717" "LBPSEMA4BRAIN779"

nrow(finalcov[Bank=="COLUMBIA" & is.na(cold_pmi_CORRECTED)])
#[1] 15


finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN007", ] #Empty data.table
finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN245", ] #Empty data.table
finalcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN598", ] #Empty data.table

#These samples do not exist in the 518 used for the liv-pm analysis, must have been excluded from analysis as outliers, need to get their PMI and add it

pmichart <- as.data.frame(read.csv("/sc/arion/work/liharl02/Postmortem_Brain_Samples_Chart_LEL_8APR2021.csv"))
#HBHJ_17_002_FC --> Miami sample that has NA for PMI in the pmi samples chart??

finalcov[IID_ISMMS=="HBHJ_17_002_FC", cold_pmi_CORRECTED]
#doesn't exist, PM samples chart says this is a "fixed" sample, not a sample that was sequenced -- this is fine

# master <- read.xlsx("/sc/arion/projects/psychgen2/lbp/files/lbp_batch1_reRAPiD_QC3/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD.xlsx") #need the iid column in master

#The psychgen2 version
#pmix <- readRDS("/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/lbp_allBatches_230-PostmortemBrainSamples_Harvard_Columbia_Miami_PMI_time-formatted-as-decimal_7JUN2021.RDS")


pmix <- readRDS("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_230-PostmortemBrainSamples_Harvard_Columbia_Miami_PMI_time-formatted-as-decimal_7JUN2021.RDS")


pmix[SAMPLE_ISMMS %in% c("LBPSEMA4BRAIN007", "LBPSEMA4BRAIN245", "LBPSEMA4BRAIN598", "LBPSEMA4BRAIN535"), ]

#  INDIVIDUAL_NAME          ISM_SEMA4     SAMPLE_ISMMS              LIMS Name
# 1:           T-138 Sample_ISM229905-2 LBPSEMA4BRAIN245 LBPSEMA4BRAIN245_T-138
# 2:           T-177 Sample_ISM255848-2 LBPSEMA4BRAIN598 LBPSEMA4BRAIN598_T-177
# 3:           T-201 Sample_ISM223644-2 LBPSEMA4BRAIN007 LBPSEMA4BRAIN007_T-201
#        Bank cold_pmi_original_format cold_pmi_CORRECTED
# 1: COLUMBIA                    03:41               3.68
# 2: COLUMBIA                    03:30               3.50
# 3: COLUMBIA                    05:50               5.83


# allpmi[SAMPLE_ISMMS %in% c("LBPSEMA4BRAIN007", "LBPSEMA4BRAIN245", "LBPSEMA4BRAIN598"), cold_pmi_CORRECTED:=pmix$cold_pmi_CORRECTED[pmix$SAMPLE_ISMMS %in% c("LBPSEMA4BRAIN007", "LBPSEMA4BRAIN245", "LBPSEMA4BRAIN598")]]

# #this is wrong!! it fucks up the order for some reason, just do it by hand
# allpmi[SAMPLE_ISMMS %in% c("LBPSEMA4BRAIN007", "LBPSEMA4BRAIN245", "LBPSEMA4BRAIN598"), ]
# #        SAMPLE_ISMMS mymet_tissue cold_pmi_CORRECTED mymet_postmortem     Bank
# # 1: LBPSEMA4BRAIN007      R_Brain               3.68       postmortem COLUMBIA
# # 2: LBPSEMA4BRAIN245      R_Brain               3.50       postmortem COLUMBIA
# # 3: LBPSEMA4BRAIN598      L_Brain               5.83       postmortem COLUMBIA


allpmi[SAMPLE_ISMMS == "LBPSEMA4BRAIN007", cold_pmi_CORRECTED:=pmix[SAMPLE_ISMMS=="LBPSEMA4BRAIN007"]$cold_pmi_CORRECTED]
allpmi[SAMPLE_ISMMS == "LBPSEMA4BRAIN245", cold_pmi_CORRECTED:=pmix[SAMPLE_ISMMS=="LBPSEMA4BRAIN245"]$cold_pmi_CORRECTED]
allpmi[SAMPLE_ISMMS == "LBPSEMA4BRAIN598", cold_pmi_CORRECTED:=pmix[SAMPLE_ISMMS=="LBPSEMA4BRAIN598"]$cold_pmi_CORRECTED]



#####################################################################################################

pmix[SAMPLE_ISMMS=="LBPSEMA4BRAIN535", ] #Empty data.table; pmix is 230 pm samples but there is a total of 246 pm samples???
#confirm that the excluded samples are the outliers that were removed from qc 

pm <- finalcov[mymet_postmortem==1]$SAMPLE_ISMMS
pmitable <- pmix$SAMPLE_ISMMS

notinpmi <-pm[!pm %in% pmitable]

#  [1] "LBPSEMA4BRAIN043" "LBPSEMA4BRAIN176" "LBPSEMA4BRAIN178" "LBPSEMA4BRAIN258"
#  [5] "LBPSEMA4BRAIN285" "LBPSEMA4BRAIN400" "LBPSEMA4BRAIN406" "LBPSEMA4BRAIN407"
#  [9] "LBPSEMA4BRAIN474" "LBPSEMA4BRAIN490" "LBPSEMA4BRAIN535" "LBPSEMA4BRAIN570"
# [13] "LBPSEMA4BRAIN639" "LBPSEMA4BRAIN711" "LBPSEMA4BRAIN717" "LBPSEMA4BRAIN779"

notinpmiids <- finalcov[SAMPLE_ISMMS %in% notinpmi]$IID_ISMMS

#  [1] "T-5373"         "T-5320"         "T-4815"         "T-108"         
#  [5] "T-5354"         "T-4924"         "T-156"          "T-4542"        
#  [9] "T-5259"         "T-52"           "HBHJ_17_001_FC" "T-4518"        
# [13] "T-5293"         "T-139"          "T-107"          "T-16"   


pmichart$Cold.PMI[pmichart$Sample.ID %in% notinpmiids]

#  [1] "NOT PROVIDED" "NOT PROVIDED" "NOT PROVIDED" "NOT PROVIDED" "NOT PROVIDED"
#  [6] "NOT PROVIDED" "NOT PROVIDED" "NaN"          "NOT PROVIDED" "NOT PROVIDED"
# [11] "NaN"          "NOT PROVIDED" "NOT PROVIDED" "NOT PROVIDED" "NOT PROVIDED"
# [16] "N/A"

#ok so pmix doesn't have any samples that don't have pmi? weird

length(which(is.na(pmichart$Cold.PMI)))
#[1] 0

#####################################################################################################
pminas2 <- allpmi[is.na(allpmi$cold_pmi_CORRECTED), ]

#         SAMPLE_ISMMS mymet_tissue cold_pmi_CORRECTED mymet_postmortem     Bank
#  1: LBPSEMA4BRAIN043      R_Brain                 NA       postmortem COLUMBIA
#  2: LBPSEMA4BRAIN176      R_Brain                 NA       postmortem COLUMBIA
#  3: LBPSEMA4BRAIN178      R_Brain                 NA       postmortem COLUMBIA
#  4: LBPSEMA4BRAIN258      R_Brain                 NA       postmortem COLUMBIA
#  5: LBPSEMA4BRAIN285      L_Brain                 NA       postmortem COLUMBIA
#  6: LBPSEMA4BRAIN400      R_Brain                 NA       postmortem COLUMBIA
#  7: LBPSEMA4BRAIN406      L_Brain                 NA       postmortem COLUMBIA
#  8: LBPSEMA4BRAIN407      R_Brain                 NA       postmortem COLUMBIA
#  9: LBPSEMA4BRAIN474      L_Brain                 NA       postmortem COLUMBIA
# 10: LBPSEMA4BRAIN490      L_Brain                 NA       postmortem COLUMBIA
# 11: LBPSEMA4BRAIN535      R_Brain                 NA       postmortem    MIAMI
# 12: LBPSEMA4BRAIN570      R_Brain                 NA       postmortem COLUMBIA
# 13: LBPSEMA4BRAIN639      L_Brain                 NA       postmortem COLUMBIA
# 14: LBPSEMA4BRAIN711      L_Brain                 NA       postmortem COLUMBIA
# 15: LBPSEMA4BRAIN717      L_Brain                 NA       postmortem COLUMBIA
# 16: LBPSEMA4BRAIN779      L_Brain                 NA       postmortem COLUMBIA


# allpmix <- merge(allpmi, lbpcov$IID_ISMMS, by="SAMPLE_ISMMS", all=TRUE)
# How do you expect to merge by SAMPLE_ISMMS when you're only selecting IID_ISMMS?

#take out the selected columns before merging, i.e., SAMPLE_ISMMS and IID_ISMMS
allpmix <- merge(allpmi, lbpcov[, c(1,6)], by="SAMPLE_ISMMS", all=TRUE)

#merge(dt1[,.(keepme1, keepme2, keepme3)], dt2[,.(keepme1, keepme2, keepme4)]) 
#will automatically merge by the column names present in dt1 and dt2 (so, keepme1 and keepme2)

allpmix[IID_ISMMS %in% notinpmiids, ]

pmichart$Cold.PMI[pmichart$Sample.ID %in% notinpmiids]

pmichart[pmichart$Sample.ID %in% notinpmiids, colnames(pmichart) %in% c("Sample.ID", "Cold.PMI", "Bank")]

#          Sample.ID     Bank     Cold.PMI
# 2           T-5373 COLUMBIA NOT PROVIDED
# 4           T-5354 COLUMBIA NOT PROVIDED
# 6           T-4924 COLUMBIA NOT PROVIDED
# 8           T-5320 COLUMBIA NOT PROVIDED
# 10          T-4542 COLUMBIA NOT PROVIDED
# 14          T-4518 COLUMBIA NOT PROVIDED
# 50            T-16 COLUMBIA NOT PROVIDED
# 51            T-52 COLUMBIA          NaN
# 54           T-107 COLUMBIA NOT PROVIDED
# 55           T-108 COLUMBIA NOT PROVIDED
# 65           T-139 COLUMBIA          NaN
# 71           T-156 COLUMBIA NOT PROVIDED
# 130         T-4815 COLUMBIA NOT PROVIDED
# 132         T-5259 COLUMBIA NOT PROVIDED
# 133         T-5293 COLUMBIA NOT PROVIDED
# 353 HBHJ_17_001_FC    MIAMI          N/A

#confirmed that these samples do not have PMI available

lbpcov3 <- merge(lbpcov2, allpmix[, c(1,3)], by="SAMPLE_ISMMS", all=TRUE)

summary(comparedf(lbpcov2, lbpcov3))


lbpcov3[is.na(Bank)]$SAMPLE_ISMMS #character(0)


saveRDS(lbpcov3, file = "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMI_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_18JAN2022.RDS")

saveRDS(allpmix, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_PMITable_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_18JAN2022.RDS")


# allpmi[is.na(Bank)]$SAMPLE_ISMMS
# #[1] "LBPSEMA4BRAIN107"

# lbpcov[SAMPLE_ISMMS=="LBPSEMA4BRAIN107"]$mymet_living
# # [1] 1
# # Levels: 0 1

# allpmi[SAMPLE_ISMMS=="LBPSEMA4BRAIN107", Bank:="LIVING"]

# length(which(is.na(allpmi$Bank)))
# #[1] 0

# allpmi2 <- allpmi[, c(1,3)]


# lbpcov2 <- merge(lbpcov, allbank2, by="SAMPLE_ISMMS")

# identical(lbpcov2$Bank, allbank2$Bank)
# #[1] TRUE


#NOT DONE; will be done after cibersortx
#Last step after adding cell type prop, PMI, and Bank: 
########################################################################################
#Mislabeled samples from Alex's identity concordance
# LBPSEMA4BRAIN658 is PT−0163 (not PT−0165)
# LBPSEMA4BRAIN327 is PT−0079 (not PT−0015)
# LBPSEMA4BLOOD445 is PT−0073 (not PT−0076)

#remove the mislabeled samples 
remove <- c("LBPSEMA4BLOOD445", "LBPSEMA4BRAIN658", "LBPSEMA4BRAIN327")

lbpcov <- lbpcov[! SAMPLE_ISMMS %in% remove, ]
dim(lbpcov)
#[1] 773 155

fc <- fc[, match(lbpcov$SAMPLE_ISMMS, colnames(fc))]

dim(fc)
#[1] 58929   773

identical(lbpcov$SAMPLE_ISMMS, colnames(fc))
#[1] TRUE

#make the rownames of lbpcov the SAMPLE_ISMMS ids in order for the formula to work
rownames(lbpcov) <- lbpcov$SAMPLE_ISMMS


