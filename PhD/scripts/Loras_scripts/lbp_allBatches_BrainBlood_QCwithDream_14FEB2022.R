#how to consider interactions in this data? since samples were processed differently based on being brain or blood, seems like there might be interaction effects bwn tech covs and samples
#since bbstatus is fully correlated with PC1 of the gene expression data, need to determine the correlation theshold that doesn't erase the signal
#which genes are NOT de bwn brain and blood? very few based on the vobject, does this make sense biologically based on which ones they are i.e., housekeeping etc -- do GO on those genes


module load R/4.0.3 #very important that you load this R version, dream will not work with the regular version that auto loads on minerva
R

rm(list=ls())
options(stringsAsFactors=F)
##################################################
#load libraries
##################################################

.libPaths(c("~/.Rlib", .libPaths())) #  need to do this before loading any packages!!!

library(edgeR)
library(BiocParallel)
library(batchtools)
library(devtools)
library(withr)
library(limma)
library(Glimma)
library(variancePartition)#, lib.loc="~/.Rlib")
library(ggplot2)
library(gridExtra)
library(grid)
library(doParallel)
registerDoParallel(20)
library(sp)
library(biomaRt)
library(gsubfn)
library(data.table)
library(sp)
library(Matrix)

#############################################################################################################
#Run the functions first
#############################################################################################################

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



#options(width=150)


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

# lbpcov <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_MYMET_QCMetrics_PlusCellTypeProportions_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_10MAY2021.RDS")

# fc <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_onlyBRAIN_532samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7MAY2021.RDS")

lbpcov <- readRDS(file = "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMICellTypes_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_11FEB2022.RDS")

fc <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")


dim(lbpcov)
#[1] 776 173

#make the rownames of lbpcov the SAMPLE_ISMMS ids in order for the de formula below to work
rownames(lbpcov) <- lbpcov$SAMPLE_ISMMS

#need to add brain-blood status covariate 

lbpcov$bbstatus <- lbpcov$mymet_tissue
lbpcov$bbstatus <- gsub("L_", "", lbpcov$bbstatus)
lbpcov$bbstatus <- gsub("R_", "", lbpcov$bbstatus)

lbpcov$bbstatus <- as.factor(lbpcov$bbstatus)
#lbpcov$bbstatus <- droplevels(lbpcov$bbstatus)


table(lbpcov$bbstatus)
# Blood Brain 
#   243   533 

table(lbpcov$mymet_tissue)
# L_Blood L_Brain R_Blood R_Brain 
#     128     284     115     249


#remove the PM samples and variables relating to the PM samples
lbpcov2 <- lbpcov[mymet_postmortem=="living", ]


table(lbpcov$mymet_postmortem)
    # living postmortem 
    #    530        246

table(lbpcov2$mymet_postmortem)
	 # living 
  #  	 530

table(lbpcov2$mymet_living)


dim(lbpcov2)
#[1] 530 174

lbpcov2$mymet_postmortem <- NULL

lbpcov2$mymet_living <- NULL


unique(lbpcov2$Bank)
#[1] LIVING
#Levels: COLUMBIA HARVARD MIAMI LIVING

unique(lbpcov2$cold_pmi_CORRECTED)
#[1] 0

unique(lbpcov2$mymet_bank)

lbpcov2$Bank <- NULL
lbpcov2$cold_pmi_CORRECTED <- NULL
lbpcov2$mymet_bank <- droplevels(lbpcov2$mymet_bank)

cellcounts <- tail(colnames(lbpcov2), n=17)
cellcounts <- cellcounts[1:16] #remove the last one, bbstatus

saveRDS(lbpcov2, file = "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusCellTypes_530LIVINGSamples_14FEB2022.RDS")

#same samples in featureCounts
dim(fc)
#[1] 58929   777

head(rownames(fc))
#[1] "ENSG00000000003.14" "ENSG00000000005.6"  "ENSG00000000419.12"

fc$Geneid <- NULL
#lbpcov <- lbpcov[match(colnames(fc), lbpcov$SAMPLE_ISMMS),]
fc <- fc[, match(lbpcov2$SAMPLE_ISMMS, colnames(fc))]

dim(fc)
#[1] 58929   530

identical(lbpcov2$SAMPLE_ISMMS, colnames(fc))
#[1] TRUE

saveRDS(fc, "/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_530LIVINGsamples_15FEB2022.RDS")

#need to account the cell count variables somehow and also make neuronal prop into a variable; later, since we won't be correcting for any of these variables
###########################################################################################
#canCorPairs Table -- correlations between all covariates 
###########################################################################################

library(variancePartition)
#form=as.formula(paste("~",paste(colnames(lbpcov)[colnames(lbpcov) !="SAMPLE_ISMMS"],collapse="+")))
#C = canCorPairs(form,lbpcov[,colnames(lbpcov)!= "SAMPLE_ISMMS", with=F])

#remove the cell count variables from the cancors because blood samples don't have brain estimates and vice versa so there's missing data for those columns
remove <- c("SAMPLE_ISMMS", cellcounts)

form=as.formula(paste("~",paste(colnames(lbpcov)[!colnames(lbpcov) %in% remove],collapse="+")))
C = canCorPairs(form,lbpcov[,!colnames(lbpcov) %in% remove, with=F])


pdf("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_canCorPairs_Covariates_14FEB2022.pdf", height=15,width=15)
plotCorrMatrix(C ,margins = c(25, 25))
dev.off()


  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_canCorPairs_Covariates_14FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC

##########################################################################################################

#No reason to do this other than that the code below all uses these variable names and it is much easier and safer to just rename the variables than to change the code below 

metadata <- lbpcov2
countMatrix <- fc
###########################################################################################################
#No Covs -- start with no covariates in the formula, except for case/control status, and the individual id (the random effect)
###########################################################################################################

form <- ~ (1|IID_ISMMS) + bbstatus 
resform <- ~ (1|IID_ISMMS) #this is the formula you will use for the residuals, it does not include the case/control status because this is the effect you want to preserve, not subtract from the model


isexpr <- rowSums(cpm(countMatrix)>=1) >= 0.1*ncol(countMatrix)
# Standard usage of limma/voom
geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

dim(geneExpr)
#[1] 23202   530

library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 20) #this is to make the functions run -- you want this number and the number below in MultiCoreParam() to multiply to 100 
#Sys.setenv(OMP_NUM_THREADS = 10)

identical(rownames(metadata), colnames(geneExpr))
rownames(metadata) <- metadata$SAMPLE_ISMMS


vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM = MulticoreParam(5))


identical(metadata$SAMPLE_ISMMS, colnames(vobjDream$E))


#########################################Not done now but will be relevant later
vp = fitExtractVarPartModel( vobjDream, form, metadata, BPPARAM = MulticoreParam(5))

pdf(file = "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_vobjDream_NoCovs_variancePartitionPlot_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_27MAY2021.pdf")
plotVarPart( sortCols(vp))
dev.off()

scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_vobjDream_NoCovs_variancePartitionPlot_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_27MAY2021.pdf ~/Desktop/pca_plots/lbp_allBatches_QC
#########################################

fitmm = dream( vobjDream, form, metadata, BPPARAM = MulticoreParam(5)) 

#########################################Can't do for first pass with no fixed effects in model
resfit = dream(vobjDream, resform, metadata, BPPARAM = MulticoreParam(5), computeResiduals = TRUE) #calculating the residuals
#with just the IID variable, get this warning: No testable fixed effects were included in the model.
res <- residuals(resfit)
#########################################Can't do for first pass 

lmgroup_DE <- fitmm #ignoring the eBayes step

summary(lmgroup_DE)
head(lmgroup_DE$coefficients)

coefcol <- "bbstatusBrain" #the name of your case/control status column; it's a two-level factor so it's dummy coded like this with blood being 0 and brain being 1 
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)

de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]
de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

length(which(de$DEG == "DEG"))
#[1] 22324 -- whoa that's a lot 

#In this case, since we have no resids because there is no resform because there are no covariates added, we are just looking at the means in the vobject (the vobject is just the normalized gene expression to start, it does not take into account any covariates)
#so, in order to not have to change this code below, set res <- vobjDream$E for this case ONLY, and moving forward do NOT confuse the residuals with normalized gene expression object

res <- vobjDream$E #just so that the code below doesn't need to have "res" changed to "vobjDream$E" but these two are NOT equivalent things, just in this case you don't have residuals yet so you are just looking at the average normalized gene expression for each gene 

brain <- metadata[bbstatus=="Brain"]$SAMPLE_ISMMS
blood <- metadata[bbstatus=="Blood"]$SAMPLE_ISMMS
tmp <- data.table(gene=names(rowMeans(res)), computedAvgExp=rowMeans(res))
de <- merge(de, tmp)
blmn <- data.table(gene=names(rowMeans(res[,blood])), mean_blood=rowMeans(res[,blood]))
brmn <- data.table(gene=names(rowMeans(res[,brain])), mean_brain=rowMeans(res[,brain]))
de <- merge(merge(de, blmn), brmn)
de[,DIFF:=mean_brain-mean_blood] #for each gene, the difference between its mean expression in the brain and blood samples 


pdf("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_Dream_NoCovs_DIFFPlot_15FEB2022.pdf")
    # ggplot(de, aes(mean_liv, mean_pm)) + geom_point(size=3, pch=21) + facet_wrap(~DEG+LFC) 
  ggplot(de, aes(DIFF, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
dev.off()

  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_Dream_NoCovs_DIFFPlot_15FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC


###############################################
#calculating the PCs and generating the covariates/PC correlation table and the PCA plots 
###############################################

#remove the columns with missing data otherwise this isn't going to work 
DATA <- metadata[, !colnames(metadata) %in% cellcounts, with=FALSE]
plotpath = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/" 
vobj <- vobjAllGenes

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

SampleByVariable=t(cov(vobj$E)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)


#The covariates-PCs correlation table
resCor=canCorAllAgainstAll_Original(metadata,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "bbstatus"])
colnames(C2)[1] <- "bbstatus_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),] #ordered by highest correlation with PC1
head(all,40)

#The giant covariates-PCs file with all of the PCA plots
pdf(paste(plotpath, "lbp_allBatches_BloodBrain_QC_pca-plot_noCovariates_15FEB2022.pdf",sep=""))
count=1
  total=ncol(DATA)
  for(col in colnames(DATA)){
    cat("on column",count,"/",total,col,"\n")
  # #=======pca-1 vs pca-2=======
    if(is.numeric(DATA[[col]])==F){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
  if(length(unique(DATA[[col]]))<7){
    a = a + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(color = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    a = a + scale_color_discrete(guide =FALSE)+ labs(color=c(substring(col,1,8)))
  }
  build <- ggplot_build(a)$data
  points <- build[[1]]
  ell <- build[[3]]


  # Find which points are inside the ellipse, and add this to the data
  dat <- data.frame(points[1:2], 
                    in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))
  outliers_3SD_PC1_PC2=points$label[which(dat$in.ell==F)]
  if(is.empty.vector(outliers_3SD_PC1_PC2)==F){
    
  }
  # show(a)
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    b = b + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    b = b + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    c = c + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    c = c + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +
    ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    d = d + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    d = d + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
}else if(is.numeric(DATA[[col]])){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=DATA[[col]][!is.na(DATA[[col]])]))+geom_point(size=3) +
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  a = a + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  b = b + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  c = c + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  d = d + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
}
# show(a)
# show(b)
# show(c)
# show(d)
multiplot_same_legend(a,b,c,d,cols=2) 
count=count+1
}

dev.off()

#   outliers_3SD_PC1_PC2
# }

  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_BloodBrain_QC_pca-plot_noCovariates_15FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC


vobjAllGenes <- vobjDream
deAllGenes <- de
#resAllGenes <- res

notdeg <- deAllGenes[deAllGenes$DEG == "NOTDEG"]$gene
length(notdeg)
#[1] 878

##########################################################################################################
#Based on this data we see that PC1 is entirely tissue so it is not possible to minimize the correlation bwn tech covs and tissue so regress out effect of tissue and correct that way?

form <- ~ (1|IID_ISMMS) + bbstatus 

#just do resids with the formula that we have since we want to subtract out bbstatus just to see what happens
resfit = dream(vobjAllGenes, form, metadata, BPPARAM = MulticoreParam(5), computeResiduals = TRUE) #calculating the residuals
#with just the IID variable, get this warning: No testable fixed effects were included in the model.
res <- residuals(resfit)

DATA <- metadata[, !colnames(metadata) %in% cellcounts, with=FALSE]
plotpath = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/" 
vobj <- vobjAllGenes

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

SampleByVariable=t(cov(res)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)


#The covariates-PCs correlation table
resCor=canCorAllAgainstAll_Original(metadata,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "bbstatus"])
colnames(C2)[1] <- "bbstatus_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),] #ordered by highest correlation with PC1
head(all,40)

head(summary(pca)$importance[3,])
#     PC1     PC2     PC3     PC4     PC5     PC6 
# 0.65888 0.86042 0.91882 0.94031 0.95486 0.96448 


pdf(paste(plotpath, "lbp_allBatches_BloodBrain_QC_pca-plot_IID_ISMMS_bbstatus_17FEB2022.pdf",sep=""))
count=1
  total=ncol(DATA)
  for(col in colnames(DATA)){
    cat("on column",count,"/",total,col,"\n")
  # #=======pca-1 vs pca-2=======
    if(is.numeric(DATA[[col]])==F){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
  if(length(unique(DATA[[col]]))<7){
    a = a + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(color = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    a = a + scale_color_discrete(guide =FALSE)+ labs(color=c(substring(col,1,8)))
  }
  build <- ggplot_build(a)$data
  points <- build[[1]]
  ell <- build[[3]]


  # Find which points are inside the ellipse, and add this to the data
  dat <- data.frame(points[1:2], 
                    in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))
  outliers_3SD_PC1_PC2=points$label[which(dat$in.ell==F)]
  if(is.empty.vector(outliers_3SD_PC1_PC2)==F){
    
  }
  # show(a)
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    b = b + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    b = b + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    c = c + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    c = c + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +
    ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    d = d + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    d = d + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
}else if(is.numeric(DATA[[col]])){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=DATA[[col]][!is.na(DATA[[col]])]))+geom_point(size=3) +
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  a = a + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  b = b + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  c = c + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  d = d + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
}
# show(a)
# show(b)
# show(c)
# show(d)
multiplot_same_legend(a,b,c,d,cols=2) 
count=count+1
}

dev.off()

#   outliers_3SD_PC1_PC2
# }

  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_BloodBrain_QC_pca-plot_IID_ISMMS_bbstatus_17FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC

##########################################################################################################

head(summary(pca)$importance[3,])
#     PC1     PC2     PC3     PC4     PC5     PC6 
# 0.91605 0.99373 0.99869 0.99914 0.99944 0.99961 

#PC1 explains 91.6% of the variance in the gene expression data and it is pretty much fully colinear with bbstatus

#New approach: try to define the genes that are common to both

# To Noam, first question is what is the right threshold for expression when considering both the brain and blood data -- percentage of samples that you need, which genes have at least 1CPM in 75% of the samples, together; use absurdly high percentage 

# Start by figuring out this set of genes that are the least variable between brain and blood

# Take set of genes expressed in X% of samples and then follow the same strategy as for liv-pm to start; look at the pc plots to see what happens when you do that

isexpr2 <- rowSums(cpm(countMatrix)>=1) >= 0.75*ncol(countMatrix)
# Standard usage of limma/voom
geneExpr2 = DGEList( countMatrix[isexpr2,] )
geneExpr2 = calcNormFactors( geneExpr2 )

dim(geneExpr2)
#[1] 13194   530 #Try it with this threshold it's way too few genes obv but just to see if it changes anything and then go from there 


form <- ~ (1|IID_ISMMS) + bbstatus
resform <- ~ (1|IID_ISMMS) #this is the formula you will use for the residuals, it does not include the case/control status because this is the effect you want to preserve, not subtract from the model


library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 20) #this is to make the functions run -- you want this number and the number below in MultiCoreParam() to multiply to 100 
#Sys.setenv(OMP_NUM_THREADS = 10)

identical(rownames(metadata), colnames(geneExpr2))
rownames(metadata) <- metadata$SAMPLE_ISMMS


vobjDream = voomWithDreamWeights( geneExpr2, form, metadata, BPPARAM = MulticoreParam(5))


identical(metadata$SAMPLE_ISMMS, colnames(vobjDream$E))

fitmm = dream( vobjDream, form, metadata, BPPARAM = MulticoreParam(5)) 

lmgroup_DE <- fitmm #ignoring the eBayes step

summary(lmgroup_DE)
head(lmgroup_DE$coefficients)

coefcol <- "bbstatusBrain" #the name of your case/control status column; it's a two-level factor so it's dummy coded like this with blood being 0 and brain being 1 
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)

de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]
de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

length(which(de$DEG == "DEG"))
#[1] 12574 out of 13194

#In this case, since we have no resids because there is no resform because there are no covariates added, we are just looking at the means in the vobject (the vobject is just the normalized gene expression to start, it does not take into account any covariates)
#so, in order to not have to change this code below, set res <- vobjDream$E for this case ONLY, and moving forward do NOT confuse the residuals with normalized gene expression object

res <- vobjDream$E #just so that the code below doesn't need to have "res" changed to "vobjDream$E" but these two are NOT equivalent things, just in this case you don't have residuals yet so you are just looking at the average normalized gene expression for each gene 

brain <- metadata[bbstatus=="Brain"]$SAMPLE_ISMMS
blood <- metadata[bbstatus=="Blood"]$SAMPLE_ISMMS
tmp <- data.table(gene=names(rowMeans(res)), computedAvgExp=rowMeans(res))
de <- merge(de, tmp)
blmn <- data.table(gene=names(rowMeans(res[,blood])), mean_blood=rowMeans(res[,blood]))
brmn <- data.table(gene=names(rowMeans(res[,brain])), mean_brain=rowMeans(res[,brain]))
de <- merge(merge(de, blmn), brmn)
de[,DIFF:=mean_brain-mean_blood] #for each gene, the difference between its mean expression in the living and postmortem samples 


pdf("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_Dream_NoCovs_DIFFPlot_16FEB2022.pdf")
    # ggplot(de, aes(mean_liv, mean_pm)) + geom_point(size=3, pch=21) + facet_wrap(~DEG+LFC) 
  ggplot(de, aes(DIFF, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
dev.off()

  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_RAPiD_BloodBrainQC_onlyLIVING_530Samples_Dream_NoCovs_DIFFPlot_15FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC


###############################################
#calculating the PCs and generating the covariates/PC correlation table and the PCA plots 
###############################################

#remove the columns with missing data otherwise this isn't going to work 
DATA <- metadata[, !colnames(metadata) %in% cellcounts, with=FALSE]
plotpath = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/" 
vobj <- vobjDream

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

SampleByVariable=t(cov(vobj$E)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)


#The covariates-PCs correlation table
resCor=canCorAllAgainstAll_Original(metadata,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "bbstatus"])
colnames(C2)[1] <- "bbstatus_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),] #ordered by highest correlation with PC1
head(all,40)

#The giant covariates-PCs file with all of the PCA plots
pdf(paste(plotpath, "lbp_allBatches_BloodBrain_QC_pca-plot_IID_ISMMS_bbstatus_17FEB2022.pdf",sep=""))
count=1
  total=ncol(DATA)
  for(col in colnames(DATA)){
    cat("on column",count,"/",total,col,"\n")
  # #=======pca-1 vs pca-2=======
    if(is.numeric(DATA[[col]])==F){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
  if(length(unique(DATA[[col]]))<7){
    a = a + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(color = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    a = a + scale_color_discrete(guide =FALSE)+ labs(color=c(substring(col,1,8)))
  }
  build <- ggplot_build(a)$data
  points <- build[[1]]
  ell <- build[[3]]


  # Find which points are inside the ellipse, and add this to the data
  dat <- data.frame(points[1:2], 
                    in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))
  outliers_3SD_PC1_PC2=points$label[which(dat$in.ell==F)]
  if(is.empty.vector(outliers_3SD_PC1_PC2)==F){
    
  }
  # show(a)
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    b = b + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    b = b + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    c = c + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    c = c + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +
    ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    d = d + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    d = d + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
}else if(is.numeric(DATA[[col]])){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=DATA[[col]][!is.na(DATA[[col]])]))+geom_point(size=3) +
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  a = a + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  b = b + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  c = c + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  d = d + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
}
# show(a)
# show(b)
# show(c)
# show(d)
multiplot_same_legend(a,b,c,d,cols=2) 
count=count+1
}

dev.off()

#   outliers_3SD_PC1_PC2
# }


scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/lbp_allBatches_BloodBrain_QC_pca-plot_noCovariates_15FEB2022.pdf ~/Desktop/pca_plots/lbp_allBatches_BloodBrain_QC

vobjDream_1cpm75p <- vobjDream
de_1cpm75p <- de
##############################################################################################
#Doing this changes nothing, the corr bwn bbstatus and PC1 gets even higher
#What if we try 10CPM in 75% of samples, how many genes is that 
##############################################################################################

isexpr3 <- rowSums(cpm(countMatrix)>=10) >= 0.75*ncol(countMatrix)
# Standard usage of limma/voom
geneExpr3 = DGEList( countMatrix[isexpr3,] )
geneExpr3 = calcNormFactors( geneExpr3 )

dim(geneExpr3)
#[1] 5901  530


geneExpr2 <- geneExpr3 

library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 20) #this is to make the functions run -- you want this number and the number below in MultiCoreParam() to multiply to 100 
#Sys.setenv(OMP_NUM_THREADS = 10)

identical(rownames(metadata), colnames(geneExpr2))

vobjDream = voomWithDreamWeights( geneExpr2, form, metadata, BPPARAM = MulticoreParam(5))


identical(metadata$SAMPLE_ISMMS, colnames(vobjDream$E))

fitmm = dream( vobjDream, form, metadata, BPPARAM = MulticoreParam(5)) 

lmgroup_DE <- fitmm #ignoring the eBayes step

# summary(lmgroup_DE)
# head(lmgroup_DE$coefficients)

coefcol <- "bbstatusBrain" #the name of your case/control status column; it's a two-level factor so it's dummy coded like this with blood being 0 and brain being 1 
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)

de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]
de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

length(which(de$DEG == "DEG"))
#[1] 5605 out of 5901


#remove the columns with missing data otherwise this isn't going to work 
DATA <- metadata[, !colnames(metadata) %in% cellcounts, with=FALSE]
plotpath = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_BloodBrain_QC/" 
vobj <- vobjDream

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

SampleByVariable=t(cov(vobj$E)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)


#The covariates-PCs correlation table
resCor=canCorAllAgainstAll_Original(metadata,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "bbstatus"])
colnames(C2)[1] <- "bbstatus_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),] #ordered by highest correlation with PC1
head(all,35)


head(summary(pca)$importance[3,])
#   PC1     PC2     PC3     PC4     PC5     PC6 
# 0.93572 0.97581 0.99042 0.99711 0.99800 0.99873 


de10CPM75p <- de
##############################################################################################
#DO GO ON THE GENES THAT ARE NOT DE
##############################################################################################

#Or a try, split into brain and blood samples, do intersect of all brain, intersect of all blood, and then intersect of the intersects to get ALL of the genes expressed in both -- this won't matter given what we've seen above

length(notdeg)
[1] 878

de <- deAllGenes 


de <- de[order(de$adj.P.Val),]
nrow(de[de$adj.P.Val<0.05,])
#[1] 22324 DE genes 


#define which genes are DE upregulated, DE downregulated, and not DE
#de.up<-de[threshold=="DEG" & logFCsign=="POS"]$gene
de.up<-de[DEG=="DEG" & LFC=="POSLFC"]$gene
de.down<-de[DEG=="DEG" & LFC=="NEGLFC"]$gene
de.not <- de[DEG=="NOTDEG"]$gene

de.background <- c(de.up, de.down)

length(de.background)
#[1] 22324


#DE results -- pos vs neg logFC separate and then do 2 GO analyses in order to get the upregulated and downregulated then do 2 pathway enrichments
#####################################################################################################
#Step 1: Make mel dataframe
#Save all the genes from my DE in first column and then 0 or 1 in the second col for if DE or not

colnames(de)
#  [1] "gene"      "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val"
#  [7] "DEG"       "threshold" "logFCsign"


de2 <- de[gene %in% de.background, DEnot:=1]
de2[!gene %in% de.background, DEnot:=0]

length(which(de2$DEnot==1))
#[1] 22324

table(de2$DEnot)
#    0    1 
# 878 22324 


de2[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1)] #removes the . and the numbers after it from the gene name for input to the getgo function below


mel <- de2[, c(1,10)]

#notdeg.go <- annotate_GOterms_DE(mel = mel.up,sample_name = sample_name, date = date, ensGene=TRUE,revigo=TRUE) #,gene.map=gene.map)
notdeg.go <- annotate_GOterms_DE(mel = mel, sample_name="NotDeg") #,gene.map=gene.map)


#THE GO FUNCTION
#############################################################################################################
annotate_GOterms_DE <- function(mel, sample_name) { #,sample_name,date, WGCNAFolder){
  # Performs GO annotation of the genes of interest
  # inputs are mel, a 2 column data.frame with number of rows corresponding
  # to your sample size, where the first column is gene name and the second
  # column is a binary vector of 0s and 1s, with 0s for genes in background
  # (not significantly differentially expressed genes for example) and 1s
  # for genes of interest (significantly differentially expressed genes for
  # example); sample_name is a string containing the name of the file that 
  # will be saved with the output of the analysis; and date, a string of the
  # date of the experiment that will also be used in the name of the output
  # file. The function has no output, output is directly saved as a table.
  options(stringsAsFactors = FALSE)
  library("R.matlab")
  library(goseq)
  library(topGO)
  library(org.Hs.eg.db)
  library(Rgraphviz)
  library(gplots)
  
  # load in the data set
  #######################################################################################
  # now carry out GO analysis
  # read in analysis results above if needed
  #pwf = nullp(a, 'hg19', 'geneSymbol')
  
  mel <-as.data.frame(mel) #this MUST be a dataframe, cannot be a datatable, it won't work
  net = mel

  
  # gene.map = getgo(net[,1],'hg19','geneSymbol')       # extracting go IDs from gene names
  genemap = getgo(net[,1],'hg38','ensGene')       # extracting go IDs from gene names
  genemap = genemap[!is.na(names(genemap))]        # removing non existing entries
  number_of_GO = length(unique(unlist(genemap)))
  
  # if (exists("res")){
  #   remove("res")
  # }
  # counter = 0
  # counter = counter + 1
  # flush.console()
  
  x = net[,1]         # extracting gene names
  a = rep(0,length(x))    # creating vector of length module names filled with 0s and name of row gene name
  names(a) = x
  a[net[,2] == 1] = 1     # changing indices in this vector at module positions to 1, vector used in definition of GO element
  
  ipsBP = new("topGOdata", description = paste("Enrichment in ", sample_name, sep = ""), # creating GO element and graph
              ontology = c("BP"),
              # ontology = c("BP", "MF", "CC", "KEGG"),
              allGenes = as.factor(a), 
              geneSel = names(a[a==1]),
              nodeSize = 10,
              annot = annFUN.gene2GO,
              gene2GO = genemap)
  
  teststat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
  resfisherBP = getSigGroups(ipsBP, teststat)   # doing Fisher statistics on GO graph
  resfinalBP = GenTable(ipsBP, classic = resfisherBP, topNodes=length(ipsBP@graph@nodes))   # selecting top 200 of resfisher
  
  ipsMF = new("topGOdata", description = paste("Enrichment in ", sample_name, sep = ""), # creating GO element and graph
              ontology = c("MF"),
              # ontology = c("BP", "MF", "CC", "KEGG"),
              allGenes = as.factor(a), 
              geneSel = names(a[a==1]),
              nodeSize = 10,
              annot = annFUN.gene2GO,
              gene2GO = genemap)
  
  teststat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
  resfisherMF = getSigGroups(ipsMF, teststat)   # doing Fisher statistics on GO graph
  resfinalMF = GenTable(ipsMF, classic = resfisherMF, topNodes=length(ipsMF@graph@nodes))   # selecting top 200 of resfisher
  
  ipsCC = new("topGOdata", description = paste("Enrichment in ", sample_name, sep = ""), # creating GO element and graph
              ontology = c("CC"),
              # ontology = c("BP", "MF", "CC", "KEGG"),
              allGenes = as.factor(a), 
              geneSel = names(a[a==1]),
              nodeSize = 10,
              annot = annFUN.gene2GO,
              gene2GO = genemap)
  
  teststat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
  resfisherCC = getSigGroups(ipsCC, teststat)   # doing Fisher statistics on GO graph
  resfinalCC = GenTable(ipsCC, classic = resfisherCC, topNodes=length(ipsCC@graph@nodes))   # selecting top 200 of resfisher
  
  # genelist=get_genes_in_signif_pathways_GO_ensembl(resfinal[,"GO.ID"],genemap)
  
  res=rbind(resfinalBP, resfinalMF, resfinalCC)
  
  
  # if (exists("res")){
  #     res = rbind(res, cbind(resfinal))
  # }else{
  #     res = cbind(resfinal) # creating res table with 1st column module color and rest resfinal
  # }
  
  resmod = cbind(res, res[,"Significant"]/res[,"Expected"],p.adjust(res[,"classic"],method="BH",n=number_of_GO))  # adding a column of ratio of significant over expected (fold enrichment)
  names(resmod)[c(7,8)] = c("fold_enrichment","BH")
  resmod = resmod[,c( "GO.ID", "Term", "Annotated", "Significant", "Expected", "fold_enrichment", "classic","BH")]
  
  resmod = resmod[order(resmod$BH),]
  row.names(resmod) = 1:nrow(resmod)
  
  # directory = paste(WGCNAFolder, "/", sample_name, "/", sep = "")
  # if(!dir.exists(file.path(directory, sep = ""))){
  #   dir.create(file.path(directory, sep = ""))
  #}
  
  #write.csv(resmod, paste(directory, sample_name, "_GO_enrichments", ".csv",sep=""))
  # write.table(resmod, output.annotation.file, sep="\t", quote=FALSE, row.names=FALSE)
  #For Revigo
   Revigo = resmod[which(resmod$BH <= 0.05),]
  return(resmod)
  #return(Revigo)
  #write.csv(Revigo, paste(directory, sample_name, "_REVIGO", ".csv",sep=""))
  # resmod
}

##########################################################################################################
#NEXT -- try without accounting for individual ID
##########################################################################################################
#OR can we do linear mixed model with tissue type as a random effect

library(ggplot2)

# Scatter plot by group
ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point()

ggplot(aes(x='x', y='y'), data=df1) + geom_point() + 
       geom_point(aes(x='x', y='y'), data=df2)


#Difference between brain and blood FOR EACH individual?























##########################################################################################################
#Everything -- this is all covariates I selected being added to the model, but you will do this iteratively starting with no covariates which you did in the code above and now you will add your first covariate -- sex, and then you'll add your first and second -- sex and rin, and then you'll add your first, second, and third -- sex, rin, X, and so forth. 
###########################################################################################################

form <- ~(1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + RNASeqMetrics_PCT_MRNA_BASES + mymet_postmortem + (1|IID_ISMMS) + (1|mymet_depletionbatch)

resform <- ~(1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch)


isexpr <- rowSums(cpm(countMatrix)>=1) >= 0.1*ncol(countMatrix)
# Standard usage of limma/voom
geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

dim(geneExpr)

library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 20)
#Sys.setenv(OMP_NUM_THREADS = 10)

#fitDream <- dream(..., BPPARAM = MulticoreParam(5))

#per Noam, there should be no big difference between these two vobjects because the formula is not used for the vobject, the vobject normalizes the distribution
vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM = MulticoreParam(5))

identical(metadata$SAMPLE_ISMMS, colnames(vobjDream$E))

vp = fitExtractVarPartModel( vobjDream, form, metadata, BPPARAM = MulticoreParam(5))


pdf(file = "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_vobjDream_allCovsPlusDepletionBatch_variancePartitionPlot_onlyBRAIN_531Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_NoBRAIN364-Outlier_14JUN2021.pdf")
plotVarPart( sortCols(vp))
dev.off()


scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_vobjDream_allCovsPlusDepletionBatch_variancePartitionPlot_onlyBRAIN_531Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_NoBRAIN364-Outlier_14JUN2021.pdf ~/Desktop/pca_plots/lbp_allBatches_QC

fitmm = dream( vobjDream, form, metadata, BPPARAM = MulticoreParam(5))

resfit = dream(vobjDream, resform, metadata, BPPARAM = MulticoreParam(5), computeResiduals = TRUE)
res <- residuals(resfit)


lmgroup_DE <- fitmm #ignoring the eBayes step

coefcol <- "mymet_postmortem"
#coefcol <- which(grepl("mymet_postmortem", colnames(group_ascov_design)))
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)

de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]
de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

length(which(de$DEG == "DEG"))
#[1] 17128 vs 17211 without deletion batch added as cov 


liv <- metadata[mymet_postmortem==0]$SAMPLE_ISMMS
pmt <- metadata[mymet_postmortem==1]$SAMPLE_ISMMS
tmp <- data.table(gene=names(rowMeans(res)), computedAvgExp=rowMeans(res))
de <- merge(de, tmp)
pmn <- data.table(gene=names(rowMeans(res[,pmt])), mean_pm=rowMeans(res[,pmt]))
lmn <- data.table(gene=names(rowMeans(res[,liv])), mean_liv=rowMeans(res[,liv]))
de <- merge(merge(de, pmn), lmn)
de[,DIFF:=mean_liv-mean_pm] #for each gene, the difference between its mean expression in the living and postmortem samples 


pdf("/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/lbp_allBatches_RAPiD_QC_Dream_AllCovs_DIFFPlot_NoBRAIN364-Outlier_14JUN2021.pdf")
    # ggplot(de, aes(mean_liv, mean_pm)) + geom_point(size=3, pch=21) + facet_wrap(~DEG+LFC) 
  ggplot(de, aes(DIFF, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
dev.off()

  scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/lbp_allBatches_RAPiD_QC_Dream_AllCovs_DIFFPlot_28MAY2021.pdf ~/Desktop/pca_plots/lbp_allBatches_QC


###############################################
DATA <- metadata
plotpath = "/sc/arion/projects/psychgen2/lbp/files/lbp_allBatches_QC/" 
vobj <- vobjDream

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA


covariance=cov(res) #here when you have covariates added to the model, you use the residuals not the vobject to calculate the PCs
SampleByVariable=t(covariance)
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)


resCor=canCorAllAgainstAll_Original(lbpcov,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "mymet_postmortem"])
colnames(C2)[1] <- "postmortem_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),]
head(all,40)

all.all <- all

pdf(paste(plotpath, "lbp_allBatches_RAPiD_QC_pca-plot_AllCovariates2_1JUN2021.pdf",sep=""))
count=1
  total=ncol(DATA)
  for(col in colnames(DATA)){
    cat("on column",count,"/",total,col,"\n")
  # #=======pca-1 vs pca-2=======
    if(is.numeric(DATA[[col]])==F){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
  if(length(unique(DATA[[col]]))<7){
    a = a + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(color = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    a = a + scale_color_discrete(guide =FALSE)+ labs(color=c(substring(col,1,8)))
  }
  build <- ggplot_build(a)$data
  points <- build[[1]]
  ell <- build[[3]]


  # Find which points are inside the ellipse, and add this to the data
  dat <- data.frame(points[1:2], 
                    in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))
  outliers_3SD_PC1_PC2=points$label[which(dat$in.ell==F)]
  if(is.empty.vector(outliers_3SD_PC1_PC2)==F){
    
  }
  # show(a)
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    b = b + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    b = b + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    c = c + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    c = c + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +
    ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    d = d + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    d = d + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
}else if(is.numeric(DATA[[col]])){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$SAMPLE_ISMMS
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=DATA[[col]][!is.na(DATA[[col]])]))+geom_point(size=3) +
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  a = a + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  b = b + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  c = c + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  d = d + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
}
# show(a)
# show(b)
# show(c)
# show(d)
multiplot_same_legend(a,b,c,d,cols=2) 
count=count+1
}

dev.off()

#   outliers_3SD_PC1_PC2
# }


outliers_3SD_PC1_PC2

