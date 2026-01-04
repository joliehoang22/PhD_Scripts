library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)
library(foreach)
library(readr)
library(sva)

########################################
###### MSBB with and without SVA #######
########################################

#1. LOAD AND FORMAT METADATA AND GENE EXPRESSION DATA
genedata_org_v30=fread("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/msbb/MSBBv30/expression/allcount_matrix_2023-08-29.txt", data.table=FALSE)
rownames(genedata_org_v30)=genedata_org_v30$Geneid
genedata_org_v30$Geneid=NULL
dim(genedata_org_v30)

all_var <- 	read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/msbb/run_dge_AD.txt", header = FALSE)
all_var =  all_var[all_var$V4=="v30",]

ids_to_keep=readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/msbb/metadata/ids_to_keep.RDS" )
info_all_base=readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/msbb/metadata/infoall_09.08.23.RDS")

main_path="/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/"
date="2023-10-23"

i = 1
if(length(Sys.glob(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS")))==0){
	genedata_org <- get(all_var[i,1])
	genedata=genedata_org[,6:ncol(genedata_org)]
	genedata <- as.data.frame(genedata)
	# dim(genedata)
	genedata=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
	saveRDS(get(all_var[i,1]),paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
} else {
	genedata=readRDS(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
}

print(all_var[i,2])
info_all=info_all_base
assign(all_var[i,2], info_all)
saveRDS(get(all_var[i,2]),paste0(main_path, all_var[i,3],"/covariate/info_all",sep="_",all_var[i,4], sep="_", date, ".RDS"))

info_all3 <- get(all_var[i,2])
values_to_match <- make.names(as.character(info_all3$sample))
columns_to_keep <- colnames(genedata) %in% values_to_match
df_without_nonmatching <- genedata[, columns_to_keep]
re_order <- match(make.names(as.character(info_all3$sample)), colnames(df_without_nonmatching))
colnames(df_without_nonmatching) <- colnames(df_without_nonmatching)[re_order]
assert_that(identical(colnames(df_without_nonmatching),make.names(as.character(info_all3$sample))), msg=paste("assertion 1 on iteration ",i)) #checks

genedata <- df_without_nonmatching

assert_that(length(setdiff(make.names(as.character(info_all3$sample)), colnames(genedata)))==0, msg=paste("assertion 2 on iteration ",i))
assert_that(length(setdiff(colnames(genedata),make.names(as.character(info_all3$sample))))==0, msg=paste("assertion 3 on iteration ",i))

info_all = info_all3
info_all$id = info_all$synapse_id
rownames(info_all)=make.names(as.character(info_all$sample))
assert_that(identical(colnames(genedata),rownames(info_all)), msg=paste("assertion 4 on iteration ",i)) #checks
		
setwd(paste(main_path,all_var[i,3],sep = ""))

#2. VOOM AND LIMMA ANALYSES PREPARATION AND APPLICATION TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
formVoom = formula(paste0("~ PMI + RACE + correct_SEX + RIN + Exonic.Rate + Batch_for_correct")) 
forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

info_all2=info_all
design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formVoom)[2],fixed=TRUE))),info_all2)
info_all2=info_all2[rownames(design),]
genedata=genedata[,rownames(design)]
assert_that(identical(colnames(genedata),make.names(as.character(info_all2$sample)))) 
info_all2$sample <- as.factor(info_all2$sample)

#Refining and Categorizing Covariate Data Based on Disease Status
#Categorize disease status using CERAD score, which is commonly used to assess Alzheimer's disease pathology.
info_all=info_all2
info_all$CERJ_defvsctl=substring(as.character(info_all$CERJ),1,1)
info_all$CERJ_defvsctl[as.numeric(substring(as.character(info_all$CERJ),1,1))==1]="NL" #No or low pathology
info_all$CERJ_defvsctl[as.numeric(substring(as.character(info_all$CERJ),1,1))==2]="AD" #Alzheimer's disease
info_all$CERJ_defvsctl=factor(info_all$CERJ_defvsctl)

disease_statuses=c("CERJ_defvsctl")

for (status in disease_statuses){
    info_alltmp=info_all2[which(is.na(info_all2[,status])==F),] ##there's sth wrong here
    resVPtmp=v[,which(is.na(info_all2[,status])==F)]
    for(i in 1:ncol(info_alltmp)){
        if(class(info_alltmp[,i])=="factor"){
            info_alltmp[,i]=factor(info_alltmp[,i])
        }
    }
}

table(info_alltmp$Batch_for_correct)
names(which(table(info_alltmp$Batch_for_correct)==1)) 
info_alltmp<- info_alltmp[!info_alltmp$Batch_for_correct %in% c("E007_C014","L43_C014"),]
info_alltmp<-droplevels(info_alltmp)

#Define statistical models for analysis
status="PlaqueMean"
formula_test=formula(paste("~ ", status," + PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))
formula_test0=formula(paste("~ PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))

design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

#Applying Voom Transformation 
isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata) 
dge <- DGEList(counts=genedata[isexpr,]) 
dge <- calcNormFactors(dge) 
dge <- dge[, colnames(dge) %in% rownames(design)] 
 
pdf("/sc/arion/projects/mscic1/results/jolie/my_voom_plot.pdf") 
v <- voom(dge, design, plot=TRUE) 
dev.off() 

fit <- lmFit(v, design) 

##WHY ARE THERE NO MAKE CONTRAST STEP LIKE IN ROSMAP???
fit2 <- eBayes(fit) 
results<-topTable(fit2,coef="PlaqueMean",number=Inf) 
#######might have some missing results line here!!!

#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_no_sva_1.11.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#4a. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "BE" METHOD
n.sv = num.sv(v$E,design,method="be")
if(n.sv>0){
svobj = sva(v$E,design,design0,n.sv=n.sv) #Step 1: extract surrogate variables
temp<-svobj$sv
colnames(temp)<-paste0('sv',1:ncol(temp)) #Step 2: name the surrogate variables for clarity
design_sv = cbind(design, temp) #Step 3: combine original design matrix with surrogate variables
} else {
	design_sv=design #Use the original design matrix if no SVs are identified
}

#Analyze the correlation between PlaqueMean and surrogate variables
cor(design_sv[,c("PlaqueMean",grep("sv",colnames(design_sv),value=TRUE))]) 

#5a. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv) 
fit2 <- eBayes(fit_sv)  
results_sv<-topTable(fit2,coef="PlaqueMean",number=Inf) 

#6a. SUMMARIZE AND SAVE DEA WITH SVA RESULTS
results_sv$X<-rownames(results_sv)
msbb_sva_be<- results_sv[, c(7,1,2,3,4,5,6)]
row.names(msbb_sva_be) <- 1:nrow(msbb_sva_be)


##SVA LEEK
#4b. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "LEEK" METHOD
n.sv = num.sv(v$E,design,method="leek")
if(n.sv>0){
svobj = sva(v$E,design,design0,n.sv=n.sv) #Step 1: extract surrogate variables
temp<-svobj$sv
colnames(temp)<-paste0('sv',1:ncol(temp)) #Step 2: name the surrogate variables for clarity
design_sv = cbind(design, temp) #Step 3: combine original design matrix with surrogate variables
} else {
	design_sv=design #Use the original design matrix if no SVs are identified
}

#Analyze the correlation between PlaqueMean and surrogate variables
cor(design_sv[,c("PlaqueMean",grep("sv",colnames(design_sv),value=TRUE))]) 

#5b. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv) 
fit2 <- eBayes(fit_sv)  
results_sv<-topTable(fit2,coef="PlaqueMean",number=Inf) 

#6b. SUMMARIZE AND SAVE DEA WITH SVA RESULTS
results_sv$X<-rownames(results_sv)
msbb_sva_leek<- results_sv[, c(7,1,2,3,4,5,6)]
row.names(msbb_sva_leek) <- 1:nrow(msbb_sva_leek)

write.table(msbb_sva_leek, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_leek_1.11.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(msbb_sva_be, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_be_1.11.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


########################################
##### ROSMAP with and without SVA ######
########################################

#1. LOAD AND FORMAT GENE EXPRESSION DATA
genedata_org_v30=fread("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmapv30/expression/allcount_matrix_2023-04-11.txt", data.table=FALSE)
rownames(genedata_org_v30)=genedata_org_v30$Geneid
genedata_org_v30$Geneid=NULL
dim(genedata_org_v30)

all_var <- 	read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/run_dge_AD.txt", header = FALSE)
all_var =  all_var[all_var$V4=="v30",]

ids_to_keep=readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/metadata/ids_to_keep.RDS")
info_all_base=readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/metadata/infoall_09.13.23.RDS")
info_all2=info_all_base
info_all2 <- info_all2[info_all2$ceradsc %in% c(1, 4), ]
info_all2$ceradsc <- droplevels(info_all2$ceradsc)

main_path="/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/"
date="2025-11-19"

i = 1
if(length(Sys.glob(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS")))==0){
	genedata_org <- get(all_var[i,1])
	genedata=genedata_org[,6:ncol(genedata_org)]
	genedata <- as.data.frame(genedata)
	# dim(genedata)
	genedata=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
	saveRDS(get(all_var[i,1]),paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
} else {
	genedata=readRDS(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
}

#covariate
print(all_var[i,2])
info_all=info_all_base
assign(all_var[i,2], info_all)
saveRDS(get(all_var[i,2]),paste0(main_path, all_var[i,3],"/covariate/info_all",sep="_",all_var[i,4], sep="_", date, ".RDS"))

##add covariate info
info_all3 <- info_all

values_to_match <- make.names(as.character(info_all3$sample)) 
columns_to_keep <- colnames(genedata) %in% values_to_match 
df_without_nonmatching <- genedata[, columns_to_keep] 
re_order <- match(make.names(as.character(info_all3$sample)), colnames(df_without_nonmatching))
df_without_nonmatching <- df_without_nonmatching[re_order]
assert_that(identical(colnames(df_without_nonmatching),make.names(as.character(info_all3$sample))), msg=paste("assertion 1 on iteration ",i))

genedata <- df_without_nonmatching

assert_that(length(setdiff(make.names(as.character(info_all3$sample)), colnames(genedata)))==0, msg=paste("assertion 2 on iteration ",i))
assert_that(length(setdiff(colnames(genedata),make.names(as.character(info_all3$sample))))==0, msg=paste("assertion 3 on iteration ",i))

info_all = info_all3
info_all$id = info_all$original_sample
rownames(info_all)=make.names(as.character(info_all$sample))
assert_that(identical(colnames(genedata),rownames(info_all)), msg=paste("assertion 4 on iteration ",i)) #checks

setwd(paste(main_path,all_var[i,3],sep = ""))

#2. VOOM AND LIMMA ANALYSES PREPARATION AND APPLICATION TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
formVoom = formula(paste0("~ Batch + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + msex + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + race")) 
forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

info_all2=info_all
design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formVoom)[2],fixed=TRUE))),info_all2)
info_all2=info_all2[rownames(design),]
genedata=genedata[,rownames(design)]
assert_that(identical(colnames(genedata),make.names(as.character(info_all2$sample)))) 

info_all2$sample <- as.factor(info_all2$sample)

#Refining and Categorizing Covariate Data Based on Disease Status
#Categorize disease status using CERAD score, which is commonly used to assess Alzheimer's disease pathology.disease_statuses=c("ceradsc")

disease_statuses=c("ceradsc")

for (status in disease_statuses){
if(status!="ceradsc_defvsctl"){
	info_alltmp=info_all2[which(is.na(info_all2[,status])==F),]
}else{
	info_alltmp=info_all2[which(is.na(info_all2[,"ceradsc"])==F),]
}
	for(i in 1:ncol(info_alltmp)){
	if(class(info_alltmp[,i])=="factor"){
	    info_alltmp[,i]=factor(info_alltmp[,i])}
	}
}

#Define statistical models for analysis
#formula_test=formula(paste("~ 0 + ceradsc + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

#for VPA use this: is ceradsc fixed or random effect? fixed!
# formula_test=formula(paste("~ (1|ceradsc) + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

##for regular DEA, use this: 
formula_test=formula(paste("~ 0 + ceradsc + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

formula_test0=formula(paste("~ 0 + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

#Applying Voom Transformation 
isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
dge <- DGEList(counts=genedata[isexpr,]) 
dge <- calcNormFactors(dge)

#pdf("/sc/arion/projects/mscic1/results/jolie/my_voom_plot.pdf")
v <- voom(dge, design, plot=TRUE)
#dev.off()

#save a Rdata so dont have to rerun every time 
#save.image(file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/ROSMAP_everything_until_v.RData")
#load("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/ROSMAP_everything_until_v.RData") #remember still have to reload libraries

#Applying limma for Differential Expression Analysis
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(ceradsc1-ceradsc4, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results<-topTable(fit2,coef=1,number=Inf)

#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
results$X<-rownames(results)
rosmap_no_sva<-results[, c(7,1,2,3,4,5,6)]
row.names(rosmap_no_sva) <- 1:nrow(rosmap_no_sva)

write.table(rosmap_no_sva, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_11.18.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#4. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "BE" AND "LEEK" METHODS
library(sva)

# Filter lowly expressed genes and calculate normalization factors
##for regular DEA, use this: 
formula_test=formula(paste("~ 0 + ceradsc + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))
formula_test0=formula(paste("~ 0 + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
design0=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test0)[2],fixed=TRUE))),info_alltmp)

isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
dge <- DGEList(counts=genedata[isexpr,]) 
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

n.sv = num.sv(v$E,design,method="be") #Change method to "leek" when appropriate
if(n.sv>0){
	svobj = sva(v$E,design,design0,n.sv=n.sv) 
	temp<-svobj$sv
	colnames(temp)<-paste0('sv',1:ncol(temp))
	design_sv = cbind(design, temp) 
	} else {
		design_sv=design 
	}

#5. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv)
contrast.matrix <- makeContrasts(ceradsc1-ceradsc4, levels=design_sv) #ceradsc4 = no symptoms of AD (doesn't mean they don't have AD); ceradsc1 = severe AD
#with ceradsc1-ceradsc4 contrast, if logfold change = 2, then the gene is 4x lower in AD than HC - I CHANGED THE CONTRAST SO IT MATCHES MSBB!!

fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)
results_sv<-topTable(fit2_sv,coef=1,number=Inf)

#6. SUMMARIZE AND SAVE DEA WITH SVA RESULTS
results_sv$X<-rownames(results_sv)
results_sv <- results_sv[, c(7,1,2,3,4,5,6)]
row.names(results_sv) <- 1:nrow(results_sv)

## TODO: JOLIE TO ADD ROSMAP SVA PATH, PLACEHOLDER BELOW
## TODO: JOLIE LOOK FOR LEEK ROSMAP
# write.table(rosmap_no_sva, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_11.18.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


