#1. LOAD AND FORMAT GENE EXPRESSION DATA
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)
library(readr)
library(assertthat)
library(sva)

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
formula_test=formula(paste("~ (1|ceradsc) + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))

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

vp.ROSMAP <- fitExtractVarPartModel(v$E, formula_test, info_alltmp)

percent_status_mean <- mean(vp.ROSMAP$ceradsc) * 100 #0.2282867%

write.table(vp.ROSMAP, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/rosmap_variance_partition.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


################################################# START HERE #######################################################
#save a Rdata so dont have to rerun every time 
save.image(file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/ROSMAP_everything_until_v.RData")
load("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/ROSMAP_everything_until_v.RData") #remember still have to reload libraries

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

# #3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
# write.table(rosmap_no_sva, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_11.18.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# library(qvalue)
# p1 = 1 - pi0est(results$P.Value)$pi0

# bp_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp_11302024.txt",data.table=FALSE)# 19017 genes
# bp_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_no_sva_bp_11302024.txt",data.table=FALSE)#

# bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp.txt",data.table=FALSE)# 19017 genes
# bp_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/mssm_sva_be_bp.txt",data.table=FALSE)#
# p1 = 1 - pi0est(bp_hbcc_no_sva$P.Value)$pi0
# p1
# table(bp_hbcc_no_sva$adj.P.Val<0.05)
# p1 = 1 - pi0est(bp_mssm_no_sva$P.Value)$pi0
# p1
# table(bp_mssm_no_sva$adj.P.Val<0.05)
# p1 = 1 - pi0est(bp_hbcc_sva_be$P.Value)$pi0
# p1
# table(bp_hbcc_sva_be$adj.P.Val<0.05)
# p1 = 1 - pi0est(bp_mssm_sva_be$P.Value)$pi0
# p1
# table(bp_mssm_sva_be$adj.P.Val<0.05)


# #######might have some missing results line here!!!

# write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# rosmap_OG <-read_tsv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt")

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

######## testing B=50 not 5
n.sv = num.sv(v$E,design,method="be") #Change method to "leek" when appropriate
if(n.sv>0){
	svobj = sva(v$E,design,design0,n.sv=n.sv,B=50) 
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
results_sv$X<-rownames(results_sv)
results_sv <- results_sv[, c(7,1,2,3,4,5,6)]
row.names(results_sv) <- 1:nrow(results_sv)

sv_matrix_b50<-design_sv[,22:69]

##b=5
######## testing B=50 not 5
n.sv = num.sv(v$E,design,method="be") #Change method to "leek" when appropriate
if(n.sv>0){
	svobj = sva(v$E,design,design0,n.sv=n.sv,B=5) 
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
results_sv2<-topTable(fit2_sv,coef=1,number=Inf)
results_sv2$X<-rownames(results_sv2)
results_sv2 <- results_sv2[, c(7,1,2,3,4,5,6)]
row.names(results_sv2) <- 1:nrow(results_sv2)

sv_matrix_b5<-design_sv[,22:69]

write.csv(sv_matrix_b5, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_matrix_with_B=5.csv", row.names = FALSE)

write.csv(sv_matrix_b50, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_matrix_with_B=50.csv", row.names = FALSE)

sva_B50<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_B=50.csv")

sv_df_b50 <- as.data.frame(sv_matrix_b50)
sv_df_b5 <- as.data.frame(sv_matrix_b5)

# Compute Pearson correlation for each SV
correlations <- sapply(colnames(sv_df_b50), function(sv) {
  cor(sv_df_b50[[sv]], sv_df_b5[[sv]], method = "pearson")
})

# Convert to a dataframe for easy viewing
correlations_df <- data.frame(SV = names(correlations), Correlation = correlations)
correlations_df

# Compute correlation matrix
cor_matrix <- cor(sv_matrix_b50, sv_matrix_b5, method = "pearson")
cor_matrix
# Print correlation matrix
print(cor_matrix)


##CURRENTLY: 
results_sv

#write.table(results_sv, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_B=50.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

write.csv(results_sv, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_B=50.csv", row.names = FALSE)


test<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp12.txt", data.table = FALSE)

sva_B50<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_B=50.csv")

sva_B5<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_48_SVs.csv")

#both dim are  19476 x 7

# Merge datasets by gene identifier (X)
merged_sva <- merge(sva_B5[, c("X", "logFC")], 
                    sva_B50[, c("X", "logFC")], 
                    by = "X", 
                    suffixes = c("_B5", "_B50"))

# Compute correlation between logFC values
logFC_correlation <- cor(merged_sva$logFC_B5, merged_sva$logFC_B50, method = "pearson")
logFC_correlation #0.9415713


#n.sv = num.sv(v$E,design,method="be") #48 SVs for "be" 0 SV for "leek"
##run n.sv = 11 and 12 to see the difference in results 
#either rerun this or i already generated results with +1 SV from [0:max SV] for ROSMAP and MSBB
#n.sv = 11
n.sv = 13

# ##this one is for kappa 
# if (n.sv > 0) {
#   svobj <- sva(v$E, design, design0, n.sv = n.sv)
#   temp <- svobj$sv
#   colnames(temp) <- paste0("sv", 1:ncol(temp))
  
#   # Initialize the design matrix and store kappa values
#   design_sv <- design
#   kappa_values <- numeric(n.sv)  # To store kappa for each step
  
#   for (i in 1:n.sv) {
#     # Incrementally add each SV to the design matrix
#     design_sv <- cbind(design_sv, temp[, i, drop = FALSE])
    
#     # Calculate kappa for the current design matrix
#     kappa_values[i] <- kappa(design_sv)
#     cat("Kappa statistic after adding SV", i, ":", kappa_values[i], "\n")
    
#     # Optional: Check for singularity at each step
#     if (kappa_values[i] > 1e7) {
#       cat("Warning: Design matrix becomes singular after adding SV", i, "\n")
#     }
#   }
  
#   # Optional: Store kappa values for further analysis or plotting
#   kappa_results <- data.frame(SV = 1:n.sv, Kappa = kappa_values)
#   print(kappa_results)
# } else {
#   design_sv <- design
# }

# Number of significant surrogate variables is:  11 
# Iteration (out of 5 ):1  2  3  4  5  Kappa statistic after adding SV 1 : 2618.019 
# Kappa statistic after adding SV 2 : 2985.542 
# Kappa statistic after adding SV 3 : 4272.128 
# Kappa statistic after adding SV 4 : 4527.639 
# Kappa statistic after adding SV 5 : 4795.038 
# Kappa statistic after adding SV 6 : 5186.369 
# Kappa statistic after adding SV 7 : 5702.244 
# Kappa statistic after adding SV 8 : 5690.722 
# Kappa statistic after adding SV 9 : 5585.465 
# Kappa statistic after adding SV 10 : 5427.84 
# Kappa statistic after adding SV 11 : 5311.215 
#    SV    Kappa
# 1   1 2618.019
# 2   2 2985.542
# 3   3 4272.128
# 4   4 4527.639
# 5   5 4795.038
# 6   6 5186.369
# 7   7 5702.244
# 8   8 5690.722
# 9   9 5585.465
# 10 10 5427.840
# 11 11 5311.215

# Number of significant surrogate variables is:  12 
# Iteration (out of 5 ):1  2  3  4  5  Kappa statistic after adding SV 1 : 2374.56 
# Kappa statistic after adding SV 2 : 3002.08 
# Kappa statistic after adding SV 3 : 4077.648 
# Kappa statistic after adding SV 4 : 4438.499 
# Kappa statistic after adding SV 5 : 4954.099 
# Kappa statistic after adding SV 6 : 5317.357 
# Kappa statistic after adding SV 7 : 5704.109 
# Kappa statistic after adding SV 8 : 5687.128 
# Kappa statistic after adding SV 9 : 5575.457 
# Kappa statistic after adding SV 10 : 5481.977 
# Kappa statistic after adding SV 11 : 5259.471 
# Kappa statistic after adding SV 12 : 4969.515 
#    SV    Kappa
# 1   1 2374.560
# 2   2 3002.080
# 3   3 4077.648
# 4   4 4438.499
# 5   5 4954.099
# 6   6 5317.357
# 7   7 5704.109
# 8   8 5687.128
# 9   9 5575.457
# 10 10 5481.977
# 11 11 5259.471
# 12 12 4969.


# plot <- ggplot(kappa_results, aes(x = SV, y = Kappa)) +
#   geom_line() +
#   geom_point() +
#   xlab("Number of SVs Added") +
#   ylab("Kappa Statistic") +
#   ggtitle("Kappa Statistic vs. Number of Surrogate Variables") +
#   scale_x_continuous(breaks = seq(min(kappa_results$SV), max(kappa_results$SV), by = 1)) + # Add x-axis ticks
#   theme_minimal()

# ggsave("/hpc/users/hoangd02/www/plots/kappa_by_nsv.pdf", plot = plot, width = 12, height = 7)

#one solution: ml R/4.3.0
install.package("nloptr") #need this to run VPA
library(nloptr)

n.sv = 13
##regular sva
if(n.sv>0){
svobj = sva(v$E,design,design0,n.sv=n.sv)
temp<-svobj$sv
colnames(temp)<-paste0('sv',1:ncol(temp))
design_sv = cbind(design, temp) 
} else{
	design_sv=design
}

##temp has all the SVs matrices and design_sv has all the covariates + all SVs 

##vpa with sva
#library(variancePartition)

#Prepare formula for variance partitioning (include SVs)
sv_terms <- paste(paste0("sv", 1:n.sv), collapse = " + ")
# updated_formula <- as.formula(paste(
#   "~ (1|ceradsc) + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race) +",
#   sv_terms
# ))

updated_formula_20241205 <- as.formula(paste(
  "~ ceradsc + Batch + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + msex + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + race+ ",
  sv_terms
))

# Step 6: Combine metadata with SVs
colnames(svobj$sv) <- paste0("sv", 1:ncol(svobj$sv))
metadata <- cbind(info_alltmp, svobj$sv)

# Step 7: Perform variance partitioning analysis
vp9 <- fitExtractVarPartModel(v$E, updated_formula_20241205, metadata)
vp10 <- fitExtractVarPartModel(v$E, updated_formula_20241205, metadata)
vp11 <- fitExtractVarPartModel(v$E, updated_formula_20241205, metadata)
vp12 <- fitExtractVarPartModel(v$E, updated_formula_20241205, metadata)

colnames(vp9)[colnames(vp9) == "ceradsc"] <- "ceradsc_vp9"
colnames(vp10)[colnames(vp10) == "ceradsc"] <- "ceradsc_vp10"
colnames(vp11)[colnames(vp11) == "ceradsc"] <- "ceradsc_vp11"
colnames(vp12)[colnames(vp12) == "ceradsc"] <- "ceradsc_vp12"

vp9_new<-vp9["ceradsc_vp9"]
vp10_new<-vp10["ceradsc_vp10"]
vp11_new<-vp11["ceradsc_vp11"]
vp12_new<-vp12["ceradsc_vp12"]

# Combine the datasets by column, ensuring row names match
combo <- cbind(vp9_new[rownames(vp10_new), , drop = FALSE], 
               vp10_new[rownames(vp9_new), , drop = FALSE],
               vp11_new[rownames(vp11_new), , drop = FALSE], 
               vp12_new[rownames(vp12_new), , drop = FALSE])

write.table(combo, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp12.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

##running vpa code above for 13
n.sv = 13

# Step 1: Run SVA if n.sv > 0
if (n.sv > 0) {
  svobj = sva(v$E, design, design0, n.sv = n.sv)
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  design_sv = cbind(design, temp) 
} else {
  design_sv = design
}

# Step 2: Prepare formula for variance partitioning (include SVs)
sv_terms <- paste(paste0("sv", 1:n.sv), collapse = " + ")
updated_formula_20241205 <- as.formula(paste(
  "~ ceradsc + Batch + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + msex + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + race + ",
  sv_terms
))

# Step 3: Combine metadata with SVs
colnames(svobj$sv) <- paste0("sv", 1:ncol(svobj$sv))
metadata <- cbind(info_alltmp, svobj$sv)

# Step 4: Perform variance partitioning analysis for n.sv = 13
vp13 <- fitExtractVarPartModel(v$E, updated_formula_20241205, metadata)
colnames(vp13)[colnames(vp13) == "ceradsc"] <- "ceradsc_vp13"

# Extract the ceradsc_vp13 column
vp13_new <- vp13["ceradsc_vp13"]

# Ensure vp13_new has a column for rownames
vp13_new <- data.frame(ID = rownames(vp13_new), vp13_new, row.names = NULL)

# Merge combo with vp13_new using the identifier column
combo <- merge(combo, vp13_new, by.x = "V1", by.y = "ID", all.x = TRUE)


write.table(combo, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp13.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


test<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp12.txt", data.table = FALSE)



test<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp12.csv")

combo<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/combo_vp9_through_vp12.csv")

library(dgof)
#Kolmogorov-Smirnov Test
ks_result <- ks.test(test$ceradsc_vp12, test$ceradsc_vp11, alternative = "greater")
ks_result
# Warning message:
# In ks.test(test$ceradsc_vp12, test$ceradsc_vp11, alternative = "greater") :
#   cannot compute correct p-values with ties

# 	Two-sample Kolmogorov-Smirnov test

# data:  test$ceradsc_vp12 and test$ceradsc_vp11
# D^+ = 0.059406, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies above that of y

#values in ceradsc_vp12 tend to be systematically higher than those in ceradsc_vp11.


tmp<-SIGN.test(test$ceradsc_vp12, test$ceradsc_vp11)
# 	Dependent-samples Sign-Test

# data:  test$ceradsc_vp12 and test$ceradsc_vp11
# S = 7693, p-value < 2.2e-16
# alternative hypothesis: true median difference is not equal to 0
# 95 percent confidence interval:
#  -4.667820e-05 -3.895539e-05
# sample estimates:
# median of x-y 
# -4.275662e-05 

#plot only FDR < 0.05 and non-significant = NA
##do this for ROSMAP -- 0 to 48 SVs 
tmp$estimate
median of x-y 
-4.275662e-05

SIGN.test(test$ceradsc_vp11, test$ceradsc_vp12)


combo=combo[order(sub,decreasing=TRUE),]
combo_subset=combo[1:500,]
melted_subset=reshape2::melt(as.matrix(combo_subset))

plot<-ggplot(melted_subset,aes(x=Var2,y=value))+ geom_violin() + geom_boxplot()
ggsave("/hpc/users/hoangd02/www/plots/vpa_sv9_to_sv12.pdf", plot = plot)

tmp=apply(combo,2,function(x){sort(x,decreasing=TRUE)[1:500]})
melted_subset=reshape2::melt(as.matrix(tmp))
plot<-ggplot(melted_subset,aes(x=Var2,y=value))+ geom_violin() + geom_boxplot()
ggsave("/hpc/users/hoangd02/www/plots/vpa_sv9_to_sv12_test.pdf", plot = plot)

##final plot as of 01/28/2025
## ceradsc_vp11 and ceradsc_vp12 - drop
##scaling
x_limits <- range(c(combo$ceradsc_vp9, combo$ceradsc_vp10, combo$ceradsc_vp11, combo$ceradsc_vp12, combo$ceradsc_vp13), na.rm = TRUE)
y_limits <- x_limits # Since it's a diagonal comparison, x and y limits should match.

ceradsc_vp11_and_vp12 <- ggplot(combo, aes(x = ceradsc_vp11, y = ceradsc_vp12)) +
  geom_point() +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Comparison of ceradsc_vp11 and ceradsc_vp12",
       x = "ceradsc_vp11",
       y = "ceradsc_vp12") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/ceradsc_vp11_and_ceradsc_vp12.pdf", plot = plot)
#https://hoangd02.u.hpc.mssm.edu/plots/ceradsc_vp11_and_ceradsc_vp12.pdf

##ceradsc_vp10 and ceradsc_vp11
ceradsc_vp10_and_vp11 <- ggplot(combo, aes(x = ceradsc_vp10, y = ceradsc_vp11)) +
  geom_point() +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Comparison of ceradsc_vp10 and ceradsc_vp11",
       x = "ceradsc_vp10",
       y = "ceradsc_vp11") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/ceradsc_vp11_and_ceradsc_vp10.pdf", plot = plot)
#https://hoangd02.u.hpc.mssm.edu/plots/ceradsc_vp11_and_ceradsc_vp10.pdf

##ceradsc_vp9 and ceradsc_vp10
ceradsc_vp9_and_vp10 <- ggplot(combo, aes(x = ceradsc_vp9, y = ceradsc_vp10)) +
  geom_point() +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Comparison of ceradsc_vp9 and ceradsc_vp10",
       x = "ceradsc_vp9",
       y = "ceradsc_vp10") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/ceradsc_vp9_and_ceradsc_vp10.pdf", plot = plot)
#https://hoangd02.u.hpc.mssm.edu/plots/ceradsc_vp9_and_ceradsc_vp10.pdf

#### RUN 12 and 13
ceradsc_vp12_and_vp13 <- ggplot(combo, aes(x = ceradsc_vp12, y = ceradsc_vp13)) +
  geom_point() +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Comparison of ceradsc_vp12 and ceradsc_vp13",
       x = "ceradsc_vp12",
       y = "ceradsc_vp13") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/ceradsc_vp12_and_ceradsc_vp13.pdf", plot = plot)
#https://hoangd02.u.hpc.mssm.edu/plots/ceradsc_vp9_and_ceradsc_vp10.pdf

library(patchwork)

# Combine plots into a grid
combined_plot <- (ceradsc_vp9_and_vp10 | ceradsc_vp10_and_vp11) /
                 (ceradsc_vp11_and_vp12 | ceradsc_vp12_and_vp13)

# Adjust spacing if necessary
#combined_plot <- combined_plot + plot_layout(guides = "collect")
ggsave("/hpc/users/hoangd02/www/plots/combo_plots_of_ceradsc_vp9_through_ceradsc_vp13.pdf", plot = combined_plot, width = 10, height = 10)




cor(combo$ceradsc_vp11,combo$ceradsc_vp12)
cor(combo$ceradsc_vp11,combo$ceradsc_vp12,method="spearman")

wilcox.test(combo$ceradsc_vp10,combo$ceradsc_vp11,paired=TRUE)

library(BSDA)
SIGN.test(
  x,
  y = NULL,
  md = 0,
  alternative = "two.sided",
  conf.level = 0.95,
  ...
)


mean_contributions <- colMeans(vp); mean_contributions 
# with sv=9: ceradsc contribution to overall vp: 0.001340641 #fixed 
# with sv=10: ceradsc contribution to overall vp: 0.001233090 ##random
# with sv=11: ceradsc contribution to overall vp: 0.001376581 #fixed 
# with sv=12: ceradsc contribution to overall vp: 0.0008602750 #fixed 

head(vp)
#extract  ceradsc column in VP for SVs 9 through 12 and plot 

##cutoff after merge ? 

plot <- plotVarPart(combo) + scale_y_continuous(limits = c(0, max(combo)))

ggsave("/hpc/users/hoangd02/www/plots/vpa_sv9_to_sv12.pdf", plot = plot)




# Step 8: Plot variance partitioning results
plot <- plotVarPart(vp)
ggsave("/hpc/users/hoangd02/www/plots/rosmap_vpa_with_sv10.pdf", plot = plot, width = 40, height = 7)

##to get the numbers
variance_summary <- colMeans(vp); variance_summary
# sv10 : 0.0149763608 
# sv11 : 0.0111720498 
# sv12 : 0.0092526232 
# sv13 : 0.0128019366 
# sv14 : 0.0065898266 

##getting SVs matrices
n_sv_values <- c(13, 11)
# Initialize a list to store the results
sv_matrices <- list()

for (n.sv in n_sv_values) {
  if (n.sv > 0) {
    # Run SVA
    svobj <- sva(v$E, design, design0, n.sv = n.sv)
    
    # Extract SV matrix and set column names
    temp <- svobj$sv
    colnames(temp) <- paste0('sv', 1:ncol(temp))
    
    # Store SV matrix in the list
    sv_matrices[[paste0("n_sv_", n.sv)]] <- temp
  } else {
    sv_matrices[[paste0("n_sv_", n.sv)]] <- design  # No SVs
  }
}

# Access the SV matrices
sv_matrix_13 <- sv_matrices[["n_sv_13"]]  # For n.sv = 11
sv_matrix_12 <- sv_matrices[["n_sv_12"]]
sv_matrix_40 <- sv_matrices[["n_sv_40"]]

cor_matrix <- cor(sv_matrix_13)

write.table(cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/cor_matrix_sv11_sv40.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

correlations <- sapply(1:min(ncol(sv_matrix_11), ncol(sv_matrix_40)), function(i) {
  cor(sv_matrix_11[, i], sv_matrix_40[, i])})
 # [1]  0.9989589  0.9981628  0.9912898  0.9892864  0.9546978  0.9269997
 # [7]  0.9613571  0.9920198  0.8982501 -0.5972662  0.1646049

# Create a dataframe for easier visualization
cor_df <- data.frame(SV = paste0("sv", 1:length(correlations)),
                     Correlation = correlations)

#5. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv)
contrast.matrix <- makeContrasts(ceradsc1-ceradsc4, levels=design_sv) #ceradsc4 = no symptoms of AD (doesn't mean they don't have AD); ceradsc1 = severe AD
#with ceradsc1-ceradsc4 contrast, if logfold change = 2, then the gene is 4x lower in AD than HC - I CHANGED THE CONTRAST SO IT MATCHES MSBB!!

fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)
results_sv<-topTable(fit2_sv,coef=1,number=Inf)
results_sv$X<-rownames(results_sv)
results_sv <- results_sv[, c(7,1,2,3,4,5,6)]
row.names(results_sv) <- 1:nrow(results_sv)

write.table(results_sv, file = "rosmap_sva_be_with_11_SV.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(results_sv, file = "rosmap_sva_be_with_12_SV.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

##running SVA with 12 SVs and remove each SV10, SV11, SV12 individually or a combination of them and plot logFC and p-values until the plots look like SV10 vs SV11 

## doing it manually
##regular sva
n.sv = 40
if(n.sv>0){
svobj = sva(v$E,design,design0,n.sv=n.sv)
temp<-svobj$sv
colnames(temp)<-paste0('sv',1:ncol(temp))
design_sv = cbind(design, temp) 
} else{
	design_sv=design
}
##temp has all the SVs matrices and design_sv has all the covariates + all SVs 

fit_sv <- lmFit(v, design_sv)
contrast.matrix <- makeContrasts(ceradsc1-ceradsc4, levels=design_sv) 
#ceradsc4 = no symptoms of AD (doesn't mean they don't have AD); ceradsc1 = severe AD
#with ceradsc1-ceradsc4 contrast, if logfold change = 2, then the gene is 4x lower in AD than HC 

fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)
results_sv<-topTable(fit2_sv,coef=1,number=Inf)
results_sv$X<-rownames(results_sv)
results_sv <- results_sv[, c(7,1,2,3,4,5,6)]
row.names(results_sv) <- 1:nrow(results_sv)

write.table(results_sv, file = "rosmap_sva_be_with_11_SV.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

##so i have SV=12 results 
sv12<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_12_SV.txt", data.table = FALSE)

##and i have SV=11 results (removed sv12)
sv11<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_11_SV.txt", data.table = FALSE)

##delete SV11 from the 12 SVs matrix 
design_sv <- design_sv[, !(colnames(design_sv) == "sv11")]

##delete SV10 from the 12 SVs matrix 
design_sv <- design_sv[, !(colnames(design_sv) == "sv10")]

##delete SV10 and SV11 and keep SV12 from the 12 SVs matrix 
design_sv <- design_sv[, !(colnames(design_sv) %in% c("sv10", "sv11"))]

##chatgpt code to not do it manually
# Required Libraries
library(limma)

# All columns to test exclusion
columns_to_test <- c("sv10", "sv11", "sv12")

# Generate all combinations of exclusions
exclusion_combinations <- unlist(lapply(1:length(columns_to_test), 
                                        function(x) combn(columns_to_test, x, simplify = FALSE)), 
                                 recursive = FALSE)

# Add cases for excluding none
exclusion_combinations <- c(list(NULL), exclusion_combinations)
# Directory to save results
output_dir <- "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct" # Specify your desired directory
dir.create(output_dir, showWarnings = FALSE)

# Loop through each combination
for (i in seq_along(exclusion_combinations)) {
  
  # Exclude specific columns
  columns_to_exclude <- exclusion_combinations[[i]]
  design_subset <- design_sv[, !(colnames(design_sv) %in% columns_to_exclude)]
  
  # Fit the linear model
  fit_sv <- lmFit(v, design_subset)
  
  # Create contrast matrix
  contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = design_subset)
  
  # Fit contrasts
  fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
  fit2_sv <- eBayes(fit2_sv)
  
  # Extract results
  results_sv <- topTable(fit2_sv, coef = 1, number = Inf)
  results_sv$X <- rownames(results_sv)
  results_sv <- results_sv[, c(7, 1, 2, 3, 4, 5, 6)]
  rownames(results_sv) <- 1:nrow(results_sv)
  
  # Create unique filename
  if (is.null(columns_to_exclude)) {
    exclusion_label <- "No_Exclusions"
  } else {
    exclusion_label <- paste(columns_to_exclude, collapse = "_and_")
  }
  
  output_file <- file.path(output_dir, paste0("DEA_results_Excluding_", exclusion_label, ".csv"))
  
  # Save results to file
  write.csv(results_sv, file = output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n") # Logging progress
}

# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_No_Exclusions.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv11.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv12.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv11.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv12.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv11_and_sv12.csv 
# Results saved to: /sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv11_and_sv12.csv 


sv11<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_11_SVs.csv")

sv11<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_11_SVs.csv")

sv40<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_40_SVs.csv")

merged <- inner_join(sv11, sv40, by = "X", suffix = c("_sv11", "_sv40"))
cor(merged$logFC_sv11, merged$logFC_sv40)


#sv27<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_27_SV.txt", data.table = FALSE)



sv12_all<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_No_Exclusions.csv")
sv12_no_sv10<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10.csv")
sv12_no_sv11<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv11.csv")
sv12_no_sv10_sv11<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv11.csv")

sv12_no_sv10_sv12<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv12.csv")


sv12_no_sv12<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv12.csv")
sv12_no_sv11_sv12<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv11_and_sv12.csv")
sv12_no_sv10_sv11_sv12<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv10_and_sv11_and_sv12.csv")

###
#sv10 = control 
#plotting sv12_all, sv12_no_sv10, sv12_no_sv11, sv12_no_sv10_sv11 against sv10

# Define a function to generate and save plots
plot_sv10_comparison <- function(sv10, sv_other, sv_other_label, output_dir) {
  # Merge data by gene ID (X)
  merged <- inner_join(sv10, sv_other, by = "X", suffix = c("_sv10", paste0("_", sv_other_label)))
  
  # Calculate -log10 P-Values
  merged$logP_sv10 <- -log10(merged$P.Value_sv10)
  merged[[paste0("logP_", sv_other_label)]] <- -log10(merged[[paste0("P.Value_", sv_other_label)]])
  
  # Calculate thresholds for adjusted p-values
  threshold_sv10 <- max(merged$P.Value_sv10[merged$adj.P.Val_sv10 <= 0.05], na.rm = TRUE)
  threshold_other <- max(merged[[paste0("P.Value_", sv_other_label)]][merged[[paste0("adj.P.Val_", sv_other_label)]] <= 0.05], na.rm = TRUE)
  
  # Transform thresholds to -log10 scale
  log_threshold_sv10 <- -log10(threshold_sv10)
  log_threshold_other <- -log10(threshold_other)
  
  # Create plot
  plot <- ggplot(merged, aes(x = logP_sv10, y = merged[[paste0("logP_", sv_other_label)]])) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_density2d(alpha = 0.7, color = "yellow") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = log_threshold_sv10, linetype = "dotted", color = "black") +
    geom_hline(yintercept = log_threshold_other, linetype = "dotted", color = "black") +
    scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    labs(
      title = paste("Comparison of Unadjusted P-Values: SV10 vs", sv_other_label),
      x = "-log10(Unadjusted P-Value) in SV10",
      y = paste0("-log10(Unadjusted P-Value) in ", sv_other_label)
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal()
  
  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0("sv10_vs_", sv_other_label, "_unadjusted_pvalues.pdf")),
    plot = plot
  )
  
  cat("Plot saved for comparison: SV10 vs", sv_other_label, "\n")
}

# Define datasets and their labels
datasets <- list(
  sv12_all = "sv12_all",
  sv12_no_sv10 = "sv12_no_sv10",
  sv12_no_sv11 = "sv12_no_sv11",
  sv12_no_sv10_sv11 = "sv12_no_sv10_sv11"
)

# Specify output directory
output_dir <- "/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/dropped_DEGs/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through datasets and generate plots
for (name in names(datasets)) {
  plot_sv10_comparison(
    sv10 = sv10,
    sv_other = get(name),  # Dynamically get the dataset
    sv_other_label = datasets[[name]],
    output_dir = output_dir
  )
}



sv11<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_11_SV.txt", data.table = FALSE)

#plotting sv12_all, sv12_no_sv10, sv12_no_sv11, sv12_no_sv10_sv11, sv12_no_sv12, sv12_no_sv11_sv12, sv12_no_sv10_sv11_sv12 against sv11
# Define a function to generate and save plots comparing datasets against sv11
plot_sv11_comparison <- function(sv11, sv_other, sv_other_label, output_dir) {
  # Merge data by gene ID (X)
  merged <- inner_join(sv11, sv_other, by = "X", suffix = c("_sv11", paste0("_", sv_other_label)))
  
  # Calculate -log10 P-Values
  merged$logP_sv11 <- -log10(merged$P.Value_sv11)
  merged[[paste0("logP_", sv_other_label)]] <- -log10(merged[[paste0("P.Value_", sv_other_label)]])
  
  # Calculate thresholds for adjusted p-values
  threshold_sv11 <- max(merged$P.Value_sv11[merged$adj.P.Val_sv11 <= 0.05], na.rm = TRUE)
  threshold_other <- max(merged[[paste0("P.Value_", sv_other_label)]][merged[[paste0("adj.P.Val_", sv_other_label)]] <= 0.05], na.rm = TRUE)
  
  # Transform thresholds to -log10 scale
  log_threshold_sv11 <- -log10(threshold_sv11)
  log_threshold_other <- -log10(threshold_other)
  
  # Create plot
  plot <- ggplot(merged, aes(x = logP_sv11, y = merged[[paste0("logP_", sv_other_label)]])) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_density2d(alpha = 0.7, color = "yellow") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = log_threshold_sv11, linetype = "dotted", color = "black") +
    geom_hline(yintercept = log_threshold_other, linetype = "dotted", color = "black") +
    scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    labs(
      title = paste("Comparison of Unadjusted P-Values: SV11 vs", sv_other_label),
      x = "-log10(Unadjusted P-Value) in SV11",
      y = paste0("-log10(Unadjusted P-Value) in ", sv_other_label)
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal()
  
  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0("sv11_vs_", sv_other_label, "_unadjusted_pvalues.pdf")),
    plot = plot
  )
  
  cat("Plot saved for comparison: SV11 vs", sv_other_label, "\n")
}

# Define datasets and their labels
datasets <- list(
  sv12_all = sv12_all,
  sv12_no_sv10 = sv12_no_sv10,
  sv12_no_sv11 = sv12_no_sv11,
  sv12_no_sv10_sv11 = sv12_no_sv10_sv11,
  sv12_no_sv10_sv12 = sv12_no_sv10_sv12,
  sv12_no_sv12 = sv12_no_sv12,
  sv12_no_sv11_sv12 = sv12_no_sv11_sv12,
  sv12_no_sv10_sv11_sv12 = sv12_no_sv10_sv11_sv12
)

# Specify output directory
output_dir <- "/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/dropped_DEGs_correct/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through datasets and generate plots
for (name in names(datasets)) {
  plot_sv11_comparison(
    sv11 = sv11,
    sv_other = datasets[[name]],  
    sv_other_label = name,
    output_dir = output_dir
  )
}






#plotting sv12_no_sv10, sv12_no_sv11, sv12_no_sv10_sv11 against sv12_all
# Define a function to generate and save plots comparing datasets against sv12_all
plot_sv12_all_comparison <- function(sv12_all, sv_other, sv_other_label, output_dir) {
  # Merge data by gene ID (X)
  merged <- inner_join(sv12_all, sv_other, by = "X", suffix = c("_sv12_all", paste0("_", sv_other_label)))
  
  # Calculate -log10 P-Values
  merged$logP_sv12_all <- -log10(merged$P.Value_sv12_all)
  merged[[paste0("logP_", sv_other_label)]] <- -log10(merged[[paste0("P.Value_", sv_other_label)]])
  
  # Calculate thresholds for adjusted p-values
  threshold_sv12_all <- max(merged$P.Value_sv12_all[merged$adj.P.Val_sv12_all <= 0.05], na.rm = TRUE)
  threshold_other <- max(merged[[paste0("P.Value_", sv_other_label)]][merged[[paste0("adj.P.Val_", sv_other_label)]] <= 0.05], na.rm = TRUE)
  
  # Transform thresholds to -log10 scale
  log_threshold_sv12_all <- -log10(threshold_sv12_all)
  log_threshold_other <- -log10(threshold_other)
  
  # Create plot
  plot <- ggplot(merged, aes(x = logP_sv12_all, y = merged[[paste0("logP_", sv_other_label)]])) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_density2d(alpha = 0.7, color = "yellow") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = log_threshold_sv12_all, linetype = "dotted", color = "black") +
    geom_hline(yintercept = log_threshold_other, linetype = "dotted", color = "black") +
    scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
    labs(
      title = paste("Comparison of Unadjusted P-Values: SV12_all vs", sv_other_label),
      x = "-log10(Unadjusted P-Value) in SV12_all",
      y = paste0("-log10(Unadjusted P-Value) in ", sv_other_label)
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal()
  
  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0("sv12_all_vs_", sv_other_label, "_unadjusted_pvalues.pdf")),
    plot = plot
  )
  
  cat("Plot saved for comparison: SV12_all vs", sv_other_label, "\n")
}

# Define datasets and their labels
datasets <- list(
  sv12_no_sv10 = "sv12_no_sv10",
  sv12_no_sv11 = "sv12_no_sv11",
  sv12_no_sv10_sv11 = "sv12_no_sv10_sv11"
)

# Specify output directory
output_dir <- "/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/dropped_DEGs/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through datasets and generate plots
for (name in names(datasets)) {
  plot_sv12_all_comparison(
    sv12_all = sv12_all,
    sv_other = get(name),  # Dynamically get the dataset
    sv_other_label = datasets[[name]],
    output_dir = output_dir
  )
}

sv12
sv12_all

sv11_og<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_11_SV.txt", data.table = FALSE)

sv11<-read.csv("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/DEA_results_correct/DEA_results_Excluding_sv12.csv")



merged <- inner_join(sv12_all, sv12, by = "X", suffix = c("_sv12_all", "_sv12"))
cor(merged$logFC_sv12_all, merged$logFC_sv12)




merged$logP_sv12_all <- -log10(merged$P.Value_sv12_all)
merged$logP_sv12 <- -log10(merged$P.Value_sv12)

# Calculate thresholds based on adjusted p-value criteria
threshold_sv12_all <- max(merged$P.Value_sv12_all[merged$adj.P.Val_sv12_all <= 0.05]) #0.005729018
threshold_sv12 <- max(merged$P.Value_sv12[merged$adj.P.Val_sv12 <= 0.05]) #0.001265276

# Transform thresholds to log10 scale for plotting
log_threshold_sv12_all <- -log10(threshold_sv12_all)
log_threshold_sv12 <- -log10(threshold_sv12)


plot <- ggplot(merged, aes(x = logP_sv12_all, y = logP_sv12)) +
  geom_point(alpha = 0.5, color = "blue") +  #Add points
  geom_density2d(alpha = 0.7, color = "yellow") +  # Add 2D density lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  #Identity line
  geom_vline(xintercept = threshold_sv12_all, linetype = "dotted", color = "black") +  #Vertical line
  geom_hline(yintercept = log_threshold_sv12, linetype = "dotted", color = "black") +  #Horizontal line
  scale_x_continuous(
    limits = c(0, 14),  # Set x-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 14),  # Set y-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set y-axis tick marks
  ) +
  labs(
    title = "Comparison of Unadjusted P-Values in SV12 and SV12_all",
    x = "-log10(Unadjusted P-Value) in logP_sv12_all",
    y = "-log10(Unadjusted P-Value) in logP_sv12"
  ) +
  coord_fixed(ratio = 1) +  #Ensure equal scaling of axes
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_sv12_all.pdf", plot = plot)




















###EXAMINING THE DIFFERENCE BTW nsv = 11 and 12 

sv10<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_10_SVs.csv", data.table = FALSE)
sv11<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_11_SV.txt", data.table = FALSE)
sv12<-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_with_12_SV.txt", data.table = FALSE)

if (kappa(design) > 1e7) {
    stop("Design matrix is singular, covariates are very correlated")
}

##number of DEGs in each dataset
sum(sv11$adj.P.Val <= 0.05, na.rm = TRUE) #2233 DEGs
sum(sv12$adj.P.Val <= 0.05, na.rm = TRUE) #493 DEGs 

range(sv11$adj.P.val, na.rm = TRUE) #1.890428e-11 9.999713e-01
range(sv12$adj.P.Val, na.rm = TRUE) #9.904912e-14 9.999446e-01

#does SVA impact p-value or logFC?
# Histogram for adj.P.Val
plot <- ggplot(sv12, aes(x = -log10(adj.P.Val))) +
  geom_histogram(binwidth = 0.1, fill = "black", color = "black") +
  scale_x_continuous(
    limits = c(0, 12),  # Set x-axis range
    breaks = seq(0, 12, by = 2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 800),  # Set y-axis range
    breaks = seq(0, 800, by = 100)  # Set y-axis tick marks
  ) +
  labs(
    title = "Distribution of adj.P.Val at SV = 12 (493 DEGs)",
    x = "-log10(Adjusted P-Value)",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_adjPvalue_distribution_histogram.pdf", plot = plot)


#including a vertical line at the maximum unadjusted p-value where the adjusted p-value is ≤ 0.05; this vertical corresponds to the largest unadjusted p-value that still results in an adjusted p-value ≤ 0.05.
#Any unadjusted p-value smaller than this threshold is statistically significant after correction (adjusted p-value ≤ 0.05).
#Conversely, any unadjusted p-value larger than this threshold will have an adjusted p-value > 0.05 and is not considered statistically significant.

# # Filter the data to find the maximum unadjusted p-value with adj.P.Val <= 0.05
# max_unadjusted_p <- sv12 %>%
#   filter(adj.P.Val <= 0.05) %>%
#   summarise(max_p = max(P.Value, na.rm = TRUE)) %>%
#   pull(max_p)

# # Create the plot using regular adjusted p-values
# plot <- ggplot(sv12, aes(x = adj.P.Val)) +
#   geom_histogram(boundary = 0, binwidth = 1/100, fill = "black", color = "black") +  # Ensure bins align with 0 and divide evenly
#   geom_vline(xintercept = max_unadjusted_p, color = "red", linetype = "dashed", size = 1) +  # Line for max unadjusted p-value
#   annotate("text", x = max_unadjusted_p, y = 750, label = paste0("Max Unadjusted P-Value: ", round(max_unadjusted_p, 3)), 
#            color = "red", angle = 90, vjust = -1) +  # Add label for the red line
#   scale_x_continuous(
#     limits = c(0, 1),  # Set x-axis range for p-values
#     breaks = seq(0, 1, by = 0.1)  # Set x-axis tick marks
#   ) +
#   scale_y_continuous(
#     limits = c(0, 800),  # Set y-axis range
#     breaks = seq(0, 800, by = 100)  # Set y-axis tick marks
#   ) +
#   labs(
#     title = "Distribution of Adjusted P-Values at SV = 12 (493 DEGs)",
#     x = "Adjusted P-Value",
#     y = "Frequency"
#   ) +
#   theme_minimal()

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_adjPvalue_distribution_histogram2.pdf", plot = plot,height=12)

range(sv11$P.val, na.rm = TRUE) #1.890428e-11 9.999713e-01
range(sv12$P.Val, na.rm = TRUE) #5.085701e-18 9.999446e-01

# Histogram for P.Val
plot <- ggplot(sv12, aes(x = P.Value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  scale_x_continuous(
    limits = c(0, 1),  # Set x-axis range for raw p-values
    breaks = seq(0, 1, by = 0.1)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 800),  # Set y-axis range
    breaks = seq(0, 800, by = 100)  # Set y-axis tick marks
  ) +
  labs(
    title = "Distribution of P.Value at SV = 12 (493 DEGs)",
    x = "P-Value",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_raw_Pvalue_distribution_histogram.pdf", plot = plot)


# Filter the data to find the maximum unadjusted p-value with adj.P.Val <= 0.05
max_unadjusted_p <- sv11 %>%
  filter(adj.P.Val <= 0.05) %>%
  summarise(max_p = max(P.Value, na.rm = TRUE)) %>%
  pull(max_p)

# Add the vertical line to the plot
plot <- ggplot(sv11, aes(x = P.Value)) +
  geom_histogram(boundary = 0, binwidth = 1/100, fill = "skyblue", color = "black") +  # Align bins at 0 and set binwidth
  geom_vline(xintercept = max_unadjusted_p, color = "red", linetype = "dashed", size = 1) +  # Add vertical line
  scale_x_continuous(
    limits = c(0, 1),  # Set x-axis range for raw p-values
    breaks = seq(0, 1, by = 0.1)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 800),  # Set y-axis range
    breaks = seq(0, 800, by = 100)  # Set y-axis tick marks
  ) +
  labs(
    title = "Distribution of P.Value at SV = 11 (2233 DEGs)",
    x = "P-Value",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_raw_Pvalue_distribution_histogram2.pdf", plot = plot)

##################################################################################
##plot sv11 and sv12 tgt for unadjusted p-values
library(ggplot2)
library(dplyr)

# Combine sv11 and sv12 by gene ID (X)
merged <- inner_join(sv11, sv12, by = "X", suffix = c("_sv11", "_sv12"))
merged10 <- inner_join(sv10, sv11, by = "X", suffix = c("_sv10", "_sv11"))

merged$logP_sv11 <- -log10(merged$P.Value_sv11)
merged$logP_sv12 <- -log10(merged$P.Value_sv12)

merged10$logP_sv10 <- -log10(merged10$P.Value_sv10)
merged10$logP_sv11 <- -log10(merged10$P.Value_sv11)

# Calculate thresholds based on adjusted p-value criteria
threshold_sv11 <- max(merged$P.Value_sv11[merged$adj.P.Val_sv11 <= 0.05]) #0.005729018
threshold_sv12 <- max(merged$P.Value_sv12[merged$adj.P.Val_sv12 <= 0.05]) #0.001265276

threshold_sv10 <- max(merged10$P.Value_sv10[merged10$adj.P.Val_sv10 <= 0.05]) #0.006282186
threshold_sv11 <- max(merged10$P.Value_sv11[merged10$adj.P.Val_sv11 <= 0.05]) # 0.005729018

# Transform thresholds to log10 scale for plotting
log_threshold_sv11 <- -log10(threshold_sv11)
log_threshold_sv12 <- -log10(threshold_sv12)

log_threshold_sv10 <- -log10(threshold_sv10)
log_threshold_sv11 <- -log10(threshold_sv11)

plot <- ggplot(merged10, aes(x = logP_sv10, y = logP_sv11)) +
  geom_point(alpha = 0.5, color = "blue") +  #Add points
  geom_density2d(alpha = 0.7, color = "yellow") +  # Add 2D density lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  #Identity line
  geom_vline(xintercept = log_threshold_sv11, linetype = "dotted", color = "black") +  #Vertical line
  geom_hline(yintercept = log_threshold_sv12, linetype = "dotted", color = "black") +  #Horizontal line
  scale_x_continuous(
    limits = c(0, 14),  # Set x-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 14),  # Set y-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set y-axis tick marks
  ) +
  labs(
    title = "Comparison of Unadjusted P-Values in SV10 and SV11",
    x = "-log10(Unadjusted P-Value) in SV10",
    y = "-log10(Unadjusted P-Value) in SV11"
  ) +
  coord_fixed(ratio = 1) +  #Ensure equal scaling of axes
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv10_and_sv11_unadjusted_pvalues.pdf", plot = plot)

#ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_and_sv12_unadjusted_pvalues3.pdf", plot = plot)


# Transform thresholds to log10 scale for plotting
log_threshold_sv11 <- -log10(threshold_sv11)
log_threshold_sv12 <- -log10(threshold_sv12)
#11 and 12 SVs
merged <- merged %>%
  mutate(
    quadrant = case_when(
      logP_sv11 >= log_threshold_sv11 & logP_sv12 >= log_threshold_sv12 ~ "Q1",  # Significant in both
      logP_sv11 < log_threshold_sv11 & logP_sv12 >= log_threshold_sv12 ~ "Q2",   # Significant in SV12 only
      logP_sv11 < log_threshold_sv11 & logP_sv12 < log_threshold_sv12 ~ "Q3",    # Not significant in either
      logP_sv11 >= log_threshold_sv11 & logP_sv12 < log_threshold_sv12 ~ "Q4"    # Significant in SV11 only
    )
  )

contingency_table <- table(merged$quadrant); contingency_table


log_threshold_sv10 <- -log10(threshold_sv10)
log_threshold_sv11 <- -log10(threshold_sv11)
#10 and 11 SVs
merged <- merged10 %>%
  mutate(
    quadrant = case_when(
      logP_sv10 >= log_threshold_sv10 & logP_sv11 >= log_threshold_sv11 ~ "Q1",  # Significant in both
      logP_sv10 < log_threshold_sv10 & logP_sv11 >= log_threshold_sv11 ~ "Q2",   # Significant in sv11 only
      logP_sv10 < log_threshold_sv10 & logP_sv11 < log_threshold_sv11 ~ "Q3",    # Not significant in either
      logP_sv10 >= log_threshold_sv10 & logP_sv11 < log_threshold_sv11 ~ "Q4"    # Significant in sv10 only
    )
  )
# Create contingency table
contingency_table <- table(merged$quadrant); contingency_table


##################################################################################
#sv11 and sv12 logFC
merged <- inner_join(sv11, sv12, by = "X", suffix = c("_sv11", "_sv12"))
merged10 <- inner_join(sv10, sv11, by = "X", suffix = c("_sv10", "_sv11"))

plot <- ggplot(merged10, aes(x = logFC_sv10, y = logFC_sv11)) +
  geom_point(alpha = 0.5, color = "blue") +  #Add points
  geom_density2d(alpha = 0.7, color = "yellow") +  # Add 2D density lines
  scale_x_continuous(
    limits = c(-0.8, 1.2),  
    breaks = seq(-0.8, 1.2, by = 0.2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(-0.8, 1.2),  
    breaks = seq(-0.8, 1.2, by = 0.2)  # Set x-axis tick marks
  ) +
  labs(
    title = "Comparison of logFC in SV10 and SV11",
    x = "logFC in SV10",
    y = "logFC in SV11"
  ) +
  coord_fixed(ratio = 1) +  #Ensure equal scaling of axes
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv10_and_sv11_logFC.pdf", plot = plot)

plot <- ggplot(merged, aes(x = logFC_sv11, y = logFC_sv12)) +
  geom_point(alpha = 0.5, color = "blue") +  #Add points
  geom_density2d(alpha = 0.7, color = "yellow") +  # Add 2D density lines
  scale_x_continuous(
    limits = c(-0.8, 1.2),  
    breaks = seq(-0.8, 1.2, by = 0.2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(-0.8, 1.2),  
    breaks = seq(-0.8, 1.2, by = 0.2)  # Set x-axis tick marks
  ) +
  labs(
    title = "Comparison of logFC in SV11 and SV12",
    x = "logFC in SV11",
    y = "logFC in SV12"
  ) +
  coord_fixed(ratio = 1) +  #Ensure equal scaling of axes
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_and_sv12_logFC.pdf", plot = plot)



##################################################################################
##plot sv11 and sv12 tgt for adjusted p-values
# Merge the two tables by gene ID (X)
merged <- inner_join(sv11, sv12, by = "X", suffix = c("_sv11", "_sv12"))

# Transform adjusted p-values to -log10 scale
plot_data <- merged %>%
  mutate(
    log_adjP_sv11 = -log10(adj.P.Val_sv11),  # Transform adj.P.Val for sv11
    log_adjP_sv12 = -log10(adj.P.Val_sv12)   # Transform adj.P.Val for sv12
  )

# Set the significance threshold (-log10(0.05))
significance_threshold <- -log10(0.05)

# Create the scatter plot
plot <- ggplot(plot_data, aes(x = log_adjP_sv11, y = log_adjP_sv12)) +
  geom_point(alpha = 0.6, color = "blue") +  # Add points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Identity line
  geom_vline(xintercept = significance_threshold, linetype = "dotted", color = "black") +  # Vertical threshold
  geom_hline(yintercept = significance_threshold, linetype = "dotted", color = "black") +  # Horizontal threshold
  scale_x_continuous(
    limits = c(0, 14),  # Set x-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set x-axis tick marks
  ) +
  scale_y_continuous(
    limits = c(0, 14),  # Set y-axis range for log-adjusted p-values
    breaks = seq(0, 14, by = 2)  # Set y-axis tick marks
  ) +
  labs(
    title = "Comparison of Adjusted P-Values in SV11 and SV12",
    x = "-log10(Adjusted P-Value) in SV11",
    y = "-log10(Adjusted P-Value) in SV12"
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_and_sv12_adjusted_pvalues.pdf", plot = plot)

###interpretation of the lines
#horizontal line is for sv12
# Any points above this line have adjusted p-values in sv12 that are statistically significant (adj.P.Val≤0.05).
# Points below this line have adjusted p-values in sv12 that are not statistically significant (adj.P.Val>0.05).

#vertical line is for sv11
# Any points to the right of this line have adjusted p-values in sv11 that are statistically significant (adj.P.Val≤0.05).
# Points to the left of this line have adjusted p-values in sv11 that are not statistically significant (adj.P.Val>0.05).

##QUADRANTS
# Top Right = Significant in Both
# Top Left = Significant in sv12 Only
# Bottom Right = Significant in sv11 Only 
# Bottom Left = Not Significant in Either

# ##adding the colors
# library(tidyr)
# # Prepare data for plotting and reshape to long format
# plot_data <- merged %>%
#   select(X, P.Value_sv11 = P.Value_sv11, P.Value_sv12 = P.Value_sv12) %>%
#   mutate(
#     logP_sv11 = -log10(P.Value_sv11),  # Transform sv11 p-values to -log10
#     logP_sv12 = -log10(P.Value_sv12)   # Transform sv12 p-values to -log10
#   ) %>%
#   pivot_longer(
#     cols = c(logP_sv11, logP_sv12),  # Reshape the log-transformed columns
#     names_to = "Dataset",
#     values_to = "logP_Value"
#   ) %>%
#   mutate(Dataset = recode(Dataset, "logP_sv11" = "SV11", "logP_sv12" = "SV12"))  # Rename for clarity

# # Set the significance threshold (-log10(0.05))
# significance_threshold <- -log10(0.05)

# # Create the scatter plot with colors
# plot <- ggplot(plot_data, aes(x = X, y = logP_Value, color = Dataset)) +
#   geom_point(alpha = 0.6) +  # Add points colored by dataset
#   geom_hline(yintercept = significance_threshold, linetype = "dotted", color = "black") +  # Horizontal threshold
#   labs(
#     title = "Comparison of Unadjusted P-Values (Log Scale)",
#     x = "Gene ID",
#     y = "-log10(P.Value)",
#     color = "Dataset"
#   ) +
#   scale_color_manual(values = c("SV11" = "blue", "SV12" = "orange")) +  # Custom colors
#   theme_minimal()

# ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_and_sv12_2.pdf", plot = plot)





range(sv11$logFC) #-0.6548974  1.1799402
range(sv12$logFC) #-0.6548974  1.1799402

# Histogram for logFC
plot2 <- ggplot(sv12, aes(x = logFC)) +
  geom_histogram(binwidth = 0.1, fill = "orange", color = "black") +
  labs(
    title = "Distribution of logFC at SV = 12 (493 DEGs)",
    x = "Log Fold Change",
    y = "Frequency"
  ) +
  coord_cartesian(xlim = c(-0.8, 1.2), ylim = c(0, 16000)) +  # Adjust x and y limits
  scale_x_continuous(
    breaks = seq(-0.8, 1.2, by = 0.2)  # Define x-axis tick positions
  ) +
  scale_y_continuous(
    breaks = seq(0, 16000, by = 2000)  # Define y-axis tick positions
  ) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_logFC_distribution_histogram.pdf", plot = plot2)


# Scatter plot of logFC vs adj.P.Val
plot3 <- ggplot(sv11, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(color = "purple", alpha = 0.7) +
  coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(0, 14)) +  # Set axis limits
  scale_x_continuous(
    breaks = seq(-0.6, 1.2, by = 0.4)  # Define x-axis tick marks
  ) +
  scale_y_continuous(
    breaks = seq(0, 14, by = 2)  # Define y-axis tick marks
  ) +
  labs(
    title = "logFC vs -log10(adj.P.Val) at SV = 11 (2233 DEGs)",
    x = "Log Fold Change (logFC)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv11_logFC_and_pvalue_scatterplot.pdf", plot = plot3)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sv12_logFC_and_pvalue_scatterplot.pdf", plot = plot3)

library(variancePartition)

expr_matrix <- as.matrix(sv11$AveExpr)
rownames(expr_matrix) <- sv11$X 
design <- ~ logFC
vp11 <- fitExtractVarPartModel(expr_matrix, design, data = sv11)


vp11 <- fitExtractVarPartModel(sv11$AveExpr ~ logFC, data = sv11)
vp12 <- fitExtractVarPartModel(sv12$AveExpr ~ logFC, data = sv12)

vp11_plot<-plotVarPart(vp11, title = "Variance Partitioning for SV = 11")
vp12_plot<-plotVarPart(vp12, title = "Variance Partitioning for SV = 12")

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/vp11_vpa_plot.pdf", plot = vp11_plot)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/vp12_vpa_plot.pdf", plot = vp12_plot)

save.image('/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_with_vpa.Rdata')



#compared to sv11, at sv12, p-values got a lot smaller 


#rosmap_sva_be_with_100_SV <- results_sv[, c(7,1,2,3,4,5,6)]
#row.names(rosmap_sva_be_with_100_SV) <- 1:nrow(rosmap_sva_be_with_100_SV)
write.table(rosmap_sva_be_with_100_SV, file = "rosmap_sva_be_with_100_SV.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

dim(rosmap_sva_be_with_100_SV)
#19476 x 7

total_significant_rosmap_with_100_SV <- sum(rosmap_sva_be_with_100_SV$adj.P.Val <= 0.05, na.rm = TRUE)
total_significant_rosmap_with_100_SV  #59 DEGs


##### RUN SVA BUT SAY NUMBERS OF SV = 100
library(sva)

# Specify the number of surrogate variables to estimate
n.sv <- 100  # Force the use of 100 SVs

# Estimate and include the surrogate variables (if n.sv > 0)
if (n.sv > 0) {
  svobj <- sva(v$E, design, design0, n.sv = n.sv)  # Use specified number of SVs
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  design_sv <- cbind(design, temp)  # Add SVs to the design matrix
} else {
  design_sv <- design  # Use the original design matrix if no SVs
}

# LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv)

# Define the contrast matrix for differential expression analysis
contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = design_sv) 
# ceradsc4 = no symptoms of AD; ceradsc1 = severe AD

# Fit contrasts and apply empirical Bayes moderation
fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)

# Extract the results and format the output
results_sv <- topTable(fit2_sv, coef = 1, number = Inf)
results_sv$X <- rownames(results_sv)  # Add gene names as a column
rosmap_sva_be_1.16 <- results_sv[, c(7, 1, 2, 3, 4, 5, 6)]  # Reorder columns
row.names(rosmap_sva_be_1.16) <- 1:nrow(rosmap_sva_be_1.16)  # Reset row indices


#############################################OLD CODE##################################################
jolie_sva_be<- results_sv[, c(7,1,2,3,4,5,6)]
row.names(jolie_sva_be) <- 1:nrow(jolie_sva_be)

jolie_sva_leek<- results_sv[, c(7,1,2,3,4,5,6)]
row.names(jolie_sva_leek) <- 1:nrow(jolie_sva_leek)

rosmap_sva_be_1.16 <- jolie_sva_be
rosmap_sva_leek_1.16 <- jolie_sva_leek

write.table(rosmap_sva_be_1.16, file = "rosmap_sva_be_1.16.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(rosmap_sva_leek_1.16, file = "rosmap_sva_leek_1.16.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

save.image('/sc/arion/projects/mscic1/results/jolie/DEA/Rosmapv30/rosmap_sva_1.16.Rdata')

#do corr btw jolie (without SVA) and jolie_sva_be and jolie_sva_leek
combined2<-merge(jolie,jolie_sva_be, by='X', all.x=TRUE, all.y=TRUE) 
cor.test(combined2$logFC_j, combined2$logFC) #this is no sva and leek #r=-0.6088539 

table(combined2$adj.P.Val.x <=0.05)
table(combined2$adj.P.Val.y <=0.05)

################################################END###############################################
#note: 
##run the following but make i = 3 which will make if v30. disease_statuses = ceradsc_defvsctl

###NEW CODE 11-11-2024 ###############################################
###add 1 SV at a time for ROSMAP + MSBB and test their correlation 
library(sva)
library(limma)

# Step 1: Estimate the number of SVs
n.sv = num.sv(v$E, design, method = "be") # 48 SVs for "be"

# Initialize a list to store results for each DEA
all_results_sv <- list()

# Loop through the number of SVs from 1 to the estimated n.sv
for (i in 1:n.sv) {
  cat("Running DEA with", i, "SV(s)\n")
  
  # Step 2: Generate SVs with the current number `i`
  svobj = sva(v$E, design, design0, n.sv = i)
  temp <- svobj$sv
  colnames(temp) <- paste0('sv', 1:ncol(temp))
  
  # Step 3: Create the design matrix with `i` SVs
  design_sv = cbind(design, temp)
  
  # Step 4: Perform DEA using limma
  fit_sv <- lmFit(v, design_sv)
  contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = design_sv)
  fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
  fit2_sv <- eBayes(fit2_sv)
  results_sv <- topTable(fit2_sv, coef = 1, number = Inf)
  
  # Add gene names to the results and reformat
  results_sv$X <- rownames(results_sv)
  results_sv_formatted <- results_sv[, c(7, 1, 2, 3, 4, 5, 6)]
  row.names(results_sv_formatted) <- 1:nrow(results_sv_formatted)
  
  # Step 5: Save the results for this iteration
  all_results_sv[[i]] <- results_sv_formatted
}

# Optional: Combine results into a single data structure or save them to files
# For example, save all results to a list of dataframes or write them as CSV files
for (i in 1:length(all_results_sv)) {
  write.csv(all_results_sv[[i]], file = paste0("rosmap_sva_be_", i, "_SVs.csv"), row.names = FALSE)
}

###TO OPEN THE RESULTS ---- 
###OPTION 1:
# Accessing results from the saved CSV files
result_sv_1 <- read.csv("rosmap_sva_be_1_SVs.csv")
result_sv_2 <- read.csv("rosmap_sva_be_2_SVs.csv")
# Continue for other files, or automate this process

###OPTION 2:
# Get all CSV files in the working directory
file_list <- list.files(pattern = "rosmap_sva_be_\\d+_SVs.csv")

# Read all results into a list
all_results_from_csv <- lapply(file_list, read.csv)

# Access individual results
result_sv_1 <- all_results_from_csv[[1]]  # Results with 1 SV
result_sv_2 <- all_results_from_csv[[2]]  # Results with 2 SVs


#start: 3:28pm, ended 3:51pm
#this is saved in "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30"

###PERFECT - JUST NEED TO DO THIS FOR MSBB!! 

#ROSMAP – randomly shuffle case and control status 100x and run DEA with and without SVA plot distribution of numbers of DEGs color by with and without SVA 
#Null distribution centers at 0 without SVA, with SVA  center > 0 
#Nsv = 0 to 49, set seed 

library(sva)
library(limma)
library(ggplot2)
library(dplyr)

#define function for DEA with SVA for varying n.sv
run_dea_sva_varying_nsv <- function(expression_data, design, design0, n_shuffle, max_nsv) {
  results_list <- list()
  
  for (i in 1:n_shuffle) {
    set.seed(i)  # Set a different seed for each shuffle
    
    # Create group labels for shuffling
    group_labels <- apply(design[, c("ceradsc1", "ceradsc4")], 1, function(row) {
      if (row["ceradsc1"] == 1) {
        return("ceradsc1")
      } else if (row["ceradsc4"] == 1) {
        return("ceradsc4")
      } else {
        return("other")
      }
    })
    
    #shuffle group labels
    shuffled_labels <- sample(group_labels)
    
    #update design matrix based on shuffled labels
    shuffled_design <- design
    shuffled_design[, "ceradsc1"] <- ifelse(shuffled_labels == "ceradsc1", 1, 0)
    shuffled_design[, "ceradsc4"] <- ifelse(shuffled_labels == "ceradsc4", 1, 0)
    
    for (n_sv in 0:max_nsv) {
      if (n_sv > 0) {
        svobj <- sva(expression_data$E, shuffled_design, design0, n.sv = n_sv)
        temp <- svobj$sv
        colnames(temp) <- paste0('sv', 1:ncol(temp))
        design_sv <- cbind(shuffled_design, temp)
      } else {
        design_sv <- shuffled_design
      }
      
      #LIMMA with SVA
      fit_sv <- lmFit(expression_data, design_sv)
      contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = design_sv)
      fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
      fit2_sv <- eBayes(fit2_sv)
      results_sv <- topTable(fit2_sv, coef = 1, number = Inf)
      
      #count DEGs
      deg_count_sv <- sum(results_sv$adj.P.Val < 0.05)
      
      #save results
      results_list[[length(results_list) + 1]] <- data.frame(
        Method = "SVA",
        Shuffle = i,
        Seed = i,  # Explicitly record the seed
        n_sv = n_sv,
        DEGs = deg_count_sv
      )
    }
  }
  
  return(do.call(rbind, results_list))
}

#define function for DEA without SVA
run_dea_no_sva <- function(expression_data, design, n_shuffle) {
  results_list <- list()
  
  for (i in 1:n_shuffle) {
    set.seed(i)  # Set a different seed for each shuffle
    
    #create group labels for shuffling
    group_labels <- apply(design[, c("ceradsc1", "ceradsc4")], 1, function(row) {
      if (row["ceradsc1"] == 1) {
        return("ceradsc1")
      } else if (row["ceradsc4"] == 1) {
        return("ceradsc4")
      } else {
        return("other")
      }
    })
    
    #shuffle group labels
    shuffled_labels <- sample(group_labels)
    
    #update design matrix based on shuffled labels
    shuffled_design <- design
    shuffled_design[, "ceradsc1"] <- ifelse(shuffled_labels == "ceradsc1", 1, 0)
    shuffled_design[, "ceradsc4"] <- ifelse(shuffled_labels == "ceradsc4", 1, 0)
    
    #LIMMA without SVA
    fit <- lmFit(expression_data, shuffled_design)
    contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = shuffled_design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    results <- topTable(fit2, coef = 1, number = Inf)
    
    #count DEGs
    deg_count <- sum(results$adj.P.Val < 0.05)
    
    #save results
    results_list[[length(results_list) + 1]] <- data.frame(
      Method = "No SVA",
      Shuffle = i,
      Seed = i,  # Explicitly record the seed
      n_sv = NA,  # Not applicable
      DEGs = deg_count
    )
  }
  
  return(do.call(rbind, results_list))
}

#run analysis
n_shuffle <- 100
max_nsv <- 49
results_with_sva <- run_dea_sva_varying_nsv(v, design, design0, n_shuffle, max_nsv)
results_without_sva <- run_dea_no_sva(v, design, n_shuffle)

#combine results
combined_results <- rbind(results_with_sva, results_without_sva)

#start: 4:03pm for a screen section rosmap2.0 or sth? killed
#start: 4:15pm for a screen section called rosmap in sklar1 - still need to save results

#start: 4:58pm 11/18 on sklar1 

#save results to CSV
write.csv(combined_results, file = "/sc/arion/projects/mscic1/results/jolie/combined_results_with_seed.csv", row.names = FALSE)

#plot results
ggplot(combined_results, aes(x = Method, y = DEGs, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~n_sv, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of DEGs with and without SVA across varying n_sv",
       x = "Method",
       y = "Number of DEGs",
       fill = "Method")














