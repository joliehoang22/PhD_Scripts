# path=/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/
# cd $path
# ml R

# bsub -q premium -P acc_mscic1 -n 15 -W 144:00 -R rusage[mem=4000] -R span[hosts=1] -o %J.stdout -eo %J.stderr Rscript ${path}ROSMAP_random_shuffling.R

149459727; 149459728; 149459768; 149464838; 149464844
#/hpc/packages/minerva-centos7/R/4.2.0/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.6: cannot open shared object file: No such file or directory

#1. LOAD AND FORMAT GENE EXPRESSION DATA
library(ggplot2)
library(dplyr)
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)
library(sva)
library(readr)
library(assertthat)

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
date="2024-11-24"

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
    
    # Print progress
    if (i %% 10 == 0 || i == 1 || i == n_shuffle) {  # Print every 10th iteration, the first, and the last
      message(sprintf("Running shuffle %d out of %d...", i, n_shuffle))
    }
    
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
    
    # Shuffle group labels
    shuffled_labels <- sample(group_labels)
    
    # Update design matrix based on shuffled labels
    shuffled_design <- design
    shuffled_design[, "ceradsc1"] <- ifelse(shuffled_labels == "ceradsc1", 1, 0)
    shuffled_design[, "ceradsc4"] <- ifelse(shuffled_labels == "ceradsc4", 1, 0)
    
    # LIMMA without SVA
    fit <- lmFit(expression_data, shuffled_design)
    contrast.matrix <- makeContrasts(ceradsc1 - ceradsc4, levels = shuffled_design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    results <- topTable(fit2, coef = 1, number = Inf)
    
    # Count DEGs
    deg_count <- sum(results$adj.P.Val < 0.05)
    
    # Save results
    results_list[[length(results_list) + 1]] <- data.frame(
      Method = "No SVA",
      Shuffle = i,
      Seed = i,  # Explicitly record the seed
      n_sv = NA,  # Not applicable
      DEGs = deg_count
    )
  }
  
  # Print completion message
  message("All shuffles completed successfully!")
  
  return(do.call(rbind, results_list))
}

#run analysis
n_shuffle <- 100
max_nsv <- 49
results_with_sva <- run_dea_sva_varying_nsv(v, design, design0, n_shuffle, max_nsv)
results_without_sva <- run_dea_no_sva(v, design, n_shuffle)

#combine results
combined_results <- rbind(results_with_sva, results_without_sva)

write.csv(combined_results, file = "/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/random_shuffling_combined_results_with_seed.csv", row.names = FALSE)

