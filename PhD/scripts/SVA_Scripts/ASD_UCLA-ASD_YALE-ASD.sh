#cd /sc/arion/projects/psychgen/resources/psychencode/NRGR/rna/files
#clinical_metadata.tsv  expected_counts_all.rds  metadata.rds  qc_metrics_all.rds  tpm_all.rds
library(GEOquery)
library(limma)
library(sva)
library(readr)
library(data.table)
library(preprocessCore)
library(stringr)
library(dplyr)
library(edgeR)
library(variancePartition)

#0. LOAD AND SUBSET UCLA-ASD AND YALE-ASD
metadata_all <- read.delim("/sc/arion/projects/mscic1/results/jolie/PsychENCODE/clinical_metadata.tsv", header = TRUE) #dont seem to have much info 
table(metadata_all$study)

metadata <- readRDS("/sc/arion/projects/mscic1/results/jolie/PsychENCODE/metadata.rds")
qc_metrics_all <- readRDS("/sc/arion/projects/mscic1/results/jolie/PsychENCODE/qc_metrics_all.rds")

##These 2 are identical and both contain gene expression info
tpm_all <- readRDS("/sc/arion/projects/mscic1/results/jolie/PsychENCODE/tpm_all.rds") #this is a list 
expected_counts_all <- readRDS("/sc/arion/projects/mscic1/results/jolie/PsychENCODE/expected_counts_all.rds") #this is a list
identical(tpm_all, expected_counts_all) #TRUE
summary(expected_counts_all)
str(expected_counts_all) #bingo

#subset metadata by study; $diagnosis, this is more comprehensive
brainGVEX <- subset(qc_metrics_all, study == "BrainGVEX") #73 BP, 1 BP (not BP) - huh?, 259 NC, 95 SZ
libd <- subset(qc_metrics_all, study == "LIBD_szControl") #318 NC, 175 SZ
ucla_asd <- subset(qc_metrics_all, study == "UCLA-ASD") #120 ASD, 133 NC 
yale_asd <- subset(qc_metrics_all, study == "Yale-ASD") #15 ASD, 30 NC 

######################################################
#######################YALE-ASD#######################
######################################################
#1. LOAD AND FORMAT METADATA AND GENE EXPRESSION DATA
table(yale_asd$diagnosis) #15 ASD and 30 HC
yale_asd$sex #has 4 NAs 
a<-yale_asd[,c('individualID','sex')]
#IDs that don't have sex: HSB270, HSB275, HSB292, HSB332

yale_asd <- yale_asd[!is.na(yale_asd$sex), ]
table(yale_asd$diagnosis) #after remove sex, 15 ASD and 26 HC

yale_asd$ageDeath #has PCW --> need to divide that number by 52.14 (365 days/7 weeks)

yale_asd$ageDeath_yr<-c(18.9, 5.6, 20.8, 8.8, 8.8, 14.4, 
						19.1, 16.7, 15.3, 15.8, 9.6, 5.4,
						0.3068661, 48, 0.3644035, 0.2848101, 0.3835827, 0.3279632,
						0.2915228, 0.4027618, 85, 49, 0.3145378, 0.1169927,
						0.2953586, 41, 37, 47, 19.1, 16.7,
						20,18.4,15.9,15.8, 8.8, 8.8, 
						5.6, 4.5, 14.4, 15.3, 4.7)
dim(yale_asd) #41 x 245

metadata <- as.data.frame(yale_asd)
rownames(metadata) <- metadata$Sample_name
metadata$diagnosis <- gsub("Autism Spectrum Disorder", "ASD", metadata$diagnosis)
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "ASD"))
metadata$sex <- as.factor(metadata$sex)
metadata$ageDeath_yr<-as.numeric(metadata$ageDeath_yr)

expected_counts_all$yale_asd #dim: 57820 x 46
gene_expression <- as.data.frame(expected_counts_all$yale_asd)

rownames(gene_expression) <-gene_expression$gene_id
gene_expression$gene_id <- NULL
gene_expression$DLPFC_HSB270 <- NULL
gene_expression$DLPFC_HSB275 <- NULL
gene_expression$DLPFC_HSB292 <- NULL
gene_expression$DLPFC_HSB332 <- NULL

dim(gene_expression) #57820  x  41

gene_expression <- gene_expression[, colnames(gene_expression) %in% rownames(metadata)]
gene_expression <- gene_expression[, match(rownames(metadata), colnames(gene_expression))] #match order

#2. VOOM AND LIMMA ANALYSES PREPARATION AND APPLICATION TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
design <- model.matrix(~ 0 + diagnosis + ageDeath_yr + sex, data = metadata); head(design)
colnames(design) <- c("Control", "ASD", "age","sex")

##2A. VPA START
# Remove genes with zero or near-zero variance
rv <- rowVars(gene_expression)
keep_gene <- is.finite(rv) & rv > 0
gene_expression <- gene_expression[keep_gene, , drop = FALSE]

formula <- ~ diagnosis + ageDeath_yr + sex 
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/ASD/yale_asd_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
## VPA ENDS 

isexpr = rowSums(cpm(gene_expression)>=1) >= 0.1*ncol(gene_expression) 
dge <- DGEList(counts=gene_expression[isexpr,]) 
dge <- calcNormFactors(dge) 
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design) 
contrast.matrix <- makeContrasts(ASD_vs_control = ASD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
results <- topTable(fit2, coef = "ASD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results) #rename rownames to 1, 2, 3, ... : nrows
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_without_sva_09232024.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#double checking if i did filter out for lowly expressed genes - i did.
# results[results$X == "ENSG00000165661.11", ]
# results2 <- read.table(
#   "/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_without_sva_09232024.txt",
#   sep = "\t",
#   header = TRUE
# )
# results2[results2$X == "ENSG00000165661.11", ]

#4. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "BE" AND "LEEK" METHODS
design <- model.matrix(~ 0 + diagnosis + ageDeath_yr + sex, data = metadata)
colnames(design) <- c("Control", "ASD", "age", "sex")
design0 <- model.matrix(~ ageDeath_yr + sex, data = metadata)  

n.sv <- num.sv(v$E, design, method = "leek") #4 SVs for "be"; 3 SVs for "leek"

if (n.sv > 0) {
  svobj <- sva(v$E, design, design0, n.sv = n.sv)
  sv_columns <- svobj$sv
  colnames(sv_columns) <- paste0('sv', 1:ncol(sv_columns))
  design_sv <- cbind(design, sv_columns)  #Add surrogate variables to the design matrix
} else {
  design_sv <- design  #No surrogate variables detected, proceed with original design
}

#5. LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv)
contrast.matrix <- makeContrasts(ASD_vs_control = ASD - Control, levels = design_sv)
fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)

#6. SUMMARIZE AND SAVE DEA WITH SVA RESULTS
results <- topTable(fit2_sv, coef = "ASD_vs_control", adjust.method = "fdr", number = Inf)
results$X <- rownames(results)
results <- results[, c(7, 1, 2, 3, 4, 5, 6)]  # Reorder columns
row.names(results) <- 1:nrow(results)
head(results, 10)
 
write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_with_sva_leek.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/Yale-ASD_DEA_results_with_sva_09232024.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

######################################################
#######################UCLA-ASD#######################
######################################################
#1. LOAD AND FORMAT METADATA AND GENE EXPRESSION DATA
table(ucla_asd$diagnosis) #120 ASD and 133 HC
#have no missing value for age and sex, proceed

metadata <- as.data.frame(ucla_asd) #253 x 244
rownames(metadata) <- metadata$Sample_name
metadata$diagnosis <- gsub("Autism Spectrum Disorder", "ASD", metadata$diagnosis)
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "ASD"))
metadata$sex <- as.factor(metadata$sex)
metadata$ageDeath<-as.numeric(metadata$ageDeath)

expected_counts_all$asd #dim: 57820 x 254
gene_expression <- as.data.frame(expected_counts_all$asd)
rownames(gene_expression) <-gene_expression$gene_id
gene_expression$gene_id <- NULL

dim(gene_expression) #57820  x  253

gene_expression <- gene_expression[, colnames(gene_expression) %in% rownames(metadata)]
gene_expression <- gene_expression[, match(rownames(metadata), colnames(gene_expression))] #match order

#2. VOOM AND LIMMA ANALYSES PREPARATION AND APPLICATION TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
design <- model.matrix(~ 0 + diagnosis + ageDeath + sex, data = metadata); head(design)
colnames(design) <- c("Control", "ASD", "age","sex")

##2A. VPA START
# Remove genes with zero or near-zero variance
rv <- rowVars(gene_expression)
keep_gene <- is.finite(rv) & rv > 0
gene_expression <- gene_expression[keep_gene, , drop = FALSE]

formula <- ~ diagnosis + ageDeath + sex
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/ASD/ucla_asd_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
## VPA ENDS 

isexpr = rowSums(cpm(gene_expression)>=1) >= 0.1*ncol(gene_expression) 
dge <- DGEList(counts=gene_expression[isexpr,]) 
dge <- calcNormFactors(dge) 
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)

contrast.matrix <- makeContrasts(ASD_vs_control = ASD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#3. SUMMARIZE AND SAVE DEA WITHOUT SVA RESULTS 
results <- topTable(fit2, coef = "ASD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results) #rename rownames to 1, 2, 3, ... : nrows
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_without_sva_09232024.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#4. RUN SURROGATE VARIABLE ANALYSIS (SVA) USING "BE" AND "LEEK" METHODS
design <- model.matrix(~ 0 + diagnosis + ageDeath + sex, data = metadata)
colnames(design) <- c("Control", "ASD", "age", "sex")
design0 <- model.matrix(~ ageDeath + sex, data = metadata) 

n.sv <- num.sv(v$E, design, method = "leek") #13 SVs for "be"; 2 SVs for "leek"

if (n.sv > 0) {
  svobj <- sva(v$E, design, design0, n.sv = n.sv)
  sv_columns <- svobj$sv
  colnames(sv_columns) <- paste0('sv', 1:ncol(sv_columns))
  design_sv <- cbind(design, sv_columns)  #Add surrogate variables to the design matrix
} else {
  design_sv <- design  #No surrogate variables detected, proceed with original design
}

#5.LIMMA ANALYSIS TO IDENTIFY DIFFERENTIALLY EXPRESSED GENES
fit_sv <- lmFit(v, design_sv)
contrast.matrix <- makeContrasts(ASD_vs_control = ASD - Control, levels = design_sv)
fit2_sv <- contrasts.fit(fit_sv, contrast.matrix)
fit2_sv <- eBayes(fit2_sv)

#6. SUMMARIZE AND SAVE DEA WITH SVA RESULTS
results <- topTable(fit2_sv, coef = "ASD_vs_control", adjust.method = "fdr", number = Inf)
results$X <- rownames(results)
results <- results[, c(7, 1, 2, 3, 4, 5, 6)]  # Reorder columns
row.names(results) <- 1:nrow(results)
head(results, 10)
 
write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_with_sva_leek.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/ASD/UCLA-ASD_DEA_results_with_sva_09232024.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
