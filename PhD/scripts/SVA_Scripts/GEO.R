###GEO: GENE EXPRESSION OMNIBUS (https://www.ncbi.nlm.nih.gov/geo/)
#https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#GEOquery is the bridge between GEO and BioConductor.

#8 total datasets
###################################### BASE CODE ####################################
library(GEOquery)

# Define the GEO accession number
gse_id <- "GSE8397"

# Download the GEO dataset
gse <- getGEO(gse_id, GSEMatrix = TRUE)

# Extract the expression set
expression_set <- gse[[1]]

# Extract metadata
metadata <- pData(expression_set)

# Extract gene expression matrix
gene_expression <- exprs(expression_set)

# Save the metadata and gene expression matrix to files
write.table(metadata, file = "metadata_GSE8397.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(gene_expression, file = "gene_expression_GSE8397.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Print the first 10 rows of metadata and gene expression matrix
head(metadata, 10)
head(gene_expression, 10)


#in gene expression matrix,
# Columns are GSM SAMPLE Identifiers (GSM184354, GSM184355, etc), and Rows are Probe Identifiers (1053_at, 117_at) correspond to specific probe sets on the microarray chip, and the rows labeled with these identifiers represent the expression levels of these probes across different samples.

########################################### END ###############################################
###############################################################################################

#SEARCH FREAD TO FIND CORRELATIONS!!!
library(GEOquery)
library(limma)
library(sva)
library(readr)
library(data.table)
library(preprocessCore)
library(stringr)
library(dplyr)
library(variancePartition)
library(BiocParallel)

##############################GSE8397, 15 donors: 8 NC, 14PD?; samples n=18 control, 29 PD, no RIN?
#adjust for brain region - double check whether they're the same people or not
#multiple samples per individual 
#lmfit - fixed effect model burns a lot degree of freedom; dream - random effect only burn degree of freedom
gse_id <- "GSE8397"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix

#remove 3 people that might be mislabeled
columns_to_remove <- c("GSM208627", "GSM208644", "GSM208664")
gene_expression <- gene_expression[, !(colnames(gene_expression) %in% columns_to_remove)]

#write.table(metadata, file = "metadata_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)
#metadata <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/metadata_GSE8397.txt", header = TRUE, sep = "\t")
#gene_expression <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/gene_expression_GSE8397.txt", header = TRUE, sep = "\t")

metadata$case_number <- str_extract(metadata$title, "case \\d+")
metadata$case_number <- str_replace(metadata$case_number, "case ", "")

# Assign value 99 to the 'case_number' column for the row 'GSM208635'
metadata["GSM208635", "case_number"] <- 99
metadata["GSM208652", "case_number"] <- 99

metadata["GSM208645", "case_number"] <- 90
metadata["GSM208668", "case_number"] <- 90

metadata["GSM208625", "case_number"] <- 100
metadata["GSM208636", "case_number"] <- 100
metadata["GSM208653", "case_number"] <- 100

metadata["GSM208637", "case_number"] <- 200
metadata["GSM208654", "case_number"] <- 200

metadata["GSM208655", "case_number"] <- 400
metadata["GSM208626", "case_number"] <- 400
metadata["GSM208638", "case_number"] <- 400

metadata["GSM208639", "case_number"] <- 700
metadata["GSM208656", "case_number"] <- 700

metadata["GSM208640", "case_number"] <- 900
metadata["GSM208657", "case_number"] <- 900

metadata["GSM208630", "case_number"] <- 11
metadata["GSM208646", "case_number"] <- 11

##clean up disease status
metadata$disease_status<-metadata$title
metadata$disease_status <- ifelse(grepl("Parkinson's disease", metadata$title), "PD", metadata$title)
metadata$disease_status <- ifelse(metadata$disease_status != "PD", "control", metadata$disease_status)
metadata$disease_status<-as.factor(metadata$disease_status)
table(metadata$disease_status)

#region
metadata$region <- str_extract(metadata$description, "(?<=Gene expression data from ).*(?= \\()")
# Specify the rows to be updated and the new region value
la_rows_to_change <- c("GSM208630", "GSM208631", "GSM208633","GSM208634","GSM208636","GSM208637",
  "GSM208638","GSM208639","GSM208640","GSM208641","GSM208645")
metadata[la_rows_to_change, "region"] <- "lateral substantia nigra"

me_rows_to_change <- c("GSM208647", "GSM208649", "GSM208652","GSM208654","GSM208655","GSM208660",
  "GSM208661")
metadata[me_rows_to_change, "region"] <- "medial substantia nigra"
metadata$region<-as.factor(metadata$region)

##clean up covariates
metadata$age<-metadata$'age:ch1'
metadata$age <- as.numeric(gsub(";.*", "", metadata$age)) #The regular expression ";.*" matches the semicolon and any characters that follow it until the end of the string. By replacing it with an empty string "", it effectively removes everything after the semicolon in each row.

metadata$gender<-metadata$'age:ch1'
metadata$gender <- as.factor(gsub("[^MF]", "", metadata$gender)) #The regular expression [^MF] matches any character that is not "M" or "F". By replacing it with an empty string "", it effectively removes all characters except "M" and "F" from each row.

metadata$case_number <- as.numeric(metadata$case_number)
metadata_sorted <- metadata %>% arrange(case_number)

metadata_sorted$case_number<-as.factor(metadata_sorted$case_number)
metadata_sorted[,c("age","gender","case_number","disease_status","region")]

rows_to_remove <- c("GSM208627", "GSM208644", "GSM208664") #removed 3 rows that might be mislabeled
# Remove the specified rows from the data frame
metadata_sorted <- metadata_sorted[!rownames(metadata_sorted) %in% rows_to_remove, ]

metadata_sorted$disease_status <- factor(metadata_sorted$disease_status)
metadata_sorted$gender <- factor(metadata_sorted$gender)
metadata_sorted$region <- factor(metadata_sorted$region)
metadata_sorted$case_number <- factor(metadata_sorted$case_number)


#check for number of samples per case number 
table(metadata_sorted$case_number)

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 0.5 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status + age + gender + region, data = metadata_sorted); head(design) #
svd(design)$d #~eigenvalue of pca <10^-10 --> multicoll

colnames(design) <- c("control", "PD", "age", "gender","msn","sfg")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + age + gender + region #in the manuscript there's age too so which version is correct?
vp<- fitExtractVarPartModel(gene_expression, formula, metadata_sorted)
# Error in .fitExtractVarPartModel(exprObj, formula, data, REML = REML,  : 
#   Response variable 5954 has a variance of 0

# Remove genes with zero or near-zero variance
rv <- rowVars(gene_expression)
keep_gene <- is.finite(rv) & rv > 0
gene_expression <- gene_expression[keep_gene, , drop = FALSE]

formula <- ~ disease_status + age + gender + region #in the manuscript there's age too so which version is correct?
vp<- fitExtractVarPartModel(gene_expression, formula, metadata_sorted)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

############## #############  VPA ENDS ############# ############# 

#run random effect #might save it and reload, dont have to do it if run dream
dupcor <- duplicateCorrelation(gene_expression,design,block=metadata_sorted$case_number)
dupcor$consensus.correlation

#fit linear model
fit <- lmFit(gene_expression,design,block=metadata_sorted$case_number,correlation=dupcor$consensus.correlation)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

#negative logFC = the gene is downreugulated in PD relative to control and positive logFC = upregulation in PD
#aveExp = average log2 expression level of the gene across all samples
#t-stat = magnitude of the difference relative to the variation in the data. Higher absolute values indicate more significant differences between conditions.
#B = log-odds of differential expression aka the log-odds that the gene is differentially expressed. Higher values indicate higher confidence that the gene is differentially expressed.
write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + age + gender + region, data = metadata_sorted) 
mod0 <- model.matrix(~ age + gender + region, data = metadata_sorted)
mod1 <- model.matrix(~ 0 + disease_status + age + gender + region, data = metadata_sorted)

sva_method <- "be" #Change to "be" or "LEEK" as needed; "be" has 8 SVs; "LEEK" has 0 Svs

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))
colnames(modSv)[1:6] <- c("control", "PD", "age","gender","msn","sfg")
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[6 + i] <- paste0("SV", i)
}

#run random effect #might save it and reload, dont have to do it if run dream
dupcor <- duplicateCorrelation(gene_expression,modSv,block=metadata_sorted$case_number)
dupcor$consensus.correlation

fit <- lmFit(gene_expression,modSv ,block=metadata_sorted$case_number,correlation=dupcor$consensus.correlation)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

save.image('/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva.Rdata') 


##############################GSE7621, n=9 control, 16 PD, no RIN or age?
gse_id <- "GSE7621" #see GSE7621_sup doc
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix; dim:  54339  x 25
sum(is.na(gene_expression)) #there are 525 NAs in the gene expression matrix
gene_expression <- na.omit(gene_expression) #drop rows with NAs; dim: 54318 x 25 (delta = 21 rows)

metadata[,c("characteristics_ch1.1","source_name_ch1")]
#write.table(metadata, file = "metadata_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$characteristics_ch1 <- ifelse(metadata$characteristics_ch1 == "Parkinson's Disease", "PD", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- ifelse(metadata$characteristics_ch1 == "Old Control", "control", metadata$characteristics_ch1)
metadata$disease_status <- as.factor(metadata$characteristics_ch1)
table(metadata$disease_status)

#covariates
metadata$gender<-metadata$characteristics_ch1.1
metadata$gender <- ifelse(metadata$gender == "female", "F", metadata$gender)
metadata$gender <- ifelse(metadata$gender == "male", "M", metadata$gender)
metadata$gender<-as.factor(metadata$gender); table(metadata$gender)

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 0.5 to avoid taking log of 0

# Quantile normalization
gene_expression<- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status +  gender, data = metadata); head(design)
colnames(design) <- c("control", "PD", "gender")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + gender #in the manuscript there's age too so which version is correct?
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

############## #############  VPA ENDS ############# ############# 
# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + gender, data = metadata)
mod0 <- model.matrix(~ gender, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + gender, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

colnames(modSv)[1:3] <- c("control", "PD", "gender") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[3 + i] <- paste0("SV", i)
}
 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)fg

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

###############################GSE24378, n=9 control, 8 PD, no gender and RIN? removed 4 NAs in age --> n=7 control and 6 PD
gse_id <- "GSE24378"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
rownames(metadata) <- metadata$geo_accession
gene_expression <- exprs(expression_set) # Extract gene expression matrix  61359 x 18
sum(is.na(gene_expression))

#write.table(metadata, file = "metadata_GSE24378.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE24378.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$disease_status <- metadata$'disease state:ch1'
metadata$disease_status <- ifelse(metadata$disease_status == "Parkinson’s disease", "PD", metadata$disease_status)
metadata$disease_status <- as.factor(metadata$disease_status)
table(metadata$disease_status)

#covariates
# Remove samples with NA in age from metadata and align gene expression matrix
metadata$age<-as.numeric(metadata$'age (y):ch1') #4 NAs
metadata <- metadata[!is.na(metadata$age), ] #removed 4 NAs without age
gene_expression <- gene_expression[, colnames(gene_expression) %in% rownames(metadata)]
gene_expression <- gene_expression[, match(rownames(metadata), colnames(gene_expression))]

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 1 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status +  age, data = metadata)
colnames(design) <- c("control", "PD", "age")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + age #in the manuscript there's age too so which version is correct?
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 

# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + age, data = metadata)
mod0 <- model.matrix(~ age, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + age, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

colnames(modSv)[1:3] <- c("control", "PD", "age") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[3 + i] <- paste0("SV", i)
}

head(modSv) 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


################################GSE20292, n=18 control, 11 PD, no RIN?
gse_id <- "GSE20292"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix

#write.table(metadata, file = "metadata_GSE20292.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE20292.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$disease_status<-metadata$'disease state:ch1'
metadata$disease_status<- ifelse(metadata$disease_status == "Parkinsons disease", "PD", metadata$disease_status)
metadata$disease_status <- ifelse(metadata$disease_status == "Control", "control", metadata$disease_status)
metadata$disease_status<-as.factor(metadata$disease_status)
table(metadata$disease_status)

#covariates
metadata$age<-as.numeric(metadata$'age:ch1')
metadata$gender<-as.factor(metadata$gender)

metadata$region <-  gsub("[0-9]+", "", metadata_sorted$title)
metadata$region <- trimws(metadata$region)
metadata["GSM508713", "region"] <- "SN Cm"
metadata["GSM508735", "region"] <- "SN Cm"

#this approach didn't change the region correctly 
#metadata$region<-as.integer(as.factor(metadata$region))
#metadata$region<-as.factor(metadata$region)

#so use this appraoch
cm_rows_to_change <- c("GSM606624", "GSM606625", "GSM606626","GSM508722","GSM508735","GSM508708","GSM521253","GSM508721","GSM508729","GSM508733","GSM508720","GSM508734","GSM508724")
metadata[cm_rows_to_change, "region"] <- 1

pm_rows_to_change <- c("GSM508732", "GSM508715","GSM508718","GSM508710","GSM508728", "GSM508716")
metadata[pm_rows_to_change, "region"] <- 2

pf_rows_to_change <- c("GSM508713", "GSM508731","GSM508711","GSM508714","GSM508712")
metadata[pf_rows_to_change, "region"] <- 3

cf_rows_to_change <- c("GSM508725", "GSM508730","GSM508723","GSM508726","GSM508717")
metadata[cf_rows_to_change, "region"] <- 4

#1 is SN Cm, 2 is SN Pm, 3 is SN Pf, 4 is SN Cf

#just to see
metadata_sorted <- metadata %>% arrange(age)
metadata_sorted[,c("age","gender","disease_status","title","region")]

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 0.5 to avoid taking log of 0

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

# Design matrix for limma
design <- model.matrix(~ 0 + disease_status + age + gender, data = metadata)
colnames(design) <- c("control", "PD", "age", "gender")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + age + gender  
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 

dupcor <- duplicateCorrelation(gene_expression,design,block=metadata$region)
dupcor$consensus.correlation

# Fit the linear model
fit <- lmFit(gene_expression,design,block=metadata$region,correlation=dupcor$consensus.correlation)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + age + gender, data = metadata)
mod0 <- model.matrix(~ age + gender , data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + age + gender, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed #5 SVs in be

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}


# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))
colnames(modSv)[1:4] <- c("control", "PD", "age","gender")
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[4 + i] <- paste0("SV", i)
}

#run random effect #might save it and reload, dont have to do it if run dream
dupcor <- duplicateCorrelation(gene_expression,modSv,block=metadata$region)
dupcor$consensus.correlation

fit <- lmFit(gene_expression,modSv ,block=metadata$region,correlation=dupcor$consensus.correlation)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


################################GSE20141, n=8 control, 10 PD, no age, gender, and RIN? :(
gse_id <- "GSE20141"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix

#write.table(metadata, file = "metadata_GSE20141.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE20141.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$disease_status <- metadata$'disease state:ch1'
metadata$disease_status <- ifelse(metadata$disease_status == "Parkinson's disease", "PD", metadata$disease_status)
metadata$disease_status <- as.factor(metadata$disease_status)
table(metadata$disease_status)

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 1 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status, data = metadata)
colnames(design) <- c("control", "PD")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status 
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

# Remove genes with zero or near-zero variance
rv <- rowVars(gene_expression)
keep_gene <- is.finite(rv) & rv > 0
gene_expression <- gene_expression[keep_gene, , drop = FALSE]

formula <- ~ disease_status 
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 

# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status, data = metadata)


sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))
colnames(modSv)[1:2] <- c("control", "PD") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[2 + i] <- paste0("SV", i)
}
 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


################################GSE20163, n=9 control, 8 PD, no gender or RIN - contact Bin Zheng at Brigham?
#sup file: ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM505nnn/GSM505996/suppl/GSM505996.CEL.gz
gse_id <- "GSE20163"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix
sum(is.na(gene_expression))

#write.table(metadata, file = "metadata_GSE20163.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE20163.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$characteristics_ch1.2 <- ifelse(metadata$characteristics_ch1.2 == "diagnosis: Parkinson's disease", "PD", metadata$characteristics_ch1.2)
metadata$characteristics_ch1.2 <- ifelse(metadata$characteristics_ch1.2 == "sample type: control", "control", metadata$characteristics_ch1.2)
metadata$disease_status<-as.factor(metadata$characteristics_ch1.2)
table(metadata$characteristics_ch1.2)

#covariates
# Remove samples with NA in age from metadata and align gene expression matrix
metadata$age<-as.numeric(metadata$'age:ch1') #5 NAs
metadata <- metadata[!is.na(metadata$age), ] #removed 5 NAs without age
gene_expression <- gene_expression[, colnames(gene_expression) %in% rownames(metadata)]
gene_expression <- gene_expression[, match(rownames(metadata), colnames(gene_expression))]

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 1 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status + age, data = metadata)
colnames(design) <- c("control", "PD", "age")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + age
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

# # Remove genes with zero or near-zero variance
# rv <- rowVars(gene_expression)
# keep_gene <- is.finite(rv) & rv > 0
# gene_expression <- gene_expression[keep_gene, , drop = FALSE]

# formula <- ~ disease_status 
# vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 

# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + age, data = metadata)
mod0 <- model.matrix(~ age, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + age, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

colnames(modSv)[1:3] <- c("control", "PD", "age") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[3 + i] <- paste0("SV", i)
}
 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


################################GSE20164, n=5 control, 6 PD
gse_id <- "GSE20164"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix

#write.table(metadata, file = "metadata_GSE20164.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE20164.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$'diagnosis:ch1' <- ifelse(metadata$'diagnosis:ch1' == "Parkinson's disease", "PD", metadata$'diagnosis:ch1')
metadata$'diagnosis:ch1' <- ifelse(metadata$'diagnosis:ch1' == "normal", "control", metadata$'diagnosis:ch1')
metadata$disease_status <- as.factor(metadata$'diagnosis:ch1')
table(metadata$disease_status)

#covariates
metadata$age <- as.numeric(metadata$'age:ch1')
metadata$gender <- as.factor(metadata$'gender:ch1')
##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 1) #add 1 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status + age + gender, data = metadata)
colnames(design) <- c("control", "PD", "age", "gender")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + age + gender
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

# # Remove genes with zero or near-zero variance
# rv <- rowVars(gene_expression)
# keep_gene <- is.finite(rv) & rv > 0
# gene_expression <- gene_expression[keep_gene, , drop = FALSE]

# formula <- ~ disease_status 
# vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 

# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + age + gender, data = metadata)
mod0 <- model.matrix(~ age + gender, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + age + gender, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

colnames(modSv)[1:4] <- c("control", "PD", "age","gender") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[4 + i] <- paste0("SV", i)
}
 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


################################GSE49036, n=8 control, 5 Lewy body disease, 15 PD, no age or gender
#genes = probe set IDs from the Affymetrix GeneChip Human Genome U133 Plus 2.0 Array
gse_id <- "GSE49036"
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata, 28 x 36
gene_expression <- exprs(expression_set) # Extract gene expression matrix, 54675 x 28

metadata[,c("tissue:ch1","disease state:ch1", "source_name_ch1","rin:ch1","title")]
#write.table(metadata, file = "metadata_GSE49036.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE49036.txt", sep = "\t", row.names = TRUE, col.names = NA)

#disease status
# Change "Parkinson's disease" to "PD" in the metadata
metadata <- metadata[metadata$'disease state:ch1' != "incidental Lewy body disease", ]
metadata$disease_status <- metadata$'disease state:ch1'
metadata$disease_status <- ifelse(metadata$disease_status == "Parkinson's disease", "PD", metadata$disease_status)
metadata$disease_status <- as.factor(metadata$disease_status)
table(metadata$disease_status)

#cannot control for braak as its info is contained in the disease status; braak: 1 is NC (stage 0), 2 is Braak α-synuclein 3-4, 3 is Braak α-synuclein 5-6

#covariates
metadata$rin <- as.numeric(metadata$'rin:ch1') #range: 60 to 81

gene_expression <- gene_expression[, rownames(metadata)]

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 1 to avoid taking log of 0, keep first col the same

# Quantile normalization
gene_expression <- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status + rin, data = metadata)
colnames(design) <- c("control", "PD", "rin")

############# ############# VPA START ############# ############# 
formula <- ~ disease_status + rin
vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

# # Remove genes with zero or near-zero variance
# rv <- rowVars(gene_expression)
# keep_gene <- is.finite(rv) & rv > 0
# gene_expression <- gene_expression[keep_gene, , drop = FALSE]

# formula <- ~ disease_status 
# vp<- fitExtractVarPartModel(gene_expression, formula, metadata)

write.table(vp, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_vpa.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
############## #############  VPA ENDS ############# ############# 


# Fit the linear model
fit <- lmFit(gene_expression, design)
#Coefficients not estimable: braak3 
#Warning message:
#Partial NA coefficients for 54675 probe(s)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#SVA
mod <- model.matrix(~ disease_status + rin, data = metadata)
mod0 <- model.matrix(~ rin, data = metadata)
mod1 <- model.matrix(~ 0 + disease_status + rin, data = metadata)

sva_method <- "be" #Change to "be" or "LEEK" as needed

if (sva_method == "be") {
  svobj <- sva(as.matrix(gene_expression), mod, mod0)
} else if (sva_method == "LEEK") {
  n.sv <- num.sv(as.matrix(gene_expression), mod, method = "leek")
  svobj <- sva(as.matrix(gene_expression), mod, mod0, n.sv = n.sv)
} else {
  stop("Invalid SVA method chosen. Use 'be' or 'LEEK'.")
}

# Include surrogate variables in the design matrix
modSv <- cbind(mod1, svobj$sv)
# Ensure column names are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

colnames(modSv)[1:3] <- c("control", "PD", "rin") ##change this for other df
for (i in 1:ncol(svobj$sv)) {
  colnames(modSv)[3 + i] <- paste0("SV", i)
}
 
# Fit the linear model with surrogate variables
fit <- lmFit(gene_expression, modSv)
  
# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = modSv)
  
# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$X<-rownames(results)
head(results, 10) 
results<- results[, c(7,1,2,3,4,5,6)]
row.names(results) <- 1:nrow(results)
head(results, 10) 

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_with_sva_be_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

######################CORRELATIONS########################
##########################################################
GSE8397_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE8397_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE7621_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE7621_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE24378_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE24378_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20292_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20292_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20141_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20141_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20163_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20163_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE20164_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE20164_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE49036_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
GSE49036_sva <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")


combo_no_sva<-merge(GSE8397_none,GSE7621_none, by='X', all.x=TRUE, all.y=TRUE) #22283x13
combo_sva<-merge(GSE8397_sva,GSE7621_sva, by='X', all.x=TRUE, all.y=TRUE) #22283x13

cor.test(combo_no_sva$logFC.x, combo_no_sva$logFC.y,method="spearman") #r=0.2091094
cor.test(combo_sva$logFC.x, combo_sva$logFC.y,method="spearman") #r=0.6351632 UNEXPECTED RESULTS!!!!


# Create lists for _none and _sva datasets
datasets <- list(
  GSE8397_none = GSE8397_none,
  GSE7621_none = GSE7621_none,
  GSE24378_none = GSE24378_none,
  GSE20292_none = GSE20292_none,
  GSE20141_none = GSE20141_none,
  GSE20163_none = GSE20163_none,
  GSE20164_none = GSE20164_none,
  GSE49036_none = GSE49036_none
)

datasets <- list(
  GSE8397_sva = GSE8397_sva,
  GSE7621_sva = GSE7621_sva,
  GSE24378_sva = GSE24378_sva,
  GSE20292_sva = GSE20292_sva,
  GSE20141_sva = GSE20141_sva,
  GSE20163_sva = GSE20163_sva,
  GSE20164_sva = GSE20164_sva,
  GSE49036_sva = GSE49036_sva
)

###single cor: 
combo_no_sva<-merge(GSE8397_none,GSE7621_none, by='X', all.x=TRUE, all.y=TRUE) #22283x13
combo_sva<-merge(GSE8397_sva,GSE7621_sva, by='X', all.x=TRUE, all.y=TRUE) #22283x13

cor.test(combo_no_sva$logFC.x, combo_no_sva$logFC.y,method="spearman") #r=0.009507782 
cor.test(combo_sva$logFC.x, combo_sva$logFC.y,method="spearman") #r=-0.009109724 UNEXPECTED RESULTS!!!!

######################################MULTIPLE CORRELATIONS####################################################
###############################################################################################################
###this matches the manual correlation! 
#run this TWICE, each time for the the datasets above
calculate_spearman_correlations <- function(datasets) {
  dataset_names <- names(datasets)
  n <- length(datasets)
  cor_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(dataset_names, dataset_names))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      data_i <- datasets[[i]]
      data_j <- datasets[[j]]
      merged_data <- merge(data_i, data_j, by = 'X', all.x = TRUE, all.y = TRUE)
      cor_value <- cor(merged_data$logFC.x, merged_data$logFC.y, method = "spearman", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
    }
  }
  return(cor_matrix)
}

spearman_cor_matrix <- calculate_spearman_correlations(datasets)
diag(spearman_cor_matrix) <- 1
spearman_cor_matrix

write.table(spearman_cor_matrix , file = "/sc/arion/projects/mscic1/results/jolie/GEO/no_sva_spearman_cor_matrix_final_final.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(spearman_cor_matrix, file = "/sc/arion/projects/mscic1/results/jolie/GEO/sva_spearman_cor_matrix_final_final.txt", sep = "\t", row.names = TRUE, col.names = NA)

###HEATMAP!
library(ggplot2)
library(reshape2)
library(seriation)

no_sva_spearman_cor_matrix <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/no_sva_spearman_cor_matrix_final_final.txt", header = TRUE, sep = "\t")
sva_spearman_cor_matrix <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/sva_spearman_cor_matrix_final_final.txt", header = TRUE, sep = "\t")


rownames(no_sva_spearman_cor_matrix) <- no_sva_spearman_cor_matrix$X
rownames(sva_spearman_cor_matrix) <- sva_spearman_cor_matrix$X

no_sva_spearman_cor_matrix$X <- NULL
sva_spearman_cor_matrix$X <- NULL

rownames(no_sva_spearman_cor_matrix) <- gsub("_none$", "", rownames(no_sva_spearman_cor_matrix))
rownames(sva_spearman_cor_matrix) <- gsub("_sva$", "", rownames(sva_spearman_cor_matrix))

#make sure they have the same col name
colnames(no_sva_spearman_cor_matrix) <- sub("_none$", "", colnames(no_sva_spearman_cor_matrix))
colnames(sva_spearman_cor_matrix) <- sub("_sva$", "", colnames(sva_spearman_cor_matrix))


# Convert correlation matrix to long format for ggplot
sva_cor_matrix_long <- melt(as.matrix(sva_spearman_cor_matrix))
no_sva_cor_matrix_long <- melt(as.matrix(no_sva_spearman_cor_matrix))

##repeat the same thing for no_sva_spearman_cor_matrix AND sva_spearman_cor_matrix
names(sva_cor_matrix_long) <- c("Dataset1", "Dataset2", "value")
names(no_sva_cor_matrix_long) <- c("Dataset1", "Dataset2", "value")


order<-c("GSE20292", "GSE8397", "GSE20164", "GSE20163","GSE24378", "GSE7621", "GSE20141", "GSE49036")
# Ensure the levels of Dataset1 and Dataset2 match the order of dataset names
#dataset_order <- names(datasets) #make sure the right datasets is loaded so your graph is not a blob lol 
sva_cor_matrix_long$Dataset1 <- factor(sva_cor_matrix_long$Dataset1, levels = order)
sva_cor_matrix_long$Dataset2 <- factor(sva_cor_matrix_long$Dataset2, levels = order)

no_sva_cor_matrix_long$Dataset1 <- factor(no_sva_cor_matrix_long$Dataset1, levels = order)
no_sva_cor_matrix_long$Dataset2 <- factor(no_sva_cor_matrix_long$Dataset2, levels = order)

# Generate heatmap
heatmap_plot <- ggplot(no_sva_cor_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid= "white",
                       limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        axis.title.x = element_text(size = 16),  # Adjust x-axis title size
        axis.title.y = element_text(size = 16)) +  # Adjust y-axis title size
  coord_fixed() +
  labs(title = "Dataset Replicability without SVA") + 
  xlab("Datasets") + ylab("Datasets"); heatmap_plot 


ggsave("/hpc/users/hoangd02/www/plots/GEO/no_sva_spearman_cor_matrix_heatmap_final.png", plot = heatmap_plot, width = 10, height = 10)

#in bash:
cd /hpc/users/hoangd02/www/plots/GEO


##MAKING THE DIFFERENCE MATRIX
library(seriation) #create cluster based on data!

rownames(no_sva_spearman_cor_matrix) <- no_sva_spearman_cor_matrix$X
rownames(sva_spearman_cor_matrix) <- sva_spearman_cor_matrix$X

no_sva_spearman_cor_matrix$X <- NULL
sva_spearman_cor_matrix$X <- NULL

rownames(no_sva_spearman_cor_matrix) <- gsub("_none$", "", rownames(no_sva_spearman_cor_matrix))
rownames(sva_spearman_cor_matrix) <- gsub("_sva$", "", rownames(sva_spearman_cor_matrix))

#make sure they have the same col name
colnames(no_sva_spearman_cor_matrix) <- sub("_none$", "", colnames(no_sva_spearman_cor_matrix))
colnames(sva_spearman_cor_matrix) <- sub("_sva$", "", colnames(sva_spearman_cor_matrix))

common_columns <- intersect(colnames(sva_spearman_cor_matrix), colnames(no_sva_spearman_cor_matrix))
sva_spearman_cor_matrix <- sva_spearman_cor_matrix[, common_columns]
no_sva_spearman_cor_matrix <- no_sva_spearman_cor_matrix[, common_columns]

sva_matrix <- as.matrix(sva_spearman_cor_matrix)
no_sva_matrix <- as.matrix(no_sva_spearman_cor_matrix)

# Create the difference matrix
difference_matrix <- no_sva_matrix - sva_matrix
rowdist <-dist(difference_matrix)
coldist <-dist(t(difference_matrix))
roworder<-seriate(rowdist)
colorder<-seriate(coldist)

#PLOT
difference_matrix_long <- melt(difference_matrix)

# Rename columns to match the expected format
names(difference_matrix_long) <- c("Dataset1", "Dataset2", "value")

# Ensure the levels of Dataset1 and Dataset2 match the order of dataset names
dataset_order <- rownames(difference_matrix)
difference_matrix_long$Dataset1 <- factor(difference_matrix_long$Dataset1, levels = rownames(difference_matrix)[unlist(roworder)])
difference_matrix_long$Dataset2 <- factor(difference_matrix_long$Dataset2, levels = colnames(difference_matrix)[unlist(colorder)])


# Generate the heatmap
heatmap_plot <- ggplot(difference_matrix_long, aes(Dataset1, Dataset2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid= "white",  # purple to yellow
                       limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        axis.title.x = element_text(size = 16),  # Adjust x-axis title size
        axis.title.y = element_text(size = 16)) +  # Adjust y-axis title size
  coord_fixed() +
  labs(title = "Difference between LogFC Correlations Without SVA and With SVA", x = "Datasets", y = "Datasets"); heatmap_plot

ggsave("/hpc/users/hoangd02/www/plots/GEO/difference_no_sva_spearman_cor_matrix_heatmap_final.png", plot = heatmap_plot, width = 10, height = 10)






















