library(biomaRt)
library(data.table)
library(ggplot2)
library(readr)
library(readxl)

public_bp<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/public_BP_signatures.csv", data.table = FALSE)
dim(public_bp) #3413 x 3

public_bp <- read_csv("Desktop/DEA for SVA/public_BP_signatures.csv")
dim(public_BP_signatures) #4342 x3

colnames(public_bp)[colnames(public_bp) == "group"] <- "dataset"

#ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast", host = "http://useast.ensembl.org")
result_public_bp <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "hgnc_symbol", #hugo gene symbol
                          values = public_bp$symbol,
                          mart = ensembl)
dim(result_public_bp) #2495 x 2 -- why is there 1000 less? @Noam
head(result_public_bp)

public_bp=merge(unique(result_public_bp[,1:2]),public_bp,by.x="hgnc_symbol",by.y="symbol",all.y=TRUE)

colnames(result_public_bp)[colnames(result_public_bp) == "ensembl_gene_id"] <- "public_ensembl_id"
result_public_bp$public_DE <- 1
head(result_public_bp) 

#check for duplicates: 
any(duplicated(result_public_bp$public_ensembl_id))
#FALSE - no dup