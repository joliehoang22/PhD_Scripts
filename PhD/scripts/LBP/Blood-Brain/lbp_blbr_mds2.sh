## MDS
library(data.table)
library(edgeR)
library(limma)
library(assertthat)
library(foreach)
library(doParallel)
library(ggplot2)

#the GTEx data v8 is downloaded here
setwd("/sc/arion/projects/mscic1/results/Noam/gtex")
#setwd("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain")
#made a screen section of MDS here 

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  #' Convert counts to transcripts per million (TPM).
  #' 
  #' Convert a numeric matrix of features (rows) and conditions (columns) with
  #' raw feature counts to transcripts per million.
  #' 
  #'    Lior Pachter. Models for transcript quantification from RNA-Seq.
  #'    arXiv:1104.3889v2 
  #'    
  #'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
  #'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
  #'    doi:10.1007/s12064-012-0162-3
  #'    
  #' @param counts A numeric matrix of raw feature counts i.e.
  #'  fragments assigned to each gene.
  #' @param featureLength A numeric vector with feature lengths.
  #' @param meanFragmentLength A numeric vector with mean fragment lengths.
  #' @return tpm A numeric matrix normalized by library size and feature length.

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
    }

#getting the raw counts for ROSMAP
genedata_org_rosmap=fread("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmapv43/expression/allcount_matrix_2023-04-11.txt", data.table=FALSE)
bed_rosmap=genedata_org_rosmap[,1:6]
#first 6 columns have Geneid (ensembl), Chr, Start, End, Strand, Length 
bed_rosmap$ensembl=unlist(lapply(strsplit(bed_rosmap$Geneid,".",fixed=TRUE),function(x){x[1]}))
#make a new column, reformatting the ensembl id ENSG00000284332.1 --> ENSG00000284332
genedata_org_rosmap=genedata_org_rosmap[7:ncol(genedata_org_rosmap)]
rownames(genedata_org_rosmap)=bed_rosmap$Geneid

#getting the raw counts for MSBB
genedata_org_msbb=fread("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/msbb/MSBBv43/expression/allcount_matrix_2023-08-29.txt", data.table=FALSE)
bed_msbb=genedata_org_msbb[,1:6]
bed_msbb$ensembl=unlist(lapply(strsplit(bed_msbb$Geneid,".",fixed=TRUE),function(x){x[1]}))
genedata_org_msbb=genedata_org_msbb[7:ncol(genedata_org_msbb)]
rownames(genedata_org_msbb)=bed_msbb$Geneid

#getting the raw counts for GTEx
geneexp_all=fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",data.table=FALSE)
rownames(geneexp_all)=geneexp_all$Name
geneexp_all$Name=NULL
geneexp_all$Description=NULL
bed_all=data.frame(Geneid=rownames(geneexp_all))
bed_all$ensembl=unlist(lapply(strsplit(bed_all$Geneid,".",fixed=TRUE),function(x){x[1]}))

#getting the raw counts for BRAIN LBP
genedata_org_lbp=fread("/sc/arion/projects/mscic1/results/Noam/LBP_imaging/expression/geneCounts_2022-07-26.txt",data.table=FALSE)
#creating a new dataframe called bed_lbp that has only 1 column called Geneid
bed_lbp=data.frame(Geneid=genedata_org_lbp$V1)
rownames(genedata_org_lbp)=bed_lbp$Geneid
genedata_org_lbp$V1=NULL #remove 1st column
bed_lbp$ensembl=unlist(lapply(strsplit(bed_lbp$Geneid,".",fixed=TRUE),function(x){x[1]}))

#getting the raw counts for BLOOD LBP
genedata_blood_lbp=fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_expression_raw_ALL_samples_not_paired_samples.txt",data.table=FALSE)
rownames(genedata_blood_lbp)=genedata_blood_lbp$V1
#creating a new dataframe called bed_lbp that has only 1 column called Geneid
bed_blood_lbp=data.frame(Geneid=genedata_blood_lbp$V1)
genedata_blood_lbp$V1 <- NULL
bed_blood_lbp$ensembl=unlist(lapply(strsplit(bed_blood_lbp$Geneid,".",fixed=TRUE),function(x){x[1]}))

length(intersect(bed_msbb$Geneid,bed_lbp$Geneid)) #43555

identical(bed_msbb,bed_rosmap) #TRUE - why is this identical?

#frequency of how many times each ensembl ID appears in bed_lbp$ensembl, etc
table(table(bed_rosmap$ensembl)) 
table(table(bed_msbb$ensembl)) 
table(table(bed_all$ensembl)) #GTEx
table(table(bed_lbp$ensembl))
#58839 ensembl IDs appear exactly once in bed_lbp$ensembl
#45 ensembl IDs appears twice in bed_lbp$ensembl - these are multiple timepoints?
table(table(bed_blood_lbp$ensembl)) #this is same as brain

##keeping only genes that appear once
bed_msbb=bed_msbb[bed_msbb$ensembl %in% names(which(table(bed_msbb$ensembl)==1)),]
bed_lbp=bed_lbp[bed_lbp$ensembl %in% names(which(table(bed_lbp$ensembl)==1)),]
bed_blood_lbp=bed_blood_lbp[bed_blood_lbp$ensembl %in% names(which(table(bed_blood_lbp$ensembl)==1)),]
bed_all=bed_all[bed_all$ensembl %in% names(which(table(bed_all$ensembl)==1)),]

table(table(bed_msbb$ensembl)) 
table(table(bed_all$ensembl)) #GTEx
table(table(bed_lbp$ensembl))

bed=merge(bed_msbb,bed_lbp,by="ensembl",suffixes=c(".msbb.rosmap",".lbpBrain.lbpBlood"))
colnames(bed_all)[colnames(bed_all)=="Geneid"]="Geneid.gtex"
bed=merge(bed,bed_all,by="ensembl")

identical(bed$Geneid.lbpBrain.lbpBlood,bed$Geneid.gtex) #FALSE
identical(bed$Geneid.msbb,bed$Geneid.gtex) #FALSE

#subsetting genedata_org_rosmap based on the row names that match the values in bed$Geneid.msbb.
dim(genedata_org_rosmap)
#62757   639
genedata_org_rosmap=genedata_org_rosmap[bed$Geneid.msbb.rosmap,] #length(bed$Geneid.msbb.rosmap) is 55437 - this is the new row
rownames(genedata_org_rosmap)=bed$ensembl
genedata_org_msbb=genedata_org_msbb[bed$Geneid.msbb.rosmap,]
rownames(genedata_org_msbb)=bed$ensembl

genedata_org_lbp=genedata_org_lbp[bed$Geneid.lbpBrain.lbpBlood,]
rownames(genedata_org_lbp)=bed$ensembl

genedata_blood_lbp=genedata_blood_lbp[bed$Geneid.lbpBrain.lbpBlood,]
rownames(genedata_blood_lbp)=bed$ensembl

geneexp_all=geneexp_all[bed$Geneid.gtex,]
rownames(geneexp_all)=bed$ensembl

assert_that(identical(rownames(genedata_org_rosmap),rownames(genedata_org_msbb))) #TRUE
assert_that(identical(rownames(genedata_org_rosmap),rownames(genedata_org_lbp))) #TRUE
assert_that(identical(rownames(genedata_org_rosmap),rownames(genedata_blood_lbp))) #TRUE
assert_that(identical(rownames(genedata_org_rosmap),rownames(geneexp_all))) #TRUE

##for GTEx
sampleMetaDic=as.data.frame(readxl::read_excel("GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx"))
subjectMetaDic=as.data.frame(readxl::read_excel("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx"))

sampleMeta=fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",data.table=FALSE)
subjectMeta=fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",data.table=FALSE)

sampleMeta$SUBJID=unlist(lapply(strsplit(sampleMeta$SAMPID,"-",fixed=TRUE),function(x){paste0(x[1],"-",x[2])}))

intersect(sampleMeta$SUBJID,subjectMeta$SUBJID)
setdiff(subjectMeta$SUBJID,sampleMeta$SUBJID) #finds the values in subjectMeta$SUBJID that are not in sampleMeta$SUBJID
setdiff(sampleMeta$SUBJID,subjectMeta$SUBJID)

sampleMeta=merge(sampleMeta,subjectMeta,by="SUBJID")
allSamples=intersect(sampleMeta$SAMPID,colnames(geneexp_all))
geneexp_all=geneexp_all[,allSamples]
sampleMeta=sampleMeta[match(colnames(geneexp_all),sampleMeta$SAMPID),]
assert_that(identical(colnames(geneexp_all),sampleMeta$SAMPID)) #TRUE

#tissues in GTEx are categorized in two different ways SMTS (global) and SMTSD (more detailed)
table(sampleMeta$SMTS)
table(sampleMeta$SMTSD)

cols=c("SUBJID", "SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH", "SMTSPAX", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", "SMAFRZE", "SMGTC", "SME2MPRT", "SMCHMPRS", "SMNTRART", "SMNUMGPS", "SMMAPRT", "SMEXNCRT", "SM550NRM", "SMGNSDTC", "SMUNMPRT", "SM350NRM", "SMRDLGTH", "SMMNCPB", "SME1MMRT", "SMSFLGTH", "SMESTLBS", "SMMPPD", "SMNTERRT", "SMRRNANM", "SMRDTTL", "SMVQCFL", "SMMNCV", "SMTRSCPT", "SMMPPDPR", "SMCGLGTH", "SMGAPPCT", "SMUNPDRD", "SMNTRNRT", "SMMPUNRT", "SMEXPEFF", "SMMPPDUN", "SME2MMRT", "SME2ANTI", "SMALTALG", "SME2SNSE", "SMMFLGTH", "SME1ANTI", "SMSPLTRD", "SMBSMMRT", "SME1SNSE", "SME1PCTS", "SMRRNART", "SME1MPRT", "SMNUM5CD", "SMDPMPRT", "SME2PCTS", "SEX", "AGE", "DTHHRDY")

#design your own matrix, make a new matrix with specified numbers of rows and columns
meta_rosmap=data.frame(matrix(data=NA,nrow=ncol(genedata_org_rosmap),ncol=length(cols)))
colnames(meta_rosmap)=cols #add in the column names
#fill in the values for the columns
meta_rosmap$SUBJID=colnames(genedata_org_rosmap) 
meta_rosmap$SAMPID=colnames(genedata_org_rosmap)
meta_rosmap$SMTS="PFC_rosmap"
meta_rosmap$SMTSD="PFC_rosmap"
meta_rosmap$readLength=101

meta_msbb=data.frame(matrix(data=NA,nrow=ncol(genedata_org_msbb),ncol=length(cols)))
colnames(meta_msbb)=cols
meta_msbb$SUBJID=colnames(genedata_org_msbb)
meta_msbb$SAMPID=colnames(genedata_org_msbb)
meta_msbb$SMTS="PFC_msbb"
meta_msbb$SMTSD="PFC_msbb"
meta_msbb$readLength=100

meta_lbp=data.frame(matrix(data=NA,nrow=ncol(genedata_org_lbp),ncol=length(cols)))
colnames(meta_lbp)=cols
meta_lbp$SUBJID=colnames(genedata_org_lbp)
meta_lbp$SAMPID=colnames(genedata_org_lbp)
meta_lbp$SMTS="PFC_lbp"
meta_lbp$SMTSD="PFC_lbp"
meta_lbp$readLength=100

meta_blood_lbp=data.frame(matrix(data=NA,nrow=ncol(genedata_blood_lbp),ncol=length(cols)))
colnames(meta_blood_lbp)=cols
meta_blood_lbp$SUBJID=colnames(genedata_blood_lbp)
meta_blood_lbp$SAMPID=colnames(genedata_blood_lbp)
meta_blood_lbp$SMTS="blood_lbp"
meta_blood_lbp$SMTSD="blood_lbp"
meta_blood_lbp$readLength=100

table(sampleMeta$SMTSD)
sampleMeta$readLength=76

#combine all the metadata tgt
meta_all=rbind(sampleMeta,meta_rosmap,meta_msbb,meta_lbp,meta_blood_lbp)

#combine all the raw counts tgt
counts_all=cbind(genedata_org_rosmap, genedata_org_msbb, genedata_org_lbp, genedata_blood_lbp, geneexp_all)

#reordering the meta_all dataframe so that its rows match the order of colnames(counts_all), based on the SAMPID column
meta_all=meta_all[match(colnames(counts_all), meta_all$SAMPID),]
assert_that(identical(colnames(counts_all),meta_all$SAMPID)) #TRUE
assert_that(identical(rownames(counts_all),bed$ensembl)) #TRUE

##FOR ALL -- SKIPPPPPPPPPPPP TO SUBSET!
tpms=counts_to_tpm(counts_all, featureLength=bed$Length, meanFragmentLength=meta_all$readLength)
#performs quantile normalization (or specifically, normal quantile transformation) on each sample (column) in the tpms matrix
transformed=apply(tpms,2,function(sample){qqnorm(sample, plot.it = FALSE)$x})
rownames(transformed)=rownames(tpms)
assert_that(identical(colnames(transformed),meta_all$SAMPID))

meta_all$SMTS
meta_all$SMTSD 

##subsetting metadata for only relevant features
meta_subset <- meta_all[meta_all$SMTS %in% c("Blood", "blood_lbp", "Brain", "PFC_lbp", "PFC_msbb", "PFC_rosmap"), ]
table(meta_subset$SMTS) #our sample size
# Blood  blood_lbp      Brain    PFC_lbp   PFC_msbb PFC_rosmap 
#   929        243       2642        279        325        639

##subsetting counts for only relevant features
counts_subset <- counts_all[, colnames(counts_all) %in% meta_subset$SAMPID]

# Verify that the column names now match the subset
all(colnames(counts_subset) %in% meta_subset$SAMPID)  #TRUE

tpms_subset=counts_to_tpm(counts_subset, featureLength=bed$Length, meanFragmentLength=meta_subset$readLength)

transformed_subset=apply(tpms_subset,2,function(sample){qqnorm(sample, plot.it = FALSE)$x})
rownames(transformed_subset)=rownames(tpms_subset)
assert_that(identical(colnames(transformed_subset),meta_subset$SAMPID))

library(parallelDist)
distObj=parDist(as.matrix(t(transformed_subset)), method = "euclidean", diag = FALSE, upper = FALSE, threads = 30) #Computing the Distance Matrix
saveRDS(distObj,file="/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/distObj_gtex_lbp_rosmap_msbb.RDS")

mds=cmdscale(distObj)
colnames(mds)=c("MDS1","MDS2")
saveRDS(mds,file="/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/mds_gtex_lbp_rosmap_msbb.RDS")

##start 12:47pm on 3/13/2025
#end ~1pm

# if(!file.exists("/sc/arion/projects/mscic1/results/Noam/tmp/distObj_gtex_lbp_rosmap_msbb.RDS")){
# 	library(parallelDist)
# 	distObj=parDist(as.matrix(t(transformed)), method = "euclidean", diag = FALSE, upper = FALSE, threads = 30) #Computing the Distance Matrix
# 	saveRDS(distObj,file="/sc/arion/projects/mscic1/results/Noam/tmp/distObj_gtex_lbp_rosmap_msbb.RDS")
# }else{
# 	distObj=readRDS("/sc/arion/projects/mscic1/results/Noam/tmp/distObj_gtex_lbp_rosmap_msbb.RDS")
# }

# if(!file.exists("/sc/arion/projects/mscic1/results/Noam/tmp/mds_gtex_lbp_rosmap_msbb.RDS")){
# 	mds=cmdscale(distObj)
# 	colnames(mds)=c("MDS1","MDS2")
# 	saveRDS(mds,file="/sc/arion/projects/mscic1/results/Noam/tmp/mds_gtex_lbp_rosmap_msbb.RDS")
# }else{
# 	mds=readRDS("/sc/arion/projects/mscic1/results/Noam/tmp/mds_gtex_lbp_rosmap_msbb.RDS")
# }

library(ggplot2)

#unique(meta_all$SMTS[meta_all$SMTS %in% c("Brain","PFC_lbp","PFC_msbb","PFC_rosmap","Heart","Lung")])

mds2=merge(mds,meta_subset,by.x="row.names",by.y="SAMPID")
#mds2=merge(mds,meta_all,by.x="row.names",by.y="SAMPID")

table(mds2$SMTS)
     # Blood  blood_lbp      Brain    PFC_lbp   PFC_msbb PFC_rosmap 
     #   929        243       2642        279        325        639 

shuffled=sample(nrow(mds2))

mds2$SMTS <- factor(mds2$SMTS, 
    levels = c("Blood", "blood_lbp", "Brain", "PFC_lbp", "PFC_msbb", "PFC_rosmap"), 
    labels = c("GTEx_blood", "LBP_blood", "GTEx_brain", "LBP_PFC", "MSBB_PFC", "ROSMAP_PFC"))

g1=ggplot(mds2[shuffled,],aes(x=MDS1, y=MDS2, color=SMTS, label = SMTS)) + geom_text(size = 3) + scale_color_discrete(name = "Tissue") + theme_bw()
g2=ggplot(mds2,aes(x=MDS1, y=MDS2, color=SMTSD)) + geom_point() + scale_color_discrete(name = "Tissue") + theme_bw()
# g3=ggplot(mds2[mds2$SMTS %in% c("Brain","PFC_lbp","PFC_msbb","PFC_rosmap","Heart","Lung"),],aes(x=MDS1, y=MDS2, color=SMTSD)) + geom_point() + scale_color_discrete(name = "Tissue")

ggsave("/hpc/users/hoangd02/www/plots/lbp_mds.png", plot = g1, width = 8, height=8) 
ggsave("/hpc/users/hoangd02/www/plots/lbp_mds2.png", plot = g2) 


# all genes with no filters (unused because weird)
# pdf("~/hpc/users/hoangd02/www/plots/lbp_mds.pdf")
# show(g1)
# show(g2)
# # show(g3)
# dev.off()












