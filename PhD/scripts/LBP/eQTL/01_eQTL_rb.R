#Calculations of rb for eqtl 
### GOALS: 
# 1. Use BrainMeta eQTLs as the "discovery" set (large, well-powered)
# 2. Extract effect sizes for the same SNP-gene pairs from:
# 3. Do LIV-only analysis → rb(BrainMeta, LIV)
# 4. Do PM-only analysis → rb(BrainMeta, PM)
# 5. Test if rb(BrainMeta, LIV) ≠ rb(BrainMeta, PM)
#### Alex's suggestion: subset to the snps from brainmeta before running qtltools since we dont care about any other snps

###which nsv to use? look at manuscript or check all nsv to see how rb changes

##to test difference in rb, calculate z = (r1-r2) / square root ( (se1)^2 + (se2)^2 ) and p-values ~
##null hypothesis: rb(BrainMeta, LIV) = rb(BrainMeta, PM)

############################################################################################
################# for a more streamlined workflow, scroll until it says STREAMLINED WORKFLOW
############################################################################################
library(dplyr)
library(data.table)

brainmeta <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/BrainMeta_cis_eqtl_summary/BrainMeta_all_chr.txt")
dim(brainmeta) #283,666,305        14
length(unique(brainmeta$Probe)) #16743 unique genes
brainmeta$snp_key <- paste0("chr", brainmeta$Chr, ":", brainmeta$BP)
brainmeta$Probe <- sub("\\..*$", "", brainmeta$Probe) #getting rid of decimal and numbers after

head(brainmeta,3)
#            SNP Chr    BP A1                                                 A2
# 1  rs766767872   1 10390  C       CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA
# 2 rs1206847539   1 10397  * CCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAA
# 3 rs1481664201   1 10412  *                                                  C
#        Freq           Probe Probe_Chr Probe_bp   Gene Orientation        b
# 1 0.0121212 ENSG00000227232         1    21987 WASH7P           - 0.142414
# 2 0.0138889 ENSG00000227232         1    21987 WASH7P           - 0.114071
# 3 0.1996530 ENSG00000227232         1    21987 WASH7P           - 0.146967
#         SE        p    snp_key
# 1 0.226814 0.530076 chr1:10390
# 2 0.217267 0.599566 chr1:10397
# 3 0.122866 0.231636 chr1:10412

nrow(brainmeta) #283,666,305

uniqueN(brainmeta[,.(SNP)]) #11,580,188 unique SNP only

uniqueN(brainmeta[,.(SNP,Chr)]) #11,580,188 unique SNP + Chr

uniqueN(brainmeta[,.(SNP,Chr,A1)]) #11,580,188 unique SNP + Chr + A1 

uniqueN(brainmeta[,.(SNP,Chr,A1,A2)]) #11,580,188 unique SNP + Chr + A1 + A2

uniqueN(brainmeta[,.(SNP,Chr,A1,A2, Probe)]) #283,666,305  unique 

unique(length(brainmeta$Probe)) #16743

uniqueN(brainmeta[,.(SNP)]) #nrow(brainmeta[,.N,list(SNP)])

table(brainmeta$Chr)
#        1        2        3        4        5        6        7        8 
# 25638104 18508554 14655251 12546663 12879989 20967906 14970831 12031584 
#        9       10       11       12       13       14       15       16 
# 11881715 12488454 14536907 14475582  6161080  8886834  9835590 12673618 
#       17       18       19       20       21       22 
# 14773253  5413319 21403285  6878896  3315988  8742902

pm_28sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/PM_nominal_nsv28.txt", header = FALSE, sep = " ", data.table=FALSE)
dim(pm_28sv) #102960666        17
head(pm_28sv)
table(pm_28sv$V2) #this is chromosome info, making sure that all the chromosomes are there
#       1       2       3       4       5       6       7       8       9      10 
# 9602159 7173246 5834485 4560757 5311693 8193559 5290814 4551768 4299624 4763210 
#      11      12      13      14      15      16      17      18      19      20 
# 5940482 5714165 2498660 3786223 3767532 5033949 5692506 2081185 1469069 3149726 
#      21      22 
# 1332875 2912979 

dim(liv_28sv)#102960666        14
table(liv_28sv$V2) 
#       1       2       3       4       5       6       7       8       9      10 
# 9602159 7173246 5834485 4560757 5311693 8193559 5290814 4551768 4299624 4763210 
#      11      12      13      14      15      16      17      18      19      20 
# 5940482 5714165 2498660 3786223 3767532 5033949 5692506 2081185 1469069 3149726 
#      21      22 
# 1332875 2912979 

liv_9sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/LIV_nominal_nsv9.txt", header = FALSE, sep = " ", data.table=FALSE)
dim(liv_9sv) #102960666        14
head(liv_9sv)

liv_28sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/LIV_nominal_nsv28.txt", header = FALSE, sep = " ", data.table=FALSE)
liv = liv_28sv #liv_9sv
pm = pm_28sv
# ---- Assign column names ----
colnames(liv) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                  "strand", "n_variants", "distance", "variant_id",
                  "var_chr", "var_start", "var_end",
                  "pval", "beta", "top_flag")        
liv$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", liv$variant_id)) # Extract variant position
liv$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", liv$variant_id)) # Extract chromosome
liv$snp_key <- paste0("chr", liv$variant_chr, ":", liv$variant_pos) # Create SNP key

colnames(pm) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                  "strand", "n_variants", "distance", "variant_id",
                  "var_chr", "var_start", "var_end",
                  "pval", "beta", "top_flag")   

pm$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", pm$variant_id)) # Extract variant position
pm$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", pm$variant_id)) # Extract chromosome
pm$snp_key <- paste0("chr", pm$variant_chr, ":", pm$variant_pos) # Create SNP key

nrow(brainmeta) #283,666,305

dim(brainmeta[!snp_key %in% pm$snp_key]) #282298934        15

#283,666,305 - 282,298,934
##hmmm so only 1,367,371 rows match in brainmeta and pm

########################check if liftover hg38 is incorporated into the new brainmeta !!!
##sanity check 1 snp --> hg19 position vs hg38 position !

#Pick a random SNP from BrainMeta:

brainmeta[1, .(SNP, Chr, BP, snp_key)]
#            SNP   Chr    BP    snp_key
#         <char> <int> <int>     <char>
# 1: rs766767872     1 10390 chr1:10390

#And compare to your hg38 BED (the one you used in -R brainmeta_positions_hg38.bed). 
#If the BED position for the same SNP is different from BP, then that’s confirmation: 
#brainmeta in R is hg19, BED/PM are hg38.


hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38_chr.bed", header = FALSE)
head(hg38,3)
#        V1    V2    V3
#    <char> <int> <int>
# 1:   chr1 10389 10390
# 2:   chr1 10396 10397
# 3:   chr1 10411 10412

setnames(hg38, c("chr", "start", "end"))

hg38[, snp_key_hg38 := paste0(chr, ":", end)]  # or V2+1; here V3 == V2+1
hg38[snp_key_hg38 == "chr1:10390", ]
#58:   chr1 10389 10390   chr1:10390
# 59:   chr1 10389 10390   chr1:10390
# 60:   chr1 10389 10390   chr1:10390

###this should be 1:1 mapping i think
save.image("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/debugging_20251130.RData")

load("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/debugging_20251130.RData")


######### Calculating SE #########
# compute z-statistic 
liv$z <- ifelse(
  liv$pval > 0 & liv$pval < 1,
  qnorm(1 - liv$pval / 2),
  ifelse(liv$pval == 0, 38, 1e-6)  # extreme p-values
)
# standard error
liv$se <- abs(liv$beta) / abs(liv$z)

# compute z-statistic 
pm$z <- ifelse(
  pm$pval > 0 & pm$pval < 1,
  qnorm(1 - pm$pval / 2),
  ifelse(pm$pval == 0, 38, 1e-6)  # extreme p-values
)
# standard error
pm$se <- abs(pm$beta) / abs(pm$z)

########### pm_28sv ###########
###############################
summary(pm$se)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.03578 0.08742 0.13492 0.18266 3.76277 
sum(is.na(pm$se)) #0
head(pm,3)
#           gene_id gene_chr gene_start gene_end strand n_variants distance
# 1 ENSG00000201033        3   27513872 27513872      +       3012  -999805
# 2 ENSG00000201033        3   27513872 27513872      +       3012  -999769
# 3 ENSG00000201033        3   27513872 27513872      +       3012  -999506
#              variant_id var_chr var_start  var_end     pval       beta top_flag
# 1 chr3:26514067[b38]A,G       3  26514067 26514067 0.700851 -0.0139834        0
# 2 chr3:26514103[b38]A,G       3  26514103 26514103 0.499979 -0.0244403        0
# 3 chr3:26514366[b38]A,G       3  26514366 26514366 0.703826 -0.0138348        0
#   variant_pos variant_chr       snp_key         z         se
# 1    26514067           3 chr3:26514067 0.3841720 0.03639880
# 2    26514103           3 chr3:26514103 0.6745228 0.03623347
# 3    26514366           3 chr3:26514366 0.3801609 0.03639196

########### liv_9sv ###########
###############################
summary(liv$se)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.04941 0.11373 0.16809 0.22447 4.18903

sum(is.na(liv$se)) #0, No NAs
head(liv,3)
#           gene_id gene_chr gene_start gene_end strand n_variants distance
# 1 ENSG00000201033        3   27513872 27513872      +       3012  -999805
# 2 ENSG00000201033        3   27513872 27513872      +       3012  -999769
# 3 ENSG00000201033        3   27513872 27513872      +       3012  -999506
#              variant_id var_chr var_start  var_end     pval       beta top_flag
# 1 chr3:26514067[b38]A,G       3  26514067 26514067 0.510722 -0.0189967        0
# 2 chr3:26514103[b38]A,G       3  26514103 26514103 0.556034 -0.0171315        0
# 3 chr3:26514366[b38]A,G       3  26514366 26514366 0.517350 -0.0186154        0
#   variant_pos variant_chr       snp_key         z         se
# 1    26514067           3 chr3:26514067 0.6577139 0.02888292
# 2    26514103           3 chr3:26514103 0.5887425 0.02909846
# 3    26514366           3 chr3:26514366 0.6474357 0.02875251

########### liv_28sv ###########
###############################
summary(liv$se)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.04272 0.10517 0.15777 0.21286 4.13157 
sum(is.na(liv$se)) #0
head(liv,3)
#           gene_id gene_chr gene_start gene_end strand n_variants distance
# 1 ENSG00000201033        3   27513872 27513872      +       3012  -999805
# 2 ENSG00000201033        3   27513872 27513872      +       3012  -999769
# 3 ENSG00000201033        3   27513872 27513872      +       3012  -999506
#              variant_id var_chr var_start  var_end     pval       beta top_flag
# 1 chr3:26514067[b38]A,G       3  26514067 26514067 0.489020 -0.0183572        0
# 2 chr3:26514103[b38]A,G       3  26514103 26514103 0.556933 -0.0157030        0
# 3 chr3:26514366[b38]A,G       3  26514366 26514366 0.520977 -0.0169544        0
#   variant_pos variant_chr       snp_key         z         se
# 1    26514067           3 chr3:26514067 0.6918684 0.02653279
# 2    26514103           3 chr3:26514103 0.5874031 0.02673292
# 3    26514366           3 chr3:26514366 0.6418401 0.02641530

############# rb calculations #############
###########################################
# step1: select the same number of eGenes from both datasets
#data1=fread(paste0(INDIRT,data1_name,"_and_",data2_name,".txt"),head=T,stringsAsFactors=F,data.table=F)
## Merge on SNP + gene (Probe == gene_id)
data1 <- merge(
  brainmeta,
  liv,
  by.x = c("snp_key", "Probe"),   # BrainMeta
  by.y = c("snp_key", "gene_id"), # LIV
  all = FALSE
)
dim(data1) #657096     32

## BrainMeta = dataset 1 (discovery)
## LIV       = dataset 2 (test)
data1$SE  <- data1$SE          # SE in BrainMeta  (already named SE)
data1$se2 <- data1$se          # SE in LIV
data1$p   <- data1$p           # p in BrainMeta   (already named p)
data1$p2  <- data1$pval        # p in LIV

index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true 
#SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
dim(data1)
#656937     34
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);

shared_index=intersect(index1,index2)
length(shared_index) #331 rows that are significant in both datasets
unique_index1=index1[!index1 %in% shared_index] #rows that are significant in dataset 1 but not dataset 2
length(unique_index1) #17030
unique_index2=index2[!index2 %in% shared_index] #rows that are significant in dataset 2 but not dataset 1
length(unique_index2) #237

##This block is making the number of unique eGenes from each dataset the same, by downsampling the larger set.
if(length(unique_index2)<length(unique_index1)){
    tmp_index=c(shared_index,unique_index2,sample(unique_index1,length(unique_index2)))
}else{
    tmp_index=c(shared_index,unique_index1,sample(unique_index2,length(unique_index1)))
}
#Take all shared eGenes, Take all unique eGenes from dataset1 or 2,
#Randomly sample from unique_index1 (dataset1) or 2 so you only keep as many as dataset2/1 has.
#Result: end up with shared eGenes and equal number of unique eGenes from dataset1 and dataset2

data2=data1[tmp_index,] 
dim(data2) #805  34, all eGenes that are significant in both datasets.

# step2: standarized beta and SE in SD unit
# z represent z statistics, p presents allele frequency, n represents sample size.\
calcu_std_b_se<-function(z,p,n){
    std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
    std_se=1/sqrt(2*p*(1-p)*(n+z^2))
    res<-data.frame(std_b_hat,std_se);
    return(res)
}

## ---- BrainMeta standardized effects (dataset 1) ----
# z = beta / SE
z1 <- data2$bm_beta / data2$SE      # <-- change 'bm_beta' to your BM beta column
p1 <- data2$bm_eaf                  # <-- change to your BM EAF/MAF column
n1 <- data2$bm_n                    # <-- change to your BM sample size column

std1 <- calcu_std_b_se(z = z1, p = p1, n = n1)

b1_std  <- std1$std_b_hat           # standardized beta (BrainMeta)
se1_std <- std1$std_se              # standardized SE   (BrainMeta)

## ---- LIV standardized effects (dataset 2) ----
z2 <- data2$liv_beta / data2$se2    # <-- change 'liv_beta' to your LIV beta column
p2 <- data2$liv_eaf                 # <-- change to your LIV EAF/MAF column
n2 <- data2$liv_n                   # <-- change to your LIV sample size column

std2 <- calcu_std_b_se(z = z2, p = p2, n = n2)

b2_std  <- std2$std_b_hat           # standardized beta (LIV)
se2_std <- std2$std_se              # standardized SE   (LIV)

# step3: calcualte rb
# b1 and se1 represent the estimate and SE for eQTLs across probes in one tissue, b2 and se2 represent the estimate and SE for eQTLs across probes in the other tissue
# theta = sample overlap * phnotypic correlation; theta could also be estimated from null SNPs; if two samples were independent, theta = 0.
# Please note that the effect allele of one SNP between two tissues should be the same.
calcu_cor_true<-function(b1,se1,b2,se2,theta){
    idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2));
    if(length(idx)>0){
        b1=b1[-idx];se1=se1[-idx]
        b2=b2[-idx];se2=se2[-idx]
        theta=theta[-idx]
    }

    var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
    var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
    if(var_b1<0){
      var_b1=var(b1,na.rm=T)
    }
    if(var_b2<0){
      var_b2=var(b2,na.rm=T)
    }
    cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(theta,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
    r=cov_b1_b2/sqrt(var_b1*var_b2)

    r_jack=c()
    n=length(b1)
    for(k in 1:n) {
       b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
       b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
       if(var_b1_jack<0){
        var_b1_jack=var(b1_jack,na.rm=T);
       }
       if(var_b2_jack<0){
        var_b2_jack=var(b2_jack,na.rm=T);
       }
       theta_jack=theta[-k];
       cov_e1_jack_e2_jack=mean(theta_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
       cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
       r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
       r_jack=c(r_jack,r_tmp)
    }
    r_mean=mean(r_jack,na.rm=T)
    idx=which(is.na(r_jack))
    if(length(idx)>0){
        se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
    }else{
    	se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
    }
    res<-cbind(r,se_r)
    return(res)
}

############ assuming that this is already in the right unit -- standard deviation 
rb_res <- calcu_cor_true(
  b1    = data2$b,
  se1   = data2$SE,
  b2    = data2$beta,
  se2   = data2$se,
  theta = 0
)

rb_res
#             r       se_r
# [1,] 0.358371 0.03055478

# r= estimated rb(BrainMeta, LIV)
# se_r = SE of rb

############# rb calculations for PM #############
##################################################
# step1: select the same number of eGenes from both datasets
#data1=fread(paste0(INDIRT,data1_name,"_and_",data2_name,".txt"),head=T,stringsAsFactors=F,data.table=F)
## Merge on SNP + gene (Probe == gene_id)
data1 <- merge(
  brainmeta,
  pm,
  by.x = c("snp_key", "Probe"),   # BrainMeta
  by.y = c("snp_key", "gene_id"), # PM
  all = FALSE
)
dim(data1) #657096     32

length(intersect(unique(brainmeta$Probe), unique(pm$gene_id))) #14393

length(unique(brainmeta$Probe)) #16743
length(unique(pm$gene_id)) #27569

length(unique(brainmeta$snp_key)) #11580188
length(unique(pm$snp_key)) #4,842,829



length(intersect(unique(brainmeta$snp_key), unique(pm$SNP))) #63,363


length(intersect(unique(brainmeta$SNP), unique(pm$snp_key))) #63,363



## BrainMeta = dataset 1 (discovery)
## LIV       = dataset 2 (test)
data1$SE  <- data1$SE          # SE in BrainMeta  (already named SE)
data1$se2 <- data1$se          # SE in LIV
data1$p   <- data1$p           # p in BrainMeta   (already named p)
data1$p2  <- data1$pval        # p in LIV

index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true 
#SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
dim(data1)
#656965     34
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);

shared_index=intersect(index1,index2)
length(shared_index) #460 rows that are significant in both datasets
unique_index1=index1[!index1 %in% shared_index] #rows that are significant in dataset 1 but not dataset 2
length(unique_index1) #16866
unique_index2=index2[!index2 %in% shared_index] #rows that are significant in dataset 2 but not dataset 1
length(unique_index2) #314

##This block is making the number of unique eGenes from each dataset the same, by downsampling the larger set.
if(length(unique_index2)<length(unique_index1)){
    tmp_index=c(shared_index,unique_index2,sample(unique_index1,length(unique_index2)))
}else{
    tmp_index=c(shared_index,unique_index1,sample(unique_index2,length(unique_index1)))
}
#Take all shared eGenes, Take all unique eGenes from dataset1 or 2,
#Randomly sample from unique_index1 (dataset1) or 2 so you only keep as many as dataset2/1 has.
#Result: end up with shared eGenes and equal number of unique eGenes from dataset1 and dataset2

data2=data1[tmp_index,] 
dim(data2) #1088   34, all eGenes that are significant in both datasets

rb_res_pm <- calcu_cor_true(
  b1    = data2$b,
  se1   = data2$SE,
  b2    = data2$beta,
  se2   = data2$se,
  theta = 0
)
rb_res_pm
#              r       se_r
# [1,] 0.4507408 0.02425202

########## a crude test to see if two rb's are different: 
r1  <- 0.358371      # rb(BM, LIV)
se1 <- 0.03055478

r2  <- 0.4507408     # rb(BM, PM)
se2 <- 0.02425202

# difference
diff <- r1 - r2

# SE of the difference
se_diff <- sqrt(se1^2 + se2^2)

# z-statistic
z <- diff / se_diff

# p-value (two-sided)
pval <- 2 * pnorm(-abs(z))

list(
  difference = diff,
  se_difference = se_diff,
  z = z,
  p_value = pval
)

# $difference
# [1] -0.0923698
# $se_difference
# [1] 0.03900968
# $z
# [1] -2.367869
# $p_value
# [1] 0.01789088

###############################################################
####################### STREAMLINED WORKFLOW ##################
###############################################################
library(dplyr)
library(data.table)

#### 0. Load and Prepare Brainmeta
brainmeta <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/BrainMeta_cis_eqtl_summary/BrainMeta_all_chr.txt", data.table=FALSE)
dim(brainmeta) #283,666,305        14
length(unique(brainmeta$Probe)) #16743 unique genes
brainmeta$snp_key <- paste0("chr", brainmeta$Chr, ":", brainmeta$BP)
brainmeta$Probe <- sub("\\..*$", "", brainmeta$Probe) #getting rid of decimal and numbers after

#### 1. Helper: preprocess a QTLtools nominal file (LIV / PM) ####
prep_qtl <- function(path) {
  df <- fread(path, header = FALSE, sep = " ", data.table = FALSE)
  colnames(df) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                    "strand", "n_variants", "distance", "variant_id",
                    "var_chr", "var_start", "var_end",
                    "pval", "beta", "top_flag")
  
  # parse variant info
  df$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", df$variant_id))
  df$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", df$variant_id))
  df$snp_key     <- paste0("chr", df$variant_chr, ":", df$variant_pos)
  
  # compute z and SE from p-value + beta
  df$z <- ifelse(
    df$pval > 0 & df$pval < 1,
    qnorm(1 - df$pval / 2),
    ifelse(df$pval == 0, 38, 1e-6)  # extreme p-values
  )
  df$se <- abs(df$beta) / abs(df$z)
  
  df
}

#### 2. Your rb function (unchanged) ####
calcu_cor_true <- function(b1, se1, b2, se2, theta) {
  idx = which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
  if (length(idx) > 0) {
    b1    = b1[-idx];    se1 = se1[-idx]
    b2    = b2[-idx];    se2 = se2[-idx]
    theta = theta[-idx]
  }

  var_b1 = var(b1, na.rm = TRUE) - mean(se1^2, na.rm = TRUE)
  var_b2 = var(b2, na.rm = TRUE) - mean(se2^2, na.rm = TRUE)

  if (var_b1 < 0) var_b1 = var(b1, na.rm = TRUE)
  if (var_b2 < 0) var_b2 = var(b2, na.rm = TRUE)

  cov_b1_b2 = cov(b1, b2, use = "complete.obs") -
    mean(theta, na.rm = TRUE) * sqrt(mean(se1^2, na.rm = TRUE) * mean(se2^2, na.rm = TRUE))

  r = cov_b1_b2 / sqrt(var_b1 * var_b2)

  # jackknife
  r_jack = c()
  n = length(b1)
  for (k in 1:n) {
    b1_jack = b1[-k]; se1_jack = se1[-k]
    b2_jack = b2[-k]; se2_jack = se2[-k]

    var_b1_jack = var(b1_jack, na.rm = TRUE) - mean(se1_jack^2, na.rm = TRUE)
    var_b2_jack = var(b2_jack, na.rm = TRUE) - mean(se2_jack^2, na.rm = TRUE)

    if (var_b1_jack < 0) var_b1_jack = var(b1_jack, na.rm = TRUE)
    if (var_b2_jack < 0) var_b2_jack = var(b2_jack, na.rm = TRUE)

    theta_jack = theta[-k]

    cov_e1_jack_e2_jack = mean(theta_jack, na.rm = TRUE) *
      sqrt(mean(se1_jack^2, na.rm = TRUE) * mean(se2_jack^2, na.rm = TRUE))

    cov_b1_b2_jack = cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack

    r_tmp = cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack = c(r_jack, r_tmp)
  }

  r_mean = mean(r_jack, na.rm = TRUE)
  idx_na = which(is.na(r_jack))

  if (length(idx_na) > 0) {
    se_r = sqrt((n - 1) / n * sum((r_jack[-idx_na] - r_mean)^2))
  } else {
    se_r = sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }

  cbind(r, se_r)
}

compute_rb_for_tissue <- function(brainmeta, tissue_df, p_thresh = 5e-8, theta = 0) {
  # 1) merge on SNP + gene
  dat <- merge(
    brainmeta,
    tissue_df,
    by.x = c("snp_key", "Probe"),
    by.y = c("snp_key", "gene_id"),
    all  = FALSE
  )
  
  # 2) require non-zero SEs in both
  dat <- dat[dat$SE != 0 & dat$se != 0, ]
  
  # 3) GW-significant indices
  idx1 <- which(dat$p    < p_thresh)   # BrainMeta p
  idx2 <- which(dat$pval < p_thresh)   # tissue p (LIV/PM)
  
  shared  <- intersect(idx1, idx2)
  uniq1   <- setdiff(idx1, shared)     # sig only in BrainMeta
  uniq2   <- setdiff(idx2, shared)     # sig only in tissue
  
  # 4) balance unique sets
  if (length(uniq1) == 0 | length(uniq2) == 0) {
    stop("No unique eGenes in one of the datasets after filtering.")
  }
  
  if (length(uniq2) < length(uniq1)) {
    keep_idx <- c(shared, uniq2, sample(uniq1, length(uniq2)))
  } else {
    keep_idx <- c(shared, uniq1, sample(uniq2, length(uniq1)))
  }
  
  dat_bal <- dat[keep_idx, ]
  
  # 5) compute rb (here using original b/SE, like you did)
  theta_vec <- rep(theta, nrow(dat_bal))
  
  rb_res <- calcu_cor_true(
    b1    = dat_bal$b,     # BrainMeta effect
    se1   = dat_bal$SE,    # BrainMeta SE
    b2    = dat_bal$beta,  # tissue effect
    se2   = dat_bal$se,    # tissue SE
    theta = theta_vec
  )
  
  list(
    rb      = as.numeric(rb_res[1]),
    se_rb   = as.numeric(rb_res[2]),
    n_pairs = nrow(dat_bal),
    n_shared = length(shared),
    n_uniq1  = length(uniq1),
    n_uniq2  = length(uniq2)
  )
}

#### Preprocess LIV & PM once ####
liv <- prep_qtl("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/LIV_nominal_nsv28.txt")
pm  <- prep_qtl("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/PM_nominal_nsv28.txt")

#### Compute rb(BrainMeta, LIV) and rb(BrainMeta, PM) ####
rb_liv <- compute_rb_for_tissue(brainmeta, liv, theta = 0)
rb_pm  <- compute_rb_for_tissue(brainmeta, pm,  theta = 0)

rb_liv
#rb_liv$rb        rb_liv$se_rb     rb_liv$n_pairs   rb_liv$n_shared  rb_liv$n_uniq1   rb_liv$n_uniq2
#0.3738545
rb_pm
#0.4581152 for  rb_pm$rb

###there are slight differences ~~ 

##previously for LIV
#             r       se_r
# [1,] 0.358371 0.03055478

##for PM 
#              r       se_r
# [1,] 0.4507408 0.02425202

####################THIS DOESNT HAVE TO BE THE TOP SNP-GENE PAIR RIGHT??? it's just whatever that has SNP-key

######### Function to compute rb for one tissue and one nsv
compute_rb_for_tissue <- function(brainmeta, tissue_df, p_thresh = 5e-8, theta = 0) {
  # 1) merge on SNP + gene
  dat <- merge(
    brainmeta,
    tissue_df,
    by.x = c("snp_key", "Probe"),
    by.y = c("snp_key", "gene_id"),
    all  = FALSE
  )
  
  # 2) require non-zero SEs in both
  dat <- dat[dat$SE != 0 & dat$se != 0, ]
  
  # 3) GW-significant indices
  idx1 <- which(dat$p    < p_thresh)   # BrainMeta p
  idx2 <- which(dat$pval < p_thresh)   # tissue p (LIV/PM)
  
  shared  <- intersect(idx1, idx2)
  uniq1   <- setdiff(idx1, shared)     # sig only in BrainMeta
  uniq2   <- setdiff(idx2, shared)     # sig only in tissue
  
  if (length(uniq1) == 0 | length(uniq2) == 0) {
    return(list(
      rb = NA_real_,
      se_rb = NA_real_,
      n_pairs = NA_integer_,
      n_shared = length(shared),
      n_uniq1 = length(uniq1),
      n_uniq2 = length(uniq2),
      note = "No unique eGenes in one dataset"
    ))
  }
  
  # 4) balance unique sets
  if (length(uniq2) < length(uniq1)) {
    keep_idx <- c(shared, uniq2, sample(uniq1, length(uniq2)))
  } else {
    keep_idx <- c(shared, uniq1, sample(uniq2, length(uniq1)))
  }
  
  dat_bal <- dat[keep_idx, ]
  
  # 5) compute rb (using original BrainMeta b/SE and tissue beta/se)
  theta_vec <- rep(theta, nrow(dat_bal))
  
  rb_res <- calcu_cor_true(
    b1    = dat_bal$b,     # BrainMeta effect
    se1   = dat_bal$SE,    # BrainMeta SE
    b2    = dat_bal$beta,  # tissue effect
    se2   = dat_bal$se,    # tissue SE
    theta = theta_vec
  )
  
  list(
    rb       = as.numeric(rb_res[1]),
    se_rb    = as.numeric(rb_res[2]),
    n_pairs  = nrow(dat_bal),
    n_shared = length(shared),
    n_uniq1  = length(uniq1),
    n_uniq2  = length(uniq2),
    note     = NA_character_
  )
}

#######Loop over nsv = 1:30 for LIV and PM
base_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA"

# nsv_vals <- 1:3
# nsv_vals <- 4:6
# nsv_vals <- 9:11
# nsv_vals <- 12:14
# nsv_vals <- 15:17
# nsv_vals <- 18:20
# nsv_vals <- 21:23
# nsv_vals <- 24:26
#nsv_vals <- 27:30

#nsv_vals <- 12:12
#nsv_vals <- 13:13
# nsv_vals <- 14:14

# nsv_vals <- 21:21
# nsv_vals <- 22:22
# nsv_vals <- 23:23

#nsv_vals <- 7:7
nsv_vals <- 8:8

results_list <- list()

for (nsv in nsv_vals) {
  message("Processing nsv = ", nsv)
  
  ## paths
  liv_path <- file.path(base_dir, sprintf("LIV_nominal_nsv%d.txt", nsv))
  pm_path  <- file.path(base_dir, sprintf("PM_nominal_nsv%d.txt",  nsv))
  
  ## read + preprocess LIV/PM
  liv_df <- prep_qtl(liv_path)
  pm_df  <- prep_qtl(pm_path)
  
  ## compute rb for each
  rb_liv <- compute_rb_for_tissue(brainmeta, liv_df, theta = 0)
  rb_pm  <- compute_rb_for_tissue(brainmeta, pm_df,  theta = 0)
  
  ## store in data.frame rows
  results_list[[length(results_list) + 1]] <- data.frame(
    tissue = "LIV",
    nsv    = nsv,
    rb     = rb_liv$rb,
    se_rb  = rb_liv$se_rb,
    n_pairs  = rb_liv$n_pairs,
    n_shared = rb_liv$n_shared,
    n_uniq1  = rb_liv$n_uniq1,
    n_uniq2  = rb_liv$n_uniq2,
    note     = rb_liv$note,
    stringsAsFactors = FALSE
  )
  
  results_list[[length(results_list) + 1]] <- data.frame(
    tissue = "PM",
    nsv    = nsv,
    rb     = rb_pm$rb,
    se_rb  = rb_pm$se_rb,
    n_pairs  = rb_pm$n_pairs,
    n_shared = rb_pm$n_shared,
    n_uniq1  = rb_pm$n_uniq1,
    n_uniq2  = rb_pm$n_uniq2,
    note     = rb_pm$note,
    stringsAsFactors = FALSE
  )
}

rb_summary <- bind_rows(results_list)
write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv8.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv23.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv21.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv14.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv12.csv",row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv9.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv12.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv15.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv18.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv21.csv",
#   row.names = FALSE)

# write.csv(rb_summary,
#   "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv24.csv", ##im in this screen section
#   row.names = FALSE)

write.csv(rb_summary,
  "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv27.csv",
  row.names = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)

nsv1 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1.csv")
nsv4 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4.csv")
nsv7 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7.csv")
nsv8 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv8.csv")
nsv9 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv9.csv")
nsv12 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv12.csv")
nsv13 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13.csv")
nsv14 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv14.csv")
nsv15 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv15.csv")
nsv18 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv18.csv")
nsv21 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv21.csv")
nsv22 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22.csv")
nsv23 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv23.csv")
nsv24 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv24.csv")
nsv27 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv27.csv")

rb_all <- rbind(nsv1, nsv4, nsv7, nsv8, nsv9, nsv12, nsv13, nsv14, nsv15, nsv18, nsv21, nsv22, nsv23, nsv24, nsv27)

head(rb_all)
  tissue nsv        rb      se_rb n_pairs n_shared n_uniq1 n_uniq2 note
1    LIV   1 0.2968506 0.03525594     578      230   17124     174   NA
2     PM   1 0.5164186 0.02952381     699      307   17043     196   NA
3    LIV   2 0.3111650 0.03345249     615      253   17104     181   NA
4     PM   2 0.4745134 0.02885562     787      375   16975     206   NA
5    LIV   3 0.3197780 0.03290860     617      261   17099     178   NA
6     PM   3 0.4874427 0.02729401     850      392   16957     229   NA

##is significant pvalue in our data needs to be considered 
##different set of SNP–gene pairs are being considered for LIV vs PM 
##using all brainmeta SNP–gene pairs

##would input be match n_uniq1 / should input for rb the same for LIV and PM 

# n_pairs is number of SNP–gene pairs used to compute rb
# n_shared is GW-significant in both BrainMeta and LIV
# n_uniq1 is unique BrainMeta eGenes (BM-sig, LIV-not)
# n_uniq2 is unique LIV eGenes (LIV-sig, BM-not)

##to compute n_pairs, BrainMeta has 17124 eQTLs and LIV has 174 eQTLs, 230 of which are shared in both BrainMeta and LIV
##keep all 230 shared genes AND keep all LIV-only signals AND Randomly sample BM-only signals so the count matches LIV’s count
## so 230 + 174 (LIV) + 174 (sampled Brainmeta to match LIV number of eqtl) = 578

table(rb_all$nsv)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 
# 27 28 29 30 
#  2  2  2  2

rb_tests <- rb_all %>%
  select(tissue, nsv, rb, se_rb) %>%
  # make LIV/PM columns
  pivot_wider(
    names_from  = tissue,
    values_from = c(rb, se_rb)
  ) %>%
  mutate(
    # difference in rb (LIV - PM)
    diff    = rb_LIV - rb_PM,
    # SE of difference (independent estimates)
    se_diff = sqrt(se_rb_LIV^2 + se_rb_PM^2),
    # z-statistic
    z       = diff / se_diff,
    # two-sided p-value
    p_two_sided = 2 * pnorm(-abs(z))
  )

rb_tests <- as.data.frame(rb_tests)
rb_tests

########## PLOTTING nsv by rb_pm - rb_liv
rb_tests <- rb_all %>%
  select(tissue, nsv, rb, se_rb) %>%
  pivot_wider(
    names_from  = tissue,
    values_from = c(rb, se_rb)
  ) %>%
  mutate(
    diff    = rb_PM - rb_LIV,  # PM minus LIV (easier to interpret direction)
    se_diff = sqrt(se_rb_LIV^2 + se_rb_PM^2),
    z       = diff / se_diff,
    p_value = 2 * pnorm(-abs(z)),
    sig = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  )

p <- ggplot(rb_tests, aes(x = nsv, y = diff)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_point(size = 3, color = "darkred") +
  
  # significance stars placed above each point
  geom_text(aes(label = sig), 
            vjust = -0.8, 
            size = 5) +
  scale_x_continuous(breaks = rb_tests$nsv) + 
  labs(
    title = "Difference in rb(BrainMeta, PM) and rb(BrainMeta, LIV) Across number of SVs",
    x = "Number of Surrogate Variables",
    y = "rb(BrainMeta, PM) minus rb(BrainMeta, LIV)",
    caption = "Stars indicate significance of difference in PM vs LIV; * = <0.05; ** = <0.01; *** = < 0.001; **** = <0.0001"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference.pdf", p, width = 14, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/eqtl_rb_difference.pdf

##plotting all 3 lines
rb_long <- rb_tests %>%
  select(nsv, rb_LIV, rb_PM, diff) %>%
  pivot_longer(
    cols = c(rb_LIV, rb_PM, diff),
    names_to = "metric",
    values_to = "value"
  )

p2 <- ggplot(rb_long, aes(x = nsv, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "rb_LIV" = "#1f77b4",   # blue
      "rb_PM"  = "#d62728",   # red
      "diff"   = "#2ca02c"    # green
    ),
    labels = c(
      "rb_LIV" = "rb(BrainMeta, LIV)",
      "rb_PM"  = "rb(BrainMeta, PM)",
      "diff"   = "Difference (PM minus LIV)"
    )
  ) +
  scale_x_continuous(breaks = rb_tests$nsv) +
  labs(
    title = "rb(BrainMeta, LIV) vs rb(BrainMeta, PM) vs Their Difference",
    x = "Number of Surrogate Variables",
    y = "rb value",
    color = "Metric"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference_with_3lines.pdf", p2, width = 14, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/eqtl_rb_difference_with_3lines.pdf

##look at the structure of SVs more, are they the same SVs for LIV vs PM

##rb(BrainMeta, PM) is higher than rb(BrainMeta, LIV) so
#PM replicates BrainMeta eQTLs better than LIV for that NSV

##With 1 surrogate variable included, PM eQTL effects are much more concordant with BrainMeta 
##than LIV eQTL effects (rb difference ~0.22, p=1.8×10⁻⁶).

##This supports the idea that postmortem brain tissue shares stronger genetic regulatory architecture 
##with BrainMeta’s meta-analyzed cross-brain dataset, compared to LIV blood samples.

###Across all surrogate variable settings (NSV=1–30), rb estimates for BrainMeta–PM were consistently higher 
##than those for BrainMeta–LIV. Differences were large (Δrb ≈ 0.14–0.22) and highly significant (p < 0.001), 
#indicating substantially stronger replication of BrainMeta eQTL effects in postmortem brain than in living brain.




liv_28sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/LIV_nominal_nsv28.txt", header = FALSE, sep = " ", data.table=FALSE)
pm_28sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/PM_nominal_nsv28.txt", header = FALSE, sep = " ", data.table=FALSE)

liv = liv_28sv 
pm = pm_28sv
######### Assign column names #########
colnames(liv) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                  "strand", "n_variants", "distance", "variant_id",
                  "var_chr", "var_start", "var_end",
                  "pval", "beta", "top_flag")        
liv$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", liv$variant_id)) # Extract variant position
liv$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", liv$variant_id)) # Extract chromosome
liv$snp_key <- paste0("chr", liv$variant_chr, ":", liv$variant_pos) # Create SNP key

colnames(pm) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                  "strand", "n_variants", "distance", "variant_id",
                  "var_chr", "var_start", "var_end",
                  "pval", "beta", "top_flag")   
pm$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", pm$variant_id)) # Extract variant position
pm$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", pm$variant_id)) # Extract chromosome
pm$snp_key <- paste0("chr", pm$variant_chr, ":", pm$variant_pos) # Create SNP key

######### Calculating SE #########
# compute z-statistic 
liv$z <- ifelse(
  liv$pval > 0 & liv$pval < 1,
  qnorm(1 - liv$pval / 2),
  ifelse(liv$pval == 0, 38, 1e-6)  # extreme p-values
)
# standard error
liv$se <- abs(liv$beta) / abs(liv$z)

# compute z-statistic 
pm$z <- ifelse(
  pm$pval > 0 & pm$pval < 1,
  qnorm(1 - pm$pval / 2),
  ifelse(pm$pval == 0, 38, 1e-6)  # extreme p-values
)
# standard error
pm$se <- abs(pm$beta) / abs(pm$z)

############# rb calculations #############
###########################################
# step1: select the same number of eGenes from both datasets
#data1=fread(paste0(INDIRT,data1_name,"_and_",data2_name,".txt"),head=T,stringsAsFactors=F,data.table=F)
## Merge on SNP + gene (Probe == gene_id)
data1 <- merge(
  brainmeta,
  liv,
  by.x = c("snp_key", "Probe"),   # BrainMeta
  by.y = c("snp_key", "gene_id"), # LIV
  all = FALSE
)

## BrainMeta = dataset 1 (discovery)
## LIV       = dataset 2 (test)
data1$SE  <- data1$SE          # SE in BrainMeta  (already named SE)
data1$se2 <- data1$se          # SE in LIV
data1$p   <- data1$p           # p in BrainMeta   (already named p)
data1$p2  <- data1$pval        # p in LIV

index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true; SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);
shared_index=intersect(index1,index2)
unique_index1=index1[!index1 %in% shared_index] #rows that are significant in dataset 1 but not dataset 2
unique_index2=index2[!index2 %in% shared_index] #rows that are significant in dataset 2 but not dataset 1

if(length(unique_index2)<length(unique_index1)){
    tmp_index=c(shared_index,unique_index2,sample(unique_index1,length(unique_index2)))
}else{
    tmp_index=c(shared_index,unique_index1,sample(unique_index2,length(unique_index1)))
}

data2=data1[tmp_index,] 

# step2: standarized beta and SE in SD unit
# z represent z statistics, p presents allele frequency, n represents sample size.\
calcu_std_b_se<-function(z,p,n){
    std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
    std_se=1/sqrt(2*p*(1-p)*(n+z^2))
    res<-data.frame(std_b_hat,std_se);
    return(res)
}

## ---- BrainMeta standardized effects (dataset 1) ----
# z = beta / SE
z1 <- data2$bm_beta / data2$SE      # <-- change 'bm_beta' to your BM beta column
p1 <- data2$bm_eaf                  # <-- change to your BM EAF/MAF column
n1 <- data2$bm_n                    # <-- change to your BM sample size column

std1 <- calcu_std_b_se(z = z1, p = p1, n = n1)

b1_std  <- std1$std_b_hat           # standardized beta (BrainMeta)
se1_std <- std1$std_se              # standardized SE   (BrainMeta)

## ---- LIV standardized effects (dataset 2) ----
z2 <- data2$liv_beta / data2$se2    # <-- change 'liv_beta' to your LIV beta column
p2 <- data2$liv_eaf                 # <-- change to your LIV EAF/MAF column
n2 <- data2$liv_n                   # <-- change to your LIV sample size column

std2 <- calcu_std_b_se(z = z2, p = p2, n = n2)

b2_std  <- std2$std_b_hat           # standardized beta (LIV)
se2_std <- std2$std_se              # standardized SE   (LIV)

# step3: calcualte rb
# b1 and se1 represent the estimate and SE for eQTLs across probes in one tissue, b2 and se2 represent the estimate and SE for eQTLs across probes in the other tissue
# theta = sample overlap * phnotypic correlation; theta could also be estimated from null SNPs; if two samples were independent, theta = 0.
# Please note that the effect allele of one SNP between two tissues should be the same.
calcu_cor_true<-function(b1,se1,b2,se2,theta){
    idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2));
    if(length(idx)>0){
        b1=b1[-idx];se1=se1[-idx]
        b2=b2[-idx];se2=se2[-idx]
        theta=theta[-idx]
    }

    var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
    var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
    if(var_b1<0){
      var_b1=var(b1,na.rm=T)
    }
    if(var_b2<0){
      var_b2=var(b2,na.rm=T)
    }
    cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(theta,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
    r=cov_b1_b2/sqrt(var_b1*var_b2)

    r_jack=c()
    n=length(b1)
    for(k in 1:n) {
       b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
       b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
       if(var_b1_jack<0){
        var_b1_jack=var(b1_jack,na.rm=T);
       }
       if(var_b2_jack<0){
        var_b2_jack=var(b2_jack,na.rm=T);
       }
       theta_jack=theta[-k];
       cov_e1_jack_e2_jack=mean(theta_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
       cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
       r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
       r_jack=c(r_jack,r_tmp)
    }
    r_mean=mean(r_jack,na.rm=T)
    idx=which(is.na(r_jack))
    if(length(idx)>0){
        se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
    }else{
    	se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
    }
    res<-cbind(r,se_r)
    return(res)
}

rb_res_liv <- calcu_cor_true(
  b1    = data2$b,
  se1   = data2$SE,
  b2    = data2$beta,
  se2   = data2$se,
  theta = 0
)

rb_res_liv 

############ FOR PM
data1 <- merge(
  brainmeta,
  pm,
  by.x = c("snp_key", "Probe"),   # BrainMeta
  by.y = c("snp_key", "gene_id"), # PM
  all = FALSE
)
dim(data1) #657096     32

## BrainMeta = dataset 1 (discovery)
## LIV       = dataset 2 (test)
data1$SE  <- data1$SE          # SE in BrainMeta  (already named SE)
data1$se2 <- data1$se          # SE in LIV
data1$p   <- data1$p           # p in BrainMeta   (already named p)
data1$p2  <- data1$pval        # p in LIV

index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true 
#SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
dim(data1)
#656965     34
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);

shared_index=intersect(index1,index2)
length(shared_index) #460 rows that are significant in both datasets
unique_index1=index1[!index1 %in% shared_index] #rows that are significant in dataset 1 but not dataset 2
length(unique_index1) #16866
unique_index2=index2[!index2 %in% shared_index] #rows that are significant in dataset 2 but not dataset 1
length(unique_index2) #314

##This block is making the number of unique eGenes from each dataset the same, by downsampling the larger set.
if(length(unique_index2)<length(unique_index1)){
    tmp_index=c(shared_index,unique_index2,sample(unique_index1,length(unique_index2)))
}else{
    tmp_index=c(shared_index,unique_index1,sample(unique_index2,length(unique_index1)))
}
#Take all shared eGenes, Take all unique eGenes from dataset1 or 2,
#Randomly sample from unique_index1 (dataset1) or 2 so you only keep as many as dataset2/1 has.
#Result: end up with shared eGenes and equal number of unique eGenes from dataset1 and dataset2

data2=data1[tmp_index,] 
dim(data2) #1088   34, all eGenes that are significant in both datasets

rb_res_pm <- calcu_cor_true(
  b1    = data2$b,
  se1   = data2$SE,
  b2    = data2$beta,
  se2   = data2$se,
  theta = 0
)
rb_res_pm