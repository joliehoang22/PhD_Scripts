###############################################################
## fixing coordinates --> liftover 
## PLAN
#BrainMeta (hg19, with SNP IDs) → hg19 BED with SNP IDs → liftOver 
#→ hg38 BED with SNP IDs → merge back to BrainMeta.

library(dplyr)
library(data.table)

#1️. Rebuild the hg19 BED from UNIQUE SNPs and include rsID
brainmeta <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/BrainMeta_cis_eqtl_summary/BrainMeta_all_chr.txt", data.table=FALSE)
dim(brainmeta) #283,666,305        14
#One row per SNP–gene pair
#Same SNP appears many times, once per gene tested within cis-window

length(unique(brainmeta$SNP)) #11,580,188 unique SNPs
uniqueN(brainmeta[,.(SNP,Chr)]) #11,580,188 unique SNP + Chr

##checking if for a given SNP, is BP always the same?
#This checks whether any SNP appears at more than one BP:
brainmeta[, .(n_pos = uniqueN(BP)), by = SNP][n_pos > 1]
#Empty data.table (0 rows and 2 cols): SNP,n_pos -- If this returns 0 rows, then every SNP has one unique BP

##pick one SNP and inspect all rows manually
brainmeta[SNP == "rs766767872", .(SNP, Chr, BP, Probe)] #60 instances, so 60 SNPs per gene

### reasons for the same SNP for MANY genes is because you're testing for each gene G,
# you gather all SNPs within 1 Mb of its TSS and test their associations
##A single SNP can be tested with hundreds of genes if those genes are all within its ±1 Mb range.

unique(length(brainmeta$Probe)) #16743 genes


uniqueN(brainmeta[,.(SNP,Chr,A1,A2, Probe)]) #283,666,305  unique 

# One row per SNP coordinate (not per SNP–gene pair)
bm_hg19_bed <- unique(brainmeta[, .(
  V1 = paste0("chr", Chr),  # chrom
  V2 = BP - 1L,             # 0-based start
  V3 = BP,                  # 1-based end
  V4 = SNP                  # SNP ID
)])

dim(bm_hg19_bed) #11,580,188        4

head(bm_hg19_bed)
#        V1    V2    V3           V4
#    <char> <int> <int>       <char>
# 1:   chr1 10389 10390  rs766767872
# 2:   chr1 10396 10397 rs1206847539

fwrite(bm_hg19_bed, "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg19_with_SNPID.bed",
       sep = "\t", col.names = FALSE)

##check 
test <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg19_with_SNPID.bed",
                  header = FALSE, sep = "\t")

###########################
##### rerun liftover ######
###########################
cd /sc/arion/projects/mscic1/results/jolie/LBP/eqtl

ls liftOver
### dont have to rerun this
# Download the liftOver binary
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
# chmod +x liftOver

# # Download the chain file
# wget -O hg19ToHg38.over.chain.gz \
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

chmod +x liftOver

./liftOver \
  brainmeta_hg19_with_SNPID.bed \
  hg19ToHg38.over.chain.gz \
  brainmeta_hg38_with_SNPID.bed \
  brainmeta_unmapped_with_SNPID.bed

##checking results 
wc -l brainmeta_hg19_with_SNPID.bed #11,580,188
wc -l brainmeta_hg38_with_SNPID.bed #11,570,881
wc -l brainmeta_unmapped_with_SNPID.bed #18,614 unmapped

hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38_with_SNPID.bed",
                  header = FALSE, sep = "\t")
setnames(hg38, c("chr38", "start38", "end38", "SNP"))

hg38[, snp_key_hg38 := paste0(chr38, ":", start38 + 1)]

head(hg38)
#     chr38 start38 end38          SNP snp_key_hg38
#    <char>   <int> <int>       <char>       <char>
# 1:   chr1   10389 10390  rs766767872   chr1:10390
# 2:   chr1   10396 10397 rs1206847539   chr1:10397

#snp_key_hg38 is combo of  chr38 and end38

####################
#################### MERGE HG38 INTO BRAINMETA LOL
####################
brainmeta_hg38 <- merge(
  brainmeta,
  hg38[, .(SNP, snp_key_hg38)],
  by = "SNP",
  all.x = FALSE,
  all.y = FALSE
)
dim(brainmeta_hg38) #283,309,353        15

head(brainmeta_hg38)
# Key: <SNP>
#             SNP   Chr        BP     A1     A2      Freq              Probe
#          <char> <int>     <int> <char> <char>     <num>             <char>
# 1: 10:134319732    10 134319732      A      G 0.0338542  ENSG00000176769.9
# 2: 10:134319732    10 134319732      A      G 0.0338542  ENSG00000189275.4
# 3: 10:134319732    10 134319732      A      G 0.0338542  ENSG00000286295.1
# 4: 10:134319732    10 134319732      A      G 0.0338542 ENSG00000175470.14
# 5: 10:134319732    10 134319732      A      G 0.0338542  ENSG00000277959.1
# 6: 10:134319732    10 134319732      A      G 0.0338542 ENSG00000176171.11
#    Probe_Chr  Probe_bp       Gene Orientation          b       SE        p
#        <int>     <int>     <char>      <char>      <num>    <num>    <num>
# 1:        10 133000319    TCERG1L           -  0.0794240 0.140376 0.571534
# 2:        10 133613634  LINC01164           -  0.0434291 0.132064 0.742270
# 3:        10 133618168 AL450307.1           + -0.0272556 0.140883 0.846596
# 4:        10 133760643    PPP2R2D           + -0.0880263 0.136674 0.519537
# 5:        10 133784872 AL162274.2           + -0.0238172 0.135634 0.860609
# 6:        10 133787738      BNIP3           - -0.1515900 0.135917 0.264719
#       snp_key_hg38
#             <char>
# 1: chr10:132506228
# 2: chr10:132506228
# 3: chr10:132506228
# 4: chr10:132506228
# 5: chr10:132506228
# 6: chr10:132506228

#so SNP contains info of hg19, 10:134319732, and snp_key_hg38 contains info of hg38.

##quick sanity check - does 10:134319732 (hg19) match with chr10:132506228 (hg38)??
test_snp <- "10:134319732"

hg38[ SNP == test_snp, .(SNP, chr38, start38, end38) ]
#             SNP  chr38   start38     end38
#          <char> <char>     <int>     <int>
# 1: 10:134319732  chr10 132506227 132506228 ###yes they do match

pm_28sv <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA/PM_nominal_nsv28.txt", header = FALSE, sep = " ", data.table=FALSE)
dim(pm_28sv) #102960666        17
pm = pm_28sv

colnames(pm) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                  "strand", "n_variants", "distance", "variant_id",
                  "var_chr", "var_start", "var_end",
                  "pval", "beta", "top_flag")   

pm$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", pm$variant_id)) # Extract variant position
pm$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", pm$variant_id)) # Extract chromosome
pm$snp_key <- paste0("chr", pm$variant_chr, ":", pm$variant_pos) # Create SNP key

brainmeta_hg38$snp_key <- brainmeta_hg38$snp_key_hg38
brainmeta_hg38$Probe <- sub("\\..*$", "", brainmeta_hg38$Probe)

dim(brainmeta_hg38[!snp_key %in% pm$snp_key])
#old: 282298934        15
#new: 162836502        16

##overlap: 282298934 -  119462432

######### Calculating SE #########
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

# Ensure both are data.tables
setDT(pm)
setDT(brainmeta_hg38)

data_pm <- merge(
  brainmeta_hg38[, .(snp_key,
                     Probe,
                     b_bm = b,
                     se_bm = SE,
                     p_bm = p)],
  pm[, .(snp_key,
         gene_id,
         beta_pm = beta,
         se_pm = se,
         p_pm = pval)],
  by.x = c("snp_key", "Probe"),
  by.y = c("snp_key", "gene_id"),
  all = FALSE
)

dim(data_pm) #55,420,762        8 --- good, this is A LOT better than the 657,096 genes in 01_eQTL_rb.R
head(data_pm,3)
# Key: <snp_key, Probe>
#            snp_key           Probe       b_bm     se_bm        p_bm    beta_pm
#             <char>          <char>      <num>     <num>       <num>      <num>
# 1: chr10:100000235 ENSG00000014919 -0.0553866 0.0336881 1.00156e-01 -0.0108094
# 2: chr10:100000235 ENSG00000075826 -0.0103621 0.0405346 7.98232e-01 -0.0485158
# 3: chr10:100000235 ENSG00000095485  0.3977350 0.0414170 7.75028e-22  0.0600147
#         se_pm       p_pm
#         <num>      <num>
# 1: 0.01518065 0.47643300
# 2: 0.02548150 0.05691520
# 3: 0.02000490 0.00269979

data1 <- data_pm
## BrainMeta = dataset 1 (discovery)
## LIV       = dataset 2 (test)
data1$SE  <- data1$se_bm       # SE in BrainMeta  (already named SE)
data1$se2 <- data1$se_pm       # SE in PM
data1$p   <- data1$p_bm           # p in BrainMeta   (already named p)
data1$p2  <- data1$p_pm        # p in PM

index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true 
#SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
dim(data1)
#55400460       12
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);

shared_index=intersect(index1,index2)
length(shared_index) #50603 rows that are significant in both datasets
unique_index1=index1[!index1 %in% shared_index] #rows that are significant in dataset 1 but not dataset 2
length(unique_index1) #1443241
unique_index2=index2[!index2 %in% shared_index] #rows that are significant in dataset 2 but not dataset 1
length(unique_index2) #3359

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
dim(data2) #57321    12, all eGenes that are significant in both datasets.

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

rb_res <- calcu_cor_true(
  b1    = data2$b_bm,
  se1   = data2$se_bm,
  b2    = data2$beta_pm,
  se2   = data2$se_pm,
  theta = 0
)

rb_res ### A LOT HIGHER THAN THE ORIGINAL ONE!!
#              r       se_r
# [1,] 0.7126613 0.00187984

########### LETS SAVE THINGS 
fwrite(
  brainmeta_hg38,
  "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz",
  sep = "\t",
  quote = FALSE,
  compress = "gzip"
)

save.image("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/new_rb_20251201.RData")
###############################################################
####################### STREAMLINED WORKFLOW ##################
###############################################################
library(dplyr)
library(data.table)

brainmeta_hg38 <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/brainmeta_hg38.txt.gz")

setDT(brainmeta_hg38)

#### 1. Helper: preprocess a QTLtools nominal file (LIV / PM) ####
prep_qtl <- function(path) {
  df <- fread(path, header = FALSE, sep = " ", data.table = FALSE)
  colnames(df) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                    "strand", "n_variants", "distance", "variant_id",
                    "var_chr", "var_start", "var_end",
                    "pval", "beta", "top_flag")

  df$gene_id   <- sub("\\..*$", "", df$gene_id)  # strip version

  df$variant_pos <- as.integer(sub(".*:(\\d+)\\[.*", "\\1", df$variant_id))
  df$variant_chr <- as.integer(sub("chr(\\d+):.*", "\\1", df$variant_id))
  df$snp_key     <- paste0("chr", df$variant_chr, ":", df$variant_pos)

  df$z <- ifelse(
    df$pval > 0 & df$pval < 1,
    qnorm(1 - df$pval / 2),
    ifelse(df$pval == 0, 38, 1e-6)
  )
  df$se <- abs(df$beta) / abs(df$z)

  df
}

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

######### Function to compute rb for one tissue and one nsv
compute_rb_for_tissue <- function(brainmeta, tissue_df, p_thresh = 5e-8, theta = 0) {
  # 1) merge on SNP + gene
  dat <- merge(
    brainmeta,
    tissue_df,
    by.x = c("snp_key", "Probe"), #this is where im using snp_key_hg38 for brainmeta!
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

#nsv_vals <- 1:30 ##job was killed/takes a long time if do everything at once!

##so gotta do it in chunks!

#nsv_vals <- 1:3
#nsv_vals <- 4:6
#nsv_vals <- 7:9
#nsv_vals <- 10:12
#nsv_vals <- 13:15 
#nsv_vals <- 16:18
#nsv_vals <- 19:21
#nsv_vals <- 22:24
#nsv_vals <- 25:27
#nsv_vals <- 28:30

nsv_vals <- 4:4
nsv_vals <- 5:5
nsv_vals <- 6:6

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
  rb_liv <- compute_rb_for_tissue(brainmeta_hg38, liv_df, theta = 0)
  rb_pm  <- compute_rb_for_tissue(brainmeta_hg38, pm_df,  theta = 0)
  
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
write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv6to6_20251201.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4to4_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv5to5_20251201.csv",row.names = FALSE)

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1to3_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4to6_20251201.csv",row.names = FALSE) ### FAILED
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7to9_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv10to12_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13to15_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv16to18_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv19to21_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22to24_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv25to27_20251201.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv28to30_20251201.csv",row.names = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)

nsv1 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1to3_20251201.csv")
nsv7 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7to9_20251201.csv")
nsv10 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv10to12_20251201.csv")
nsv13 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13to15_20251201.csv")
nsv16 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv16to18_20251201.csv")
nsv19 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv19to21_20251201.csv")
nsv22 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22to24_20251201.csv")
nsv25 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv25to27_20251201.csv")
nsv28 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv28to30_20251201.csv")

nsv4 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4to4_20251201.csv")
nsv5 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv5to5_20251201.csv")
nsv6 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv6to6_20251201.csv")

rb_all <- rbind(nsv1,nsv4, nsv5, nsv6, nsv7, nsv10, nsv13, nsv16, nsv19, nsv22, nsv25, nsv28)

head(rb_all)

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

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference_20251201.pdf", p, width = 14, height = 8)

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
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
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

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference_with_3lines_20251201.pdf", p2, width = 14, height = 8)

######################
###################### DEC 5th 2025
######################
cd /sc/arion/projects/mscic1/results/jolie/LBP/eqtl
ml R

args <- commandArgs(trailingOnly = TRUE)
nsv <- as.integer(args[1])
nsv_vals <- nsv:nsv

bsub -q premium -P acc_mscic1 -n 10 -W 144:00 -R "rusage[mem=20000]" \
     -J "rb_array[1-30]" \
     "Rscript 00_right_eqtl_rb_for_job_submission.r \$LSB_JOBINDEX"

Job <216596432> is submitted to queue <premium>. #exited reach memory limit

bjobs -a -J rb_array

# Job <216597045> is submitted to queue <premium>.
# bjobs -l 216596432[29]

######################
###################### DEC 11th 2025
######################
cd /sc/arion/projects/mscic1/results/jolie/LBP/eqtl
ml R

bsub -q premium -P acc_mscic1 -n 1 -W 144:00 -R "rusage[mem=20000]" \
  -J "rb_array[1-30]" \
  "Rscript 0_eqtl_rb_for_job_submission.r \$LSB_JOBINDEX"

Job <217723225> is submitted to queue <premium>.

###

bsub -q premium -P acc_mscic1 -n 1 -W 144:00 -R "rusage[mem=80000]" \
  -J "rb_array[1-30]" \
  "Rscript 0_leadsnp_eqtl_rb.r \$LSB_JOBINDEX"

Job <217741007> is submitted to queue <premium>. #limit is 100GB
#Job <217736973> is submitted to queue <premium>.

### read the csvs and merge them together 
#!/usr/bin/env Rscript

library(dplyr)

# directory where all per-NSV rb files are saved
rb_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl"

# pattern matching all per-nsv CSVs
# tweak this if your date stamp / prefix changed
rb_files <- list.files(
  path    = rb_dir,
  pattern = "^rb_brainmeta_liv_pm_by_nsv[0-9]+to[0-9]+_20251201\\.csv$",
  full.names = TRUE
)

cat("Found", length(rb_files), "rb files:\n")
print(rb_files)

if (length(rb_files) == 0) {
  stop("No rb per-NSV files found. Check 'rb_dir' or pattern.")
}

# read and combine
rb_list <- lapply(rb_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  # optional: add filename as a column for debugging
  df$source_file <- basename(f)
  df
})

rb_all <- bind_rows(rb_list)

cat("Combined rows:", nrow(rb_all), "\n")
cat("Unique NSVs:", paste(sort(unique(rb_all$nsv)), collapse = ", "), "\n")
cat("Tissues:", paste(unique(rb_all$tissue), collapse = ", "), "\n")

# save combined summary
out_path <- file.path(rb_dir, "rb_brainmeta_liv_pm_all_nsv_20251201_merged.csv")
write.csv(rb_all, out_path, row.names = FALSE)

cat("Saved merged rb summary to:\n", out_path, "\n")



######################
###################### DEC 2nd 2025
######################
## new to-do:
#1. keep the original code and pick one snp-gene pair per gene
#2. update the code to use all SNP–gene pairs that are significant in the discovery dataset 
#(BrainMeta) only for all snp-gene pair for all genes 
#3. do the same as two but one snp-gene pair per gene

dim(brainmeta_hg38) #283,309,353        16

# choose the SNP–gene pair with the smallest p-value per gene (Probe)
brainmeta_lead <- brainmeta_hg38[
  , .SD[which.min(p)],  # row with min p for that gene
  by = Probe            # one row per gene
]

length(unique(brainmeta_hg38$Probe)) #16743

##expecting brainmeta_lead to have 16743 rows given there are 16743 genes
dim(brainmeta_lead) #16743    16

#base_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/results_qtltools_nominal/MATURERNA"

#nsv_vals <- 1:30 ##job was killed/takes a long time if do everything at once!
##so gotta do it in chunks!

#nsv_vals <- 1:3
#nsv_vals <- 4:6 ##
#nsv_vals <- 7:9
#nsv_vals <- 10:12
#nsv_vals <- 13:15 
#nsv_vals <- 16:18
#nsv_vals <- 19:21
nsv_vals <- 22:24
#nsv_vals <- 25:27
#nsv_vals <- 28:30

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
  rb_liv <- compute_rb_for_tissue(brainmeta_lead, liv_df, theta = 0) ## CHANGED brainmeta_lead
  rb_pm  <- compute_rb_for_tissue(brainmeta_lead, pm_df,  theta = 0) ## CHANGED brainmeta_lead
  
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
write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22to24_pergene_20251202.csv",row.names = FALSE) FAILED

#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv28to30_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv25to27_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22to24_pergene_20251202.csv",row.names = FALSE) FAILED
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv19to21_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv16to18_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13to15_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv10to12_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7to9_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4to6_pergene_20251202.csv",row.names = FALSE)
#write.csv(rb_summary,"/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1to3_pergene_20251202.csv",row.names = FALSE)

nsv1 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv1to3_pergene_20251202.csv")
nsv4 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv4to6_pergene_20251202.csv")
nsv7 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv7to9_pergene_20251202.csv")
nsv10 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv10to12_pergene_20251202.csv")
nsv13 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv13to15_pergene_20251202.csv")
nsv16 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv16to18_pergene_20251202.csv")
nsv19 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv19to21_pergene_20251202.csv")
nsv22 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv22to24_pergene_20251202.csv")
nsv25 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv25to27_pergene_20251202.csv")
nsv28 <- read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/rb_brainmeta_liv_pm_by_nsv28to30_pergene_20251202.csv")

rb_all <- rbind(nsv1,nsv4, nsv5, nsv6, nsv7, nsv10, nsv13, nsv16, nsv19, nsv22, nsv25, nsv28)

head(rb_all)

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
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
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

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference_with_3lines_pergene_20251202.pdf", p2, width = 14, height = 8)

###############
############### created a new file: 03_eqtl_rb_for_job_submission.r - this is using the old approach
# significant snp-gene pairs in BOTH brainmeta and liv/pm but just filter to one SNP-gene pair PER gene in the shared dataset
###############
CDIR="/sc/arion/projects/mscic1/results/jolie/LBP/eqtl"
RSC="$CDIR/03_eqtl_rb_for_job_submission.r"

for tissue in LIV PM; do
  for start in 1 4 7 10 13 16 19 22 25 28; do
    end=$((start + 2))
    out="${CDIR}/rb_lead_${tissue}_nsv${start}to${end}_20251202.csv"

    echo "Submitting job for ${tissue}, nsv ${start}-${end}"

    bsub -q premium \
         -P acc_mscic1 \
         -n 30 \
         -W 144:00 \
         -R "rusage[mem=4000]" \
         -R "span[hosts=1]" \
         -o %J.stdout \
         -eo %J.stderr \
         "module load R/4.2.0; Rscript $RSC ${tissue} ${start} ${end} ${out}"
  done
done

##jobs: 
#215866004

###great, next is Task 2:
#Use all SNP–gene pairs that are significant in BrainMeta (p < 5e-8), regardless of whether they are significant in LIV or PM.
#No filtering to one SNP per gene
#rb is computed using BrainMeta-significant SNP–gene pairs only
#This is the “BrainMeta-powered” version your PI wanted (avoids low LIV/PM power problem)

###############
############### created a new file: 04_eqtl_rb_job_ALL_gene_snp_pair.r
###############
CDIR="/sc/arion/projects/mscic1/results/jolie/LBP/eqtl"
RSC="$CDIR/04_eqtl_rb_job_ALL_gene_snp_pair.r"

for mode in bm_sig bm_all; do
  for start in 1 4 7 10 13 16 19 22 25 28; do
    end=$((start + 2))
    out="${CDIR}/rb_task2_${mode}_nsv${start}to${end}_20251202.csv"

    echo "Submitting job: mode=${mode}, nsv ${start}-${end}"

    bsub -q premium \
         -P acc_mscic1 \
         -n 30 \
         -W 144:00 \
         -R "rusage[mem=10000]" \
         -R "span[hosts=1]" \
         -o %J.stdout \
         -eo %J.stderr \
         "module load R/4.2.0; Rscript $RSC ${mode} ${start} ${end} ${out}"
  done
done

###############
############### created a new file: 05_eqtl_rb_job_PERGENE_shared.r
###############
CDIR="/sc/arion/projects/mscic1/results/jolie/LBP/eqtl"
RSC="$CDIR/05_eqtl_rb_job_PERGENE_shared.r"

for mode in bm_sig bm_all; do
  for start in 1 4 7 10 13 16 19 22 25 28; do
    end=$((start + 2))
    out="${CDIR}/rb_task3_pergene_${mode}_nsv${start}to${end}_20251202.csv"

    echo "Submitting job: mode=${mode}, nsv ${start}-${end}"

    bsub -q premium \
         -P acc_mscic1 \
         -n 30 \
         -W 144:00 \
         -R "rusage[mem=8000]" \
         -R "span[hosts=1]" \
         -o %J.stdout \
         -eo %J.stderr \
         "module load R/4.2.0; Rscript $RSC ${mode} ${start} ${end} ${out}"
  done
done

###############################
########### OUTPUTS ###########
###############################
## TASK 1
rb_lead_LIV_nsv1to3_20251202.csv #yes
rb_lead_LIV_nsv4to6_20251202.csv
...
rb_lead_PM_nsv28to30_20251202.csv #yes

## TASK 2 - might be too big tbh
rb_task2_bm_sig_nsv1to3_20251202.csv
rb_task2_bm_all_nsv1to3_20251202.csv
...
rb_task2_bm_sig_nsv28to30_20251202.csv
rb_task2_bm_all_nsv28to30_20251202.csv

## TASK 3
rb_task3_pergene_bm_sig_nsv1to3_20251202.csv #yes
rb_task3_pergene_bm_all_nsv1to3_20251202.csv #yes
...
rb_task3_pergene_bm_sig_nsv28to30_20251202.csv
rb_task3_pergene_bm_all_nsv28to30_20251202.csv

#########################################
########### PLOT TASK 1 AND 3 ###########
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)

# directory containing the files
cdir <- "/sc/arion/projects/mscic1/results/jolie/LBP/eqtl"

# read all rb_lead_*_20251202.csv files
rb_all <- list.files(
  path = cdir,
  pattern = "^rb_lead_.*_20251202\\.csv$",
  full.names = TRUE
) %>%
  map_df(~ read_csv(.x) %>% mutate(source_file = basename(.x)))

dim(rb_all) #60 10 - expecting 60

head(rb_all)

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
head(rb_tests)

sub <- filter(rb_tests,nsv =="1") ##difference is negative! so rb_PM > rb_LIV

rb_tests$p_two_sided ##no longer significant

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
  scale_y_continuous(limits = c(-0.1, 0.8), breaks = seq(-0.1, 0.8, by = 0.1)) +
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

ggsave("/hpc/users/hoangd02/www/plots/lbp/eqtl_rb_difference_with_3lines_pergene_20251204.pdf", p2, width = 14, height = 8)
