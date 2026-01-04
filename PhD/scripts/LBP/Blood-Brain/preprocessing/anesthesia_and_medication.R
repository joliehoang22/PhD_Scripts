## code derived from LBP_blood_brain_QC_20250729.R

## Missingness in anesthesia and medications
dop_daily_dose <- readRDS("/sc/arion/projects/psychgen/lbp/data/emr/lbp_lel2021_clinical/dopamine_daily_dosing.RDS")
dim(dop_daily_dose) #166 7
length(unique(dop_daily_dose$iid)) #108

dop_daily_dose <- dop_daily_dose[!duplicated(iid)] #getting rid of duplicates, keeping only the first row for each unique iids (the values duplicate for duplicated IDs)
dim(dop_daily_dose) #108   7

dop_daily_dose$IID_ISMMS <- dop_daily_dose$iid
dop_daily_dose_sub <- dop_daily_dose[,c("IID_ISMMS", "carbidopa","levodopa","entacapone")]
head(dop_daily_dose_sub)
#   IID_ISMMS carbidopa levodopa entacapone
#      <char>     <num>    <num>      <num>
#1:   PT-0002     25.00      100          0
#2:   PT-0010    312.50     1250          0
#3:   PT-0012    187.50      750       1000

blood_only_metadata <- merge(blood_only_metadata, dop_daily_dose_sub, by = "IID_ISMMS", all.x = TRUE)
dim(blood_only_metadata) #225 173

summary(blood_only_metadata$carbidopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   25.0   100.0   175.0   194.6   250.0   793.8      84 
summary(blood_only_metadata$levodopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    100     400     750     832    1050    3175      84 
summary(blood_only_metadata$entacapone)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.0     0.0     0.0   106.4     0.0  6000.0      84 

brain_only_metadata <- merge(brain_only_metadata, dop_daily_dose_sub, by = "IID_ISMMS", all.x = TRUE)
dim(brain_only_metadata) #225 173

summary(brain_only_metadata$carbidopa)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   25.0   100.0   175.0   194.6   250.0   793.8      84 
summary(brain_only_metadata$levodopa)
summary(brain_only_metadata$entacapone)
## same as above - for blood (the samples are paired so it should be the same!)

##need to take a closer look at the NAs..... and not sure if merge IID_ISMMS is the best way
