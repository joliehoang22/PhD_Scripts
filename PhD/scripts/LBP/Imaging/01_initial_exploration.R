## code derived from LBP_blood_brain_QC_20250729.R

##A bit of imaging -- do NOT change anything in the path! copy and paste brains in your folder if needed
#Path for Raw nifty
/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output

#Path for Freesurfer outputs
/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output_fs2

imaging_map <- fread("/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/mapping_files/scan_id_date_path_map.csv")

length(unique(imaging_map$subject_id)) #179

subset(imaging_map, subject_id == "PT-0018")
#FSPGR and FSPGR_contrast

subset(imaging_map, subject_id == "PT-0019")
#FSPGR and FSPGR_contrast


##potential confounders: look for scanner differences 
#T1-weighted: 
#GE: FSPGR | FSPGR_contrast = after a contrast agent (typically Gadolinium) used to enhanced signals in tumors, inflammation, etc.
#Siemens: MPRAGE | Ax T1 MPRAGE = Axial slice only?

# Extract unique subject IDs from the dictionary
IID_ISMMS <- unique(dictionary$IID_ISMMS)
# Subset imaging_map based on matching subject_id
imaging_map_subset <- imaging_map[subject_id %in% IID_ISMMS]

dim(imaging_map_subset) #309x9

table(imaging_map_subset$SeriesDescription2)
#   Ax T1 MPRAGE           FSPGR  FSPGR_contrast          MPRAGE MPRAGE_contrast 
#              1              56              39              95             118 

subset(imaging_map_subset, SeriesDescription2 == "Ax T1 MPRAGE") #PT-0064
subset(imaging_map_subset, subject_id == "PT-0064") 
#Ax T1 MPRAGE, FSPGR, and MPRAGE_contrast

##there are scans pre-surgery and post-surgery (which has the electrode) - subset to only scans pre-surgery?
#make two columns to compare acquisition data to surgery dates first 
imaging_map_subset$surgery1_timing <- ifelse(
  imaging_map_subset$AcquisitionDate < imaging_map_subset$surgery1,
  "pre-surgery", "post-surgery"
)

imaging_map_subset$surgery2_timing <- ifelse(
  imaging_map_subset$AcquisitionDate < imaging_map_subset$surgery2,
  "pre-surgery", "post-surgery"
)
#is this actually meaningful? surgery2 is always after surgery1 and pre-surgery of surgery2 doesn't mean true pre-surgery (aka before surgery1)

presurgery_scans <- filter(imaging_map_subset, surgery1_timing == "pre-surgery")
dim(presurgery_scans) #256  11

table(presurgery_scans$surgery1_timing) #TRUE
table(presurgery_scans$surgery2_timing) #TRUE

length(unique(presurgery_scans$subject_id)) #151 people, so some have multiple scans | also i think for dictionary i have 152 people?
table(presurgery_scans$subject_id)
table(dictionary$IID_ISMMS)

setdiff(unique(dictionary$IID_ISMMS), unique(presurgery_scans$subject_id)) 
##PT-0025 is in dictionary but not presurgery_scans

table(presurgery_scans$SeriesDescription2) #total: 256 scans ~ 1 week? double check with someone
#          FSPGR  FSPGR_contrast          MPRAGE MPRAGE_contrast 
#              4              39              95             118

length(unique(presurgery_scans$scan_id)) #256, good this match the dimension of the df so this means that all pt has unique scan number 
head(presurgery_scans)
#      V1 subject_id scan_id AcquisitionDate ContentTime SeriesDescription2
#1:    18    PT-0019  scan31      2014-03-28      140642     FSPGR_contrast
presurgery_scans$V1<-NULL
write.csv(presurgery_scans, "/sc/arion/projects/mscic1/results/jolie/LBP/imaging/presurgery_scans_subset_with_blood_brain_ge.csv")

test<-read.csv("/sc/arion/projects/mscic1/results/jolie/LBP/imaging/presurgery_scans_subset_with_blood_brain_ge.csv")
test$X <- NULL
