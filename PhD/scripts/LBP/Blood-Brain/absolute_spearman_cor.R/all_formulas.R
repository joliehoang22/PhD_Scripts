### just a list of formulas

###### ###### ###### ###### ######
###### ###### BRAIN ###### #######
###### ###### ###### ###### ######
brain_full_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full_removed_outliers_no_residID_225samples.txt", data.table=FALSE) #this is the correct version
row.names(brain_full_no_residID)<- brain_full_no_residID$V1
brain_full_no_residID$V1 <- NULL
dim(brain_full_no_residID) #21356   225

brain_full<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_brain_form_full.txt", data.table=FALSE)
row.names(brain_full) <- brain_full$V1
brain_full$V1 <- NULL

###### ###### ###### ###### ######
###### ###### BLOOD ###### #######
###### ###### ###### ###### ######
############# FORM 3 #############
blood_form3<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_resid_ID_rin_STAR_Insertion_average_length.txt",data.table=FALSE)
rownames(blood_form3) <- blood_form3$V1; blood_form3$V1 <- NULL

blood_form3_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form3_no_indivdualID.txt",data.table=FALSE) 
rownames(blood_form3_no_residID) <- blood_form3_no_residID$V1; blood_form3_no_residID$V1 <- NULL

############# FORM 5 #############
blood_form5<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5.txt",data.table=FALSE)
rownames(blood_form5) <- blood_form5$V1; blood_form5$V1 <- NULL

blood_form5_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/resid_expr_blood_form5_no_indivdualID.txt",data.table=FALSE)
rownames(blood_form5_no_residID) <- blood_form5_no_residID$V1; blood_form5_no_residID$V1 <- NULL

