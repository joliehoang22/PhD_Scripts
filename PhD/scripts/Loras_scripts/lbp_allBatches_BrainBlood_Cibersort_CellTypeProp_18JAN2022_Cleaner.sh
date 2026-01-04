#cd /sc/arion/projects/mscic1/results/Noam/cibersortx
cd /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx
#you need to cd into the directory above, it doesn't work otherwise

#the 4 refs that noam used for the blood data: 1) input_LM22, 2) input_NSCLC_PBMC, 3) input_wilk, 4) input_SCP424

# module load singularity
module load singularity/3.6.4
module unload R

TOKEN=f609e4fd5544b7021f061f9d4c54ce79 

export SINGULARITY_DOCKER_USERNAME=lora.liharska@icahn.mssm.edu
export SINGULARITY_DOCKER_PASSWORD=Rewnaser1


# export SINGULARITY_DOCKER_USERNAME=YOUR_USERNAME
# export SINGULARITY_DOCKER_PASSWORD=YOUR_PASSWORD

#TOKEN_laptop=YOURLAPTOPTOKEN
#TOKEN=YOURMINERVATOKEN

# R_LIBS=/hpc/packages/minerva-centos7/rpackages/3.5.3/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.8:/hpc/packages/minerva-centos7/R/3.5.3/lib64/R/library
	
#mkdir fractions
#singularity pull --docker-login /sc/arion/projects/mscic1/results/Noam/cibersortx/fractions docker://cibersortx/fractions:latest 

singularity pull --docker-login /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/fractions2 docker://cibersortx/fractions:latest

cp /sc/arion/projects/mscic1/results/Noam/cibersortx/fractions /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/

#mkdir outdir_cibersortx

##############################################################
mkdir data

cp /sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/data/

##############################################################

# singularity shell docker://fractions.sif

# singularity run docker://fractions.sif -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/outdir /sc/arion/projects/mscic1/results/Noam/cibersortx/fractions --token e4dc4243b3f3383afc4b4acdedec29b6
# /hpc/packages/minerva-centos7/rpackages/3.5.3/site-library


# fwrite(tpm,file="expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
# fwrite(tpm_geneNames,file="expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_geneName.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#####################################################################################################################
#this is the test to make sure it works
#####################################################################################################################

#first -B is my expression data; i created the folder "data"
#the second -B is where you save the output, can be in the same folder

# singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture Fig2b-WholeBlood_RNAseq.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100 


singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/data/ -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/data/ --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/data/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100 

#####################################################################################################################
#when you have single cell data -- single_cell above is TRUE --> this pertains to when the REFERENCE is single cell

#-- single_cell -- 2 diff reference types you can give cibersort: 

#1. LM22-ref-sample.txt --> the one that they provide, noam will give it to me LM22; matrix where each row is a gene and each col is for one cell type the expression for the gene

#2. LM22 classes.txt --> tells you for each type of immune cell, first col is the label of the cell and the rest of the cols is telling you for each sample from the other file, where they belong; the mapping for the other file, what kind of cell each column label represent

#You use LM22-classes.txt when you don't have single cell data; would use the mapping corresponding to your reference

# --mixture is your own bulk expression data

# --perm is to assign a p value of confidence to the estimation

# all of the above is for non single cell data
#in the single cell reference, there's no mapping file

#for my bulk data, run for all of the references, one is the single cell and the other is not single cell, it's facs sorted etc

#then do the correlation to the cdc blood panel for the samples to see which reference is the best reference

###########################
# --sigmatrix --> need to see what this is in it, it's a matrix that's much smaller and you have each gene and each cell type that exists in that reference in the single cell data, related to how highly expressed the gene is in the specific cell type

# starting on line 70 to see what is the thing for sigmatrix

# noam has never used sigmatrix
# sigmatrix is not required because it will make it internally but if you give it one it will be faster

#######################################################
#will need to fix the gene names based on the reference gene names

#for the Lake dataset for brain, need to reformat it per the cibersort reference, check noam's references
#do this AFTER running the provided reference

/sc/arion/projects/psychgen/lbp/files/lake_for_cibersort_3.Rdata
#######################################################
find /sc/arion/projects/mscic1/results/Noam/cibersortx -wholename "*LM22-ref-sample.txt*" -print
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22_MainCovid/LM22-ref-sample.txt
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22_MainCovid/old/LM22-ref-sample.txt
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22/LM22-ref-sample.txt
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22/._LM22-ref-sample.txt
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_KD/input_LM22/LM22-ref-sample.txt
# /sc/arion/projects/mscic1/results/Noam/cibersortx/input_KD/input_LM22/._LM22-ref-sample.txt


diff /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22_MainCovid/LM22-ref-sample.txt /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22/LM22-ref-sample.txt
#If diff shows no output, that means the two files are the same

cp -r /sc/arion/projects/mscic1/results/Noam/cibersortx/input_LM22/ /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/data

#####################################################################################################################
#Need to make sure gene names are the same for reference and mixture -- do this in R
#####################################################################################################################
#SINGLE CELL REFERENCE
# when you have single cell --data single_cell is TRUE
# this pertains to when the REFERENCE is single cell
# in the single cell reference, there's no mapping file
#############################################################################################
Name=NSCLC_PBMC_MainCovid
#cd /sc/arion/projects/mscic1/results/Noam/cibersortx/
cd /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx
#you need to cd into the directory above, it doesn't work otherwise

# mkdir input_${Name}
# mkdir outdir_${Name}
# mkdir input_${Name}/old/
# mkdir outdir_${Name}/old/
# mv input_${Name}/* input_${Name}/old/
# mv outdir_${Name}/* outdir_${Name}/old/
# cp /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_* input_${Name}/
# cp /sc/arion/projects/mscic1/results/Noam/MainCovid/expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt input_${Name}/

# singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/input_${Name}/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/outdir_${Name}/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100


cp /sc/arion/projects/mscic1/results/Noam/cibersortx/input_NSCLC_PBMC_MainCovid/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata
cp /sc/arion/projects/mscic1/results/Noam/cibersortx/input_NSCLC_PBMC_MainCovid/Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output3/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture lbp_allBatches_RAPiD_featureCounts_TPM_AllGenesWithGeneSymbol_forCibersort_BLOODOnly_243samples_NoBLOOD582_24JAN2022.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100


##################################################################################################
#cd /sc/arion/projects/mscic1/results/Noam/cibersortx
cd /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx
#you need to cd into the directory above, it doesn't work otherwise

# module load singularity
module load singularity/3.6.4
module unload R

#export SINGULARITY_DOCKER_USERNAME=YOUR_USERNAME
#export SINGULARITY_DOCKER_PASSWORD=YOUR_PASSWORD
#TOKEN_laptop=YOURLAPTOPTOKEN
#TOKEN=YOURMINERVATOKEN

# R_LIBS=/hpc/packages/minerva-centos7/rpackages/3.5.3/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.8:/hpc/packages/minerva-centos7/R/3.5.3/lib64/R/library
	
#mkdir fractions
#singularity pull --docker-login /sc/arion/projects/mscic1/results/Noam/cibersortx/fractions docker://cibersortx/fractions:latest 

#mkdir outdir_cibersortx

# singularity shell docker://fractions.sif

# singularity run docker://fractions.sif -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/outdir /sc/arion/projects/mscic1/results/Noam/cibersortx/fractions --token f609e4fd5544b7021f061f9d4c54ce79
# /hpc/packages/minerva-centos7/rpackages/3.5.3/site-library

# fwrite(tpm,file="expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
# fwrite(tpm_geneNames,file="expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_geneName.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture Fig2b-WholeBlood_RNAseq.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100 



##################################################################################################
#HERE there is no --single_cell argument set to TRUE so this means that the reference data (LM22-ref-sample.txt) is not single cell and thus requires the mapping file --phenoclasses LM22-classes.txt; notice how for the reference files where --single_cell argument is set to TRUE, there is no --phenoclasses argument
##################################################################################################
Name=LM22_MainCovid
cd /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx

module load singularity/3.6.4
module unload R
#cd /sc/arion/projects/mscic1/results/Noam/cibersortx/
#you need to cd into the directory above, it doesn't work otherwise


# mkdir input_${Name}
# mkdir outdir_${Name}
# mkdir input_${Name}/old/
# mkdir outdir_${Name}/old/
# mv input_${Name}/* input_${Name}/old/
# mv outdir_${Name}/* outdir_${Name}/old/
# cp input_LM22/LM22* input_${Name}/
# cp /sc/arion/projects/mscic1/results/Noam/MainCovid/expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt input_${Name}/

# singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/input_${Name}/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/outdir_${Name}/:/src/outdir fractions /src/CIBERSORTxFractions --username lora.liharska@icahn.mssm.edu --token f609e4fd5544b7021f061f9d4c54ce79 --refsample LM22-ref-sample.txt --phenoclasses LM22-classes.txt --mixture MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt  --QN TRUE --perm 100


#first -B alldata is being defined as the source data; cibersort creating local environment with specific structure, need to tell it where the elements are -- source folder created called src/data in which the data from all data will go; src is local evn so it doesn't exist
#second -B creates local env output directory; outdir is the parameter in the local env of cibersort and setting the path to it; fractions is the name of the software; setting /src/CIBERSORTxFractions to be the executable fractions software, syntax is different, no ":""
#--refsample -- don't put the full path, because it will look in alldata for it; everything is in reference to the -B now

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output2/:/src/outdir fractions /src/CIBERSORTxFractions --username lora.liharska@icahn.mssm.edu --token f609e4fd5544b7021f061f9d4c54ce79 --refsample LM22-ref-sample.txt --phenoclasses LM22-classes.txt --mixture lbp_allBatches_RAPiD_featureCounts_TPM_AllGenesWithGeneSymbol_forCibersort_BLOODOnly_243samples_NoBLOOD582_24JAN2022.txt --QN TRUE --perm 100

##################################################################################################
#make a new example directory and run this 

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output3/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture Fig2b-WholeBlood_RNAseq.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100 
##################################################################################################
Name=NSCLC_PBMC_MainCovid
#cd /sc/arion/projects/mscic1/results/Noam/cibersortx/
cd /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx
#you need to cd into the directory above, it doesn't work otherwise

# mkdir input_${Name}
# mkdir outdir_${Name}
# mkdir input_${Name}/old/
# mkdir outdir_${Name}/old/
# mv input_${Name}/* input_${Name}/old/
# mv outdir_${Name}/* outdir_${Name}/old/
# cp /sc/arion/projects/mscic1/results/Noam/cibersortx/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_* input_${Name}/
# cp /sc/arion/projects/mscic1/results/Noam/MainCovid/expression/MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt input_${Name}/

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/mscic1/results/Noam/cibersortx/input_${Name}/:/src/data -B /sc/arion/projects/mscic1/results/Noam/cibersortx/outdir_${Name}/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture MainBatch_all_tpm_matrix_noOutliers_BatchCTL_summed_2020-12-31_geneName.txt --fraction 0 --sigmatrix Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt --perm 100

##############################################################################################
#input_SCP424 -- non-single cell

cp /sc/arion/projects/mscic1/results/Noam/cibersortx/input_SCP424/SCP424_full_counts.umi_Sum1584toinf_Count630toinf_ds1584.tsv /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata

cp /sc/arion/projects/mscic1/results/Noam/cibersortx/input_SCP424/SCP424_full_counts.umi_Sum1584toinf_Count630toinf_ds1584_cellTypeLabels.tsv /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata


singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output4/:/src/outdir fractions /src/CIBERSORTxFractions --username lora.liharska@icahn.mssm.edu --token f609e4fd5544b7021f061f9d4c54ce79 --refsample SCP424_full_counts.umi_Sum1584toinf_Count630toinf_ds1584_EnsemblNoDots.txt --phenoclasses SCP424_full_counts.umi_Sum1584toinf_Count630toinf_ds1584_cellTypeLabels.tsv --mixture lbp_allBatches_RAPiD_featureCounts_TPM_AllGenesWithENSEMBL_forCibersort_BLOODOnly_243samples_NoBLOOD582_26JAN2022.txt --QN TRUE --perm 100


singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output4/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample SCP424_full_counts.umi_Sum1584toinf_Count630toinf_ds1584_cellTypeLabels.tsv --mixture lbp_allBatches_RAPiD_featureCounts_TPM_AllGenesWithENSEMBL_forCibersort_BLOODOnly_243samples_NoBLOOD582_26JAN2022.txt --fraction 0 --perm 100



##############################################################################################
#Single Cell 2
#input_wilk

cp /sc/arion/projects/mscic1/results/Noam/cibersortx/input_wilk/wilk_refs_forNoam_fine.tsv /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output5/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample wilk_refs_forNoam_fine.tsv --mixture lbp_allBatches_RAPiD_featureCounts_TPM_AllGenesWithGeneSymbol_forCibersort_BLOODOnly_243samples_NoBLOOD582_24JAN2022.txt --fraction 0 --perm 100



##############################################################################################
#Single Cell 2
#input_lake

singularity exec --env R_LIBS_USER="" --env R_LIBS="" -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata/:/src/data -B /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_lake/:/src/outdir fractions /src/CIBERSORTxFractions --token f609e4fd5544b7021f061f9d4c54ce79 --username lora.liharska@icahn.mssm.edu --single_cell TRUE --refsample lake_singleCell_reference_01FEB2022.txt --mixture lbp_allBatches_RAPiD_featureCounts_TPM_All58884ENSEMBLGenes_NoDups_forCibersort_BRAINOnly_533samples_01FEB2022.txt --fraction 0 --perm 100
