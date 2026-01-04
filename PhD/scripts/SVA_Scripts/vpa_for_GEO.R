#################################################################################################
####################### This figure includes VPA of real data for PD GEO ########################

#code derived from Extended_Figure5.sh and GEO.R

#there are 8 GEOs: GSE20292, GSE8397, GSE20164, GSE20163, GSE24378, GSE7621, GSE20141, GSE49036, each with its own formula listed below and in Hoang_SVA manuscript_20251014.docx in Dropbox. The goal is to make a variance partition (VP) plot for each of the GEO and combine them with the existing plot (Extended_Figure5.sh/ExtendedDataFigure5.pdfhell)

#this is a supplemental figure so it doesn't have to be perfect but it should still be good enough! Ppease be consistent with naming and ordering (if possible). in the Methods section in the text, I put "Diagnosis" in lieu of "Status," so I'd have to update that before submission. Please feel free to slack me anytime for anything!

## Instructions
#1. Load in the data (I already list them out for you)
#2. Add appropriate column names (see below; also see ExtendedDataFigure5.pdf) -- see code in Extended_Figure5.sh
#3. Plot -- see code in Extended_Figure5.sh
#4. Combine the plots -- see code in Extended_Figure5.sh

### This is the correct column names and order:
#### Status <= disease_status (ideally first column)
#### Sex <= gender
#### Residuals (for sure last column) -- note for the formula, Residuals are implied.

## Start:
### Dataset 1: GSE20292 ~ Formula: "Status"+"Age"+"Sex" -------------------------------------- DONE
vp_GSE20292 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 2: GSE8397 ~ Formula: "Status"+"Age"+"Sex"+"Brain_Region" -------------------------------------- DONE
vp_GSE8397 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 3: GSE20164 ~ Formula: "Status"+"Age"+"Sex" -------------------------------------- DONE
vp_GSE20164 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 4: GSE20163 ~ Formula: "Status"+"Age" -------------------------------------- DONE
vp_GSE20163 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 5: GSE24378 ~ Formula: "Status"+"Age" -------------------------------------- DONE
vp_GSE24378 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
 
### Dataset 6: GSE7621 ~ Formula: "Status"+"Sex" -------------------------------------- DONE
vp_GSE7621 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 7: GSE20141 ~ Formula: "Status" -------------------------------------- DONE
vp_GSE20141 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### Dataset 8: GSE49036 ~ Formula: "Status"+"RIN" 
vp_GSE49036 <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_vpa.txt", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### YAY, almost done! Last few steps,
#1. Patchwork with combo_all from Extended_Figure5.sh, add paneling A, B, C, etc. for the new 8 PD datasets
#2. Save it as a pdf and download it
#3. Put it in Dropbox > ExtendedDataFigure5_ABG.pdf (you can leave it in the main SVA folder and not other subfolder)

## THANK YOU SO MUCH!!!





















