##this made the scan paths and folders so we could run freesurfer 
##run in parallel
DIR=...
cd ${DIR}
mkdir ${DIR}/output 
while read LINE
do 
CHECK=`echo ${LINE} | awk -F"|" '{print $1}'`
if [[ ${CHECK} != '"path"' ]]
then 
  NEWPATH=${DIR}/output/`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
  OLDPATH=`echo ${LINE} | awk -F "|" '{print $1}'`
  if [[ ! -d ${NEWPATH} ]]
  then mkdir ${NEWPATH}
  fi
  eval ln -f -s ${OLDPATH} ${DIR}/input/`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
fi 
done < map_input_output_imaging_sMRI.txt

while read LINE
do 
CHECK=`echo ${LINE} | awk -F"|" '{print $1}'`
if [[ ${CHECK} != '"path"' ]]
then 
  SCAN=`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
  I=${DIR}/input/${SCAN}
  O=${DIR}/output/${SCAN}
  dcm2niix -o ${O} -f "%f" -p y -z y ${I}
  STATUS=$?   
  if [[ ${STATUS} -eq 0 ]]
  then touch ${O}/${SCAN}.success
  else touch ${O}/${SCAN}.fail;exit 1
  fi
fi 
done < map_input_output_imaging_sMRI.txt

ls ${DIR}/output/scan*/* | grep ".success" | wc -l
ls ${DIR}/output/scan*/* | grep ".fail" | wc -l
ls | wc -l
ls ../input/ | wc -l

##this is the job you're submitting 
module load freesurfer/8.0.0-1 #load freesurfer/6.0.0 (4hrs/brain v8 compared to 8 hrs/brain in v6)
mkdir ${DIR}/output_fs/
cd ${DIR}/output/
for d in `ls -d *`; do
        cd ${d}
        bsub -J ${d} -n 4 -P acc_mscic1 -n 30 -W 144:00 -q premium -R "rusage[mem=10000]" -oo "fs.stdout" -eo "fs.stderr" recon-all -i ${d}.nii.gz -s ${d} -sd ${DIR}/output_fs/ -all -qcache;
        cd ..
done

grep Success ${DIR}/output/*.sdtout | wc -l

bsub -q premium -P acc_mscic1 -n 30 -W 144:00 -R rusage[mem=4000] -R span[hosts=1] -o %J.stdout -eo %J.stderr Rscript ${path}	blood-brain_distance_cor_form3.R

################
################ TRY 1 | NOV 17 2025 | 399 SCANS
################
## WORKFLOW: 
#Convert DICOM → NIfTI → Give FreeSurfer the NIfTI
#DICOM = raw images directly out of MRI, .dcm
#NIfTI = clean, standardized research format; One file per MRI volume (e.g., T1); .nii.gz; .nii
##The standard in neuroimaging analysis pipelines (FSL, AFNI, FreeSurfer).
cd /sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output
/sc/arion/projects/mscic1/results/jolie/LBP/imaging

mkdir -p /sc/arion/projects/mscic1/results/jolie/LBP/imaging
cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging

mkdir -p input # where we’ll put symlinks to existing scans
mkdir -p freesurfer_subjects #where FreeSurfer will write everything
mkdir -p logs #to store recon-all logs

##testing out symlink
ln -s \
  /sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan12 \
  input/scan12

ls input/scan12 #success
# fs.stderr  fs.stdout  scan12.json  scan12.nii.gz  scan12.success

cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging
module avail freesurfer
module load freesurfer/8.1.0   # or whatever the module is called
source $FREESURFER_HOME/SetUpFreeSurfer.sh

##tell FS where to put subjects
export SUBJECTS_DIR=/sc/arion/projects/mscic1/results/jolie/LBP/imaging/freesurfer_subjects
export FS_LICENSE=/sc/arion/projects/mscic1/results/jolie/LBP/imaging/config/license.txt

##verify: 
echo $SUBJECTS_DIR #would tell you where FS puts subjects
echo $FS_LICENSE

##run FS on scan12 only
cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging

recon-all \
  -s scan12 \
  -i input/scan12/scan12.nii.gz \
  -all \
  > logs/scan12_recon-all.log 2>&1

tail -n 20 logs/scan12_recon-all.log
ls $SUBJECTS_DIR/scan12

################################################################################################################
############################## RUNNING EVERYTHING IN PARALLEL VIA SUBMITTING JOBS ##############################
################################################################################################################
cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging
nano run_reconall_subject.sh

chmod +x run_reconall_subject.sh

###running one job
cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging

bsub -q premium \
     -P acc_mscic1 \
     -n 20 \
     -W 144:00 \
     -R "rusage[mem=4000]" \
     -R "span[hosts=1]" \
     -o %J.scan12.stdout \
     -eo %J.scan12.stderr \
     bash run_reconall_subject.sh scan12

#Job <212671012> is submitted to queue <premium>. 1st attempt
#Job <212673697> is submitted to queue <premium>. 2nd, stupid, remove scan 12
#Job <212675489> is submitted to queue <premium>.

# SynthStrip: Skull-Stripping for Any Brain Image
# A Hoopes, JS Mora, AV Dalca, B Fischl, M Hoffmann
# NeuroImage 206 (2022), 119474
# https://doi.org/10.1016/j.neuroimage.2022.119474

##watch job
bjobs        # see if it's PEND/RUN/DONE
tail -n 20 freesurfer_subjects/scan12/scripts/recon-all.log

tail -n 20 212671012.scan12.stderr   # replace %J with the actual job ID once you see it
tail -n 20 212673697.scan12.stderr

##########run for all subjects
cd /sc/arion/projects/mscic1/results/jolie/LBP/imaging

##make sure all symlinks exist
for s in /sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan*; do
  bname=$(basename "$s")
  ln -s "$s" "input/$bname" 2>/dev/null || true
done

ls -d input/scan* | head

##submit all scans as separate jobs, avoiding scan 12
for sdir in input/scan*; do
  subj=$(basename "$sdir")
  # Skip scan12 since it is already running
  if [[ "$subj" == "scan12" ]]; then
    echo "Skipping ${subj} (already running)"
    continue
  fi

  echo "Submitting recon-all for ${subj}"

  bsub -q premium \
       -P acc_mscic1 \
       -n 20 \
       -W 144:00 \
       -R "rusage[mem=4000]" \
       -R "span[hosts=1]" \
       -o logs/%J.${subj}.stdout \
       -eo logs/%J.${subj}.stderr \
       bash run_reconall_subject.sh "${subj}"
done

bjobs 
##started ~10pm 11/17/2025
##check 9am 11/28 - 5 jobs left

##FS outputs will be in:
/sc/arion/projects/mscic1/results/jolie/LBP/imaging/freesurfer_subjects/scan12

##check for 
aparc+aseg.mgz

##some scans have it and some dont
##a volume segmentation file created by the FreeSurfer software, which combines the automatic subcortical segmentation (aseg) 
#with the automatic cortical parcellation (aparc) into a single file


















head map_input_output_imaging_sMRI.txt

##make freesurfer script
nano run_freesurfer_batch.sh
#!/usr/bin/env bash
set -euo pipefail

# Usage: bash run_freesurfer_batch.sh /path/to/base_dir /path/to/map_file
DIR=${1:-"."}
MAP_FILE=${2:-"${DIR}/map_input_output_imaging_sMRI.txt"}

cd "$DIR"

# Directories for symlinks, outputs, and FS subjects
mkdir -p input output logs
export SUBJECTS_DIR="${DIR}/freesurfer_subjects"
mkdir -p "$SUBJECTS_DIR"

echo "Base directory  : $DIR"
echo "Mapping file    : $MAP_FILE"
echo "SUBJECTS_DIR    : $SUBJECTS_DIR"
echo

# Main loop over the map file
while IFS='|' read -r raw_path scan_id rest; do
  # Strip quotes and whitespace
  raw_path=$(echo "$raw_path" | tr -d '"[:space:]')
  scan_id=$(echo "$scan_id" | tr -d '"[:space:]')

  # Skip header or empty lines
  if [[ -z "$raw_path" || "$raw_path" == "path" ]]; then
    continue
  fi

  echo "=== Processing scan: $scan_id ==="
  echo "  Raw DICOM path : $raw_path"

  # Define input symlink and output directory
  in_link="${DIR}/input/${scan_id}"
  out_dir="${DIR}/output/${scan_id}"

  mkdir -p "$out_dir"

  # Create or update symlink
  ln -sfn "$raw_path" "$in_link"
  echo "  Symlink        : $in_link -> $raw_path"

  # ---- DICOM → NIfTI with dcm2niix ----
  echo "  Running dcm2niix..."
  dcm2niix \
    -o "$out_dir" \
    -f "%f" \
    -p y \
    -z y \
    "$in_link" \
    > "logs/${scan_id}_dcm2niix.log" 2>&1 || {
      echo "  [ERROR] dcm2niix failed for $scan_id (see logs/${scan_id}_dcm2niix.log)"
      touch "${out_dir}/${scan_id}.fail"
      continue
    }

  touch "${out_dir}/${scan_id}.success"
  echo "  dcm2niix done. NIfTIs in: $out_dir"

  # ---- Run FreeSurfer recon-all ----
  # You can tweak this if you only want some scans or need more flags.
  nifti_inputs=("$out_dir"/*.nii*)
  if [[ ${#nifti_inputs[@]} -eq 0 ]]; then
    echo "  [WARNING] No NIfTI files found in $out_dir – skipping recon-all."
    continue
  fi

  echo "  Running recon-all for subject: $scan_id"
  recon-all \
    -s "$scan_id" \
    -i "${nifti_inputs[@]}" \
    -all \
    > "logs/${scan_id}_recon-all.log" 2>&1 || {
      echo "  [ERROR] recon-all failed for $scan_id (see logs/${scan_id}_recon-all.log)"
      # We don't exit here so other subjects can still run.
    }

  echo "=== Done with $scan_id ==="
  echo

done < "$MAP_FILE"

echo "--------------------------------------"
echo "Summary:"
echo "Successful dcm2niix: $(find output -name '*.success' | wc -l)"
echo "Failed dcm2niix    : $(find output -name '*.fail'    | wc -l)"
echo "FreeSurfer subjects: $(ls -1 \"$SUBJECTS_DIR\" | wc -l)"
echo "--------------------------------------"
