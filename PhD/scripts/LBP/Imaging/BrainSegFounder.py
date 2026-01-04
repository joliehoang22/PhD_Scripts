#BrainSegFounder

#https://github.com/lab-smile/BrainSegFounder

#The model was designed to segment brain tumors or stroke lesions. 
# After pre‑training (Stage 1 and 2), it attaches a decoder (UNet‑style) to the 
# pretrained encoder for supervised fine‑tuning.

#Model architecture: Swin‑UNETR

#Stage 1 pre‑training (self-supervised)
#Input: Multimodal healthy brain MRIs (T1‑w + T2‑w)
#Output: Learned encoder weights (no labels)

#Stage 2 pre‑training (SSL adaptation)
#Input: Disease‑specific datasets (BraTS or ATLAS images)
#Adapted encoder (still no segmentation outputs)

#Stage 3 Supervised fine‑tuning
#Input: BraTS: four-channel MRI (T1/T1c/T2/FLAIR); ATLAS: single T1‑ce
#Output: Voxel‑level segmentations (tumors or lesions)

## Pretrain
# This directory contains scripts and pretrained models for Stage 1 pretraining. 
# Our models are pretrained on the UK Biobank dataset using Self-Supervised Learning (SSL) heads 
# and a 3-way loss function. This stage sets the foundation for advanced feature extraction crucial 
# for downstream tasks.

#Combining these three types of SSL objectives helps the model:
#Learn local features (reconstruction loss),
#Learn global representations (contrastive/similarity loss),
#And maintain representation stability (consistency/ regularization loss).

cd /sc/arion/projects/mscic1/results/jolie/LBP/

#python -m venv brainsegfounder_env
source brainsegfounder_env/bin/activate
#otherwise: "conda create -n brainsegfounder_env" then "conda activate brainsegfounder_env"

module load git

#from their git
git clone https://www.github.com/lab-smile/BrainSegFounder.git
cd BrainSegFounder

pip install -r https://raw.githubusercontent.com/Project-MONAI/MONAI/dev/requirements-dev.txt
pip install -U git+https://github.com/npnl/bidsio
pip install numpy torch torchvision torchaudio monai nibabel scikit-learn pandas scipy tqdm tensorboard

#ERROR: Could not install packages due to an OSError: [Errno 122] Disk quota exceeded

#look at 
cd pretrain

#i am /hpc/users/hoangd02/BrainSegFounder/pretrain

###next time START HERE: 
module purge                    # unload any random modules
module load python/3.10          # or whichever stable version Minerva has

source /sc/arion/projects/mscic1/users/hoangd02/venv/brainseg_env/bin/activate
module load git

cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain

# downloaded UK Biobank images from Google Drive link - in README.md
# JSON files containing the folds used for our data and PyTorch pretrained models 
# can be downloaded from this Google Drive link.

## it's important to check what the parser arguments are in the script: main_T1T2.py
python main_T1T2.py \
  --split_json /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/mini_fold.json \
  --logdir /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke \
  --num_steps 10 \
  --batch_size 1 \
  --num_workers 0 \
  --modality T1T2

#main_T1T2.py = the BrainSegFounder training script for two-modality input (T1 and T2 MRI)
#it contains the code to load your dataset, build the neural network, and run the pretraining process.

#data_list = an argument to the script that tells the code where to find your dataset description file
#mini_fold.json has training and validation keys (training has 1 image key: [T1_path, T2_path]; val is empty)

#--output_dir runs/stage1_smoke = tells the script where to save results
##inside runs/stage1_smoke/ there is model checkpoints, training logs, etc.

#--max_steps 10 =  limit trainign to 10 udate steps. each step load a batch of images, pass it through the model, compute loss, and update the model's weights 

#--batch_size 1 = number of sample processed in 1 step. 1 means 1 subject's T1+T2 pair at a time 

#--num_workers 0 = number of worker processes for loading data in parallel; 0 menas don't use extra processes

## CAN EITHER RUN JOBS VIA GPU or CPU
######################create a script first - START
cat > smoke_t1t2.sh <<'SH'
#!/bin/bash
# 1) start from a clean module state
module purge

# 2) load the same Python you used to make your venv (adjust if yours is different)
module load python/3.10

# 3) activate your venv (so python & pip & packages come from here)
source /sc/arion/projects/mscic1/users/hoangd02/venv/brainseg_env/bin/activate

# 4) run the training script with your mini split JSON
python /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain/main_T1T2.py \
  --split_json /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/mini_fold.json \
  --logdir     /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke \
  --num_steps  10 \
  --batch_size 1 \
  --num_workers 0 \
  --modality   T1T2
SH

# make it executable
chmod +x smoke_t1t2.sh
######################create a script first - END

##SUBMINT
bsub -P acc_mscic1 -q gpu -W 00:10 \
     -n 4 -R "span[hosts=1]" -R a100 -gpu "num=1" \
     -R rusage[mem=8000] \
     -o %J.stdout -eo %J.stderr \
     /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain/smoke_t1t2.sh


Job <198961143> is submitted to queue <gpu>.
bjobs -u $USER -q gpu

more 198961143.stdout

[2354983] single-GPU training
[0] current gpu: 0
[0] Writing Tensorboard logs to /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke/run_brai
nseg_t1t2_smoke_GPU001_D2_H3__08-13-2025-11:11:13
load json keys:  dict_keys(['training', 'validation'])
Training on 1 T1T2 images.
Validation on 0 T1T2 images.
Using generic dataset

#################################
########### TRY 2 ###############
cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder

cat > mini_split_with_val.json <<'JSON'
{
  "training": [
    {
      "image": [
        "/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz",
        "/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz"
      ]
    }
  ],
  "validation": [
    {
      "image": [
        "/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz",
        "/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz"
      ]
    }
  ]
}
JSON

#################################
########### SCRIPT ##############
cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain

cat > smoke_t1t2.sh <<'SH'
#!/bin/bash
# --- clean env ---
module purge
module load python/3.10

# --- activate your venv (with CUDA PyTorch installed) ---
source /sc/arion/projects/mscic1/users/hoangd02/venv/brainseg_env/bin/activate

# --- LSF doesn't define SLURM vars; the script expects SLURM_JOB_NAME ---
export SLURM_JOB_NAME="brainseg_t1t2_smoke"

# --- run a 10-step smoke test on GPU ---
python /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain/main_T1T2.py \
  --split_json /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/mini_split_with_val.json \
  --logdir     /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke \
  --num_steps  10 \
  --batch_size 1 \
  --num_workers 0 \
  --modality   T1T2
SH

chmod +x smoke_t1t2.sh


# submit from the same pretrain folder (or use the absolute path to the script)
bsub -P acc_mscic1 -q gpu -W 00:10 \
     -n 4 -R "span[hosts=1]" -R a100 -gpu "num=1" \
     -R rusage[mem=8000] \
     -o %J.stdout -eo %J.stderr \
     ./smoke_t1t2.sh

Job <198961626> is submitted to queue <gpu>.

bjobs -u $USER -q gpu

## troubleshooting 
#checking the dimensions of the image 
python - <<'PY'
import nibabel as nib
p="/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz"
img=nib.load(p)
print("shape:", img.shape)
PY

#shape: (256, 256, 156) #3D, 256 voxels in X, 256 in Y, 156 slices in Z.

##replace AddChanneld with EnsureChannelFirstd by making a backup called utils/data_utils.py.bak
cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain
bsub -P acc_mscic1 -q gpu -W 00:10 \
     -n 4 -R "span[hosts=1]" -R a100 -gpu "num=1" \
     -R rusage[mem=8000] \
     -o %J.stdout -eo %J.stderr \
     ./smoke_t1t2.sh

Job <198963159> is submitted to queue <gpu>.

################ update script again ################ 
cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain

cat > smoke_t1t2.sh <<'SH'
#!/bin/bash
# --- clean env ---
module purge
module load python/3.10

# --- activate your venv (with CUDA PyTorch installed) ---
source /sc/arion/projects/mscic1/users/hoangd02/venv/brainseg_env/bin/activate

# Make sure Python imports your local repo first (optional but nice)
export PYTHONPATH="/sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain:${PYTHONPATH}"

# LSF doesn't define SLURM vars; the script expects SLURM_JOB_NAME
export SLURM_JOB_NAME="brainseg_t1t2_smoke"

python /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain/main_T1T2.py \
  --split_json /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/mini_split_with_val.json \
  --logdir     /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke \
  --num_steps  10 \
  --batch_size 1 \
  --sw_batch_size 1 \
  --num_workers 0 \
  --modality   T1T2 \
  --in_channels 2
SH

chmod +x smoke_t1t2.sh

bsub -P acc_mscic1 -q gpu -W 00:10 \
     -n 4 -R "span[hosts=1]" -R a100 -gpu "num=1" \
     -R rusage[mem=8000] \
     -o %J.stdout -eo %J.stderr \
     ./smoke_t1t2.sh

Job <198964139> is submitted to queue <gpu>.

bjobs -u $USER -q gpu

Job <198965143> is submitted to queue <gpu>.

################################################
################################################
cd /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain

cat > smoke_t1t2.sh <<'SH'
#!/bin/bash
# Clean environment (avoid stray modules)
module purge

# Activate your existing venv (this already has torch+cu121, monai, einops, etc.)
source /sc/arion/projects/mscic1/users/hoangd02/venv/brainseg_env/bin/activate

# Make sure Python can import your local repo first (optional but helpful)
export PYTHONPATH="/sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain:${PYTHONPATH}"

# The script expects a SLURM var in its logdir naming; set a safe value for LSF
export SLURM_JOB_NAME="brainseg_t1t2_smoke"

# ---- Initialize a 1-process distributed group so dist.all_reduce() works ----
# Required env vars for torch.distributed.init_process_group(init_method="env://")
export MASTER_ADDR="127.0.0.1"
export MASTER_PORT="29501"     # any free TCP port is fine
export RANK="0"
export WORLD_SIZE="1"
export LOCAL_RANK="0"

python /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/pretrain/main_T1T2.py \
  --split_json   /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/mini_split_with_val.json \
  --logdir       /sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke \
  --num_steps    10 \
  --batch_size   1 \
  --sw_batch_size 1 \
  --num_workers  0 \
  --modality     T1T2 \
  --in_channels  2 \
  --distributed \
  --local_rank   0
SH

chmod +x smoke_t1t2.sh

bsub -P acc_mscic1 -q gpu -W 00:10 \
     -n 4 -R "span[hosts=1]" -R a100 -gpu "num=1" \
     -R rusage[mem=8000] \
     -o %J.stdout -eo %J.stderr \
     ./smoke_t1t2.sh


Job <198966976> is submitted to queue <gpu>.

bjobs -u $USER -q gpu

#### VIEWING OUTPUTS (tensorboard log directory and model checkpoints)
LOGDIR="/sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke/run_brainseg_t1t2_smoke_GPU001_D2_H3__08-13-2025-12:24:50"

# See what's inside
ls -lah "$LOGDIR"

tensorboard --logdir "/sc/arion/projects/mscic1/results/jolie/LBP/BrainSegFounder/runs/stage1_smoke" \
            --port 6006 --bind_all



/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/



##########################################
################ ARCHIVED ################
##########################################

#check version 
python -V
python -c "import monai; print(monai.__version__)"
1.4.0 #this is not compatible with AddChanneld
pip uninstall monai
pip install monai==1.1.0

#To run with only one modalitiy of images on a single GPU for 100 epochs, run:
python main.py --json "/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.json" --epochs 100 --in_channels 1

#This will: 
# Initialize a Swin‑UNETR encoder with a first layer accepting one channel.
# Use the 3-way SSL loss function during self-supervised pretraining.
# Train on your UKB-derived or similarly structured single-modal MRI dataset for 100 epochs on one GPU.


###Main Code in ~/pretrain/main.py 
# Copyright 2020 - 2022 MONAI Consortium
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Import necessary libraries and modules
import argparse                              # For parsing command-line arguments
import os                                    # For interacting with the operating system
from time import time                        # To measure time intervals
from datetime import timedelta               # To display elapsed time in readable format
import numpy as np                           # For numerical computations
import torch                                 # PyTorch core library
import torch.distributed as dist             # For distributed (multi-GPU/machine) training
import torch.optim as optim                  # For model optimizers like Adam, SGD
from losses.loss import Loss                 # Custom loss function from the codebase
from models.ssl_head import SSLHead          # Model architecture for self-supervised learning
from optimizers.lr_scheduler import WarmupCosineSchedule  # Custom learning rate scheduler
from torch.cuda.amp import GradScaler, autocast  # For mixed-precision training
from torch.nn.parallel import DistributedDataParallel     # For parallelizing training across GPUs
from torch.utils.tensorboard import SummaryWriter         # For logging training metrics to TensorBoard
from utils.data_utils import get_loader      # Custom data loader function
from utils.ops import aug_rand, rot_rand     # Data augmentation functions

# Main training function
def main():

    # Function to save model checkpoints to a given directory
    def save_ckp(state, checkpoint_dir):
        torch.save(state, checkpoint_dir)

    # Training function for one epoch
    def train(args, global_step, train_loader, val_best, scaler, count_epoch):
        model.train()                                # Set model to training mode
        loss_train = []                              # List to record overall loss values
        loss_train_recon = []                        # List to record reconstruction loss values

        for step, batch in enumerate(train_loader):  # Loop over training batches
            t1 = time()                              # Start timer for current step
            x = batch["image"].to(args.device)       # Move batch of images to GPU or CPU

            # Generate two randomly rotated versions of input x
            x1, rot1 = rot_rand(args, x)             # First rotated view and its label
            x2, rot2 = rot_rand(args, x)             # Second rotated view and its label

            # Apply random augmentations to each rotated view
            x1_augment = aug_rand(args, x1)
            x2_augment = aug_rand(args, x2)
            x1_augment = x1_augment                   # Redundant (can be removed)
            x2_augment = x2_augment                   # Redundant (can be removed)

            # Debug print: only on first global step and main process (rank 0)
            if global_step <= 1 and args.rank == 0:
                print("x:", x.size())
                print("x1 : ", x1.size(), " rot1 : ", rot1.size())
                print("x2 : ", x2.size(), " rot2 : ", rot2.size())
                print("x1_augment:", x1_augment.size())
                print("x2_augment:", x2_augment.size())

            # Run model forward pass with mixed precision if enabled
            with autocast(enabled=args.amp):
                rot1_p, contrastive1_p, rec_x1 = model(x1_augment)  # Outputs for first view
                rot2_p, contrastive2_p, rec_x2 = model(x2_augment)  # Outputs for second view

                # Combine outputs and targets from both views
                rot_p = torch.cat([rot1_p, rot2_p], dim=0)          # Predicted rotation labels
                rots = torch.cat([rot1, rot2], dim=0)               # True rotation labels
                imgs_recon = torch.cat([rec_x1, rec_x2], dim=0)     # Reconstructed images
                imgs = torch.cat([x1, x2], dim=0)                   # Original input images

                # Calculate overall loss and task-specific losses
                loss, losses_tasks = loss_function(rot_p, rots, contrastive1_p, contrastive2_p, imgs_recon, imgs)

            # Record current loss values
            loss_train.append(loss.item())
            loss_train_recon.append(losses_tasks[2].item())  # Typically the reconstruction loss

            # Backward pass and optimizer step using AMP if enabled
            if args.amp:
                scaler.scale(loss).backward()        # Backward pass with scaled loss
                scaler.step(optimizer)               # Optimizer step
                scaler.update()                      # Update scaler for next iteration
            else:
                loss.backward()                      # Standard backward pass
                if args.grad_clip:                   # Optional gradient clipping
                    torch.nn.utils.clip_grad_norm_(model.parameters(), args.max_grad_norm)
                optimizer.step()                     # Optimizer step

            if args.lrdecay:
                scheduler.step()                     # Update learning rate if decay is used

            optimizer.zero_grad()                    # Reset gradients to zero

            # if args.distributed:
            #     if args.rank == 0:
            #         print("Step:{}/{}, Loss:{:.4f}, Time:{:.4f}".format(global_step, args.num_steps, loss, time() - t1))
            # else:
            #     print("Step:{}/{}, Loss:{:.4f}, Time:{:.4f}".format(global_step, args.num_steps, loss, time() - t1))
            # if args.rank == 0:
            if args.rank == 0:
                # Print training status info if on the main process (rank 0)
                print(f"[{args.rank}] " + f"train: " +
                        f"epoch {count_epoch}/{args.epochs - 1}, " +    # Show current epoch out of total
                        f"step_within_epoch {step}/{len(train_loader) - 1}, " +  # Show current batch step
                        f"global_step {global_step}/{args.num_steps - 1}, " +   # Global step counter
                        f"loss: {loss.item():.4f}, " +                           # Current loss value
                        f"time: {(time() - t1):.2f}s")                          # Time taken for this step

            global_step += 1  # Increment global step counter

            # Set validation condition
            val_cond = False                 # Initialize validation flag
            if global_step % args.eval_num == 0:  # Check if it's time to run validation
                val_cond = True

            if val_cond:  # If it's time to validate
                val_loss, val_loss_recon, img_list = validation(args, test_loader, count_epoch, global_step)

            if args.rank == 0 and val_cond:  # Only main process writes logs or saves model
                # Log validation and training loss to TensorBoard
                writer.add_scalar("Validation/loss_recon", scalar_value=val_loss_recon, global_step=global_step)
                writer.add_scalar("train/loss_total", scalar_value=np.mean(loss_train), global_step=global_step)
                writer.add_scalar("train/loss_recon", scalar_value=np.mean(loss_train_recon), global_step=global_step)

                # Log visual validation results
                writer.add_image("Validation/x1_gt", img_list[0], global_step, dataformats="HW")
                writer.add_image("Validation/x1_aug", img_list[1], global_step, dataformats="HW")
                writer.add_image("Validation/x1_recon", img_list[2], global_step, dataformats="HW")

                # Save model if validation loss improved
                if val_loss_recon < val_best:
                    val_best = val_loss_recon   # Update best validation loss
                    checkpoint = {
                        "global_step": global_step,
                        "epoch": count_epoch,
                        "state_dict": model.state_dict(),       # Save model weights
                        "optimizer": optimizer.state_dict(),    # Save optimizer state
                    }
                    save_ckp(checkpoint, args.logdir + "model_bestValRMSE.pt")  # Save to file

                    print(f"[{args.rank}] " + "train: Model was saved! " +
                        f"Best Recon. Val Loss {val_best:.4f} " +  
                        f"Recon. Val Loss {val_loss_recon:.4f}"
                    )                     
                else:
                    # Print message when model is not saved
                    print(f"[{args.rank}] " + "train: Model was not saved! " +
                        f"Best Recon. Val Loss {val_best:.4f} " +  
                        f"Recon. Val Loss {val_loss_recon:.4f}"
                    )                      

        return global_step, loss, val_best  # Return updated training state

    # Validation function
    def validation(args, test_loader, count_epoch, global_step):
        model.eval()                      # Set model to evaluation mode
        loss_val = []                     # Store total validation loss
        loss_val_recon = []               # Store reconstruction loss

        with torch.no_grad():             # Disable gradient calculation for validation
            for step, batch in enumerate(test_loader):  # Loop over validation batches
                val_inputs = batch["image"].to(args.device)  # Move batch to GPU or CPU

                # Generate augmented rotated views
                x1, rot1 = rot_rand(args, val_inputs)
                x2, rot2 = rot_rand(args, val_inputs)
                x1_augment = aug_rand(args, x1)
                x2_augment = aug_rand(args, x2)

                # Model forward pass (no gradients needed)
                with autocast(enabled=args.amp):
                    rot1_p, contrastive1_p, rec_x1 = model(x1_augment)
                    rot2_p, contrastive2_p, rec_x2 = model(x2_augment)
                    rot_p = torch.cat([rot1_p, rot2_p], dim=0)           # Combined predictions
                    rots = torch.cat([rot1, rot2], dim=0)               # Combined true labels
                    imgs_recon = torch.cat([rec_x1, rec_x2], dim=0)     # Combined reconstructions
                    imgs = torch.cat([x1, x2], dim=0)                   # Original inputs
                    loss, losses_tasks = loss_function(rot_p, rots, contrastive1_p, contrastive2_p, imgs_recon, imgs)

                loss_recon = losses_tasks[2]  # Get reconstruction loss
                loss_val.append(loss)         # Store total loss
                loss_val_recon.append(loss_recon)  # Store recon loss

                # Convert tensors to numpy images for logging
                x_gt = x1.detach().cpu().numpy()
                x_gt = (x_gt - np.min(x_gt)) / (np.max(x_gt) - np.min(x_gt))  # Normalize to [0, 1]
                xgt = x_gt[0][0][:, :, 48] * 255.0                             # Slice and scale to [0, 255]
                xgt = xgt.astype(np.uint8)

                x1_augment = x1_augment.detach().cpu().numpy()
                x1_augment = (x1_augment - np.min(x1_augment)) / (np.max(x1_augment) - np.min(x1_augment))
                x_aug = x1_augment[0][0][:, :, 48] * 255.0
                x_aug = x_aug.astype(np.uint8)

                rec_x1 = rec_x1.detach().cpu().numpy()
                rec_x1 = (rec_x1 - np.min(rec_x1)) / (np.max(rec_x1) - np.min(rec_x1))
                recon = rec_x1[0][0][:, :, 48] * 255.0
                recon = recon.astype(np.uint8)

                img_list = [xgt, x_aug, recon]  # Store images for logging

                # Print validation step information
                print(f"[{args.rank}] " + "validation: " +
                      f"epoch {count_epoch}/{args.epochs - 1}, " +  
                      f"global_step {global_step}/{args.num_steps - 1}, " +  
                      f"Validation step {step}/{len(test_loader)}, " +  
                      f"Loss {loss.item():.4f}, " +  
                      f"Loss Reconstruction {loss_recon.item():.4f}"
                )

        # YY
        if args.device != torch.device("cpu"):
            torch.cuda.synchronize(args.device)  # Make sure all CUDA ops are finished (for accurate timing/logging)

        # Combine all validation losses from different processes (if distributed)
        loss_val_mean = torch.sum(torch.stack(loss_val), dim=0)          # Sum all batch losses
        loss_val_recon_mean = torch.sum(torch.stack(loss_val_recon), dim=0)  # Sum all reconstruction losses

        # Use distributed reduction to sum losses across all GPUs/processes
        dist.all_reduce(loss_val_mean, op=torch.distributed.ReduceOp.SUM)         # Aggregate loss across all ranks
        dist.all_reduce(loss_val_recon_mean, op=torch.distributed.ReduceOp.SUM)   # Aggregate recon loss across all ranks

        if args.rank == 0:
            # Print the combined loss only on the main process
            print(f"Rank {args.rank}, loss_val_mean={loss_val_mean.item()}, loss_val_recon_mean={loss_val_recon_mean.item()}")

        return loss_val_mean.item(), loss_val_recon_mean.item(), img_list  # Return final loss and images

        # return np.mean(loss_val), np.mean(loss_val_recon), img_list  # (old version that used local mean only)

    parser = argparse.ArgumentParser(description="PyTorch Training")  # Setup argument parser with description
    # Add all training options below — can be set via CLI or defaults will be used
    parser.add_argument("--resume", default=None, type=str, help="resume training")    
    parser.add_argument("--logdir", default="/mnt/runs", type=str, help="directory to save the tensorboard logs")
    parser.add_argument("--workdir", default="/mnt", type=str, help="root of working directory")
    parser.add_argument("--epochs", default=100, type=int, help="number of training epochs")
    parser.add_argument("--num_steps", default=100000, type=int, help="number of training iterations")
    parser.add_argument("--eval_num", default=10, type=int, help="evaluation frequency")
    parser.add_argument("--warmup_steps", default=500, type=int, help="warmup steps")
    parser.add_argument('--num_workers', default=4, type=int, help='number of workers')
    parser.add_argument("--in_channels", default=1, type=int, help="number of input channels")
    parser.add_argument("--feature_size", default=48, type=int, help="embedding size")
    parser.add_argument("--dropout_path_rate", default=0.0, type=float, help="drop path rate")
    parser.add_argument("--use_checkpoint", action="store_true", help="use gradient checkpointing to save memory")
    parser.add_argument("--spatial_dims", default=3, type=int, help="spatial dimension of input data")
    parser.add_argument("--a_min", default=-1000, type=float, help="a_min in ScaleIntensityRanged")
    parser.add_argument("--a_max", default=1000, type=float, help="a_max in ScaleIntensityRanged")
    parser.add_argument("--b_min", default=0.0, type=float, help="b_min in ScaleIntensityRanged")
    parser.add_argument("--b_max", default=1.0, type=float, help="b_max in ScaleIntensityRanged")
    parser.add_argument("--space_x", default=1.5, type=float, help="spacing in x direction")
    parser.add_argument("--space_y", default=1.5, type=float, help="spacing in y direction")
    parser.add_argument("--space_z", default=2.0, type=float, help="spacing in z direction")
    parser.add_argument("--roi_x", default=96, type=int, help="roi size in x direction")
    parser.add_argument("--roi_y", default=96, type=int, help="roi size in y direction")
    parser.add_argument("--roi_z", default=96, type=int, help="roi size in z direction")
    parser.add_argument("--batch_size", default=2, type=int, help="number of batch size")
    parser.add_argument("--sw_batch_size", default=2, type=int, help="number of sliding window batch size")
    parser.add_argument("--lr", default=4e-4, type=float, help="learning rate")
    parser.add_argument("--decay", default=0.1, type=float, help="decay rate")
    parser.add_argument("--momentum", default=0.9, type=float, help="momentum")
    parser.add_argument("--lrdecay", action="store_true", help="enable learning rate decay")
    parser.add_argument("--max_grad_norm", default=1.0, type=float, help="maximum gradient norm")
    parser.add_argument("--loss_type", default="SSL", type=str)
    parser.add_argument("--opt", default="adamw", type=str, help="optimization algorithm")
    parser.add_argument("--lr_schedule", default="warmup_cosine", type=str)
    parser.add_argument("--grad_clip", action="store_true", help="gradient clip")
    parser.add_argument("--noamp", action="store_true", help="do NOT use amp for training")
    parser.add_argument("--smartcache_dataset", action="store_true", help="use monai smartcache Dataset")
    parser.add_argument("--cache_dataset", action="store_true", help="use monai cache Dataset")
    # parse the command-line argument --local_rank, provided by torch.distributed.launch
    parser.add_argument("--local_rank", type=int, help='provided by torch.distributed.launch')
    parser.add_argument('--distributed', action='store_true', help='start distributed training')

    args = parser.parse_args()           # Parse command-line args
    args.amp = not args.noamp            # Enable AMP by default unless `--noamp` is set

    torch.backends.cudnn.benchmark = True             # Enable cuDNN autotuning for performance
    torch.autograd.set_detect_anomaly(True)           # Enable gradient anomaly detection (helpful for debugging)
    # args.distributed = False
    # if "WORLD_SIZE" in os.environ:
    #     args.distributed = int(os.environ["WORLD_SIZE"]) > 1
    # args.device = "cuda:0"
    # args.world_size = 1
    # args.rank = 0

    # for debugging purpose
    if args.distributed:
        # Enable debugging environment variables
        os.environ["TORCH_DISTRIBUTED_DEBUG"] = "DETAIL"
        os.environ["NCCL_ASYNC_ERROR_HANDLING"] = "1" 

        # Print current environment setup
        env_dict = {
            key: os.environ[key]
            for key in ("MASTER_ADDR", "MASTER_PORT", "RANK", "WORLD_SIZE")
        }
        print(f"[{os.getpid()}] Initializing process group with: {env_dict}")        

        # Initialize distributed backend (NCCL = for multi-GPU)
        dist.init_process_group(backend="nccl", init_method="env://", timeout=timedelta(minutes=10))
        args.world_size = dist.get_world_size()    # Total number of processes (GPUs)
        args.rank = dist.get_rank()                # Rank of this process (e.g., 0, 1, ...)
        args.device = torch.device(f"cuda:{args.local_rank}")  # Assign GPU for this process
    else:
        # Single GPU training setup
        print(f"[{os.getpid()}] single-GPU training")
        args.rank = 0
        args.device = torch.device(f"cuda:{torch.cuda.current_device()}")  # Use current GPU

    assert args.rank >= 0     # Sanity check

    torch.cuda.set_device(args.device)  # Set the device for current process
    print(f"[{args.rank}] current gpu: {torch.cuda.current_device()}")

    ##Loggin step
    from datetime import datetime  # For timestamping logs

    # Create logging directory with timestamp and GPU count
    args.logdir = args.logdir + f"/run_{args.world_size:03d}G_" + datetime.now().strftime("%m-%d-%Y-%H:%M:%S") 

    if args.rank == 0:
        # Only main process creates directory and initializes TensorBoard writer
        os.makedirs(args.logdir, exist_ok=True)
        writer = SummaryWriter(log_dir=args.logdir)
        print(f"[{args.rank}] " + f"Writing Tensorboard logs to {args.logdir}")
    else:
        writer = None  # Other processes don’t need a writer

    model = SSLHead(args)                      # Initialize your model using the custom SSLHead class
    model.to(args.device)                      # Move the model to the specified device (GPU or CPU)

    # Choose optimizer based on user-specified argument
    if args.opt == "adam":
        optimizer = optim.Adam(params=model.parameters(), lr=args.lr, weight_decay=args.decay)

    elif args.opt == "adamw":
        optimizer = optim.AdamW(params=model.parameters(), lr=args.lr, weight_decay=args.decay)

    elif args.opt == "sgd":
        optimizer = optim.SGD(params=model.parameters(), lr=args.lr, momentum=args.momentum, weight_decay=args.decay)

    # Resume training from a checkpoint if provided
    if args.resume is not None:
        try:
            model_dict = torch.load(args.resume)           # Load checkpoint file
            state_dict = model_dict["state_dict"]          # Extract the model weights

            # Fix keys if model was saved using DistributedDataParallel (adds "module.")
            if "module." in list(state_dict.keys())[0]:    
                if args.rank == 0:
                    print(f"[{args.rank}] " + "Tag 'module.' found in state dict - fixing!")
                for key in list(state_dict.keys()):
                    # Replace "module." with "swinViT." to match current model keys
                    state_dict[key.replace("module.", "swinViT.")] = state_dict.pop(key)

            # Load weights into the model — strict=False ignores mismatched keys (e.g., decoder weights)
            model.load_state_dict(state_dict, strict=False)

            # If available, also restore epoch and optimizer state
            if 'epoch' in model_dict:
                model.epoch = model_dict["epoch"]
            if 'optimizer' in model_dict:
                model.optimizer = model_dict["optimizer"]

            if args.rank == 0:
                print(f"[{args.rank}] " + "Using pretrained self-supervised Swin UNETR backbone weights !")

        except ValueError:
            raise ValueError("Self-supervised pre-trained weights not available for " + str(args.model_name))

    # Set up learning rate scheduler if enabled
    if args.lrdecay:
        if args.lr_schedule == "warmup_cosine":
            # Warmup + cosine decay schedule
            scheduler = WarmupCosineSchedule(
                optimizer,
                warmup_steps=args.warmup_steps,
                t_total=args.num_steps
            )

        elif args.lr_schedule == "poly":
            # Polynomial decay learning rate scheduler
            def lambdas(epoch):
                return (1 - float(epoch) / float(args.epochs)) ** 0.9

            scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lambdas)

    # Initialize loss function with batch size (used in some contrastive losses)
    loss_function = Loss(args.batch_size * args.sw_batch_size, args)

    # Convert BatchNorm to SyncBatchNorm for multi-GPU (ensures consistent stats)
    if args.distributed:
        model = torch.nn.SyncBatchNorm.convert_sync_batchnorm(model)   # Optional but recommended
        model = DistributedDataParallel(model, device_ids=[args.local_rank])  # Wrap model for multi-GPU training

    # Load training and validation data loaders using custom get_loader()
    train_loader, test_loader = get_loader(args)

    # Calculate total number of trainable parameters in the model
    pytorch_total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    # Print parameter count (only on main process)
    if args.rank == 0:
        print(f"[{args.rank}] " + f"Total parameters count: {pytorch_total_params}")


    global_step = 0               # Initialize the global step counter
    best_val = 1e8                # Initialize best validation loss with a large number (to be improved)

    # Set up automatic mixed precision (AMP) scaler if AMP is enabled
    if args.amp:
        scaler = GradScaler()     # Enables mixed precision training (saves memory, speeds up)
    else:
        scaler = None             # Disable AMP

    count_epoch = 0               # Start counting from epoch 0

    # Start training loop
    while global_step < args.num_steps:
        # Train one epoch (or partial epoch depending on batch size and num_steps)
        global_step, loss, best_val = train(args, global_step, train_loader, best_val, scaler, count_epoch)
        count_epoch = count_epoch + 1   # Move to next epoch

    # After training finishes, prepare final checkpoint to save
    checkpoint = {
        "global_step": global_step,             # Save final global step
        "epoch": count_epoch,                   # Save number of completed epochs
        "state_dict": model.state_dict(),       # Save model weights
        "optimizer": optimizer.state_dict(),    # Save optimizer state
    }

    # Save final model only from main process
    if args.rank == 0:
        print(f"[{args.rank}] " + f"Training Finished! Best val: {best_val}")
        save_ckp(checkpoint, args.logdir + "model_final_epoch.pt")  # Save full checkpoint
        print(f"[{args.rank}] " + "Saved model_final_epoch.pt")

    # If using distributed training, clean up the process group
    if args.distributed:
        dist.destroy_process_group()  # Release resources for distributed training

# Standard Python entry point
if __name__ == "__main__":
    main()  # Call the main() function to start training

























