#UniverSeg
#https://universeg.csail.mit.edu/
#https://github.com/JJGO/UniverSeg
#https://openaccess.thecvf.com/content/ICCV2023/html/Butoi_UniverSeg_Universal_Medical_Image_Segmentation_ICCV_2023_paper.html

## in terminal
python -m venv universeg_env
source universeg_env/bin/activate

module load git

pip install git+https://github.com/JJGO/UniverSeg.git

python
## in python
from universeg import universeg

model = universeg(pretrained=True)

#raw brains are in 
/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/

#UniverSeg expects:
#2D grayscale images (shape: 1×128×128).
#Images normalized to [0, 1].
#Tensors (so we need to go from .nii.gz → NumPy array → torch tensor).

import nibabel as nib
import numpy as np
import torch
import matplotlib.pyplot as plt
from skimage.transform import resize

# Load the MRI
nii = nib.load("/sc/arion/projects/psychgen/lbp/data/neuroimaging/LBP_NEUROIMAGING_PIPELINE_RUN_17JUN2022/output/scan31/scan31.nii.gz")
data = nii.get_fdata()
##....shape=(256, 256, 156))

# Pick a middle slice (example: axial)
mid_slice = data.shape[2] // 2 #78
img = data[:, :, mid_slice]

#array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
#       [ 0.,  0.,  0., ...,  0.,  0.,  0.],
#       [ 0., 15., 22., ...,  0.,  0.,  0.],
#       ...,
#       [ 0.,  0.,  2., ...,  0.,  0.,  0.],
#       [ 0.,  0.,  0., ...,  0.,  0.,  0.],
#       [ 0.,  0.,  0., ...,  0.,  0.,  0.]], shape=(256, 256))

# Normalize to [0, 1]
img = (img - np.min(img)) / (np.max(img) - np.min(img))

# Resize to 128x128
img_resized = resize(img, (128, 128), anti_aliasing=True)
#.... shape=(128, 128))

# Convert to tensor: shape (1, 1, 128, 128)
img_tensor = torch.tensor(img_resized, dtype=torch.float32).unsqueeze(0).unsqueeze(0)

# Save image for later use
torch.save(img_tensor, "/sc/arion/projects/mscic1/results/jolie/LBP/scan31_tensor.pt")

# Optional: visualize
plt.imshow(img_resized, cmap="gray")
plt.title("Middle Slice Resized to 128x128")
plt.axis("off")
plt.savefig("/sc/arion/projects/mscic1/results/jolie/LBP/scan31_preview.png")
plt.show()

prediction = model(
    target_image,       # e.g., shape (B, 1, 128, 128)
    support_images,     # shape (B, S, 1, 128, 128)
    support_labels      # shape (B, S, 1, 128, 128)
)
# prediction -> (B, 1, 128, 128)

python preprocess_scan31.py

