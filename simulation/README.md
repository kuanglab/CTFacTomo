1. Run the tensor model on 1D spatial gene expression data (tomo-seq), protein-protein interaction (PPI) network and spatial chain graphs, and 3D binary mask (where 0s indicate no tissue covered).

```python
# Load required packages
import torch
import numpy as np
import tensorly as tl
from TFacTomo import reconstruct
# Set tensorly backend as pytorch
tl.set_backend('pytorch')

# Load 1D gene expression data along different spatial axes
X_x = np.load("data/mouse_olfactory_mucosa/normalized_fitted_lml_data.npy"); X_x = torch.from_numpy(X_x).to(torch.float)
X_y = np.load("data/mouse_olfactory_mucosa/normalized_fitted_dv_data.npy"); X_y = torch.from_numpy(X_y).to(torch.float)
X_z = np.load("data/mouse_olfactory_mucosa/normalized_fitted_ap_data.npy"); X_z = torch.from_numpy(X_z).to(torch.float)
# Load knowledge graphs along gene and spatial axes
W_g = np.load("data/mouse_olfactory_mucosa/mus_musculus_ppi_adjacency_matrix.npy"); W_g = torch.from_numpy(W_g).to(torch.float)
W_x = np.load("data/mouse_olfactory_mucosa/W_x.npy"); W_x = torch.from_numpy(W_x).to(torch.float)
W_y = np.load("data/mouse_olfactory_mucosa/W_y.npy"); W_y = torch.from_numpy(W_y).to(torch.float)
W_z = np.load("data/mouse_olfactory_mucosa/W_z.npy"); W_z = torch.from_numpy(W_z).to(torch.float)
# Load the 3D binary mask 
M = np.load("data/mouse_olfactory_mucosa/mouse_olfactory_mucosa_mask.npy"); M = torch.from_numpy(M).to(torch.float)
