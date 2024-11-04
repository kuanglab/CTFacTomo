# TFacTomo
Reconstructing 3D Spatial Structures of RNA-Tomography Transcriptomes by Tensor Factorization

![](https://github.com/kuanglab/TFacTomo/blob/main/figures/TFacTomo_Workflow.png)

System and package requirements
------------------------------------------------------------

The code was developed and tested on a node in the cluster using a Linux system machine (Ubuntu 20.04.5 LTS). The node has 16 cores (AMD EPYC 7763 CPU), 60GB memory, and one NVIDIA A100 GPU with 40GB memory.

The Python packages can be downloaded and run with the following library versions.
```
[python 3.8.12]
[numpy 1.21.5]
[pytorch 1.10.2]
[tensorly 0.6.0]
```

The R packages can be downloaded and run with the following library versions (visualization only).
```
[R 4.2.2]
[shiny 1.7.1]
[reticulate 1.25]
```

Reproducing Results on Mouse Olfactory Mucosa
------------------------------------------------------------
#### Reconstruction

*You can skip this step if you just want to use the tensor models already in our repository.*

*They are stored in vis_OM/tensor_model/\*.npy*

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

# Reconstruct 4D expression tensor by setting hyperparameters rank=500, alpha=1e2, beta=1, and lambda=1, where lambda is a
# hyperparameter controlling l2 regularization on tensor factors, we kept using 1 in our experiments reported in the paper.
A = reconstruct([X_x, X_y, X_z], 1-M, [W_g, W_x, W_y, W_z], 500, 1e2, 1, 1, stop_crit=1e-4, reduction="sum", max_epoch=1000)
T_hat = tl.cp_tensor.cp_to_tensor((None, A))
torch.save(A, "results/mouse_olfactory_mucosa_CPD_factors.pt")
```

2. Save CANDECOMP/PARAFAC (CP) decomposition factors (A_g, A_x, A_y, A_z) for visualization.

```python
A = torch.load("results/mouse_olfactory_mucosa_CPD_factors.pt")
np.save("A_g.npy", A[0])
np.save("A_x.npy", A[1])
np.save("A_y.npy", A[2])
np.save("A_z.npy", A[3])
```

You are done with the reconstruction part of this tutorial. Time to move onto visualization!

#### Visualization

*The provided visualization tool is a modified version of the one in the original Mouse OM paper[1]*

```r
# Load required packages
library(shiny)
library(reticulate)
source("vis_OM/server.R")
source("vis_OM/ui.R")

# Load in environment variables using OM_TFacTomo.RData
load("vis_OM/OM_TFacTomo.RData")

# Load in the gene symbols and CP factors (A_g, A_x, A_y, A_z)
np <- import("numpy") # using the reticulate package to load Python objects into R
genes <- np$loadtxt("genes.txt") 
A_g <- np$load("vis_OM/tensor_model/A_g.npy")
A_x <- np$load("vis_OM/tensor_model/A_x.npy")
A_y <- np$load("vis_OM/tensor_model/A_y.npy")
A_z <- np$load("vis_OM/tensor_model/A_z.npy")

# Launch R shiny App to see the 3D expression reconstruction interactively  
shinyApp(ui, server)
```

#### References
_________
1. Ruiz Tejada Segura, M.L. et al. (2022) A 3D transcriptomics atlas of the mouse nose sheds light on the anatomical logic of smell. Cell Reports, 38, 110547.


(Under construction ...)
