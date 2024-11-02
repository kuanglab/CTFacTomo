# TFacTomo
Reconstructing 3D Spatial Structures of RNA-Tomography Transcriptomes by Tensor Factorization

![](https://github.com/kuanglab/TFacTomo/blob/main/figures/TFacTomo_Workflow.png)

Reproducing Results on Mouse Olfactory Mucosa
------------------------------------------------------------
#### Reconstruction

* You can skip this step if you just want to use the tensor models already in our repository.

1. Run the tensor model on 1D expression data, PPI network and spatial graphs, and 3D binary mask (where 1s represent no tissue).

```python
if __name__ == "__main__": 
  X_x = torch.from_numpy(np.load("data/mouse_olfactory_bulb/normalized_fitted_lml_data.npy")).to(torch.float)
  X_y = torch.from_numpy(np.load("data/mouse_olfactory_bulb/normalized_fitted_dv_data.npy")).to(torch.float)
  X_z = torch.from_numpy(np.load("data/mouse_olfactory_bulb/normalized_fitted_ap_data.npy")).to(torch.float)
  W_g = torch.from_numpy(np.load("data/mouse_olfactory_bulb/expanded_ppi_adjacency_list_diagonal_filled.npy")) .to(torch.float)
  W_x = torch.from_numpy(np.load("data/mouse_olfactory_bulb/W_x.npy")).to(torch.float)
  W_y = torch.from_numpy(np.load("data/mouse_olfactory_bulb/W_y.npy")).to(torch.float)
  W_z = torch.from_numpy(np.load("data/mouse_olfactory_bulb/W_z.npy")).to(torch.float)
  M = torch.from_numpy(np.swapaxes(np.load("data/mouse_olfactory_bulb/mouse_olfactory_bulb_mask.npy"), 0, 1)).to(torch.float)
  
  A = reconstruct([X_x, X_y, X_z], 1-M, [W_g, W_x, W_y, W_z], 500, 100, 1, 1, stop_crit=0.0001, reduction="sum",  max_epoch=500)
  torch.save(A, "mouse_olfactory_bulb_results.pt")
```

2. Save your resulting 2D matrices (Ag, Ax, Ay, Az) as numpy arrays.

```python
matrices = torch.load("mouse_om_tensor.pt")
np.save("Ag.npy", matrices[0])
np.save("Ax.npy", matrices[1])
np.save("Ay.npy", matrices[2])
np.save("Az.npy", matrices[3])
```

You are done with the reconstruction part of this tutorial. Time to move onto visualization!

#### Visualization

*The provided visualization tool is a modified version of the one in the original Mouse OM paper [1]*

1. Install the reticulate package in R
2. Load in the library 
3. Load in the genes
4. Load in the 2D Matrices (Ag, Ax, Ay, Az)

```r
# Assuming your working directory is the vis_OM folder
install.packages("reticulate")
library(reticulate)
np <- import("numpy")
genes <- np$loadtxt("ggenes.txt") 
A_g <- np$load("Ag.npy")
A_x <- np$load("Ax.npy")
A_y <- np$load("Ay.npy")
A_z <- np$load("Az.npy")
```

References
_________
1. Ruiz Tejada Segura,M.L. et al. (2022) A 3D transcriptomics atlas of the mouse nose sheds light on the anatomical logic of smell. Cell Reports, 38, 110547.


(Under construction ...)
