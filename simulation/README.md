Run the tensor model on projected spatial gene expression data (simulation), protein-protein interaction (PPI) network, spatial chain graphs, and 3D binary mask (where 0s indicate no tissue covered).

#### Drosophila
```python
# Load required packages
import torch
import numpy as np
import tensorly as tl
from IPF import IPF 
from CTFacTomo import reconstruct

# Set tensorly backend as pytorch
tl.set_backend('pytorch')

# Load 1D gene expression data along different spatial axes
data = torch.load('data/Drosophila_L2_simulated_data_63x53x21.pt')
M = data['M']; 
X_x = data['X_x']; X_x = X_x.to(torch.float32); X_y = data['X_y']; X_y = X_y.to(torch.float32); X_z = data['X_z']; X_z = X_z.to(torch.float32)
W_g = data['W_g']; W_g = W_g.to(torch.float32); W_x = data['W_x']; W_x = W_x.to(torch.float32); 
W_y = data['W_y']; W_y = W_y.to(torch.float32); W_z = data['W_z']; W_z = W_z.to(torch.float32);

X = [X_x, X_y, X_z]
W = [W_g, W_x, W_y, W_z]

# Reconstruct 4D expression tensor by setting hyperparameters rank=180, alpha=1e3, beta=1e-3, and lambda=1
A = reconstruct(X, 1-M, W, 180, 1e3, 1e-3, 1, reduction='mean', stop_crit=1e-10, max_epoch=501, verbose=True, freq=10)
T_hat = tl.cp_tensor.cp_to_tensor((None, A))

# Reconstruct 4D expression tensor by using IPF
T_hat = IPF(X, M, max_epoch=501, reduction='mean', normalized=False, verbose=True, freq=10)
```
#### Human Heart and Mouse Brain 
```Matlab
% Please download [MATLAB Tensor Toolbox](https://gitlab.com/tensors/tensor_toolbox) package and add folders to the search path
addpath('tensor_toolbox');

% Human Heart
% Reconstruct 4D expression tensor by setting hyperparameters rank=99, alpha=1e3, beta=1e-3
load('Human_simulated_data_20x20x9.mat');
opts.rank = 99; opts.alpha = 1e4; opts.beta = 1e-2; opts.lambda = 1;
opts.stopcrit = 1e-2; opts.maxiters = 200; opts.miniters = 50;
A = reconstruct(T, X, M, W, opts);
T_hat = tensor(ktensor(A));
% Reconstruct 4D expression tensor by using IPF
opts.maxiters = 200;
T_hat = IPF(X{2}, X{3}, X{4}, M, opts.maxiters);

% Mouse Brain
% Reconstruct 4D expression tensor by setting hyperparameters rank=99, alpha=1e3, beta=1e-3
load('Mouse_simulated_data_30x30x30.mat');
opts.rank = 120; opts.alpha = 1e2; opts.beta = 1e-3; opts.lambda = 1;
opts.stopcrit = 1e-2; opts.maxiters = 200; opts.miniters = 50;
A = reconstruct(T, X, M, W, opts);
T_hat = tensor(ktensor(A));
% Reconstruct 4D expression tensor by using IPF
opts.maxiters = 200;
T_hat = IPF(X{2}, X{3}, X{4}, M, opts.maxiters);
```

#### For reconstructing 4D expression tensor with Tomographer, please refer to the original package [Tomographer](https://github.com/lamanno-epfl/tomographer).
