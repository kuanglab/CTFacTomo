import os
import gc
import torch
import numpy as np
import tensorly as tl
from functools import reduce

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
tl.set_backend('pytorch')


def reconstruct(X, M, max_epoch=201, reduction='sum', normalized=True, verbose=True, freq=10):
    
    X_x = X[0]; X_y = X[1]; X_z = X[2]
    n_g, n_x = X_x.size(); n_y = X_y.size()[1]; n_z = X_z.size()[1]
    
    if reduction == 'sum' and normalized:
        
        X_x = X_x / X_x.sum(dim=0) * M.squeeze().sum(dim=(1,2))
        X_y = X_y / X_y.sum(dim=0) * M.squeeze().sum(dim=(0,2))
        X_z = X_z / X_z.sum(dim=0) * M.squeeze().sum(dim=(0,1))
        
        X_x[torch.where(torch.isnan(X_x))] = 0
        X_y[torch.where(torch.isnan(X_y))] = 0
        X_z[torch.where(torch.isnan(X_z))] = 0
        
        mean = torch.mean(torch.stack((X_x.sum(), X_y.sum(), X_z.sum())))
        
        X_x = X_x / X_x.sum(dim=1)[:, None] * mean
        X_y = X_y / X_y.sum(dim=1)[:, None] * mean
        X_z = X_z / X_z.sum(dim=1)[:, None] * mean
        
    T = M.repeat(n_g, 1, 1, 1)
    
    for epoch in range(1, max_epoch):
        
        if reduction == 'sum':
            X_x_ = T.sum(dim=(2, 3))
        else:
            X_x_ = T.mean(dim=(2, 3))
        T = T * torch.permute((X_x/X_x_).repeat(n_y, n_z, 1, 1), (2, 3, 0, 1))
        T[torch.where(torch.isnan(T))] = 0; T[torch.where(torch.isinf(T))] = 0
        
        if reduction == 'sum':
            X_y_ = T.sum(dim=(1, 3))
        else:
            X_y_ = T.mean(dim=(1, 3))
        T = T * torch.permute((X_y/X_y_).repeat(n_x, n_z, 1, 1), (2, 0, 3, 1))
        T[torch.where(torch.isnan(T))] = 0; T[torch.where(torch.isinf(T))] = 0
        
        if reduction == 'sum':
            X_z_ = T.sum(dim=(1, 2))
        else:
            X_z_ = T.mean(dim=(1, 2))
        T = T * torch.permute((X_z/X_z_).repeat(n_x, n_y, 1, 1), (2, 0, 1, 3))
        T[torch.where(torch.isnan(T))] = 0; T[torch.where(torch.isinf(T))] = 0
        
        if verbose and epoch % freq == 0:
            if reduction == 'sum':
                res = torch.mean((T.sum(dim=(2, 3)) - X_x) ** 2) + torch.mean((T.sum(dim=(1, 3)) - X_y) ** 2) + torch.mean((T.sum(dim=(1, 2)) - X_z) ** 2)
            else:
                res = torch.mean((T.mean(dim=(2, 3)) - X_x) ** 2) + torch.mean((T.mean(dim=(1, 3)) - X_y) ** 2) + torch.mean((T.mean(dim=(1, 2)) - X_z) ** 2)
#         res = torch.mean((T-T_)**2)
            print(f'Training epoch: {epoch:02d}, Residual: {res:.10f}')
    
    return T