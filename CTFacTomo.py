import os
import gc
import torch
import numpy as np
import tensorly as tl
from functools import reduce

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
tl.set_backend('pytorch')

# Auxiliary function to collapse all components
def collapse(A):
    
    a = [A[i].sum(axis=0).unsqueeze(0) for i in range(len(A))]
    
    return a

# Compute residuals among old and new components
def compute_res(A, A_hat):
    res_num = 0
    res_denom = 0
    for i in range(len(A)):
        res_num += torch.sum((A[i] - A_hat[i]) ** 2)
        res_denom += torch.sum(A_hat[i] ** 2)
    return float(torch.sqrt(res_num / res_denom))

# Define the model
def reconstruct(X, M, W, rank, alpha, beta, l, stop_crit=1e-3, 
         max_epoch=501, reduction='sum', verbose=True, freq=10):
    
    n = np.array([W_.shape[0] for W_ in W])
    n_modes = len(n)
    eps = 1e-10
    
    # Normalize adjacency matrix
    D_ = []
    W_ = W.copy()
    for i in range(n_modes):
        W_[i] = W_[i] - torch.diag(torch.diag(W_[i]))
        d = torch.sum(W_[i], 0)
        nz_index = torch.where(d != 0)
        d[nz_index] = d[nz_index] ** (-0.5)
        d = d.expand([len(d), -1])
        W_[i] = d.T * W_[i] * d
        D_.append(torch.diag(torch.sum(W_[i], 0)).to_sparse().type(torch.float32))
        W_[i] = W_[i].to_sparse().type(torch.float32)
    
    # Initialization
    A = []
    ATA = []; aaT = []
    WA = []; DA = []
    ATWA = []; ATDA = []
    
    for i in range(n_modes):
        A.append(torch.rand(n[i], rank))
        ATA.append(A[i].T @ A[i])
        aaT.append(torch.outer(A[i].sum(axis=0), A[i].sum(axis=0)))
        WA.append(W_[i] @ A[i])
        DA.append(D_[i] @ A[i])
        ATWA.append(A[i].T @ WA[i])
        ATDA.append(A[i].T @ DA[i])
    
    # Convert expression matrix to tensor
    X_ = X.copy()
    X_ = [None] + X_
    for i in range(1, n_modes):
        X_[i] = X_[i][(..., ) + (None, )*2]
        X_[i] = X_[i].swapaxes(1, i)
    
    for epoch in range(1, max_epoch):
        
        A_hat = A.copy()
        
        # Update gene component and spatial components in x, y ,z order
        for i in range(n_modes):
            
            a = collapse(A)
            
            # gene component
            if i == 0:
                
                J1n = 0
                for j in range(1, n_modes):
                    a_ = a.copy()
                    a_[0] = A[0]; a_[j] = A[j]
                    theta = tl.cp_tensor.unfolding_dot_khatri_rao(X_[j], (None, a_), 0)
                    if reduction == 'mean':
                        J1n = J1n + theta / np.prod(np.delete(n, [0, j]))
                    else: # sum
                        J1n = J1n + theta 
                J1p = 0
                for j in range(1, n_modes):
                    phi = aaT.copy()
                    phi[j] = ATA[j]; del phi[0]
                    if reduction == 'mean':
                        J1p = J1p + A[0] @ reduce(torch.mul, phi) / (np.prod(np.delete(n, [0, j])) ** 2)
                    else: # sum
                        J1p = J1p + A[0] @ reduce(torch.mul, phi)
        
                A_ = A.copy()
                A_[0] = a[0]
                J2p = tl.cp_tensor.unfolding_dot_khatri_rao(tl.cp_tensor.cp_to_tensor((None, A_), M),
                                                            (None, A_), 0)
                if reduction == 'mean':
                    J2p = J2p.repeat(n[0], 1) / (n[0] ** 2)
                else:
                    J2p = J2p.repeat(n[0], 1)
            
            # spatial component
            else:
                
                J1n = 0
                for j in range(1, n_modes):
                    a_ = a.copy()
                    a_[0] = A[0]; a_[j] = A[j]
                    theta = tl.cp_tensor.unfolding_dot_khatri_rao(X_[j], (None, a_), i)
                    if j != i:
                        theta = theta.repeat(n[i], 1)
                    if reduction == 'mean':
                        J1n = J1n + theta / np.prod(np.delete(n, [0, j]))
                    else: # sum
                        J1n = J1n + theta
                
                J1p = 0
                phi = aaT.copy()
                phi[0] = ATA[0]; del phi[i]
                if reduction == 'mean':
                    J1p = J1p + A[i] @ reduce(torch.mul, phi) / (np.prod(np.delete(n, [0, i])) ** 2)
                else: # sum
                    J1p = J1p + A[i] @ reduce(torch.mul, phi)

                phi_sum = 0
                for j in range(1, n_modes):
                    if j != i:
                        phi = aaT.copy()
                        phi[0] = ATA[0]; phi[j] = ATA[j]; del phi[i]
                        if reduction == 'mean':
                            phi_sum = phi_sum + reduce(torch.mul, phi) / (np.prod(np.delete(n, [0, j])) ** 2)
                        else: # sum
                            phi_sum = phi_sum + reduce(torch.mul, phi)
                        
                phi_sum = A[i].sum(axis=0).unsqueeze(0) @ phi_sum
                phi_sum = phi_sum.repeat(n[i], 1)
                J1p = J1p + phi_sum
                
                
                A_ = A.copy()
                A_[0] = a[0]
                J2p = tl.cp_tensor.unfolding_dot_khatri_rao(tl.cp_tensor.cp_to_tensor((None, A_), M),
                                                            (None, A_), i)
                
                if reduction == 'mean':
                    J2p = J2p / (n[0] ** 2)
            
            
            ATA_ = ATA.copy()
            del ATA_[i]
            ATA_prod = reduce(torch.mul, ATA_)
            psi_sum = 0
            for j in range(n_modes):
                psi = ATA.copy()
                psi[j] = ATWA[j]
                del psi[i]
                psi_sum = psi_sum + reduce(torch.mul, psi)
            
            J3n =  WA[i] @ ATA_prod + A[i] @ psi_sum
            
            
            ATA_ = ATA.copy()
            del ATA_[i]
            ATA_prod = reduce(torch.mul, ATA_)
            psi_sum = 0
            for j in range(n_modes):
                psi = ATA.copy()
                psi[j] = ATDA[j]
                del psi[i]
                psi_sum = psi_sum + reduce(torch.mul, psi)
            
            J3p =  DA[i] @ ATA_prod + A[i] @ psi_sum
                
            J4p = A[i]
            
            num =  J1n + beta * J3n
            denom = J1p + alpha * J2p + beta * J3p + l * J4p
            A[i] = A[i] * (num/(denom+eps))
            
            ATA[i] = A[i].T @ A[i]
            aaT[i] = torch.outer(A[i].sum(axis=0), A[i].sum(axis=0))
            WA[i] = W_[i] @ A[i]
            DA[i] = D_[i] @ A[i]
            ATWA[i] = A[i].T @ WA[i]
            ATDA[i] = A[i].T @ DA[i]
            
        res = compute_res(A, A_hat)
        if verbose and epoch % freq == 0:
            print(f'Training epoch: {epoch:02d}, Residual: {res:.4f}')
        
        # Terminate the training either when the stopping criterion is met or the maximum number of iterations is reached
        if res < stop_crit:
            break        
    
    return A