"""
BYM2 model helper functions for computing scaling factors.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg
from libpysal.weights import Queen, W
from scipy import stats
from typing import Optional, Union, List

def nbs_to_adjlist(nbs: W) -> np.ndarray:
    nbs_adj =  nbs.to_adjlist(remove_symmetric=True)
    j1 = nbs_adj['focal'] + 1
    j2 = nbs_adj['neighbor'] + 1
    edge_pairs = np.vstack([j1, j2])
    return (edge_pairs)

def compute_jittered_precision_matrix(nbs: W) -> sp.spmatrix:
    """
    Compute the precision matrix Q = D - W with added numerical jitter.
    - param nbs: Neighborhood structure (W object from libpysal)
    - return: Sparse precision matrix with jitter for numerical stability
    """
    # Get adjacency matrix in CSR format
    W_sparse = sp.csr_matrix(nbs.sparse)
    # Compute diagonal matrix D (degree matrix)
    D = sp.diags(W_sparse.sum(axis=1).A1)
    # Compute precision matrix Q = D - W
    Q = D - W_sparse
    # Add jitter for numerical stability (similar to R's max(diag(Q)) * sqrt(machine precision))
    jitter = np.max(Q.diagonal()) * np.sqrt(np.finfo(float).eps)
    Q_pert = Q + sp.eye(nbs.n) * jitter
    return Q_pert

def q_inv_dense(Q: sp.spmatrix, A: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute the inverse (or pseudo-inverse) of a precision matrix Q.
    - param Q: Sparse precision matrix
    - param A: Optional constraint matrix (default: None)
    - return: Dense matrix representing Q^(-1) or constrained inverse
    """
    Q_dense = Q.todense()
    try:
        Q_inv = np.linalg.inv(Q_dense)
    except np.linalg.LinAlgError:
        Q_inv = np.linalg.pinv(Q_dense)  # use pseudo-inverse for better stability
    if A is None:
        return Q_inv
    # Apply constraint if provided (ensuring spatial structure is respected)
    W = Q_inv @ A.T
    Q_inv_const = Q_inv - W @ np.linalg.inv(A @ W) @ W.T
    return Q_inv_const


def get_scaling_factor(nbs: W) -> np.float64:
    """
    Compute the geometric mean of the variances from the spatial covariance matrix.
    - param nbs: Neighborhood structure (W object from libpysal)
    - return: Geometric mean of the variances
    """
    Q = compute_jittered_precision_matrix(nbs)
    
    # Compute inverse of precision matrix
    A = np.ones((1, Q.shape[0]))  # Ensure same constraint matrix as in R
    Q_inv = q_inv_dense(Q, A)
    
    # Compute geometric mean of positive variances (diagonal elements of Q_inv)
    variances = np.diag(Q_inv)
    valid_vars = variances[variances > 0]  # Exclude any zero or negative values
    scaling = np.exp(np.mean(np.log(valid_vars)))
    return scaling


def get_scaling_factors(N_components: int, gdf: gpd.GeoDataFrame ) -> List[np.float64]:
    scaling_factors = np.ones(N_components)
    for i in range(N_components):
        comp_gdf = gdf[gdf['comp_id'] == i].reset_index(drop=True)
        comp_nbs = Queen.from_dataframe(comp_gdf, geom_col='geometry')
        component_w = W(comp_nbs.neighbors, comp_nbs.weights)
        scaling_factors[i] = get_scaling_factor(component_w)

    return scaling_factors
