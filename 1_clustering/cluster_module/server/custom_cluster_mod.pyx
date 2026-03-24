# -*- coding: utf-8 -*-
# cython: language_level=3

# --- Cythonコンパイラへの指示子 ---
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

"""
custom_cluster_mod.pyx

This module implements a custom hierarchical clustering algorithm based on
the Lance-Williams update formula for Ward's method. The parallel distance
calculation is encapsulated in a helper function to resolve compilation errors.
"""

import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from cython.parallel import prange, parallel
from scipy.spatial.distance import squareform

# -----------------------------------------------------------------------------
# 1a. Helper function to calculate distances for one pair (for the initial matrix)
# -----------------------------------------------------------------------------
cdef void _calculate_pair_dist_for_init_matrix(
    double[:] data_i_orig, double[:] data_j_orig,
    double[:] data_i_std, double[:] data_j_std,
    int[:, :] vec_indices_view,
    int num_vec_groups, int num_vec_components,
    double[:] output_distances
) noexcept nogil:
    """
    Calculates all distance components for a single pair of data points.
    """
    cdef int k, component_idx, col_idx
    cdef double dot_product, norm_i, norm_j, dist_euc

    for k in range(num_vec_groups):
        dot_product = 0.0; norm_i = 0.0; norm_j = 0.0
        for component_idx in range(num_vec_components):
            col_idx = vec_indices_view[k, component_idx]
            dot_product += data_i_orig[col_idx] * data_j_orig[col_idx]
            norm_i += data_i_orig[col_idx]**2
            norm_j += data_j_orig[col_idx]**2
        if norm_i == 0 or norm_j == 0:
            output_distances[k] = 1.0
        else:
            output_distances[k] = 1.0 - (dot_product / (sqrt(norm_i) * sqrt(norm_j)))

    dist_euc = 0.0
    for k in range(data_i_std.shape[0]):
        dist_euc += (data_i_std[k] - data_j_std[k])**2
    output_distances[num_vec_groups] = sqrt(dist_euc)

# -----------------------------------------------------------------------------
# 1b. Initial Distance Matrix Calculation
# -----------------------------------------------------------------------------
cdef np.ndarray[np.double_t, ndim=2] _calculate_initial_distance_matrix(
    np.ndarray[np.double_t, ndim=2] data_orig_np,
    np.ndarray[np.double_t, ndim=2] data_scalar_std_np,
    list vector_cols_indices,
    dict weights):
    """
    Computes the initial n x n distance matrix using the custom metric.
    """
    cdef int i, j, k
    cdef int N = data_orig_np.shape[0]
    cdef int num_vec_groups = len(vector_cols_indices)
    cdef int num_dist_components = num_vec_groups + 1
    cdef int num_pairs = N * (N - 1) // 2
    cdef double[:, :] data_orig_view = data_orig_np
    cdef double[:, :] data_scalar_std_view = data_scalar_std_np
    cdef np.ndarray[np.double_t, ndim=2] all_distances_components = np.zeros((num_dist_components, num_pairs), dtype=np.double)
    cdef double[:, :] all_distances_components_view = all_distances_components
    cdef np.ndarray[np.int32_t, ndim=2] vec_indices_np
    if vector_cols_indices: vec_indices_np = np.array(vector_cols_indices, dtype=np.int32)
    else: vec_indices_np = np.empty((0, 0), dtype=np.int32)
    cdef int[:, :] vec_indices_view = vec_indices_np
    cdef int num_vec_components = vec_indices_np.shape[1]
    cdef int i_inner, j_inner, pair_idx_inner
    
    print("    [Cython Mod] Initial distance matrix: calculating components...")
    with nogil, parallel():
        
        for i_inner in prange(N):
            for j_inner in range(i_inner + 1, N):
                pair_idx_inner = (N * i_inner) - (i_inner * (i_inner + 1) // 2) + (j_inner - i_inner - 1)
                # Helper function call
                _calculate_pair_dist_for_init_matrix(
                    data_orig_view[i_inner], data_orig_view[j_inner],
                    data_scalar_std_view[i_inner], data_scalar_std_view[j_inner],
                    vec_indices_view, num_vec_groups, num_vec_components,
                    all_distances_components_view[:, pair_idx_inner])

    print("    [Cython Mod] Initial distance matrix: standardizing and integrating...")
    cdef np.ndarray[np.double_t, ndim=2] std_dist_components = np.zeros_like(all_distances_components)
    cdef double mean_val, std_val
    for k in range(num_dist_components):
        mean_val, std_val = np.mean(all_distances_components[k, :]), np.std(all_distances_components[k, :])
        if std_val > 1e-9:
            std_dist_components[k, :] = (all_distances_components[k, :] - mean_val) / std_val
            
    cdef double global_min_std_dist = np.min(std_dist_components)
    cdef double shift_value = -global_min_std_dist if global_min_std_dist < 0 else 0.0
    std_dist_components += shift_value
    
    cdef np.ndarray[np.double_t, ndim=1] total_distances = np.zeros(num_pairs, dtype=np.double)
    weight_keys = sorted(weights.keys())
    for k in range(len(weight_keys)):
        key = weight_keys[k]
        total_distances += weights[key] * std_dist_components[k, :]
        
    return np.asarray(squareform(total_distances), dtype=np.double)

# -----------------------------------------------------------------------------
# 2. Main Clustering Function using Lance-Williams
# -----------------------------------------------------------------------------
cpdef np.ndarray[np.double_t, ndim=2] lance_williams_ward_clustering(
    np.ndarray[np.double_t, ndim=2] data_orig_np,
    np.ndarray[np.double_t, ndim=2] data_scalar_std_np,
    list vector_cols_indices,
    dict weights):
    """
    Performs hierarchical clustering using a stored distance matrix and the
    Lance-Williams update formula for Ward's method.
    """
    cdef int n_samples = data_orig_np.shape[0]
    if n_samples < 2:
        return np.empty((0, 4), dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=2] D = _calculate_initial_distance_matrix(
        data_orig_np, data_scalar_std_np, vector_cols_indices, weights)

    cdef np.ndarray[np.double_t, ndim=2] linkage_matrix = np.zeros((n_samples - 1, 4), dtype=np.double)
    cdef long[:] cluster_sizes = np.ones(n_samples * 2 - 1, dtype=np.long)
    cdef list active_indices = list(range(n_samples))
    cdef int next_cluster_id = n_samples
    cdef int k, i, j, p, row_idx
    cdef double min_dist
    cdef int min_i = -1, min_j = -1
    cdef int cluster_id_i, cluster_id_j, cluster_id_p
    cdef long size_i, size_j, size_p, new_size
    cdef double total_n, alpha_i, alpha_j, beta, new_dist

    print("    [Cython Mod] Lance-Williams併合ループを開始します...")
    for row_idx in range(n_samples - 1):
        min_dist = np.inf
        for i in range(len(active_indices) - 1):
            for j in range(i + 1, len(active_indices)):
                cluster_id_i = active_indices[i]
                cluster_id_j = active_indices[j]
                if D[cluster_id_i, cluster_id_j] < min_dist:
                    min_dist = D[cluster_id_i, cluster_id_j]
                    min_i, min_j = i, j
        
        cluster_id_i = active_indices[min_i]
        cluster_id_j = active_indices[min_j]
        
        size_i = cluster_sizes[cluster_id_i]
        size_j = cluster_sizes[cluster_id_j]
        new_size = size_i + size_j
        
        linkage_matrix[row_idx, 0] = min(cluster_id_i, cluster_id_j)
        linkage_matrix[row_idx, 1] = max(cluster_id_i, cluster_id_j)
        linkage_matrix[row_idx, 2] = min_dist
        linkage_matrix[row_idx, 3] = new_size

        cluster_sizes[next_cluster_id] = new_size

        for p in range(len(active_indices)):
            if p == min_i or p == min_j: continue
            
            cluster_id_p = active_indices[p]
            size_p = cluster_sizes[cluster_id_p]
            
            total_n = <double>(size_i + size_j + size_p)
            alpha_i = (size_i + size_p) / total_n
            alpha_j = (size_j + size_p) / total_n
            beta = -size_p / total_n
            
            # Use distance squared for Ward's method
            new_dist = alpha_i * D[cluster_id_i, cluster_id_p]**2 + \
                       alpha_j * D[cluster_id_j, cluster_id_p]**2 + \
                       beta * min_dist**2
            
            # The distance for Ward should be the square root of the increase in variance
            if new_dist < 0: new_dist = 0 # Avoid floating point issues
            D[next_cluster_id, cluster_id_p] = sqrt(new_dist)
            D[cluster_id_p, next_cluster_id] = sqrt(new_dist)

        active_indices.pop(max(min_i, min_j))
        active_indices.pop(min(min_i, min_j))
        active_indices.append(next_cluster_id)
        
        next_cluster_id += 1

    print("    [Cython Mod] 併合ループが完了しました。")
    return linkage_matrix

