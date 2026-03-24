# -*- coding: utf-8 -*-
# cython: language_level=3

# --- Cythonコンパイラへの指示子 ---
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

"""
calculate_distances_cython.pyx

This module is dedicated to calculating the full, all-pairs custom distance
matrix at high speed using parallel processing. It is intended for use
with evaluation metrics like the silhouette score.
"""

import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from cython.parallel import prange, parallel
from scipy.spatial.distance import squareform

# -----------------------------------------------------------------------------
# ★★★ 追加点: 距離計算をカプセル化するヘルパー関数 ★★★
# -----------------------------------------------------------------------------
cdef void _calculate_pair_distances_for_matrix(
    double[:] data_i_orig, double[:] data_j_orig,
    double[:] data_i_std, double[:] data_j_std,
    int[:, :] vec_indices_view,
    int num_vec_groups, int num_vec_components,
    double[:] output_distances
) noexcept nogil:
    """
    一つのデータペアに対する全ての距離を計算するCレベルのヘルパー関数。
    """
    cdef int k, component_idx, col_idx
    cdef double dot_product, norm_i, norm_j, dist_euc

    for k in range(num_vec_groups):
        dot_product = 0.0; norm_i = 0.0; norm_j = 0.0
        for component_idx in range(num_vec_components):
            col_idx = vec_indices_view[k, component_idx]
            dot_product += data_i_orig[col_idx] * data_j_orig[col_idx]
            norm_i += data_i_orig[col_idx] * data_i_orig[col_idx]
            norm_j += data_j_orig[col_idx] * data_j_orig[col_idx]

        if norm_i == 0 or norm_j == 0:
            output_distances[k] = 1.0
        else:
            output_distances[k] = 1.0 - (dot_product / (sqrt(norm_i) * sqrt(norm_j)))

    dist_euc = 0.0
    for k in range(data_i_std.shape[0]):
        dist_euc += (data_i_std[k] - data_j_std[k])**2
    output_distances[num_vec_groups] = sqrt(dist_euc)

# -----------------------------------------------------------------------------
# Core Distance Matrix Calculation Function
# -----------------------------------------------------------------------------
cpdef np.ndarray[np.double_t, ndim=2] calculate_full_distance_matrix(
    np.ndarray[np.double_t, ndim=2] data_orig_np,
    np.ndarray[np.double_t, ndim=2] data_scalar_std_np,
    list vector_cols_indices,
    dict weights):
    """
    Computes the full all-pairs custom distance matrix at high speed.
    """
    # --- 関数全体で使用する変数の宣言 ---
    cdef int i, j, k, i_inner, j_inner, pair_idx_inner
    cdef int N = data_orig_np.shape[0]
    cdef int num_vec_groups = len(vector_cols_indices)
    cdef int num_dist_components = num_vec_groups + 1
    cdef int num_pairs = N * (N - 1) // 2
    cdef double[:, :] data_orig_view = data_orig_np
    cdef double[:, :] data_scalar_std_view = data_scalar_std_np
    cdef np.ndarray[np.double_t, ndim=2] all_distances_components = np.zeros((num_dist_components, num_pairs), dtype=np.double)
    cdef double[:, :] all_distances_components_view = all_distances_components
    cdef np.ndarray[np.int32_t, ndim=2] vec_indices_np
    cdef int[:, :] vec_indices_view
    cdef int num_vec_components
    cdef np.ndarray[np.double_t, ndim=2] std_dist_components
    cdef double mean_val, std_val, global_min_std_dist, shift_value
    cdef np.ndarray[np.double_t, ndim=1] total_distances
    
    # --- 初期化処理 ---
    if vector_cols_indices:
        vec_indices_np = np.array(vector_cols_indices, dtype=np.int32)
    else:
        vec_indices_np = np.empty((0, 0), dtype=np.int32)
    vec_indices_view = vec_indices_np
    num_vec_components = vec_indices_np.shape[1]
    
    print("    [Cython] 全ペアの距離成分を並列計算中...")
    with nogil, parallel():
        for i_inner in prange(N):
            for j_inner in range(i_inner + 1, N):
                pair_idx_inner = (N * i_inner) - (i_inner * (i_inner + 1) // 2) + (j_inner - i_inner - 1)
                
                _calculate_pair_distances_for_matrix(
                    data_orig_view[i_inner], data_orig_view[j_inner],
                    data_scalar_std_view[i_inner], data_scalar_std_view[j_inner],
                    vec_indices_view,
                    num_vec_groups, num_vec_components,
                    all_distances_components_view[:, pair_idx_inner]
                )

    print("    [Cython] 標準化と統合距離の計算中...")
    std_dist_components = np.zeros_like(all_distances_components)
    for k in range(num_dist_components):
        mean_val, std_val = np.mean(all_distances_components[k, :]), np.std(all_distances_components[k, :])
        if std_val > 1e-9:
            std_dist_components[k, :] = (all_distances_components[k, :] - mean_val) / std_val
            
    global_min_std_dist = np.min(std_dist_components)
    shift_value = -global_min_std_dist if global_min_std_dist < 0 else 0.0
    std_dist_components += shift_value
    
    total_distances = np.zeros(num_pairs, dtype=np.double)
    weight_keys = sorted(weights.keys())
    for k in range(len(weight_keys)):
        key = weight_keys[k]
        total_distances += weights[key] * std_dist_components[k, :]

    return np.asarray(squareform(total_distances), dtype=np.double)

