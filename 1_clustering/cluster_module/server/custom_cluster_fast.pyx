# -*- coding: utf-8 -*-
# cython: language_level=3

# --- Cythonコンパイラへの指示子 ---
# 各指示子の後ろのコメントは、必ず別の行に記述する必要があります。

# 配列の範囲外アクセスチェックを無効化し、速度を向上させる
# cython: boundscheck=False
# 配列のマイナスインデックス（例: a[-1]）を無効化し、速度を向上させる
# cython: wraparound=False
# ゼロ除算のチェックを無効化し、速度を向上させる
# cython: cdivision=True

"""
custom_cluster_fast.pyx

This is the final version incorporating the custom Ward's method with a
cumulative distance logic to ensure monotonicity. It also fixes a critical
bug where the linkage matrix was being created before data was populated.
"""

import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from cython.parallel import prange, parallel

# -----------------------------------------------------------------------------
# 距離計算をカプセル化するヘルパー関数
# -----------------------------------------------------------------------------
cdef void _calculate_distances_for_pair(
    double[:] centroid_i_orig, double[:] centroid_j_orig,
    double[:] centroid_i_std, double[:] centroid_j_std,
    int[:, :] vec_indices_view,
    int num_vec_groups, int num_vec_components,
    double[:] output_distances
) noexcept nogil:
    """
    一つのクラスタペアに対する全ての距離を計算するCレベルのヘルパー関数。
    """
    cdef int k, component_idx, col_idx
    cdef double dot_product, norm_i, norm_j, dist_euc

    for k in range(num_vec_groups):
        dot_product = 0.0; norm_i = 0.0; norm_j = 0.0
        for component_idx in range(num_vec_components):
            col_idx = vec_indices_view[k, component_idx]
            dot_product += centroid_i_orig[col_idx] * centroid_j_orig[col_idx]
            norm_i += centroid_i_orig[col_idx] * centroid_i_orig[col_idx]
            norm_j += centroid_j_orig[col_idx] * centroid_j_orig[col_idx]

        if norm_i == 0 or norm_j == 0:
            output_distances[k] = 1.0
        else:
            output_distances[k] = 1.0 - (dot_product / (sqrt(norm_i) * sqrt(norm_j)))

    dist_euc = 0.0
    for k in range(centroid_i_std.shape[0]):
        dist_euc += (centroid_i_std[k] - centroid_j_std[k])**2
    output_distances[num_vec_groups] = sqrt(dist_euc)


# -----------------------------------------------------------------------------
# Core Clustering Function (callable from Python)
# -----------------------------------------------------------------------------
cpdef tuple custom_hierarchical_clustering(
    np.ndarray[np.double_t, ndim=2] data_orig_np,
    np.ndarray[np.double_t, ndim=2] data_scalar_std_np,
    list vector_cols_indices,
    dict weights):
    """
    Performs custom hierarchical clustering using a cumulative Ward-like method.
    """
    cdef int i, j, k, best_i, best_j, pair_idx, min_pair_idx, pair_idx_inner
    cdef int N = data_orig_np.shape[0]
    cdef int num_clusters = N
    cdef int next_cluster_id = N
    cdef list clusters = []
    cdef list linkage_matrix_list = []
    cdef list std_history_list = []
    # ★★★ 修正点: 変数宣言を関数の先頭に移動 ★★★
    cdef list current_stds
    
    cdef int num_vec_groups = len(vector_cols_indices)
    cdef int num_dist_components = num_vec_groups + 1
    cdef int num_pairs
    cdef double[:, :] all_distances_components_view, all_centroids_orig_view, all_centroids_std_view
    cdef np.ndarray[np.double_t, ndim=2] all_distances_components, std_dist_components, all_centroids_orig, all_centroids_std
    cdef np.ndarray[np.double_t, ndim=1] total_distances
    cdef np.ndarray[np.double_t, ndim=1] ward_costs
    cdef int nA, nB
    
    cdef np.ndarray[np.int32_t, ndim=2] vec_indices_np
    cdef int[:, :] vec_indices_view
    cdef int num_vec_components
    cdef double min_dist, mean_val, std_val, global_min_std_dist, shift_value
    cdef double max_dist_so_far = 0.0
    cdef double final_merge_dist

    if vector_cols_indices:
        vec_indices_np = np.array(vector_cols_indices, dtype=np.int32)
    else:
        vec_indices_np = np.empty((0, 0), dtype=np.int32)
    vec_indices_view = vec_indices_np
    num_vec_components = vec_indices_np.shape[1]

    print("    [Cython] クラスタリストを初期化しています...")
    for i in range(N):
        clusters.append({
            'id': i, 'points_indices': [i],
            'centroid_orig': data_orig_np[i, :], 'centroid_std': data_scalar_std_np[i, :]
        })

    print("    [Cython] メインの併合ループを開始します...")
    while num_clusters > 1:
        if num_clusters % 50 == 0 and num_clusters != N:
            print(f"      [Cython] 残りクラスタ数: {num_clusters}")

        num_pairs = num_clusters * (num_clusters - 1) // 2
        all_distances_components = np.zeros((num_dist_components, num_pairs), dtype=np.double)

        all_centroids_orig = np.array([c['centroid_orig'] for c in clusters], dtype=np.double)
        all_centroids_std = np.array([c['centroid_std'] for c in clusters], dtype=np.double)

        all_distances_components_view = all_distances_components
        all_centroids_orig_view = all_centroids_orig
        all_centroids_std_view = all_centroids_std

        with nogil, parallel():
            for i in prange(num_clusters):
                for j in range(i + 1, num_clusters):
                    pair_idx_inner = (num_clusters * i) - (i * (i + 1) // 2) + (j - i - 1)
                    _calculate_distances_for_pair(
                        all_centroids_orig_view[i], all_centroids_orig_view[j],
                        all_centroids_std_view[i], all_centroids_std_view[j],
                        vec_indices_view, num_vec_groups, num_vec_components,
                        all_distances_components_view[:, pair_idx_inner])

        std_dist_components = np.zeros_like(all_distances_components)
        # ★★★ 修正点: ループの各周回でリストを初期化 ★★★
        current_stds = []
        for k in range(num_dist_components):
            mean_val = np.mean(all_distances_components[k, :])
            std_val = np.std(all_distances_components[k, :])
            current_stds.append(std_val)
            if std_val > 1e-9:
                std_dist_components[k, :] = (all_distances_components[k, :] - mean_val) / std_val
        std_history_list.append(current_stds)
        
        global_min_std_dist = np.min(std_dist_components)
        shift_value = 0.0
        if global_min_std_dist < 0:
            shift_value = -global_min_std_dist
        
        std_dist_components += shift_value
        
        total_distances = np.zeros(num_pairs, dtype=np.double)
        weight_keys = sorted(weights.keys())
        for k in range(len(weight_keys)):
            key = weight_keys[k]
            total_distances += weights[key] * std_dist_components[k, :]

        ward_costs = np.zeros(num_pairs, dtype=np.double)
        pair_idx = 0
        for i in range(num_clusters):
            for j in range(i + 1, num_clusters):
                nA = len(clusters[i]['points_indices'])
                nB = len(clusters[j]['points_indices'])
                ward_costs[pair_idx] = (float(nA * nB) / float(nA + nB)) * total_distances[pair_idx]
                pair_idx += 1
        
        min_dist = np.min(ward_costs)
        min_pair_idx = np.argmin(ward_costs)
        
        final_merge_dist = min_dist + max_dist_so_far
            
        best_i = -1; best_j = -1
        pair_idx = 0
        for i in range(num_clusters):
            for j in range(i + 1, num_clusters):
                if pair_idx == min_pair_idx: best_i, best_j = i, j; break
                pair_idx += 1
            if best_i != -1: break

        if best_i > best_j: best_i, best_j = best_j, best_i

        cluster_i = clusters[best_i]
        cluster_j = clusters[best_j]
        new_points_indices = cluster_i['points_indices'] + cluster_j['points_indices']

        linkage_matrix_list.append([
            <double>cluster_i['id'], <double>cluster_j['id'],
            final_merge_dist, <double>len(new_points_indices)])
        
        max_dist_so_far = final_merge_dist
        
        new_centroid_orig = data_orig_np[new_points_indices].mean(axis=0)
        new_centroid_std = data_scalar_std_np[new_points_indices].mean(axis=0)

        new_cluster = {
            'id': next_cluster_id, 'points_indices': new_points_indices,
            'centroid_orig': new_centroid_orig, 'centroid_std': new_centroid_std}
        next_cluster_id += 1

        clusters.pop(best_j)
        clusters.pop(best_i)
        clusters.append(new_cluster)
        num_clusters = len(clusters)

    print("    [Cython] 併合ループが完了しました。")
    
    cdef np.ndarray linkage_matrix_np = np.array(linkage_matrix_list, dtype=np.double)
    cdef np.ndarray std_history_np = np.array(std_history_list, dtype=np.double)
    
    return (linkage_matrix_np, std_history_np)

