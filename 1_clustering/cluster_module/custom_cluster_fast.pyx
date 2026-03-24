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

This is the high-performance core of the custom hierarchical clustering algorithm.

MODIFICATION:
- Added logic to shift the final distance distribution to be non-negative,
  preventing the "Linkage 'Z' contains negative distances" error in SciPy.
"""

import numpy as np
cimport numpy as np
from libc.math cimport sqrt

# -----------------------------------------------------------------------------
# Core Clustering Function (callable from Python)
# -----------------------------------------------------------------------------
cpdef np.ndarray[np.double_t, ndim=2] custom_hierarchical_clustering(
    np.ndarray[np.double_t, ndim=2] data_orig_np,
    np.ndarray[np.double_t, ndim=2] data_scalar_std_np,
    list vector_cols_indices,
    dict weights):
    """
    Performs custom hierarchical clustering with a non-negative distance guarantee.
    """
    cdef int N = data_orig_np.shape[0]
    cdef int num_clusters = N
    cdef int next_cluster_id = N
    cdef list clusters = []
    cdef list linkage_matrix_list = []
    cdef int i, j, k, best_i, best_j, pair_idx, min_pair_idx
    cdef int num_vec_groups = len(vector_cols_indices)
    cdef int num_dist_components = num_vec_groups + 1
    cdef int num_pairs

    cdef np.ndarray[np.double_t, ndim=1] centroid_i_orig, centroid_j_orig
    cdef np.ndarray[np.double_t, ndim=1] centroid_i_std, centroid_j_std
    
    cdef np.ndarray[np.double_t, ndim=2] all_distances_components
    cdef np.ndarray[np.double_t, ndim=2] std_dist_components
    cdef np.ndarray[np.double_t, ndim=1] total_distances # ★★★ 追加 ★★★

    cdef double dist_cos, dist_euc, norm_i, norm_j, dot_product, diff
    cdef double min_dist, total_dist, mean_val, std_val
    cdef double min_total_dist_in_loop, shift_value # ★★★ 追加 ★★★

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
        # 全ての種類の距離（コサイン距離、ユークリッド距離など）を格納するために使用。事前に一括でメモリを確保する
        all_distances_components = np.zeros((num_dist_components, num_pairs), dtype=np.double)

        pair_idx = 0
        for i in range(num_clusters):
            for j in range(i + 1, num_clusters):
                centroid_i_orig = clusters[i]['centroid_orig']
                centroid_j_orig = clusters[j]['centroid_orig']
                centroid_i_std = clusters[i]['centroid_std']
                centroid_j_std = clusters[j]['centroid_std']

                # コサイン類似度計算
                # k=0のとき最初のベクトル(n1)、k=1のとき2番目のベクトル(n2)を処理するループ
                for k in range(num_vec_groups):
                    vec_indices = vector_cols_indices[k]
                    dot_product = 0.0; norm_i = 0.0; norm_j = 0.0

                    # 各ベクトルの成分（3次元なら3回）ループする
                    for col_idx in vec_indices:
                        # 内積: ベクトルiとjの同じ次元の値を掛け合わせ、足し込んでいきます
                        dot_product += centroid_i_orig[col_idx] * centroid_j_orig[col_idx]
                        # ノルム(大きさ)の2乗: ベクトルiの各値を2乗して足し込んでいきます
                        norm_i += centroid_i_orig[col_idx] * centroid_i_orig[col_idx]
                        norm_j += centroid_j_orig[col_idx] * centroid_j_orig[col_idx]
                    
                    # ゼロ除算を避けるためのチェック
                    if norm_i == 0 or norm_j == 0: dist_cos = 1.0
                     # コサイン距離算出
                    else: dist_cos = 1.0 - (dot_product / (sqrt(norm_i) * sqrt(norm_j)))
                    # 計算結果を、事前に確保した配列の(k, pair_idx)の位置に格納
                    all_distances_components[k, pair_idx] = dist_cos

                # ユークリッド距離計算
                # ユークリッド距離の計算結果を保持する変数を0で初期化します
                dist_euc = 0.0
                # 標準化された重心ベクトルの全変数でループ
                for k in range(centroid_i_std.shape[0]):
                    diff = centroid_i_std[k] - centroid_j_std[k]
                    # 差を2乗して、変数に足し込んでいきます
                    dist_euc += diff * diff
                # 最後に全体の平方根を取り、ユークリッド距離を確定
                all_distances_components[num_vec_groups, pair_idx] = sqrt(dist_euc)
                pair_idx += 1

        # 各距離の標準化
        # 標準化後の値を格納するための、元と同じ形の配列を0で作成
        std_dist_components = np.zeros_like(all_distances_components)
        # 各距離の種類ごと(配列の行ごと)にループ
        for k in range(num_dist_components):
            # k番目の距離について、全ペアの距離の平均値を計算
            mean_val = np.mean(all_distances_components[k, :])
            # k番目の距離について、全ペアの距離の標準偏差を計算
            std_val = np.std(all_distances_components[k, :])
            if std_val > 1e-9:
                std_dist_components[k, :] = (all_distances_components[k, :] - mean_val) / std_val
        
        # C. 統合距離の計算
        # 全ペアの最終的な統合距離を格納する1次元配列を作成
        total_distances = np.zeros(num_pairs, dtype=np.double)
        weight_keys = sorted(weights.keys())

        # 各距離の種類ごとにループ
        for k in range(len(weight_keys)):
            key = weight_keys[k]
            # k番目の標準化済み距離の全ペア分に、対応する重みを掛け合わせ、その結果をtotal_distancesに足し込んでいきます
            total_distances += weights[key] * std_dist_components[k, :]
        
        # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
        # ★★★ 負の距離を解消するためのシフト処理 ★★★
        # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
        min_total_dist_in_loop = np.min(total_distances)
        shift_value = 0.0
        # 最小値がマイナスだった場合のみシフト処理
        if min_total_dist_in_loop < 0:
            shift_value = -min_total_dist_in_loop # 最小値の絶対値をシフト量とする
        
        # シフトした後の距離で最小値を探す
        min_dist = np.min(total_distances) + shift_value
        min_pair_idx = np.argmin(total_distances)

        # 最小ペアのインデックス (best_i, best_j) を特定
        best_i = -1; best_j = -1
        pair_idx = 0
        for i in range(num_clusters):
            for j in range(i + 1, num_clusters):
                if pair_idx == min_pair_idx:
                    best_i = i
                    best_j = j
                    break
                pair_idx += 1
            if best_i != -1: break
        
        if best_i > best_j: best_i, best_j = best_j, best_i

        cluster_i = clusters[best_i]
        cluster_j = clusters[best_j]
        new_points_indices = cluster_i['points_indices'] + cluster_j['points_indices']

        # リンケージマトリックスには、シフト後の非負の距離を記録する
        linkage_matrix_list.append([
            <double>cluster_i['id'], <double>cluster_j['id'],
            min_dist, <double>len(new_points_indices)
        ])

        # 新しいクラスタの重心を計算
        new_centroid_orig = data_orig_np[new_points_indices].mean(axis=0)
        new_centroid_std = data_scalar_std_np[new_points_indices].mean(axis=0)

        # 新しいクラスタの情報を格納
        new_cluster = {
            'id': next_cluster_id, 'points_indices': new_points_indices,
            'centroid_orig': new_centroid_orig, 'centroid_std': new_centroid_std
        }
        next_cluster_id += 1

        # 併合元となったクラスタのうち、インデックス番号が大きい方を先に削除
        clusters.pop(best_j)
        # 次に、インデックス番号が小さい方を削除
        clusters.pop(best_i)
        # リストの末尾に、新しく作成したクラスタを追加
        clusters.append(new_cluster)
        
        # 現在のクラスタの総数を更新
        num_clusters = len(clusters)

    print("    [Cython] 併合ループが完了しました。")
    return np.array(linkage_matrix_list, dtype=np.double)
