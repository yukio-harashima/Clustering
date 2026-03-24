#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# relabel_clusters.py
#
# This script re-assigns cluster labels using a pre-computed linkage matrix.
# It allows for rapid exploration of different numbers of clusters without
# re-running the expensive clustering algorithm.
# It outputs both the merge ancestry (uncut) and the new final cluster label (cut).

import pandas as pd
import numpy as np
import os
import sys
from scipy.cluster.hierarchy import fcluster
import shutil

# =============================================================================
# --- スクリプト設定セクション ---
# ユーザーはこのセクションの値を直接編集してください。
# =============================================================================

# 1. 元の分析結果が保存されているディレクトリ
#    linkage_matrix.dat と filtered_indices.dat が含まれるディレクトリを指定します。
ANALYSIS_OUTPUT_DIRECTORY = "~/work/clustering/results_20250730-100041/filt-arbP75_distCUSTw10-10-10_nClust8"

# 2. 新しい目標クラスタ数
#    試したいクラスタの数を指定します。
NEW_TARGET_NUM_CLUSTERS = 10

# 3. 元データがあるディレクトリ
#    snap_y.dat と n_vector.dat がある作業ディレクトリを指定します。
#    auto_clustering_cython_server.py の WORKING_DIRECTORY と同じ値を設定します。
ORIGINAL_DATA_DIRECTORY = "~/work/clustering/results_20250730-100041/"

# =============================================================================
# --- 以下、スクリプト本体 ---
# =============================================================================

def generate_ancestry(linkage_matrix, num_points):
    """
    linkage_matrixから各データポイントのマージ履歴（祖先クラスタIDのリスト）を生成する。
    """
    cluster_contents = {i: {i} for i in range(num_points)}
    ancestry_map = {i: [] for i in range(num_points)}
    for i, row in enumerate(linkage_matrix):
        new_cluster_id = num_points + i
        cluster1_id, cluster2_id = int(row[0]), int(row[1])
        points1 = cluster_contents.get(cluster1_id, set())
        points2 = cluster_contents.get(cluster2_id, set())
        merged_points = points1.union(points2)
        for point_idx in merged_points:
            ancestry_map[point_idx].append(float(new_cluster_id))
        cluster_contents[new_cluster_id] = merged_points
        if cluster1_id >= num_points: del cluster_contents[cluster1_id]
        if cluster2_id >= num_points: del cluster_contents[cluster2_id]
    return ancestry_map

def stage_load_and_merge_data(wd, snap_y_file="snap_y.dat", n_vector_file="n_vector.dat"):
    # '~' をユーザーのホームディレクトリに展開
    wd = os.path.expanduser(wd)
    ps_snap = os.path.join(wd, snap_y_file)
    if not os.path.exists(ps_snap): raise FileNotFoundError(f"{ps_snap} が見つかりません。")
    df_snap = pd.read_table(ps_snap, sep=r'\s+', header=None,
                            names=['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                                   'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                   'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                   'trendp', 'trendt', 'trendb', 'plungp',
                                   'plungt', 'plungb', 'NDC'])
    ps_nvector = os.path.join(wd, n_vector_file)
    if not os.path.exists(ps_nvector): raise FileNotFoundError(f"{ps_nvector} が見つかりません。")
    df_nvector = pd.read_table(ps_nvector, sep=r'\s+', header=None,
                               names=['n', 'm', 'tw','t', 'x_coord', 'y_coord',
                                      'n1_(1)','n1_(2)','n1_(3)',
                                      'n2_(1)','n2_(2)','n2_(3)',
                                      'sliprate_nv'])
    df_nvector.drop('sliprate_nv', axis=1, inplace=True, errors='ignore')
    return pd.merge(df_snap, df_nvector, on=['n', 'm', 'tw'], how='inner')

def main():
    try:
        print("クラスタ再割り当てスクリプトを開始します。")
        # '~' をユーザーのホームディレクトリに展開
        analysis_dir = os.path.expanduser(ANALYSIS_OUTPUT_DIRECTORY)
        linkage_matrix_path = os.path.join(analysis_dir, "linkage_matrix.dat")
        filtered_indices_path = os.path.join(analysis_dir, "filtered_indices.dat")

        if not os.path.exists(linkage_matrix_path) or not os.path.exists(filtered_indices_path):
            print(f"エラー: {analysis_dir} に linkage_matrix.dat または filtered_indices.dat が見つかりません。")
            sys.exit(1)

        print("\n--- 1. 保存済みデータの読み込み ---")
        linkage_matrix = np.loadtxt(linkage_matrix_path)
        filtered_indices = pd.read_csv(filtered_indices_path, header=None, index_col=0).index
        df_merged = stage_load_and_merge_data(ORIGINAL_DATA_DIRECTORY)
        
        print(f"\n--- 2. 新しい目標クラスタ数 ({NEW_TARGET_NUM_CLUSTERS}) で再割り当て ---")

        # 30列目: 新しいカットされたクラスタ番号
        new_labels_cut = fcluster(linkage_matrix, NEW_TARGET_NUM_CLUSTERS, criterion='maxclust')
        df_new_labels_cut = pd.DataFrame({'label_cut': new_labels_cut}, index=filtered_indices)

        # 29列目: 階層クラスタ番号（マージ履歴）
        print("  データポイントごとのマージ履歴を生成中...")
        ancestry_map = generate_ancestry(linkage_matrix, len(filtered_indices))
        
        # ダブルクォーテーションを避けるため、スペースなしの文字列を生成
        ancestry_series_data = {
            original_idx: '[' + ','.join(map(str, ancestry_map[clustering_idx])) + ']'
            for clustering_idx, original_idx in enumerate(filtered_indices)
        }
        df_ancestry = pd.DataFrame.from_dict(ancestry_series_data, orient='index', columns=['ancestry'])

        # 結果をマージして保存
        print("\n--- 3. 結果の保存 ---")
        df_merged = df_merged.join(df_ancestry)
        df_merged = df_merged.join(df_new_labels_cut)
        
        # 警告を回避し、対象外データを指定の値で埋める
        df_merged['ancestry'] = df_merged['ancestry'].fillna('[0]')
        df_merged['label_cut'] = df_merged['label_cut'].fillna(0)
        df_merged['label_cut'] = df_merged['label_cut'].astype(int)

        # --- 新しい結果を保存するサブディレクトリを作成 ---
        new_output_dirname = f"relabeled_nClust{NEW_TARGET_NUM_CLUSTERS}"
        new_output_dir = os.path.join(analysis_dir, new_output_dirname)
        os.makedirs(new_output_dir, exist_ok=True)
        
        output_path = os.path.join(new_output_dir, "clusteringSnap2.dat")
        original_snap_y_cols = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                               'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                               'trendp', 'trendt', 'trendb', 'plungp',
                               'plungt', 'plungb', 'NDC']
        output_cols = original_snap_y_cols + ['ancestry', 'label_cut']
        df_to_save = df_merged[output_cols]
        
        df_to_save.to_csv(output_path, sep=' ', index=False, header=False)




        print(f"新しいクラスタリング結果を {output_path} に保存しました。")
        print("\n=== 再割り当て処理が正常に完了しました ===")

    except Exception as e:
        print(f"\nエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

