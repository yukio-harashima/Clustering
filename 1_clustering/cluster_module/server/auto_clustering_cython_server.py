#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# auto_clustering_cython_server.py
#
# This script is a non-interactive, server-oriented version for custom
# hierarchical clustering. It uses a pre-compiled, parallelized Cython module.
# Final version incorporating custom Ward's method with cumulative distance,
# extended output format, and file copying.

import pandas as pd
import numpy as np
import os
import sys
import shutil
from datetime import datetime

# Cythonでコンパイルされた高速なクラスタリング関数をインポート
try:
    import custom_cluster_fast
except ImportError:
    print("\n" + "="*80)
    print("エラー: Cythonモジュール 'custom_cluster_fast' が見つかりません。")
    print("スクリプトを実行する前に、ターミナルで以下のコマンドを実行してコンパイルしてください:")
    print("python3 setup.py build_ext --inplace")
    print("="*80)
    sys.exit(1)

from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import fcluster

# =============================================================================
# --- スクリプト設定セクション ---
# =============================================================================

# 1. 作業ディレクトリ (元データが存在する場所)
WORKING_DIRECTORY = "~/work/clustering/results_20250730-100041/"

# 2. 目標クラスタ数
TARGET_NUM_CLUSTERS = 8

# 3. 分析パラメータ設定
ANALYSIS_CONFIG = {
    "filter_method": "arbitrary_percentage",
    "filter_quantile": 0.80,
    "vector_cols_list": [['n1_(1)', 'n1_(2)', 'n1_(3)'], ['n2_(1)', 'n2_(2)', 'n2_(3)']],
    "scalar_cols": ['n', 'm', 'tw', 'sliprate', 'lat', 'lon', 'depth'],
    "weights": {'cos1': 1, 'cos2': 1, 'euc': 5},
}

# =============================================================================
# --- 設定の検証と前処理 ---
# =============================================================================
_original_path = WORKING_DIRECTORY
WORKING_DIRECTORY = os.path.expanduser(WORKING_DIRECTORY)
if not os.path.isdir(WORKING_DIRECTORY):
    print(f"エラー: 指定された作業ディレクトリが見つかりません: {_original_path}")
    print(f"       展開後のパス: {WORKING_DIRECTORY}")
    sys.exit(1)

# =============================================================================
# --- 以下、スクリプト本体 ---
# =============================================================================

def ensure_dir(directory_path):
    os.makedirs(directory_path, exist_ok=True)

def generate_output_dirname(config, num_clusters):
    parts = []
    q_val = config.get('filter_quantile', 0)
    parts.append(f"filt-arbP{int(q_val*100)}")
    weights = config.get('weights', {})
    w_str = f"w{int(weights.get('cos1', 0)*10)}-{int(weights.get('cos2', 0)*10)}-{int(weights.get('euc', 0)*10)}"
    parts.append(f"distCUST{w_str}")
    parts.append(f"nClust{num_clusters}")
    return "_".join(parts)

def save_metadata(output_dir, config, num_clusters):
    filepath = os.path.join(output_dir, "metadata.txt")
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Execution Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Number of CPU threads (OMP_NUM_THREADS): {os.environ.get('OMP_NUM_THREADS', 'Default')}\n")
        f.write("--- Configuration ---\n")
        if config:
            for key, value in config.items():
                f.write(f"{key}: {value}\n")
        f.write("\n--- Clustering ---\n")
        f.write(f"Target number of clusters: {num_clusters}\n")
    print(f"メタデータを {filepath} に保存しました。")

def stage_load_and_merge_data(wd, snap_y_file="snap_y.dat", n_vector_file="n_vector.dat"):
    print("\n--- 1. データの読み込み ---")
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
    merged_df = pd.merge(df_snap, df_nvector, on=['n', 'm', 'tw'], how='inner')
    return merged_df

def stage_filter_data(df, method, quantile_val=0.75):
    print(f"\n--- 2. データフィルタリング ---")
    print(f"手法: {method}, Quantile: {quantile_val}")
    if method == 'arbitrary_percentage':
        threshold = df['sliprate'].quantile(quantile_val)
        filtered_df = df[df['sliprate'] > threshold].copy()
        print(f"  フィルタリング前: {df.shape[0]}行, 後: {filtered_df.shape[0]}行")
    else:
        raise NotImplementedError(f"フィルタリング手法 '{method}' は未実装です。")
    if filtered_df.empty:
        raise ValueError("フィルタリングの結果、データが空になりました。")
    return filtered_df

def generate_ancestry(linkage_matrix, num_points):
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

def stage_custom_clustering_and_save(df_merged, df_filtered, config, output_dir, num_clusters):
    print("\n--- 3. カスタム階層クラスタリング実行 ---")
    vector_cols_list, scalar_cols, weights = config['vector_cols_list'], config['scalar_cols'], config['weights']
    data_orig_filtered = df_merged.loc[df_filtered.index]
    all_clustering_cols = [col for group in vector_cols_list for col in group] + scalar_cols
    if data_orig_filtered[all_clustering_cols].isnull().values.any():
        print("  警告: NaN(欠損値)が含まれています。該当する行を削除します。")
        data_orig_filtered = data_orig_filtered.dropna(subset=all_clustering_cols)
        if data_orig_filtered.empty: raise ValueError("NaNを削除した結果、データが空になりました。")
        df_filtered = df_filtered.loc[data_orig_filtered.index]

    scaler = StandardScaler()
    data_scalar_std_np = scaler.fit_transform(data_orig_filtered[scalar_cols])
    all_cols = data_orig_filtered.columns.tolist()
    vector_cols_indices = [[all_cols.index(c) for c in v_group] for v_group in vector_cols_list]
    data_orig_np = data_orig_filtered.to_numpy(dtype=np.double)
    
    print("  並列化されたCythonモジュールを呼び出します...")
    linkage_matrix, std_history = custom_cluster_fast.custom_hierarchical_clustering(
        data_orig_np, data_scalar_std_np, vector_cols_indices, weights)
    print("  クラスタリング計算が完了しました。")
    
    std_history_path = os.path.join(output_dir, "std_values.dat")
    header_parts = [f"cos{i+1}_std" for i in range(len(config['vector_cols_list']))] + ["euc_std"]
    np.savetxt(std_history_path, std_history, header=" ".join(header_parts), comments='')
    print(f"  Standard Deviations History を {std_history_path} に保存しました。")
    
    linkage_matrix_path = os.path.join(output_dir, "linkage_matrix.dat")
    np.savetxt(linkage_matrix_path, linkage_matrix)
    print(f"  Linkage Matrix を {linkage_matrix_path} に保存しました。")
    filtered_indices_path = os.path.join(output_dir, "filtered_indices.dat")
    df_filtered.index.to_series().to_csv(filtered_indices_path, index=False, header=False)
    print(f"  Filtered Indices を {filtered_indices_path} に保存しました。")
    
    print("\n--- 4. ラベル割り当てと結果保存 ---")
    labels_cut = fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    df_labels_cut = pd.DataFrame({'label_cut': labels_cut}, index=df_filtered.index)

    print("  データポイントごとのマージ履歴を生成中...")
    ancestry_map = generate_ancestry(linkage_matrix, len(df_filtered))
    ancestry_series_data = {
        original_idx: '[' + ','.join(map(str, ancestry_map[clustering_idx])) + ']'
        for clustering_idx, original_idx in enumerate(df_filtered.index)
    }
    df_ancestry = pd.DataFrame.from_dict(ancestry_series_data, orient='index', columns=['ancestry'])

    df_merged = df_merged.join(df_ancestry)
    df_merged = df_merged.join(df_labels_cut)
    
    df_merged['ancestry'] = df_merged['ancestry'].fillna('[0]')
    df_merged['label_cut'] = df_merged['label_cut'].fillna(0)
    df_merged['label_cut'] = df_merged['label_cut'].astype(int)
    
    output_path = os.path.join(output_dir, "clusteringSnap2.dat")
    original_snap_y_cols = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                           'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                           'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                           'trendp', 'trendt', 'trendb', 'plungp',
                           'plungt', 'plungb', 'NDC']
    output_cols = original_snap_y_cols + ['ancestry', 'label_cut']
    df_to_save = df_merged[output_cols]
    
    df_to_save.to_csv(output_path, sep=' ', index=False, header=False)
    print(f"クラスタリング結果を {output_path} に保存しました。")
    
    save_metadata(output_dir, config=config, num_clusters=num_clusters)

def copy_related_files(source_dir, dest_dir, filenames):
    print("\n--- 5. 関連ファイルのコピー ---")
    for filename in filenames:
        source_path = os.path.join(source_dir, filename)
        dest_path = os.path.join(dest_dir, filename)
        if os.path.exists(source_path):
            try:
                shutil.copy2(source_path, dest_path)
                print(f"  コピーしました: {filename}")
            except Exception as e:
                print(f"  エラー: {filename} のコピーに失敗しました。理由: {e}")
        else:
            print(f"  スキップ: {filename} は存在しません。")

def main():
    try:
        print("サーバー実行用クラスタリングスクリプトを開始します。")
        print(f"設定: 作業ディレクトリ = {WORKING_DIRECTORY}")
        print(f"設定: 目標クラスタ数 = {TARGET_NUM_CLUSTERS}")
        
        print("\n--- 分析設定 ---")
        for key, value in ANALYSIS_CONFIG.items():
            print(f"  {key}: {value}")
        
        df_merged = stage_load_and_merge_data(WORKING_DIRECTORY)
        df_filtered = stage_filter_data(df_merged, ANALYSIS_CONFIG['filter_method'], ANALYSIS_CONFIG['filter_quantile'])
        
        output_dirname = generate_output_dirname(ANALYSIS_CONFIG, TARGET_NUM_CLUSTERS)
        output_dir = os.path.join(WORKING_DIRECTORY, output_dirname)
        if os.path.exists(output_dir):
            timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
            output_dir = f"{output_dir}_{timestamp}"
        ensure_dir(output_dir)
        print(f"\n結果は次のディレクトリに保存されます:\n{output_dir}")
        
        stage_custom_clustering_and_save(df_merged, df_filtered, ANALYSIS_CONFIG, output_dir, TARGET_NUM_CLUSTERS)
        
        files_to_copy = [
            'snap_yr.dat', 'snap2_y.dat', 'fort.40', 'faultline.dat',
            'mrf.dat', 'snap2.dat', 'rigid_amp.info'
        ]
        copy_related_files(WORKING_DIRECTORY, output_dir, files_to_copy)
        
        print("\n=== すべての処理が正常に完了しました ===")

    except Exception as e:
        print("\n" + "="*80)
        print(f"エラーが発生したため、処理を中断しました: {e}")
        import traceback
        traceback.print_exc()
        print("="*80)

if __name__ == "__main__":
    main()

