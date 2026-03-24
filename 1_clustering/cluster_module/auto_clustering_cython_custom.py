#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# auto_clustering_cython_custom.py
# This script uses a Cython-accelerated module for custom hierarchical clustering.
#
# MODIFICATION:
# - Added a NaN check before calling the Cython module to improve stability.

print("現在テストコードです。 スクリプトが開始されました。ログ確認用。")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import re
from datetime import datetime
from tqdm import tqdm

# Cythonでコンパイルされた高速なクラスタリング関数をインポート
try:
    import custom_cluster_fast
except ImportError:
    print("\n" + "="*80)
    print("エラー: Cythonモジュール 'custom_cluster_fast' が見つかりません。")
    print("スクリプトを実行する前に、ターミナルで以下のコマンドを実行してコンパイルしてください:")
    print("python setup.py build_ext --inplace")
    print("="*80)
    exit()

from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, fcluster

# --- ヘルパー関数 ---
def ensure_dir(directory_path):
    """指定されたパスのディレクトリが存在しない場合、作成する"""
    os.makedirs(directory_path, exist_ok=True)

def generate_output_dirname(config):
    """設定辞書から分かりやすいディレクトリ名を生成する"""
    parts = []
    
    filter_method = config.get('filter_method', 'none')
    if filter_method == 'arbitrary_percentage':
        q_val = config.get('filter_quantile', 0)
        parts.append(f"filt-arbP{q_val*100:.0f}")
    
    weights = config.get('weights', {})
    w_str = f"w{weights.get('cos1', 0)*10:.0f}-{weights.get('cos2', 0)*10:.0f}-{weights.get('euc', 0)*10:.0f}"
    parts.append(f"distCUST{w_str}")
    
    return "_".join(parts)

def save_metadata(output_dir, config=None, comments=None):
    """実行設定やコメントをメタデータファイルに保存する"""
    filepath = os.path.join(output_dir, "metadata.txt")
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Execution Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("--- Configuration ---\n")
        if config:
            for key, value in config.items():
                f.write(f"{key}: {value}\n")
        if comments:
            f.write("\n--- Comments ---\n")
            f.write(comments)
    print(f"メタデータを {filepath} に保存しました。")

# --- データ処理ステージ関数 ---
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
    print(f"{snap_y_file} を読み込み。形状: {df_snap.shape}")

    ps_nvector = os.path.join(wd, n_vector_file)
    if not os.path.exists(ps_nvector): raise FileNotFoundError(f"{ps_nvector} が見つかりません。")
    df_nvector = pd.read_table(ps_nvector, sep=r'\s+', header=None,
                               names=['n', 'm', 'tw','t', 'x_coord', 'y_coord',
                                      'n1_(1)','n1_(2)','n1_(3)',
                                      'n2_(1)','n2_(2)','n2_(3)',
                                      'sliprate_nv'])
    print(f"{n_vector_file} を読み込み。形状: {df_nvector.shape}")
    
    df_nvector.drop('sliprate_nv', axis=1, inplace=True, errors='ignore')
    merged_df = pd.merge(df_snap, df_nvector, on=['n', 'm', 'tw'], how='inner')
    print(f"マージ後のデータ形状: {merged_df.shape}")
    
    return merged_df

def stage_filter_data(df, method, quantile_val=0.75):
    print(f"\n--- 2. データフィルタリング ---")
    print(f"手法: {method}")
    
    if method == 'arbitrary_percentage':
        threshold = df['sliprate'].quantile(quantile_val)
        filtered_df = df[df['sliprate'] > threshold].copy()
        print(f"  sliprate > {threshold:.4f} (quantile={quantile_val}) でフィルタリング。")
        print(f"  フィルタリング前: {df.shape[0]}行, 後: {filtered_df.shape[0]}行")
    else:
        raise NotImplementedError(f"フィルタリング手法 '{method}' は未実装です。")
        
    if filtered_df.empty:
        raise ValueError("フィルタリングの結果、データが空になりました。")
        
    return filtered_df

def stage_custom_clustering_and_save(df_merged, df_filtered, config, output_dir):
    print("\n--- 3. カスタム階層クラスタリング実行 ---")
    
    vector_cols_list = config['vector_cols_list']
    scalar_cols = config['scalar_cols']
    weights = config['weights']
    
    # 元のデータから、フィルタリング後も残っている行を抽出
    data_orig_filtered = df_merged.loc[df_filtered.index]
    
    # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    # ★★★ 安定性向上のためのNaNチェックを追加 ★★★
    # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    all_clustering_cols = [col for group in vector_cols_list for col in group] + scalar_cols
    if data_orig_filtered[all_clustering_cols].isnull().values.any():
        print("  警告: クラスタリング対象のデータにNaN(欠損値)が含まれています。該当する行を削除します。")
        nan_rows = data_orig_filtered[all_clustering_cols].isnull().any(axis=1)
        print(f"    削除される行数: {nan_rows.sum()}")
        data_orig_filtered = data_orig_filtered.dropna(subset=all_clustering_cols)
        if data_orig_filtered.empty:
            raise ValueError("NaNを削除した結果、データが空になりました。")
        # フィルタリング後のDataFrameも更新
        df_filtered = df_filtered.loc[data_orig_filtered.index]

    # スカラー変数を標準化
    scaler = StandardScaler()
    data_scalar_std_np = scaler.fit_transform(data_orig_filtered[scalar_cols])
    
    # Cython関数に渡す列のインデックスを取得
    all_cols = data_orig_filtered.columns.tolist()
    vector_cols_indices = [[all_cols.index(c) for c in v_group] for v_group in vector_cols_list]
    
    data_orig_np = data_orig_filtered.to_numpy(dtype=np.double)
    
    linkage_matrix = custom_cluster_fast.custom_hierarchical_clustering(
        data_orig_np,
        data_scalar_std_np,
        vector_cols_indices,
        weights
    )
    
    print("\n--- 4. 閾値決定と結果保存 ---")
    
    plt.figure(figsize=(15, 10))
    dendrogram(linkage_matrix)
    plt.title("Custom Hierarchical Clustering Dendrogram")
    plt.xlabel("Sample Index")
    plt.ylabel("Distance (Standardized Weighted Sum)")
    plt.show()
    
    threshold = -1.0
    while threshold <= 0:
        try:
            threshold_str = input("  デンドログラムを参考に、クラスタを分割する距離の閾値を入力してください: ")
            threshold = float(threshold_str)
            if threshold <= 0:
                print("  エラー: 閾値は正の数で入力してください。")
        except ValueError:
            print("  エラー: 無効な数値です。")

    labels = fcluster(linkage_matrix, threshold, criterion='distance')
    
    # 結果をフィルタリング後のデータフレームに追加
    df_filtered_with_labels = df_filtered.copy()
    df_filtered_with_labels['label'] = labels
    
    # 元のデータフレームに結果をマージ（ラベルがないものは0とする）
    df_merged['label'] = df_filtered_with_labels['label']
    df_merged['label'].fillna(0, inplace=True)
    df_merged['label'] = df_merged['label'].astype(int)
    
    output_path = os.path.join(output_dir, "clusteringSnap_custom.dat")
    # 出力する列を元のsnap_y.datの形式に近づける
    # 元のsnap_yの列 + label
    original_snap_y_cols = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                           'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                           'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                           'trendp', 'trendt', 'trendb', 'plungp',
                           'plungt', 'plungb', 'NDC']
    output_cols = original_snap_y_cols + ['label']
    df_to_save = df_merged[output_cols]
    
    df_to_save.to_csv(output_path, sep=' ', index=False, header=False)
    print(f"クラスタリング結果を {output_path} に保存しました。")
    
    save_metadata(output_dir, config=config, comments=f"Manual threshold set to: {threshold}")

def main():
    try:
        print("現在テストコードです。 スクリプトが開始されました。ログ確認用。")
        default_working_directory = "./"

        working_directory = input(f"1. 作業ディレクトリ (デフォルト: {default_working_directory}): ") or default_working_directory
        if not os.path.isdir(working_directory):
            print(f"エラー: ディレクトリ '{working_directory}' が存在しません。終了します。")
            return
            
        fixed_config = {
            "filter_method": "arbitrary_percentage",
            "filter_quantile": 0.80,
            "vector_cols_list": [['n1_(1)', 'n1_(2)', 'n1_(3)'], ['n2_(1)', 'n2_(2)', 'n2_(3)']],
            "scalar_cols": ['n', 'm', 'tw', 'sliprate'],
            "weights": {'cos1': 1, 'cos2': 1, 'euc': 1},
        }
        
        print("\n--- 固定設定による分析を実行します ---")
        for key, value in fixed_config.items():
            print(f"  {key}: {value}")
        
        df_merged = stage_load_and_merge_data(working_directory)
        df_filtered = stage_filter_data(df_merged, fixed_config['filter_method'], fixed_config['filter_quantile'])
        
        output_dirname = generate_output_dirname(fixed_config)
        output_dir = os.path.join(working_directory, output_dirname)
        if os.path.exists(output_dir):
            timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
            output_dir = f"{output_dir}_{timestamp}"
        ensure_dir(output_dir)
        print(f"\n結果は次のディレクトリに保存されます:\n{output_dir}")
        
        stage_custom_clustering_and_save(df_merged, df_filtered, fixed_config, output_dir)
        
        print("\n=== すべての処理が正常に完了しました ===")

    except Exception as e:
        print("\n" + "="*80)
        print(f"エラーが発生したため、処理を中断しました: {e}")
        import traceback
        traceback.print_exc()
        print("="*80)

if __name__ == "__main__":
    main()

