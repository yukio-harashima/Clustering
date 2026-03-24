#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# calculate_silhouette.py
#
# This script calculates the silhouette score for a range of cluster numbers.
# It uses the dedicated 'calculate_distances_cython' module for high performance
# and automatically restores the analysis configuration from metadata.txt.

import pandas as pd
import numpy as np
import os
import sys
import ast
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
# ★★★ 変更点: 新しい距離計算モジュール 'calculate_distances_cython' をインポート ★★★
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
try:
    import calculate_distances_cython
except ImportError:
    print("\n" + "="*80)
    print("エラー: Cythonモジュール 'calculate_distances_cython' が見つかりません。")
    print("スクリプトを実行する前に、ターミナルで以下のコマンドを実行してコンパイルしてください:")
    print("python3 setup_distance.py build_ext --inplace")
    print("="*80)
    sys.exit(1)

# =============================================================================
# --- スクリプト設定セクション ---
# =============================================================================

# 1. 分析結果が保存されているディレクトリ
#    metadata.txt, linkage_matrix.dat などが含まれるディレクトリを指定します。
ANALYSIS_OUTPUT_DIRECTORY = "~/work/clustering/results_20250730-100041/filt-arbP75_distCUSTw10-10-10_nClust8_20250928-204431"

# 2. 評価するクラスター数の範囲
CLUSTER_RANGE_TO_EVALUATE = range(3, 16)

# 3. 元データがあるディレクトリ
ORIGINAL_DATA_DIRECTORY = "~/work/clustering/results_20250730-100041/"

# =============================================================================
# --- 以下、スクリプト本体 ---
# =============================================================================

def load_config_from_metadata(result_dir):
    """
    metadata.txtを読み込み、ANALYSIS_CONFIGを辞書として復元する。
    """
    print("--- 1a. metadata.txt から設定を復元 ---")
    metadata_path = os.path.join(result_dir, "metadata.txt")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"metadata.txt が見つかりません: {metadata_path}")

    config = {}
    is_config_section = False
    with open(metadata_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith("--- Configuration ---"):
                is_config_section = True
                continue
            if is_config_section and line.startswith("---"):
                break
            
            if is_config_section and ':' in line:
                key, value_str = line.split(':', 1)
                key, value_str = key.strip(), value_str.strip()
                try:
                    config[key] = ast.literal_eval(value_str)
                except (ValueError, SyntaxError):
                    config[key] = value_str
    
    if not config:
        raise ValueError("metadata.txtから設定を復元できませんでした。")
    print("  設定の復元が完了しました。")
    return config

def load_and_prepare_data(data_dir, result_dir):
    """
    元データとLinkage Matrixを読み込み、スコア計算に必要なデータを準備する。
    """
    print("--- 1b. データの読み込みと準備 ---")
    data_dir = os.path.expanduser(data_dir)
    snap_y_path = os.path.join(data_dir, "snap_y.dat")
    n_vector_path = os.path.join(data_dir, "n_vector.dat")
    df_snap = pd.read_table(snap_y_path, sep=r'\s+', header=None, names=['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth', 'trendp', 'trendt', 'trendb', 'plungp', 'plungt', 'plungb', 'NDC'])
    df_nvector = pd.read_table(n_vector_path, sep=r'\s+', header=None, names=['n', 'm', 'tw','t', 'x_coord', 'y_coord', 'n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)', 'sliprate_nv'])
    df_nvector.drop('sliprate_nv', axis=1, inplace=True, errors='ignore')
    df_merged = pd.merge(df_snap, df_nvector, on=['n', 'm', 'tw'], how='inner')

    result_dir = os.path.expanduser(result_dir)
    linkage_path = os.path.join(result_dir, "linkage_matrix.dat")
    indices_path = os.path.join(result_dir, "filtered_indices.dat")
    if not os.path.exists(linkage_path) or not os.path.exists(indices_path):
        raise FileNotFoundError(f"{result_dir} に必要なファイルが見つかりません。")

    linkage_matrix = np.loadtxt(linkage_path)
    filtered_indices = pd.read_csv(indices_path, header=None, index_col=0).index
    df_clustered = df_merged.loc[filtered_indices].copy()

    if len(df_clustered) < 2:
        raise ValueError("クラスタリングされたデータが2点未満のため、スコアを計算できません。")
    print(f"評価対象のデータ点数: {len(df_clustered)}")
    return df_clustered, linkage_matrix

def main():
    try:
        analysis_dir = os.path.expanduser(ANALYSIS_OUTPUT_DIRECTORY)
        analysis_config = load_config_from_metadata(analysis_dir)

        df_clustered, linkage_matrix = load_and_prepare_data(
            ORIGINAL_DATA_DIRECTORY, ANALYSIS_OUTPUT_DIRECTORY
        )
        
        print("--- 2. Cythonによるカスタム距離行列の計算 ---")
        data_orig_np = df_clustered.to_numpy(dtype=np.double)
        scaler = StandardScaler()
        data_scalar_std_np = scaler.fit_transform(df_clustered[analysis_config['scalar_cols']])
        all_cols = df_clustered.columns.tolist()
        vector_cols_indices = [[all_cols.index(c) for c in v_group] for v_group in analysis_config['vector_cols_list']]

        # ★★★ 変更点: 'calculate_distances_cython' モジュールを呼び出す ★★★
        distance_matrix = calculate_distances_cython.calculate_full_distance_matrix(
            data_orig_np, data_scalar_std_np, vector_cols_indices, analysis_config['weights']
        )
        print("  距離行列の計算完了。")

        print("\n--- 3. 各クラスター数でのシルエットスコアの計算 ---")
        results = []
        for n_clusters in CLUSTER_RANGE_TO_EVALUATE:
            if len(df_clustered) <= n_clusters:
                print(f"  クラスター数 {n_clusters}: スキップ (データ点数が不足)")
                continue
            labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            if len(np.unique(labels)) < 2:
                print(f"  クラスター数 {n_clusters}: スキップ (生成されたクラスタが1種類のみ)")
                continue
            score = silhouette_score(distance_matrix, labels, metric='precomputed')
            results.append((n_clusters, score))
            print(f"  クラスター数 {n_clusters}: シルエットスコア = {score:.4f}")

        if not results:
            print("\n有効なシルエットスコアが計算されませんでした。")
            return
            
        best_result = max(results, key=lambda item: item[1])
        
        print("\n" + "="*50)
        print("                 --- 結果概要 ---")
        print("="*50)
        print("  クラスター数 | シルエットスコア")
        print("  ----------------|-----------------")
        for n_clusters, score in results:
            print(f"    {n_clusters:<12} | {score:.4f}")
        print("="*50)
        print(f"  最適なクラスタ数 (最高スコア): {best_result[0]}")
        print(f"  最高シルエットスコア: {best_result[1]:.4f}")
        print("="*50)

        score_path = os.path.join(analysis_dir, "silhouette_score_summary.txt")
        with open(score_path, 'w') as f:
            f.write("Silhouette Score Analysis\n")
            f.write(f"Analyzed Directory: {analysis_dir}\n\n")
            f.write("="*30 + "\n")
            f.write("Num Clusters | Silhouette Score\n")
            f.write("-------------|-----------------\n")
            for n_clusters, score in results:
                f.write(f"{n_clusters:<12} | {score:<.4f}\n")
            f.write("="*30 + "\n\n")
            f.write(f"Best number of clusters found: {best_result[0]}\n")
            f.write(f"Best silhouette score: {best_result[1]:.4f}\n")
        print(f"\n結果の概要を {score_path} に保存しました。")

    except Exception as e:
        print(f"\nエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

