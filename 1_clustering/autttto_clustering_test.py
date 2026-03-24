#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# auto_clustering_fixed_manual_threshold.py
# PDTIによって出力されたデータsnap.datを使ってクラスタリングを行うスクリプト
#
# MODIFICATION:
# - This script is configured for a fixed analysis pipeline for testing.
# - Analysis is fixed to use "custom_combined" distance metric.
# - Threshold determination is fixed to manual user input ("threshold_arbitrary").
# - To change parameters, edit the 'fixed_config' dictionary in the main() function.

print("現在テストコードです。 スクリプトが開始されました。ログ確認用。")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import itertools
import re
from datetime import datetime
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import silhouette_score

# --- ヘルパー関数 ---
def ensure_dir(directory_path):
    """指定されたパスのディレクトリが存在しない場合、作成する"""
    os.makedirs(directory_path, exist_ok=True)

def generate_output_dirname_from_config(config, arbitrary_threshold_val=None):
    """設定辞書から短縮された分かりやすいディレクトリ名を生成する (通常モード用)"""
    parts = []
    filter_method = config.get('filter_method', 'none')
    if filter_method == 'arbitrary_percentage':
        q_val = config.get('filter_quantile', 0)
        parts.append(f"filt-arbP{q_val*100:.0f}")
    elif filter_method == 'skewness_kurtosis':
        parts.append("filt-skKur")
    elif filter_method == 'advanced_filtering':
        ext_perc = config.get('target_extraction_percentage',0)
        parts.append(f"filt-adv{ext_perc*100:.0f}")
    else:
        parts.append(f"filt-{filter_method[:3]}")

    var_method = config.get('variable_selection_method', 'none')
    if var_method == 'custom_columns':
        parts.append("var-CUSTOM")
    elif "matrix_components_" in var_method:
        base_name = "var-mat"
        suffix_key = var_method.replace("matrix_components_", "")
        if suffix_key == "only": suffix_display = "ONLY"
        elif suffix_key == "nm": suffix_display = "NM"
        elif suffix_key == "nmtw": suffix_display = "NMTW"
        elif suffix_key == "nmtw_sliprate": suffix_display = "NMTWSL"
        else: suffix_display = suffix_key.upper()[:4]
        parts.append(f"{base_name}{suffix_display}")
    elif "nv6_components_" in var_method:
        base_name = "var-nv6"
        suffix_key = var_method.replace("nv6_components_", "")
        if suffix_key == "only": suffix_display = "ONLY"
        elif suffix_key == "nm": suffix_display = "NM"
        elif suffix_key == "nmtw": suffix_display = "NMTW"
        elif suffix_key == "nmtw_sliprate": suffix_display = "NMTWSL"
        else: suffix_display = suffix_key.upper()[:4]
        parts.append(f"{base_name}{suffix_display}")
    else:
        parts.append(f"var-{var_method[:6]}")
    
    pca_method = config.get('pca_method', 'none')
    if pca_method == 'pca_80': parts.append("pca80")
    elif pca_method == 'pca_90': parts.append("pca90")
    elif pca_method == 'pca_none': parts.append("pcaNo")
    elif pca_method == 'pca_n_components':
        n_comp = config.get('pca_n_components_value', 'Err')
        parts.append(f"pcaN{n_comp}")
    else: parts.append(f"pca-{pca_method[:3]}")

    std_method = config.get('standardization_method', 'none')
    if std_method == 'standardization_on': parts.append("stdOn")
    elif std_method == 'standardization_off': parts.append("stdOff")
    else: parts.append(f"std-{std_method[:3]}")

    dist_metric = config.get('distance_metric', 'none')
    if dist_metric == 'cosine': parts.append("distCos")
    elif dist_metric == 'euclidean': parts.append("distEuc")
    elif dist_metric == 'custom_combined':
        weights = config.get('distance_weights', {'d1': 0, 'd2': 0, 'd3': 0})
        w_str = f"w{weights['d1']*10:.0f}-{weights['d2']*10:.0f}-{weights['d3']*10:.0f}"
        parts.append(f"distCUST{w_str}")
    else: parts.append(f"dist-{dist_metric[:3]}")
    
    thresh_method = config.get('threshold_method', 'none')
    if thresh_method == 'threshold_arbitrary':
        if config.get('run_all_combinations_mode') and config.get('target_cluster_count_for_arbitrary'):
            k_val = config.get('target_cluster_count_for_arbitrary')
            parts.append(f"thrArbK{k_val}")
        elif arbitrary_threshold_val is not None:
            parts.append(f"thrArb{arbitrary_threshold_val:.2f}".replace(".","p"))
        else:
            parts.append("thrArb")
    elif thresh_method == 'threshold_silhouette': parts.append("thrSil")
    elif thresh_method == 'threshold_elbow': parts.append("thrElb")
    elif thresh_method == 'threshold_hybrid': parts.append("thrHyb")
    else: parts.append(f"thr-{thresh_method[:3]}")
    
    return "_".join(parts)

def save_metadata(output_dir, filename="metadata.txt", config=None, comments=None):
    """実行設定やコメントをメタデータファイルに保存する"""
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Execution Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("--- Configuration ---\n")
        if config:
            if config.get('variable_selection_method') == 'custom_columns':
                    f.write(f"custom_columns: {config.get('custom_columns')}\n")
            if config.get('distance_metric') == 'custom_combined':
                f.write(f"distance_weights: {config.get('distance_weights')}\n")
            for key, value in config.items():
                if key not in ['custom_columns', 'distance_weights']:
                    f.write(f"{key}: {value}\n")
        if comments:
            f.write("\n--- Comments ---\n")
            f.write(comments)
    print(f"メタデータを {filepath} に保存しました。")

def read_run_pdti_ta_params(working_directory):
    """run_PDTI_TA.bash からパラメータを読み取る"""
    bash_file_path = os.path.join(working_directory, "const", "run_PDTI_TA.bash")
    params = {}
    if not os.path.exists(bash_file_path):
        raise FileNotFoundError(f"{bash_file_path} が見つかりません。")

    print(f"{bash_file_path} からパラメータを読み取ります...")
    try:
        with open(bash_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            match_xx = re.search(r"^\s*XX\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_xx: params['n_interval'] = float(match_xx.group(1))
            else: raise ValueError("run_PDTI_TA.bash から XX (n_interval) が見つかりません。")
            match_yy = re.search(r"^\s*YY\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_yy: params['m_interval'] = float(match_yy.group(1))
            else: raise ValueError("run_PDTI_TA.bash から YY (m_interval) が見つかりません。")
            match_tr = re.search(r"^\s*TR\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_tr: params['delta_tw'] = float(match_tr.group(1))
            else: raise ValueError("run_PDTI_TA.bash から TR (delta_tw) が見つかりません。")
        print(f"  読み取り成功: n_interval={params['n_interval']}, m_interval={params['m_interval']}, delta_tw={params['delta_tw']}")
        return params
    except Exception as e:
        print(f"エラー: {bash_file_path} の読み取りまたは解析中にエラーが発生しました: {e}")
        raise

# --- データ処理ステージ関数 ---
def stage_load_data(wd, snap_y_file="snap_y.dat", n_vector_file="n_vector.dat"):
    print(f"\n--- 1. データの読み込み ---")
    print(f"作業ディレクトリ: {wd}")
    ps_snap = os.path.join(wd, snap_y_file)
    if not os.path.exists(ps_snap): raise FileNotFoundError(f"データファイルが見つかりません: {ps_snap}")
    df_snap_original = pd.read_table(ps_snap, sep=r'\s+', header=None)
    df_snap_original.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                                'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                'trendp', 'trendt', 'trendb', 'plungp',
                                'plungt', 'plungb', 'NDC']
    print(f"{snap_y_file} を読み込みました。形状: {df_snap_original.shape}")
    ps_nvector = os.path.join(wd, n_vector_file)
    if not os.path.exists(ps_nvector): raise FileNotFoundError(f"データファイルが見つかりません: {ps_nvector}")
    snap4_original = pd.read_table(ps_nvector, sep=r'\s+', header=None)
    snap4_original.columns = ['n', 'm', 'tw','t', 'x_coord', 'y_coord',
                                'n1_(1)','n1_(2)','n1_(3)',
                                'n2_(1)','n2_(2)','n2_(3)',
                                'sliprate']
    snap4_original.rename(columns={'n': 'x_grid', 'm': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    print(f"{n_vector_file} を読み込みました。形状: {snap4_original.shape}")
    return df_snap_original, snap4_original

def stage_advanced_filtering(df_input_orig, n_interval, m_interval, delta_tw, target_extraction_percentage):
    """主要な時空間点をフィルタリングする新しいステージ"""
    print("  高度なフィルタリング処理を開始します...")
    df_input = df_input_orig.copy()
    df_input.rename(columns={'x_grid': 'x', 'y_grid': 'y', 'tw_grid': 'tw'}, inplace=True)

    df_input['S'] = df_input['sliprate']
    df_input['S_xy'] = df_input.groupby(['x', 'y'])['sliprate'].transform('sum')
    df_input['Ctg'] = df_input['sliprate'] / df_input['S_xy']
    df_input['Ctg'].fillna(0, inplace=True)
    df_input['S_tw'] = df_input.groupby('tw')['sliprate'].transform('sum')
    df_input['Cgt'] = df_input['sliprate'] / df_input['S_tw']
    df_input['Cgt'].fillna(0, inplace=True)

    df_input.sort_values(by=['x', 'y', 'tw'], inplace=True)
    df_input['sliprate_prev_tw'] = df_input.groupby(['x', 'y'])['sliprate'].shift(1)
    dt_prime_values = []
    for _, row in df_input.iterrows():
        if pd.isna(row['sliprate_prev_tw']):
            dt_prime_values.append(row['sliprate'] / delta_tw if delta_tw != 0 else 0)
        else:
            dt_prime_values.append(abs(row['sliprate'] - row['sliprate_prev_tw']) / delta_tw if delta_tw != 0 else 0)
    df_input['Dt_prime'] = dt_prime_values
    df_input['Dt_prime'].fillna(0, inplace=True)

    df_input['Ds'] = 0.0
    unique_tws = df_input['tw'].unique()
    unique_tws.sort()
    for current_tw in tqdm(unique_tws, desc="  Calculating Ds (spatial gradient) for each tw", leave=False):
        df_tw = df_input[df_input['tw'] == current_tw].copy()
        if df_tw.empty: continue
        try:
            sliprate_grid = df_tw.pivot_table(index='y', columns='x', values='sliprate', fill_value=0)
        except Exception as e_pivot:
            print(f"    警告: tw={current_tw} でpivot_table作成エラー: {e_pivot}。この時間窓のDs計算をスキップします。")
            continue
        if sliprate_grid.empty or sliprate_grid.shape[0] < 2 or sliprate_grid.shape[1] < 2:
            print(f"    警告: tw={current_tw} でsliprate_gridが小さすぎるためDs計算をスキップ。形状: {sliprate_grid.shape}")
            continue
        if n_interval == 0 or m_interval == 0:
            print(f"    警告: tw={current_tw} で n_interval または m_interval が0のためDs計算をスキップ。")
            continue
        gy, gx = np.gradient(sliprate_grid.values, m_interval, n_interval)
        ds_grid_values = np.sqrt(gx**2 + gy**2)
        for y_idx, y_val in enumerate(sliprate_grid.index):
            for x_idx, x_val in enumerate(sliprate_grid.columns):
                mask = (df_input['x'] == x_val) & (df_input['y'] == y_val) & (df_input['tw'] == current_tw)
                if mask.sum() > 0 :
                    df_input.loc[mask, 'Ds'] = ds_grid_values[y_idx, x_idx]

    metrics_to_normalize = ['S', 'Ctg', 'Cgt', 'Dt_prime', 'Ds']
    scaler = MinMaxScaler()
    for metric in metrics_to_normalize:
        col_name_norm = f"norm_{metric}"
        if df_input[metric].empty or df_input[metric].min() == df_input[metric].max():
            df_input[col_name_norm] = 0.5 if not df_input[metric].empty else 0.0
        else:
            df_input[col_name_norm] = scaler.fit_transform(df_input[[metric]])
        df_input[col_name_norm].fillna(0.5, inplace=True)

    df_input['S_total'] = (df_input['norm_S'] + df_input['norm_Ctg'] +
                           df_input['norm_Cgt'] + df_input['norm_Dt_prime'] +
                           df_input['norm_Ds']) / 5.0

    threshold_s_total = df_input['S_total'].quantile(1.0 - target_extraction_percentage)
    print(f"  S_total の閾値 ({1.0 - target_extraction_percentage:.2f} クォンタイル): {threshold_s_total:.4f}")
    df_input['is_major'] = df_input['S_total'] >= threshold_s_total

    df_input.rename(columns={'x': 'x_grid', 'y': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    filtered_df = df_input[df_input['is_major'] == True].copy()
    return df_input, filtered_df

def stage_filter_data(snap4_df, method, quantile_val=0.75, advanced_filter_params=None):
    """データフィルタリング"""
    print(f"\n--- 2. データフィルタリング ---")
    print(f"手法: {method}")
    df_processed_adv_filter_for_saving = None

    if method == 'arbitrary_percentage':
        threshold_sliprate = snap4_df['sliprate'].quantile(quantile_val)
        filtered_df = snap4_df[snap4_df['sliprate'] > threshold_sliprate].copy()
        print(f"  sliprate > {threshold_sliprate:.4f} (quantile={quantile_val}) でフィルタリングしました。")
        print(f"  フィルタリング前: {snap4_df.shape[0]}行, フィルタリング後: {filtered_df.shape[0]}行")
    elif method == 'skewness_kurtosis':
        print("  警告: 'skewness_kurtosis' フィルタリングは未実装です。フィルタリングは行われません。")
        filtered_df = snap4_df.copy()
    elif method == 'advanced_filtering':
        if advanced_filter_params is None:
            raise ValueError("高度なフィルタリングには advanced_filter_params が必要です。")
        
        df_processed_adv_filter_for_saving, filtered_df = stage_advanced_filtering(
            snap4_df.copy(),
            advanced_filter_params['n_interval'],
            advanced_filter_params['m_interval'],
            advanced_filter_params['delta_tw'],
            advanced_filter_params['target_extraction_percentage']
        )
        print(f"  高度フィルタリング適用。フィルタリング前: {snap4_df.shape[0]}行, フィルタリング後: {filtered_df.shape[0]}行")
    else:
        print(f"  警告: 未知のフィルタリング手法 '{method}' です。フィルタリングは行われません。")
        filtered_df = snap4_df.copy()
    
    if filtered_df.empty and method != 'advanced_filtering':
        raise ValueError("フィルタリングの結果、データが空になりました。フィルタリング条件を確認してください。")
    
    return filtered_df, df_processed_adv_filter_for_saving


def stage_select_variables(df_snap_original, snap4_filtered, method, custom_columns=None):
    """使用変数決定 (注意箇所1を考慮)"""
    print(f"\n--- 3. 使用変数決定 ---")
    
    if isinstance(method, str):
        method = method.strip()
    else:
        raise TypeError(f"変数選択手法(method)は文字列である必要がありますが、{type(method)}型が渡されました。")

    print(f"手法: {method}")
    
    base_matrix_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    base_nv6_cols = ['n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)']
    common_cols_suffix = {
        "only": [], "nm": ['n', 'm'], "nmtw": ['n', 'm', 'tw'],
        "nmtw_sliprate": ['n', 'm', 'tw', 'sliprate']
    }
    data_for_clustering = pd.DataFrame()

    if method == 'custom_columns':
        if not custom_columns:
            raise ValueError("任意列指定が選択されましたが、列が指定されていません。")
        print(f"  任意指定された列: {custom_columns}")
        
        snap4_temp_keys = snap4_filtered.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'})
        snap_y_cols_to_merge = [col for col in df_snap_original.columns if col not in snap4_temp_keys.columns or col in ['n','m','tw']]
        
        merged_df = pd.merge(snap4_temp_keys,
                             df_snap_original[snap_y_cols_to_merge],
                             on=['n', 'm', 'tw'], how='left')
        
        missing_cols = [col for col in custom_columns if col not in merged_df.columns]
        if missing_cols:
            raise ValueError(f"指定された列の一部がデータに存在しません: {missing_cols}")

        data_for_clustering = merged_df[custom_columns].copy()
        data_for_clustering.index = snap4_filtered.index

    elif method.startswith("matrix_components_"):
        suffix_key = method.replace("matrix_components_", "")
        if suffix_key in common_cols_suffix:
            additional_cols_keys = common_cols_suffix[suffix_key]
            
            if suffix_key == "only":
                 if not snap4_filtered.index.is_unique:
                    print("警告: snap4_filtered のインデックスが一意ではありません。重複を削除します。")
                    snap4_filtered_unique_idx = snap4_filtered[~snap4_filtered.index.duplicated(keep='first')]
                 else:
                    snap4_filtered_unique_idx = snap4_filtered

                 common_indices = df_snap_original.index.intersection(snap4_filtered_unique_idx.index)
                 df_snap_subset = df_snap_original.loc[common_indices]
                 
                 final_cluster_cols = [col for col in base_matrix_cols if col in df_snap_subset.columns]
                 if not final_cluster_cols:
                         raise ValueError(f"matrix_components_only で、df_snap_subsetに使用可能な列がありません。要求列: {base_matrix_cols}")
                 data_for_clustering = df_snap_subset[final_cluster_cols].copy()
                 data_for_clustering = data_for_clustering.reindex(snap4_filtered_unique_idx.index)
            else:
                columns_to_cluster_from_snap = base_matrix_cols + [col for col in additional_cols_keys if col in df_snap_original.columns]
                snap4_temp_keys = snap4_filtered[['x_grid', 'y_grid', 'tw_grid']].copy()
                snap4_temp_keys.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'}, inplace=True)
                
                merge_on_keys = ['n', 'm', 'tw']
                
                if not all(key in df_snap_original.columns for key in merge_on_keys):
                    missing_keys = [key for key in merge_on_keys if key not in df_snap_original.columns]
                    raise ValueError(f"df_snap_original にマージキー {missing_keys} が不足しています。")

                merged_df = pd.merge(snap4_temp_keys,
                                     df_snap_original[columns_to_cluster_from_snap],
                                     on=merge_on_keys, how='left')
                
                final_cluster_cols = base_matrix_cols[:]
                if 'sliprate' in additional_cols_keys and 'sliprate' in merged_df.columns:
                    final_cluster_cols.append('sliprate')
                if 'n' in additional_cols_keys and 'n' in merged_df.columns: final_cluster_cols.append('n')
                if 'm' in additional_cols_keys and 'm' in merged_df.columns: final_cluster_cols.append('m')
                if 'tw' in additional_cols_keys and 'tw' in merged_df.columns: final_cluster_cols.append('tw')
                
                final_cluster_cols = [col for col in final_cluster_cols if col in merged_df.columns]
                if merged_df.shape[0] != snap4_filtered.shape[0] or not final_cluster_cols or merged_df[final_cluster_cols].isnull().any().any():
                    raise ValueError("行列成分型の変数選択で、マージに失敗したか、必要な列が欠損しています。")
                data_for_clustering = merged_df[final_cluster_cols].copy()
                data_for_clustering.index = snap4_filtered.index
        else:
            raise ValueError(f"未知の行列成分型変数選択手法のsuffix: {suffix_key}")

    elif method.startswith("nv6_components_"):
        suffix_key = method.replace("nv6_components_", "")
        if suffix_key in common_cols_suffix:
            additional_cols_keys = common_cols_suffix[suffix_key]
            columns_to_cluster = base_nv6_cols[:]
            if 'n' in additional_cols_keys and 'x_grid' in snap4_filtered.columns: columns_to_cluster.append('x_grid')
            if 'm' in additional_cols_keys and 'y_grid' in snap4_filtered.columns: columns_to_cluster.append('y_grid')
            if 'tw' in additional_cols_keys and 'tw_grid' in snap4_filtered.columns: columns_to_cluster.append('tw_grid')
            if 'sliprate' in additional_cols_keys and 'sliprate' in snap4_filtered.columns: columns_to_cluster.append('sliprate')
            columns_to_cluster = [col for col in columns_to_cluster if col in snap4_filtered.columns]
            if not columns_to_cluster:
                raise ValueError(f"nv6成分型の変数選択で、使用可能な列がsnap4_filteredに見つかりませんでした。")
            data_for_clustering = snap4_filtered[columns_to_cluster].copy()
        else:
            raise ValueError(f"未知のnv6成分型変数選択手法のsuffix: {suffix_key}")
    else:
        raise ValueError(f"未知の変数選択手法カテゴリ: {method}")

    if data_for_clustering.isnull().values.any():
        print("警告: 選択された変数に欠損値が含まれています。欠損値を含む行を削除します。")
        data_for_clustering.dropna(inplace=True)
        if data_for_clustering.empty:
            raise ValueError("変数選択後、欠損値処理によりデータが空になりました。")
    print(f"  使用する変数: {data_for_clustering.columns.tolist()}")
    print(f"  クラスタリング用データの形状: {data_for_clustering.shape}")
    return data_for_clustering, data_for_clustering.columns.tolist()


def stage_standardize_data(df_for_clustering, columns_to_cluster, apply_standardization):
    print(f"\n--- 4. データ標準化 ---")
    scaler_model = None
    target_cols_for_std = [col for col in columns_to_cluster if col in df_for_clustering.columns]
    
    numeric_cols_for_std = df_for_clustering[target_cols_for_std].select_dtypes(include=np.number).columns.tolist()

    if not numeric_cols_for_std:
        print("  警告: 標準化対象の数値列がデータフレームに存在しません。標準化をスキップします。")
        return df_for_clustering, None

    if apply_standardization:
        print(f"  手法: あり (StandardScalerを列 {numeric_cols_for_std} に適用)")
        scaler_model = StandardScaler()
        df_copy = df_for_clustering.copy()
        df_copy[numeric_cols_for_std] = scaler_model.fit_transform(df_copy[numeric_cols_for_std])
        return df_copy, scaler_model
    else:
        print(f"  手法: なし")
        return df_for_clustering, scaler_model

def stage_pca(data_scaled, method, output_dir=None, config=None):
    """主成分分析 (PCA)"""
    if config is None: config = {}
    print(f"\n--- 5. 主成分分析 (PCA) ---")
    print(f"手法: {method}")
    pca_model = None
    data_after_pca = data_scaled

    if method == 'pca_none':
        print("  PCAは実行しません。")
        return data_scaled, None

    data_for_pca_fit = data_scaled
    if isinstance(data_scaled, pd.DataFrame):
        numeric_cols = data_scaled.select_dtypes(include=np.number).columns
        if len(numeric_cols) == 0:
            print("  警告: 数値データ列がないためPCAは実行できません。PCAをスキップします。")
            return data_scaled, None
        if len(numeric_cols) < data_scaled.shape[1]:
            print(f"  警告: PCAは数値列 {list(numeric_cols)} のみに対して実行されます。")
        data_for_pca_fit = data_scaled[numeric_cols]
    
    if data_for_pca_fit.shape[1] < 1:
        print("  警告: データ列がないためPCAは実行できません。PCAをスキップします。")
        return data_scaled, None
    
    max_components = data_for_pca_fit.shape[1]

    if method in ['pca_80', 'pca_90']:
        pca_temp = PCA()
        pca_temp.fit(data_for_pca_fit)
        cumulative_variance_temp = np.cumsum(pca_temp.explained_variance_ratio_)
        current_pca_threshold = 0.80 if method == 'pca_80' else 0.90
        n_components = np.argmax(cumulative_variance_temp >= current_pca_threshold) + 1
        print(f"  累積寄与率 {current_pca_threshold*100:.0f}% を満たす主成分数: {n_components}")
    
    elif method == 'pca_n_components':
        n_components = config.get('pca_n_components_value')
        if n_components is None or n_components < 1:
            raise ValueError("主成分数が指定されていません、または無効な値です。")
        print(f"  ユーザー指定の主成分数: {n_components}")
        if n_components > max_components:
            print(f"  警告: 指定された主成分数 ({n_components}) が利用可能な特徴量数 ({max_components}) を超えています。")
            print(f"  主成分数を {max_components} に調整します。")
            n_components = max_components
    else:
        raise ValueError(f"未知のPCA手法: {method}")

    n_components = max(1, min(n_components, max_components))
    pca_model = PCA(n_components=n_components)
    data_after_pca_transformed = pca_model.fit_transform(data_for_pca_fit)
    
    if isinstance(data_scaled, pd.DataFrame):
        data_after_pca = pd.DataFrame(data_after_pca_transformed, index=data_scaled.index,
                                      columns=[f'PC{i+1}' for i in range(n_components)])
    else:
        data_after_pca = data_after_pca_transformed
    print(f"  PCA適用後のデータ形状: {data_after_pca.shape}")

    if output_dir and pca_model:
        explained_var_ratio = pca_model.explained_variance_ratio_
        cumulative_var_ratio = np.cumsum(explained_var_ratio)
        pc_numbers = [f'PC{i+1}' for i in range(len(explained_var_ratio))]
        
        pca_scores_df = pd.DataFrame({
            'PrincipalComponent': pc_numbers,
            'ExplainedVarianceRatio': explained_var_ratio,
            'CumulativeExplainedVarianceRatio': cumulative_var_ratio
        })
        pca_score_path = os.path.join(output_dir, "PCAscore.csv")
        pca_scores_df.to_csv(pca_score_path, index=False)
        print(f"  PCAの寄与率情報を {pca_score_path} に保存しました。")

    return data_after_pca, pca_model

def calculate_custom_combined_distance(data_df, weights):
    """
    Calculates a custom weighted distance matrix based on different metrics
    for different variable groups.
    """
    print("  カスタム統合距離の計算を開始します...")
    
    # 1. Define variable groups
    cols_d1_cosine = ['n1_(1)', 'n1_(2)', 'n1_(3)']
    cols_d2_cosine = ['n2_(1)', 'n2_(2)', 'n2_(3)']
    # Use original column names for spatio-temporal data before potential renaming
    cols_d3_euclidean = ['n', 'm', 'tw', 'sliprate']
    # Check for renamed grid columns as fallback
    if not all(c in data_df.columns for c in cols_d3_euclidean):
        cols_d3_euclidean_alt = ['x_grid', 'y_grid', 'tw_grid', 'sliprate']
        if all(c in data_df.columns for c in cols_d3_euclidean_alt):
           cols_d3_euclidean = cols_d3_euclidean_alt
        else:
             # Find which columns are missing for a better error message
            missing_cols = [c for c in ['n/x_grid', 'm/y_grid', 'tw/tw_grid', 'sliprate'] if c.split('/')[0] not in data_df.columns and c.split('/')[1] not in data_df.columns]
            raise ValueError(f"カスタム距離計算に必要な列が見つかりません: {missing_cols}")


    # Check if all necessary columns are present
    all_required_cols = cols_d1_cosine + cols_d2_cosine + cols_d3_euclidean
    missing_cols = [col for col in all_required_cols if col not in data_df.columns]
    if missing_cols:
        raise ValueError(f"カスタム距離計算に必要な列がデータにありません: {missing_cols}")

    # Initialize a total distance matrix
    total_dist_matrix = np.zeros((data_df.shape[0], data_df.shape[0]))
    
    # 2. Calculate, Normalize, and Add each distance component
    
    # D1: Cosine Distance for n1 vectors
    if weights['d1'] > 0:
        print(f"    - D1 (Cosine) for {cols_d1_cosine} with weight {weights['d1']:.2f}")
        data_d1 = data_df[cols_d1_cosine].values
        dist_d1_condensed = pdist(data_d1, metric='cosine')
        dist_d1_square = squareform(dist_d1_condensed)
        
        # Normalize the distance matrix
        scaler_d1 = MinMaxScaler()
        dist_d1_normalized = scaler_d1.fit_transform(dist_d1_square)
        total_dist_matrix += weights['d1'] * dist_d1_normalized

    # D2: Cosine Distance for n2 vectors
    if weights['d2'] > 0:
        print(f"    - D2 (Cosine) for {cols_d2_cosine} with weight {weights['d2']:.2f}")
        data_d2 = data_df[cols_d2_cosine].values
        dist_d2_condensed = pdist(data_d2, metric='cosine')
        dist_d2_square = squareform(dist_d2_condensed)
        
        scaler_d2 = MinMaxScaler()
        dist_d2_normalized = scaler_d2.fit_transform(dist_d2_square)
        total_dist_matrix += weights['d2'] * dist_d2_normalized
        
    # D3: Euclidean Distance for spatio-temporal variables
    if weights['d3'] > 0:
        print(f"    - D3 (Euclidean) for {cols_d3_euclidean} with weight {weights['d3']:.2f}")
        data_d3 = data_df[cols_d3_euclidean].values
        
        # Standardize spatio-temporal data before Euclidean distance
        scaler_std = StandardScaler()
        data_d3_std = scaler_std.fit_transform(data_d3)
        
        dist_d3_condensed = pdist(data_d3_std, metric='euclidean')
        dist_d3_square = squareform(dist_d3_condensed)
        
        scaler_d3 = MinMaxScaler()
        dist_d3_normalized = scaler_d3.fit_transform(dist_d3_square)
        total_dist_matrix += weights['d3'] * dist_d3_normalized

    # 3. Convert final combined square matrix back to condensed form for linkage
    print("  カスタム統合距離の計算が完了しました。")
    return squareform(total_dist_matrix, checks=False)

def stage_linkage_clustering(data_for_linkage, distance_metric, config=None):
    print(f"\n--- 6. 階層的クラスタリング (ウォード法) ---")
    print(f"  距離計算法: {distance_metric}")
    
    if distance_metric == 'custom_combined':
        if not isinstance(data_for_linkage, pd.DataFrame):
            raise TypeError("カスタム距離計算にはDataFrameが必要です。")
        
        weights = config.get('distance_weights')
        if weights is None:
            raise ValueError("カスタム距離計算法が選択されましたが、重みが設定されていません。")
            
        distance_matrix_condensed = calculate_custom_combined_distance(
            data_for_linkage,
            weights
        )
    else:
        # Original logic for standard metrics
        values_for_pdist = data_for_linkage
        if isinstance(data_for_linkage, pd.DataFrame):
            numeric_cols = data_for_linkage.select_dtypes(include=np.number).columns
            if numeric_cols.empty: raise ValueError("クラスタリング対象データに数値列がありません。")
            if len(numeric_cols) < data_for_linkage.shape[1]:
                print(f"  警告: linkage計算には数値列 {list(numeric_cols)} のみ使用されます。")
            values_for_pdist = data_for_linkage[numeric_cols].values
        elif not isinstance(data_for_linkage, np.ndarray):
            raise TypeError("クラスタリング対象データの型が不正です。DataFrame または ndarray。")
        
        if values_for_pdist.shape[0] < 2:
            raise ValueError(f"データ点が少なすぎるため ({values_for_pdist.shape[0]}点)、階層的クラスタリングを実行できません。")
        if values_for_pdist.ndim == 1: values_for_pdist = values_for_pdist.reshape(-1, 1)
        
        distance_matrix_condensed = pdist(values_for_pdist, metric=distance_metric)
    
    try:
        linkage_matrix = linkage(distance_matrix_condensed, method='average')
    except Exception as e:
        print(f"エラー: linkage計算中に問題が発生しました: {e}")
        # Add more detailed debug info
        print(f"Condensed matrix shape: {distance_matrix_condensed.shape}, Contains NaN: {np.isnan(distance_matrix_condensed).any()}, Contains Inf: {np.isinf(distance_matrix_condensed).any()}")
        raise e
    
    print(f"  linkage matrix の形状: {linkage_matrix.shape}")
    return linkage_matrix

def stage_determine_cluster_threshold(linkage_matrix, data_for_silhouette, method, distance_metric_for_silhouette, output_dir, config=None, columns_used=None):
    """クラスターの距離しきい値決定"""
    if config is None: config = {}
    print(f"\n--- 7. クラスター距離しきい値決定 ---")
    print(f"  手法: {method}")
    
    fig_dendro_full, ax_dendro_full = plt.subplots(figsize=(16, 10))
    dendrogram(linkage_matrix, ax=ax_dendro_full, color_threshold=0)
    ax_dendro_full.set_title("Hierarchical Clustering Dendrogram (Full)")
    ax_dendro_full.set_xlabel("Sample index or (Cluster size)")
    ax_dendro_full.set_ylabel("Distance")
    dendrogram_full_path = os.path.join(output_dir, "dendrogram_full_for_thresholding.png")
    plt.savefig(dendrogram_full_path)
    print(f"  完全なデンドログラムを {dendrogram_full_path} に保存しました。")
    
    if not config.get('run_all_combinations_mode', False):
        plt.show()
    plt.close(fig_dendro_full)

    final_distance_threshold = None
    user_input_threshold = None

    if method == 'threshold_arbitrary':
        target_k_for_arbitrary = config.get('target_cluster_count_for_arbitrary')
        if config.get('run_all_combinations_mode') and target_k_for_arbitrary is not None:
            print(f"  全通りモード: 目標クラスタ数 {target_k_for_arbitrary} から距離しきい値を計算します。")
            if linkage_matrix.shape[0] + 1 < target_k_for_arbitrary or target_k_for_arbitrary < 2:
                raise ValueError(f"目標クラスタ数 {target_k_for_arbitrary} は不適切です。サンプル数: {linkage_matrix.shape[0]+1}")
            
            sorted_distances = np.sort(linkage_matrix[:, 2])
            if target_k_for_arbitrary > len(sorted_distances):
                final_distance_threshold = sorted_distances[-1] if len(sorted_distances) > 0 else 0.1
            else:
                final_distance_threshold = sorted_distances[-target_k_for_arbitrary]
            
            if final_distance_threshold == 0:
                non_zero_distances = sorted_distances[sorted_distances > 0]
                if len(non_zero_distances) > 0:
                    if -target_k_for_arbitrary + 1 < 0 and sorted_distances[-target_k_for_arbitrary + 1] > 0:
                        final_distance_threshold = sorted_distances[-target_k_for_arbitrary + 1] / 2.0
                        if final_distance_threshold == 0: final_distance_threshold = 1e-6
                    else:
                        final_distance_threshold = non_zero_distances[0] * 0.5 if non_zero_distances[0] * 0.5 > 0 else 1e-6
                else:
                    final_distance_threshold = 1e-6
            user_input_threshold = final_distance_threshold
            print(f"  計算された距離しきい値 (目標クラスタ数 {target_k_for_arbitrary}): {final_distance_threshold:.6f}")
        else:
            print(f"  デンドログラム ({dendrogram_full_path}) を参照し、距離しきい値を入力してください。")
            while True:
                try:
                    user_threshold_str = input("  距離しきい値を入力: ")
                    final_distance_threshold = float(user_threshold_str)
                    if final_distance_threshold <= 0:
                        print("  しきい値は正の数である必要があります。")
                        continue
                    user_input_threshold = final_distance_threshold
                    break
                except ValueError:
                    print("  無効な数値です。再度入力してください。")
            print(f"  ユーザー指定の距離しきい値: {final_distance_threshold}")

    elif method == 'threshold_silhouette':
        print("  シルエットスコアに基づいて最適なクラスタ数を探索し、対応する距離しきい値を決定します。")
        actual_data_for_sil_values = data_for_silhouette
        
        use_precomputed_metric = (distance_metric_for_silhouette == 'custom_combined')
        if use_precomputed_metric:
             print("  シルエットスコア計算に、事前計算済みのカスタム統合距離行列を使用します。")
             if not isinstance(data_for_silhouette, pd.DataFrame):
                 raise TypeError("シルエットスコアのカスタム距離計算にはDataFrameが必要です。")
             weights = config.get('distance_weights')
             condensed_matrix = calculate_custom_combined_distance(data_for_silhouette, weights)
             actual_data_for_sil_values = squareform(condensed_matrix)
             distance_metric_for_silhouette = 'precomputed' # Set metric for silhouette_score
        else:
             print(f"  シルエットスコア計算に、'{distance_metric_for_silhouette}' 距離を使用します。")
             if isinstance(data_for_silhouette, pd.DataFrame):
                 numeric_cols_sil = data_for_silhouette.select_dtypes(include=np.number).columns
                 if numeric_cols_sil.empty: raise ValueError("シルエットスコア計算用のデータに数値列がありません。")
                 actual_data_for_sil_values = data_for_silhouette[numeric_cols_sil].values

        if not use_precomputed_metric and actual_data_for_sil_values.ndim == 1:
            actual_data_for_sil_values = actual_data_for_sil_values.reshape(-1,1)
        if actual_data_for_sil_values.shape[0] < 2:
            raise ValueError("シルエットスコア分析に必要なデータ点数が不足しています（2点未満）。")
        
        upper_k_limit = min(11, actual_data_for_sil_values.shape[0])
        if upper_k_limit <=2 :
            raise ValueError(f"データ点数が非常に少ないため({actual_data_for_sil_values.shape[0]})、シルエット分析で複数のkを試せません。")
        cluster_range_list = list(range(3, upper_k_limit))
        sil_scores = []
        for k_val in tqdm(cluster_range_list, desc="Silhouette Analysis", leave=False):
            try:
                labels = fcluster(linkage_matrix, k_val, criterion='maxclust')
                if len(np.unique(labels)) < 2 :
                    sil_scores.append(-1.0)
                    continue
                score = silhouette_score(actual_data_for_sil_values, labels, metric=distance_metric_for_silhouette)
                sil_scores.append(score)
            except ValueError as e:
                print(f"  警告: k={k_val} でシルエットスコア計算エラー: {e}")
                sil_scores.append(-2.0)
        
        if not sil_scores or max(sil_scores) < -0.9 :
            raise ValueError("シルエットスコアから最適なクラスタ数を決定できませんでした。有効なスコアが得られません。")

        if output_dir:
            silhouette_scores_df = pd.DataFrame({
                'NumberOfClusters': cluster_range_list,
                'SilhouetteScore': sil_scores
            })
            slietscore_path = os.path.join(output_dir, "slietscore.csv")
            silhouette_scores_df.to_csv(slietscore_path, index=False)
            print(f"  シルエットスコアのデータを {slietscore_path} に保存しました。")

        fig_sil, ax_sil = plt.subplots(figsize=(8, 4))
        ax_sil.plot(cluster_range_list, sil_scores, marker='o')
        ax_sil.set_xlabel("Number of clusters")
        ax_sil.set_ylabel("Silhouette Score")
        ax_sil.set_title("Silhouette analysis")
        ax_sil.grid(True)

        if columns_used:
            display_text = "columns_used: " + ", ".join(columns_used)
            fig_sil.text(0.5, 0.01, display_text, ha='center', va='bottom', fontsize=8, wrap=True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

        silhouette_plot_path = os.path.join(output_dir, "silhouette_scores.png")
        plt.savefig(silhouette_plot_path)
        print(f"  シルエットスコアのプロットを {silhouette_plot_path} に保存しました。")
        if not config.get('run_all_combinations_mode', False): plt.show()
        plt.close(fig_sil)

        optimal_k_idx = np.argmax(sil_scores)
        optimal_k = cluster_range_list[optimal_k_idx]
        print(f"  最適なクラスタ数（シルエットスコア最大）: {optimal_k} (スコア: {sil_scores[optimal_k_idx]:.4f})")
        
        sorted_distances = np.sort(linkage_matrix[:, 2])
        if optimal_k > len(sorted_distances):
            final_distance_threshold = sorted_distances[-1] if len(sorted_distances) > 0 else 0.1
        else:
            final_distance_threshold = sorted_distances[-optimal_k]
        print(f"  推奨される距離しきい値 (シルエットスコアベース): {final_distance_threshold:.6f}")

    elif method == 'threshold_elbow':
        raise NotImplementedError("エルボー法は未実装です。")
    elif method == 'threshold_hybrid':
        raise NotImplementedError("ハイブリッド法は未実装です。")
    else:
        raise ValueError(f"未知のしきい値決定手法: {method}")
    
    return final_distance_threshold, user_input_threshold

def _generate_clset_for_original_data(clustered_data_indices, linkage_df_with_clno):
    clset_map = {}
    num_original_samples = len(clustered_data_indices)
    if num_original_samples == 0: return clset_map
    clmax = linkage_df_with_clno['clNo'].max()
    original_idx_to_linkage_idx = {orig_idx: i for i, orig_idx in enumerate(clustered_data_indices)}

    for original_data_idx in tqdm(clustered_data_indices, desc="Generating clset strings", leave=False):
        if original_data_idx not in original_idx_to_linkage_idx: continue
        linkage_inner_idx = original_idx_to_linkage_idx[original_data_idx]
        clsb_ids = []
        temp_current_id = linkage_inner_idx
        processed_clnos_path = set()
        while temp_current_id != clmax and temp_current_id not in processed_clnos_path:
            linkage_df_with_clno_int = linkage_df_with_clno.astype({'ind1': int, 'ind2': int, 'clNo': int})
            matching_rows = linkage_df_with_clno_int[
                (linkage_df_with_clno_int['ind1'] == temp_current_id) | \
                (linkage_df_with_clno_int['ind2'] == temp_current_id)
            ]
            if matching_rows.empty:
                if temp_current_id == clmax: break
                break
            next_clno = matching_rows.iloc[0]['clNo']
            clsb_ids.append(next_clno)
            processed_clnos_path.add(temp_current_id)
            temp_current_id = next_clno
            if len(clsb_ids) > 2 * num_original_samples:
                print(f"警告: clset生成でループが長すぎます。ID {original_data_idx}, clsb {clsb_ids}")
                break
        clset_map[original_data_idx] = f"[{', '.join(map(str, clsb_ids))}]" if clsb_ids else "[0]"
    return clset_map

def stage_assign_clusters_and_save_results(
    df_snap_original_full,
    df_n_vector_clustered_data,
    linkage_result_matrix,
    final_cluster_threshold,
    output_dir,
    run_config
):
    print(f"\n--- 8. クラスター割り当てと結果保存 ---")
    print(f"  使用する距離しきい値: {final_cluster_threshold:.6f}")

    num_original_samples_in_clustering = df_n_vector_clustered_data.shape[0]
    if num_original_samples_in_clustering == 0:
        print("  警告: クラスタリング対象データが0件のため、クラスタ割り当てと結果保存をスキップします。")
        df_for_output = df_snap_original_full.copy()
        df_for_output['clset'] = "[0]"
        df_for_output['label'] = 0
    else:
        cl_original_style_df = pd.DataFrame(linkage_result_matrix, columns=['ind1', 'ind2', 'dist', 'NoD'])
        cl_original_style_df['clNo'] = np.arange(
            num_original_samples_in_clustering,
            num_original_samples_in_clustering + len(linkage_result_matrix)
        )
        cl_original_style_df[['ind1', 'ind2']] = cl_original_style_df[['ind1', 'ind2']].astype(int)
        clset_map = _generate_clset_for_original_data(
            df_n_vector_clustered_data.index.tolist(),
            cl_original_style_df.copy()
        )
        df_for_output = df_snap_original_full.copy()
        df_for_output['clset'] = df_for_output.index.map(clset_map).fillna("[0]")
        raw_fcluster_labels = fcluster(linkage_result_matrix, final_cluster_threshold, criterion='distance')
        unique_raw_labels = sorted(list(np.unique(raw_fcluster_labels)))
        label_mapping = {raw_label: i+1 for i, raw_label in enumerate(unique_raw_labels)}
        final_mapped_labels = np.array([label_mapping[rl] for rl in raw_fcluster_labels])
        labels_df_temp = pd.DataFrame({'label': final_mapped_labels}, index=df_n_vector_clustered_data.index)
        df_for_output['label'] = pd.NA
        common_indices_for_label = df_for_output.index.intersection(labels_df_temp.index)
        df_for_output.loc[common_indices_for_label, 'label'] = labels_df_temp.loc[common_indices_for_label, 'label']
        df_for_output['label'] = df_for_output['label'].fillna(0).astype(int)

    print(f"  割り当てられた最終クラスタラベル (元コード風) のユニーク数: {df_for_output['label'].nunique()}")
    print(f"  各クラスタのメンバー数 (0は未分類):\n{df_for_output['label'].value_counts().sort_index()}")

    output_path_snap = os.path.join(output_dir, "clusteringSnap.dat")
    df_for_output.to_csv(output_path_snap, sep=' ', index=False, header=False)
    print(f"  元コード形式の clusteringSnap.dat を {output_path_snap} に保存しました。")
    output_path_snap2 = os.path.join(output_dir, "clusteringSnap2.dat")
    df_for_output.to_csv(output_path_snap2, sep='\t', index=False, header=False)
    print(f"  元コード形式の clusteringSnap2.dat を {output_path_snap2} に保存しました。")
    output_path_info = os.path.join(output_dir, "clusteringinfo.dat")
    df_for_output['label'].to_csv(output_path_info, sep='\t', index=False, header=False)
    print(f"  元コード形式の clusteringinfo.dat を {output_path_info} に保存しました。")
    
    if num_original_samples_in_clustering > 0:
        fig_final_dendro, ax_final_dendro = plt.subplots(figsize=(18, 12))
        dendrogram(
            linkage_result_matrix, ax=ax_final_dendro,
            color_threshold=final_cluster_threshold, truncate_mode=None,
        )
        title_str = f"Hierarchical Clustering Dendrogram (Threshold: {final_cluster_threshold:.4f})"
        ax_final_dendro.set_title(title_str)
        ax_final_dendro.set_xlabel("Sample index or (Cluster size)")
        ax_final_dendro.set_ylabel("Distance")
        ax_final_dendro.axhline(y=final_cluster_threshold, color='r', linestyle='--')
        info_text = "Run Configuration:\n"
        info_text += f"  Filter: {run_config.get('filter_method')} ({run_config.get('filter_quantile', 'N/A') if run_config.get('filter_method') == 'arbitrary_percentage' else run_config.get('target_extraction_percentage','N/A') if run_config.get('filter_method') == 'advanced_filtering' else 'N/A'})\n"
        info_text += f"  Variables: {run_config.get('variable_selection_method')}\n"
        if run_config.get('variable_selection_method') == 'custom_columns':
            info_text += f"    Cols: {run_config.get('custom_columns')}\n"
        info_text += f"  PCA: {run_config.get('pca_method')}\n"
        info_text += f"  Standardize: {run_config.get('standardization_method')}\n"
        info_text += f"  Distance: {run_config.get('distance_metric')}\n"
        if run_config.get('distance_metric') == 'custom_combined':
            info_text += f"    Weights: {run_config.get('distance_weights')}\n"
        info_text += f"  Threshold Method: {run_config.get('threshold_method')}\n"
        if run_config.get('threshold_method') == 'threshold_arbitrary':
            if run_config.get('target_cluster_count_for_arbitrary'):
                info_text += f"  Arbitrary Target K: {run_config.get('target_cluster_count_for_arbitrary')}\n"
            elif hasattr(run_config, 'arbitrary_threshold_value_for_plot'):
                info_text += f"  Arbitrary Threshold Val: {run_config.arbitrary_threshold_value_for_plot:.4f}\n"
        fig_final_dendro.text(0.01, 0.98, info_text, transform=fig_final_dendro.transFigure,
                              fontsize=8, verticalalignment='top',
                              bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        dendrogram_final_path = os.path.join(output_dir, "dendrogram_final.png")
        plt.savefig(dendrogram_final_path)
        print(f"  最終デンドログラム (情報付き) を {dendrogram_final_path} に保存しました。")
        plt.close(fig_final_dendro)
    else:
        print("  クラスタリング対象データがなかったため、最終デンドログラムは保存されません。")

    files_to_copy = ["snap_yr.dat", "snap2_y.dat", "fort.40",
                     "faultline.dat", "mrf.dat", "snap2.dat","rigid_amp.info"]
    print(f"\n  指定されたファイルのコピーを開始します...")
    for filename in files_to_copy:
        source_path = os.path.join(run_config['working_directory'], filename)
        destination_path = os.path.join(output_dir, filename)
        if os.path.exists(source_path):
            try:
                shutil.copy2(source_path, destination_path)
                print(f"    コピー成功: {filename} -> {destination_path}")
            except Exception as e:
                print(f"    警告: {filename} のコピー中にエラーが発生しました: {e}")
        else:
            print(f"    警告: コピー元ファイルが見つかりません: {source_path}。スキップします。")

# --- メイン処理関数 ---
def run_single_pipeline(config, df_snap_orig, snap4_orig, pipeline_run_count=None):
    """単一の設定でクラスタリングパイプラインを実行する"""
    print(f"\n=== パイプライン開始 (実行カウント: {pipeline_run_count if pipeline_run_count is not None else 'N/A'}) ===")
    arbitrary_threshold_value_for_filename = None
    current_run_output_dir = ""
    df_processed_adv_filter_details = None
    
    try:
        advanced_filter_params_dict = None
        if config['filter_method'] == 'advanced_filtering':
            bash_params = read_run_pdti_ta_params(config['working_directory'])
            advanced_filter_params_dict = {
                'n_interval': bash_params['n_interval'],
                'm_interval': bash_params['m_interval'],
                'delta_tw': bash_params['delta_tw'],
                'target_extraction_percentage': config['target_extraction_percentage']
            }
        
        snap4_filtered, df_processed_adv_filter_details = stage_filter_data(
            snap4_orig.copy(),
            config['filter_method'],
            quantile_val=config.get('filter_quantile', 0.75),
            advanced_filter_params=advanced_filter_params_dict
        )
        
        data_for_clustering, actual_cols_used = stage_select_variables(
            df_snap_orig.copy(), snap4_filtered.copy(),
            config['variable_selection_method'],
            custom_columns=config.get('custom_columns')
        )

        if data_for_clustering.empty:
            print("警告: 変数選択後にクラスタリング対象データが空になりました。この組み合わせの処理をスキップします。")
            if config.get('output_directory_path_override'):
                current_run_output_dir = config['output_directory_path_override']
                ensure_dir(current_run_output_dir)
                save_metadata(current_run_output_dir, config=config, comments="Skipped: No data after variable selection.")
            return

        if config['distance_metric'] == 'custom_combined':
            print("カスタム距離モードのため、このステージでの一括標準化をスキップします。")
            data_std = data_for_clustering.copy()
        else:
             data_std, scaler = stage_standardize_data(
                data_for_clustering.copy(), actual_cols_used,
                apply_standardization=(config['standardization_method'] == 'standardization_on')
            )

        run_specific_dirname_temp = generate_output_dirname_from_config(config)
        current_run_output_dir = os.path.join(config['working_directory'], run_specific_dirname_temp)
        
        # In fixed mode, we might not want interactive overwrite handling.
        # Let's just create a unique directory if it exists.
        if os.path.exists(current_run_output_dir):
            print(f"'{current_run_output_dir}' は既に存在するため、タイムスタンプを付加します。")
            current_run_output_dir = f"{current_run_output_dir}_{datetime.now().strftime('%Y%m%d-%H%M%S')}"

        ensure_dir(current_run_output_dir)
        print(f"結果出力先: {current_run_output_dir}")

        data_pca, pca_model = stage_pca(
            data_std.copy(),
            config['pca_method'],
            output_dir=current_run_output_dir,
            config=config
        )
        
        if config['distance_metric'] == 'custom_combined' and config['pca_method'] != 'pca_none':
            print("警告: カスタム距離モードはPCA適用前の変数を使用します。PCAはこのパイプラインでは距離計算に影響しません。")
            data_for_linkage_final = data_for_clustering.copy()
            if 'x_grid' in data_for_linkage_final.columns:
                data_for_linkage_final.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'}, inplace=True)
        else:
            data_for_linkage_final = data_pca
            
        linkage_mat = stage_linkage_clustering(
            data_for_linkage_final,
            config['distance_metric'],
            config=config
        )
        
        if config['distance_metric'] == 'custom_combined':
             sil_data_input = data_for_linkage_final
        elif config['pca_method'] != 'pca_none':
             sil_data_input = data_pca
        else:
             sil_data_input = data_std

        dist_thresh, user_arb_thresh_val = stage_determine_cluster_threshold(
            linkage_mat, sil_data_input,
            method=config['threshold_method'],
            distance_metric_for_silhouette=config['distance_metric'],
            output_dir=current_run_output_dir,
            config=config,
            columns_used=actual_cols_used
        )
        if dist_thresh is None:
            print("しきい値決定がキャンセルされたか失敗したため、パイプラインを終了します。")
            return

        if config['threshold_method'] == 'threshold_arbitrary' and user_arb_thresh_val is not None:
            # Re-generate directory name with the user-provided threshold value for clarity
            final_dirname = generate_output_dirname_from_config(config, user_arb_thresh_val)
            final_output_dir = os.path.join(config['working_directory'], final_dirname)
            if final_output_dir != current_run_output_dir:
                if os.path.exists(final_output_dir):
                    final_output_dir = f"{final_output_dir}_{datetime.now().strftime('%Y%m%d-%H%M%S')}"
                print(f"手動閾値指定のため、ディレクトリ名を {final_output_dir} に変更します。")
                os.rename(current_run_output_dir, final_output_dir)
                current_run_output_dir = final_output_dir
        
        if df_processed_adv_filter_details is not None and current_run_output_dir:
            datascore_cols = ['x_grid', 'y_grid', 'tw_grid', 'sliprate',
                              'S', 'Ctg', 'Cgt', 'Dt_prime', 'Ds',
                              'norm_S', 'norm_Ctg', 'norm_Cgt', 'norm_Dt_prime', 'norm_Ds',
                              'S_total', 'is_major']
            cols_to_save = [col for col in datascore_cols if col in df_processed_adv_filter_details.columns]
            
            datascore_path = os.path.join(current_run_output_dir, "datascore.csv")
            df_processed_adv_filter_details[cols_to_save].to_csv(datascore_path, index=False)
            print(f"  高度フィルタリングのスコア情報を {datascore_path} に保存しました。")


        stage_assign_clusters_and_save_results(
            df_snap_original_full=df_snap_orig.copy(),
            df_n_vector_clustered_data=data_for_clustering.copy(),
            linkage_result_matrix=linkage_mat,
            final_cluster_threshold=dist_thresh,
            output_dir=current_run_output_dir,
            run_config=config
        )
        save_metadata(current_run_output_dir, config=config)
        run_identifier = os.path.basename(current_run_output_dir)
        print(f"=== パイプライン正常終了: {run_identifier} ===")

    except NotImplementedError as e:
        print(f"エラー: 未実装の機能が呼び出されました: {e}")
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: NotImplementedError - {e}")
    except FileNotFoundError as e:
        print(f"エラー: ファイルが見つかりません: {e}")
    except ValueError as e:
        print(f"エラー: データ処理中に値に関する問題が発生しました: {e}")
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: ValueError - {e}")
    except Exception as e:
        import traceback
        print(f"予期せぬエラーが発生しました: {e}")
        print(traceback.format_exc())
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: Unexpected error - {e}\n{traceback.format_exc()}")

def main():
    """メイン実行関数"""
    print("メイン関数が開始されました。")
    default_working_directory = "./"

    working_directory = input(f"1. 作業ディレクトリ (デフォルトの場合Enter) (デフォルト: {default_working_directory}): ") or default_working_directory
    while not os.path.isdir(working_directory):
        print(f"エラー: 指定された作業ディレクトリ '{working_directory}' は存在しないか、ディレクトリではありません。")
        working_directory = input(f"    有効な作業ディレクトリを入力してください (例: {default_working_directory}): ")
        if not working_directory:
            working_directory = default_working_directory
            if not os.path.isdir(working_directory):
                print("デフォルトの作業ディレクトリも無効です。プログラムを終了します。"); exit()
    
    print(f"使用する作業ディレクトリ: {working_directory}")
    
    try:
        df_snap, df_nvector = stage_load_data(working_directory)
    except FileNotFoundError as e:
        print(f"エラー: データファイルの読み込みに失敗しました: {e}")
        return
    except Exception as e:
        print(f"エラー: 予期せぬ問題でデータが読み込めませんでした: {e}")
        return
    
    # --- FIXED CONFIGURATION FOR TESTING ---
    # To change the analysis method, edit this dictionary.
    fixed_config = {
        # --- フィルタリング設定 ---
        "filter_method": "advanced_filtering",
        "target_extraction_percentage": 0.25,
        
        # --- 変数選択設定 ---
        # カスタム距離モードでは、必要な変数がすべて含まれるこの手法に固定されます。
        "variable_selection_method": "nv6_components_nmtw_sliprate",
        
        # --- 標準化設定 ---
        # カスタム距離モードでは無視され、内部で適切に処理されます。
        "standardization_method": "standardization_off",
        
        # --- PCA設定 ---
        # カスタム距離モードでは、変数グループを維持するためPCAは適用できません。
        "pca_method": "pca_none",
        
        # --- 距離計算設定 ---
        # 'cosine', 'euclidean', または 'custom_combined' を選択
        "distance_metric": "custom_combined",
        # 'custom_combined' を使用する場合の重み (合計1.0を推奨)
        "distance_weights": {'d1': 0.5, 'd2': 0.5, 'd3': 0.0}, # D1:n1, D2:n2, D3:時空間
        
        # --- クラスタ閾値決定設定 ---
        # 'threshold_silhouette' (自動) または 'threshold_arbitrary' (手動)
        "threshold_method": "threshold_arbitrary",
        
        # --- 作業ディレクトリ ---
        "working_directory": working_directory,
    }

    print("\n--- 固定設定による分析を実行します ---")
    for key, value in fixed_config.items():
        print(f"  {key}: {value}")
    
    run_single_pipeline(fixed_config, df_snap_orig=df_snap, snap4_orig=df_nvector)


if __name__ == "__main__":
    main()
