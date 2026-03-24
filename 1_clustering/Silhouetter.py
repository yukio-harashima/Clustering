#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Silhouetter.py
#
# 概要:
# auttto_clustering.pyの機能から、シルエットスコアの計算に特化し、
# 最適な手法とクラスタ数を探索するための高速なスクリプト。

print("Silhouetter.py スクリプトが開始されました。") 

import pandas as pd
import numpy as np
import os
import itertools 
import re 
from datetime import datetime
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score

# --- ヘルパー関数 ---
def read_run_pdti_ta_params(working_directory):
    """run_PDTI_TA.bash からパラメータを読み取る"""
    bash_file_path = os.path.join(working_directory, "const", "run_PDTI_TA.bash")
    params = {}
    if not os.path.exists(bash_file_path):
        raise FileNotFoundError(f"パラメータファイルが見つかりません: {bash_file_path}")

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
        return params
    except Exception as e:
        print(f"エラー: {bash_file_path} の読み取りまたは解析中にエラーが発生しました: {e}")
        raise

# --- データ処理ステージ関数 (auttto_clustering.pyより抜粋・軽量化) ---
def stage_load_data(wd, snap_y_file="snap_y.dat", n_vector_file="n_vector.dat"):
    ps_snap = os.path.join(wd, snap_y_file)
    if not os.path.exists(ps_snap): raise FileNotFoundError(f"データファイルが見つかりません: {ps_snap}")
    df_snap_original = pd.read_table(ps_snap, sep=r'\s+', header=None)
    df_snap_original.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                                'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                'trendp', 'trendt', 'trendb', 'plungp', 
                                'plungt', 'plungb', 'NDC']
    ps_nvector = os.path.join(wd, n_vector_file)
    if not os.path.exists(ps_nvector): raise FileNotFoundError(f"データファイルが見つかりません: {ps_nvector}")
    snap4_original = pd.read_table(ps_nvector, sep=r'\s+', header=None)
    snap4_original.columns = ['n', 'm', 'tw','t', 'x_coord', 'y_coord',
                              'n1_(1)','n1_(2)','n1_(3)',
                              'n2_(1)','n2_(2)','n2_(3)',
                              'sliprate']
    snap4_original.rename(columns={'n': 'x_grid', 'm': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    return df_snap_original, snap4_original

def stage_advanced_filtering(df_input_orig, n_interval, m_interval, delta_tw, target_extraction_percentage):
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
    for current_tw in unique_tws:
        df_tw = df_input[df_input['tw'] == current_tw]
        if df_tw.empty: continue
        try:
            sliprate_grid = df_tw.pivot_table(index='y', columns='x', values='sliprate', fill_value=0) 
        except Exception: continue
        if sliprate_grid.empty or sliprate_grid.shape[0] < 2 or sliprate_grid.shape[1] < 2: continue
        if n_interval == 0 or m_interval == 0: continue
        gy, gx = np.gradient(sliprate_grid.values, m_interval, n_interval)
        ds_grid_values = np.sqrt(gx**2 + gy**2)
        for y_idx, y_val in enumerate(sliprate_grid.index):
            for x_idx, x_val in enumerate(sliprate_grid.columns):
                mask = (df_input['x'] == x_val) & (df_input['y'] == y_val) & (df_input['tw'] == current_tw)
                if mask.sum() > 0 : df_input.loc[mask, 'Ds'] = ds_grid_values[y_idx, x_idx]
    metrics_to_normalize = ['S', 'Ctg', 'Cgt', 'Dt_prime', 'Ds']
    scaler = StandardScaler() 
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
    df_input['is_major'] = df_input['S_total'] >= threshold_s_total
    df_input.rename(columns={'x': 'x_grid', 'y': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    filtered_df = df_input[df_input['is_major'] == True].copy()
    return filtered_df

def stage_filter_data(snap4_df, method, quantile_val=0.75, advanced_filter_params=None):
    if method == 'arbitrary_percentage':
        threshold_sliprate = snap4_df['sliprate'].quantile(quantile_val)
        filtered_df = snap4_df[snap4_df['sliprate'] > threshold_sliprate].copy()
    elif method == 'advanced_filtering':
        if advanced_filter_params is None: raise ValueError("advanced_filter_paramsが必要です。")
        filtered_df = stage_advanced_filtering(
            snap4_df.copy(), 
            advanced_filter_params['n_interval'],
            advanced_filter_params['m_interval'],
            advanced_filter_params['delta_tw'],
            advanced_filter_params['target_extraction_percentage']
        )
    else:
        filtered_df = snap4_df.copy()
    if filtered_df.empty: raise ValueError("フィルタリングでデータが空になりました。")
    return filtered_df

def stage_select_variables(df_snap_original, snap4_filtered, method, custom_columns=None):
    if isinstance(method, str): method = method.strip() 
    else: raise TypeError(f"変数選択手法(method)は文字列である必要がありますが、{type(method)}型が渡されました。")
    base_matrix_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    base_nv6_cols = ['n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)']
    common_cols_suffix = {"only": [], "nm": ['n', 'm'], "nmtw": ['n', 'm', 'tw'], "nmtw_sliprate": ['n', 'm', 'tw', 'sliprate']}
    data_for_clustering = pd.DataFrame()
    
    if method == 'custom_columns':
        if not custom_columns: raise ValueError("任意列指定が選択されましたが、列が指定されていません。")
        snap4_temp_keys = snap4_filtered.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'})
        snap_y_cols_to_merge = [col for col in df_snap_original.columns if col not in snap4_temp_keys.columns or col in ['n','m','tw']]
        merged_df = pd.merge(snap4_temp_keys, df_snap_original[snap_y_cols_to_merge], on=['n', 'm', 'tw'], how='left')
        missing_cols = [col for col in custom_columns if col not in merged_df.columns]
        if missing_cols: raise ValueError(f"指定された列の一部がデータに存在しません: {missing_cols}")
        data_for_clustering = merged_df[custom_columns].copy()
        data_for_clustering.index = snap4_filtered.index 
    elif method.startswith("matrix_components_"):
        suffix_key = method.replace("matrix_components_", "")
        if suffix_key in common_cols_suffix:
            additional_cols_keys = common_cols_suffix[suffix_key] 
            if suffix_key == "only":
                snap4_filtered_unique_idx = snap4_filtered[~snap4_filtered.index.duplicated(keep='first')] if not snap4_filtered.index.is_unique else snap4_filtered
                common_indices = df_snap_original.index.intersection(snap4_filtered_unique_idx.index)
                df_snap_subset = df_snap_original.loc[common_indices]
                final_cluster_cols = [col for col in base_matrix_cols if col in df_snap_subset.columns]
                if not final_cluster_cols: raise ValueError(f"matrix_components_only で、df_snap_subsetに使用可能な列がありません。")
                data_for_clustering = df_snap_subset[final_cluster_cols].copy().reindex(snap4_filtered_unique_idx.index)
            else: 
                columns_to_cluster_from_snap = base_matrix_cols + [col for col in additional_cols_keys if col in df_snap_original.columns]
                snap4_temp_keys = snap4_filtered[['x_grid', 'y_grid', 'tw_grid']].copy()
                snap4_temp_keys.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'}, inplace=True)
                merge_on_keys = ['n', 'm', 'tw'] 
                if not all(key in df_snap_original.columns for key in merge_on_keys):
                    raise ValueError(f"df_snap_original にマージキーが不足しています。")
                merged_df = pd.merge(snap4_temp_keys, df_snap_original[columns_to_cluster_from_snap], on=merge_on_keys, how='left')
                final_cluster_cols = base_matrix_cols[:] 
                if 'sliprate' in additional_cols_keys and 'sliprate' in merged_df.columns: final_cluster_cols.append('sliprate')
                if 'n' in additional_cols_keys and 'n' in merged_df.columns: final_cluster_cols.append('n')
                if 'm' in additional_cols_keys and 'm' in merged_df.columns: final_cluster_cols.append('m')
                if 'tw' in additional_cols_keys and 'tw' in merged_df.columns: final_cluster_cols.append('tw')
                final_cluster_cols = [col for col in final_cluster_cols if col in merged_df.columns] 
                if merged_df.shape[0] != snap4_filtered.shape[0] or not final_cluster_cols or merged_df[final_cluster_cols].isnull().any().any():
                     raise ValueError("行列成分型の変数選択で、マージに失敗したか、必要な列が欠損しています。")
                data_for_clustering = merged_df[final_cluster_cols].copy()
                data_for_clustering.index = snap4_filtered.index 
        else: raise ValueError(f"未知の行列成分型変数選択手法のsuffix: {suffix_key}")
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
            if not columns_to_cluster: raise ValueError(f"nv6成分型の変数選択で、使用可能な列がsnap4_filteredに見つかりませんでした。")
            data_for_clustering = snap4_filtered[columns_to_cluster].copy()
        else: raise ValueError(f"未知のnv6成分型変数選択手法のsuffix: {suffix_key}")
    else: raise ValueError(f"未知の変数選択手法カテゴリ: {method}")
    if data_for_clustering.isnull().values.any(): data_for_clustering.dropna(inplace=True)
    if data_for_clustering.empty: raise ValueError("変数選択後、データが空になりました。")
    return data_for_clustering

def stage_standardize_data(df_for_clustering, apply_standardization):
    if apply_standardization:
        numeric_cols = df_for_clustering.select_dtypes(include=np.number).columns.tolist()
        if numeric_cols:
            scaler = StandardScaler()
            df_for_clustering[numeric_cols] = scaler.fit_transform(df_for_clustering[numeric_cols])
    return df_for_clustering

def stage_pca(data_scaled, method, config=None):
    if config is None: config = {}
    if method == 'pca_none': return data_scaled
    data_for_pca_fit = data_scaled.select_dtypes(include=np.number)
    if data_for_pca_fit.shape[1] < 1: return data_scaled
    max_components = data_for_pca_fit.shape[1]
    if method in ['pca_80', 'pca_90']:
        pca_temp = PCA() 
        pca_temp.fit(data_for_pca_fit)
        cumulative_variance_temp = np.cumsum(pca_temp.explained_variance_ratio_)
        current_pca_threshold = 0.80 if method == 'pca_80' else 0.90
        n_components = np.argmax(cumulative_variance_temp >= current_pca_threshold) + 1
    elif method == 'pca_n_components':
        n_components = config.get('pca_n_components_value')
        if n_components is None or n_components < 1: raise ValueError("主成分数が無効です。")
        if n_components > max_components: n_components = max_components
    else: raise ValueError(f"未知のPCA手法: {method}")
    n_components = max(1, min(n_components, max_components))
    pca_model = PCA(n_components=n_components)
    data_after_pca_transformed = pca_model.fit_transform(data_for_pca_fit)
    return pd.DataFrame(data_after_pca_transformed, index=data_scaled.index, columns=[f'PC{i+1}' for i in range(n_components)])

def calculate_best_silhouette(linkage_matrix, data_for_silhouette, distance_metric):
    """シルエットスコアを計算し、最適なクラスタ数とスコアを返す"""
    # 修正: クラスタ数の計算範囲を3-10に限定
    min_k = 3
    max_k = min(11, data_for_silhouette.shape[0] - 1) # 10までなので上限は11
    if max_k <= min_k: return -1.0, -1 # 計算不可
    
    cluster_range = range(min_k, max_k)
    sil_scores = []
    for k_val in cluster_range:
        labels = fcluster(linkage_matrix, k_val, criterion='maxclust')
        if len(np.unique(labels)) < 2 : 
            sil_scores.append(-1.0) 
            continue
        score = silhouette_score(data_for_silhouette, labels, metric=distance_metric)
        sil_scores.append(score)
    if not sil_scores: return -1.0, -1
    best_score = max(sil_scores)
    best_k = list(cluster_range)[np.argmax(sil_scores)]
    return best_score, best_k

def run_all_silhouette_combinations(base_working_dir, df_snap_orig, snap4_orig):
    """全通りの組み合わせでシルエットスコアを計算する"""
    print("\n--- 全手法のシルエットスコア計算開始 ---")
    
    # 修正: パラメータ範囲の定義を指示通りに変更
    filter_methods_options = ['arbitrary_percentage', 'advanced_filtering'] 
    # 90%から0%まで10%刻み -> quantileは0.9, 0.8, ..., 0.0
    filter_quantiles_options = np.arange(0.0, 1.0, 0.1)[::-1].tolist()
    # 抽出率10%から100%まで10%刻み -> target_extraction_percentageは0.1, 0.2, ..., 1.0
    target_extraction_percentage_options = np.arange(0.1, 1.01, 0.1).tolist()
    
    variable_selection_options = [
        "matrix_components_only", "matrix_components_nm", "matrix_components_nmtw", "matrix_components_nmtw_sliprate",
        "nv6_components_only", "nv6_components_nm", "nv6_components_nmtw", "nv6_components_nmtw_sliprate",
        {'method': 'custom_columns', 'columns': ['n1_(1)','n1_(2)','n1_(3)']},
        {'method': 'custom_columns', 'columns': ['n1_(1)','n1_(2)','n1_(3)', 'n', 'm']},
        {'method': 'custom_columns', 'columns': ['n1_(1)','n1_(2)','n1_(3)', 'n', 'm', 'tw']},
        {'method': 'custom_columns', 'columns': ['n1_(1)','n1_(2)','n1_(3)', 'n', 'm', 'tw', 'sliprate']},
    ]
    standardization_options = ['standardization_on', 'standardization_off']
    pca_options = ['pca_80', 'pca_90', 'pca_n_components', 'pca_none']
    pca_n_components_values = [3, 4] 
    distance_metric_options = ['cosine', 'euclidean']

    param_combinations = []
    
    for filt_m in filter_methods_options:
        current_filter_specific_params = [{}]
        if filt_m == 'arbitrary_percentage':
            current_filter_specific_params = [{'filter_quantile': round(q,1)} for q in filter_quantiles_options]
        elif filt_m == 'advanced_filtering':
            current_filter_specific_params = [{'target_extraction_percentage': round(p,1)} for p in target_extraction_percentage_options]
        
        for filt_specific_param_set in current_filter_specific_params:
            for var_sel, std_m, pca_m, dist_m in itertools.product(
                variable_selection_options,
                standardization_options,
                pca_options,
                distance_metric_options
            ):
                current_pca_n_values = [None]
                if pca_m == 'pca_n_components':
                    current_pca_n_values = pca_n_components_values

                for pca_n_val in current_pca_n_values:
                    config = {
                        "working_directory": base_working_dir,
                        "filter_method": filt_m,
                        "standardization_method": std_m,
                        "pca_method": pca_m,
                        "distance_metric": dist_m,
                    }
                    if isinstance(var_sel, dict):
                        config["variable_selection_method"] = var_sel['method']
                        config["custom_columns"] = var_sel['columns']
                    else:
                        config["variable_selection_method"] = var_sel

                    if filt_specific_param_set: config.update(filt_specific_param_set)
                    if pca_n_val is not None: config["pca_n_components_value"] = pca_n_val
                    
                    param_combinations.append(config)

    total_combinations = len(param_combinations)
    print(f"合計 {total_combinations} 通りの組み合わせを計算します。")
    
    results = []
    
    with tqdm(total=total_combinations, desc="シルエットスコア計算進捗") as pbar:
        for i, config_item in enumerate(param_combinations):
            try:
                advanced_filter_params = None
                if config_item['filter_method'] == 'advanced_filtering':
                    bash_params = read_run_pdti_ta_params(config_item['working_directory'])
                    advanced_filter_params = {
                        'n_interval': bash_params['n_interval'],
                        'm_interval': bash_params['m_interval'],
                        'delta_tw': bash_params['delta_tw'],
                        'target_extraction_percentage': config_item['target_extraction_percentage']
                    }
                
                snap4_filtered = stage_filter_data(
                    snap4_orig.copy(), 
                    config_item['filter_method'],
                    quantile_val=config_item.get('filter_quantile', 0.9),
                    advanced_filter_params=advanced_filter_params
                )
                
                data_for_clustering = stage_select_variables(
                    df_snap_orig.copy(), snap4_filtered.copy(), 
                    config_item['variable_selection_method'],
                    custom_columns=config_item.get('custom_columns')
                )

                if data_for_clustering.shape[0] < 20: 
                    pbar.update(1)
                    continue

                data_std = stage_standardize_data(
                    data_for_clustering.copy(), 
                    apply_standardization=(config_item['standardization_method'] == 'standardization_on')
                )
                
                data_pca = stage_pca(
                    data_std.copy(), 
                    config_item['pca_method'],
                    config=config_item
                )
                
                linkage_mat = linkage(pdist(data_pca, metric=config_item['distance_metric']), method='ward')
                
                best_score, best_k = calculate_best_silhouette(linkage_mat, data_pca, config_item['distance_metric'])
                
                if best_k != -1: 
                    results.append({
                        'score': best_score,
                        'k': best_k,
                        'config': config_item
                    })
            except (ValueError, FileNotFoundError):
                pass 
            except Exception:
                pass
            pbar.update(1)
            
    print("\n--- 計算完了。結果をファイルに出力します ---")
    results.sort(key=lambda x: x['score'], reverse=True)
    
    output_path = os.path.join(base_working_dir, "silhouette_scores_summary.txt")
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(f"シルエットスコア計算結果 ({datetime.now().strftime('%Y-%m-%d %H:%M:%S')})\n")
        f.write(f"総組み合わせ数: {total_combinations}, 有効結果数: {len(results)}\n")
        f.write("="*80 + "\n")
        
        for i, res in enumerate(results):
            f.write(f"--- Rank {i+1} ---\n")
            f.write(f"Max Silhouette Score: {res['score']:.6f}\n")
            f.write(f"Optimal Number of Clusters (k): {res['k']}\n")
            f.write("Parameters:\n")
            for key, value in res['config'].items():
                if key != 'working_directory':
                    f.write(f"  - {key}: {value}\n")
            f.write("-" * 40 + "\n\n")
            
    print(f"結果を {output_path} に保存しました。")


def main():
    """メイン実行関数"""
    default_working_directory = "./" 

    working_directory = input(f"作業ディレクトリを入力してください (デフォルト: {default_working_directory}): ") or default_working_directory
    while not os.path.isdir(working_directory):
        print(f"エラー: 指定された作業ディレクトリ '{working_directory}' は存在しません。")
        working_directory = input(f"   有効な作業ディレクトリを入力してください: ")
        if not working_directory: 
            working_directory = default_working_directory
            if not os.path.isdir(working_directory):
                 print("デフォルトの作業ディレクトリも無効です。プログラムを終了します。"); exit()
    
    print(f"使用する作業ディレクトリ: {working_directory}")
    
    try:
        df_snap, df_nvector = stage_load_data(working_directory)
        run_all_silhouette_combinations(working_directory, df_snap_orig=df_snap, snap4_orig=df_nvector)
    except FileNotFoundError as e:
        print(f"エラー: データファイルの読み込みに失敗しました: {e}")
    except Exception as e:
        print(f"エラー: 予期せぬ問題が発生しました: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

