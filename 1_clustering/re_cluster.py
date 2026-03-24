#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# re_cluster.py
# サーバー環境でSHCの計算結果を高速に再クラスタリングするための専用スクリプト。
# run_clustering.txt から設定を読み込んで動作する。

print("SHC Re-Clustering Tool started.")

import pandas as pd
import numpy as np
import matplotlib
# サーバー環境でGUIバックエンドがなくても動作するように設定
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import shutil
import re
from datetime import datetime
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import silhouette_score


# --- Helper Functions ---
def ensure_dir(directory_path):
    """Creates a directory if it does not exist."""
    os.makedirs(directory_path, exist_ok=True)

def generate_output_dirname_from_config(config):
    """Generates a concise, descriptive directory name from the configuration dictionary."""
    parts = []
    filter_method = config.get('filter_method', 'none')
    if filter_method == 'arbitrary_percentage':
        q_val = config.get('filter_value', 0)
        parts.append(f"filt-arbP{q_val*100:.0f}")
    elif filter_method == 'advanced_filtering':
        ext_perc = config.get('filter_value',0)
        parts.append(f"filt-adv{ext_perc*100:.0f}")
    else:
        parts.append(f"filt-{filter_method[:3]}")

    var_method = config.get('variable_selection_method', 'none')
    parts.append(f"var-{var_method[:10]}")

    pca_method = config.get('pca_method', 'none')
    parts.append(f"pca-{pca_method.replace('pca_','')}")

    std_method = "stdOn" if config.get('standardization') else "stdOff"
    parts.append(std_method)

    clustering_method = config.get('clustering_method', 'shc')
    if clustering_method == 'shc':
        metric = config.get('shc_metric', 'euc')
        alpha_val = config.get('shc_alpha', 0)
        lambda_val = config.get('shc_lambda', 0)
        n_clusters = config.get('n_clusters', 0)
        parts.append(f"SHC_{metric[:3]}_a{str(alpha_val).replace('.','p')}_l{str(lambda_val).replace('.','p')}_k{n_clusters}")

    return "_".join(parts)


def save_metadata(output_dir, filename="metadata.txt", config=None, comments=None):
    """Saves run settings and comments to a metadata file."""
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Execution Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("--- Configuration ---\n")
        if config:
            for key, value in config.items():
                f.write(f"{key}: {value}\n")
        if comments:
            f.write("\n--- Comments ---\n")
            f.write(comments)
    print(f"Metadata saved to {filepath}.")


def read_run_pdti_ta_params(working_directory):
    """Reads parameters from run_PDTI_TA.bash."""
    bash_file_path = os.path.join(working_directory, "const", "run_PDTI_TA.bash")
    params = {}
    if not os.path.exists(bash_file_path):
        raise FileNotFoundError(f"run_PDTI_TA.bash not found at {bash_file_path}")

    print(f"Reading parameters from {bash_file_path}...")
    try:
        with open(bash_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            match_xx = re.search(r"^\s*XX\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_xx: params['n_interval'] = float(match_xx.group(1))
            else: raise ValueError("Could not find XX (n_interval) in run_PDTI_TA.bash.")
            match_yy = re.search(r"^\s*YY\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_yy: params['m_interval'] = float(match_yy.group(1))
            else: raise ValueError("Could not find YY (m_interval) in run_PDTI_TA.bash.")
            match_tr = re.search(r"^\s*TR\s*=\s*([0-9.]+)", content, re.MULTILINE)
            if match_tr: params['delta_tw'] = float(match_tr.group(1))
            else: raise ValueError("Could not find TR (delta_tw) in run_PDTI_TA.bash.")
        print(f"  Successfully read: n_interval={params['n_interval']}, m_interval={params['m_interval']}, delta_tw={params['delta_tw']}")
        return params
    except Exception as e:
        print(f"Error: An error occurred while reading or parsing {bash_file_path}: {e}")
        raise

# --- Data Processing Stage Functions ---
def stage_load_data(wd, snap_y_file="snap_y.dat", n_vector_file="n_vector.dat"):
    print(f"\n--- 1. Loading Data ---")
    print(f"Working directory: {wd}")
    ps_snap = os.path.join(wd, snap_y_file)
    if not os.path.exists(ps_snap): raise FileNotFoundError(f"Data file not found: {ps_snap}")
    df_snap_original = pd.read_table(ps_snap, sep=r'\s+', header=None)
    df_snap_original.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                                'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                'trendp', 'trendt', 'trendb', 'plungp',
                                'plungt', 'plungb', 'NDC']
    print(f"Loaded {snap_y_file}. Shape: {df_snap_original.shape}")
    ps_nvector = os.path.join(wd, n_vector_file)
    if not os.path.exists(ps_nvector): raise FileNotFoundError(f"Data file not found: {ps_nvector}")
    snap4_original = pd.read_table(ps_nvector, sep=r'\s+', header=None)
    snap4_original.columns = ['n', 'm', 'tw','t', 'x_coord', 'y_coord',
                                'n1_(1)','n1_(2)','n1_(3)',
                                'n2_(1)','n2_(2)','n2_(3)',
                                'sliprate']
    snap4_original.rename(columns={'n': 'x_grid', 'm': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    print(f"Loaded {n_vector_file}. Shape: {snap4_original.shape}")
    return df_snap_original, snap4_original

def stage_advanced_filtering(df_input_orig, n_interval, m_interval, delta_tw, target_extraction_percentage):
    """New stage for filtering major spatio-temporal points."""
    print("  Starting advanced filtering process...")
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
            print(f"    Warning: pivot_table creation failed for tw={current_tw}: {e_pivot}. Skipping Ds calculation for this time window.")
            continue
        if sliprate_grid.empty or sliprate_grid.shape[0] < 2 or sliprate_grid.shape[1] < 2:
            print(f"    Warning: Skipping Ds calculation for tw={current_tw} due to insufficient sliprate_grid size. Shape: {sliprate_grid.shape}")
            continue
        if n_interval == 0 or m_interval == 0:
            print(f"    Warning: Skipping Ds calculation for tw={current_tw} because n_interval or m_interval is 0.")
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
    print(f"  S_total threshold ({1.0 - target_extraction_percentage:.2f} quantile): {threshold_s_total:.4f}")
    df_input['is_major'] = df_input['S_total'] >= threshold_s_total

    df_input.rename(columns={'x': 'x_grid', 'y': 'y_grid', 'tw': 'tw_grid'}, inplace=True)
    filtered_df = df_input[df_input['is_major'] == True].copy()
    return df_input, filtered_df

def stage_filter_data(snap4_df, config):
    """Data Filtering"""
    method = config['filter_method']
    filter_value = config['filter_value']
    print(f"\n--- 2. Data Filtering ---")
    print(f"Method: {method}, Value: {filter_value}")
    
    if method == 'arbitrary_percentage':
        threshold_sliprate = snap4_df['sliprate'].quantile(filter_value)
        filtered_df = snap4_df[snap4_df['sliprate'] > threshold_sliprate].copy()
        print(f"  Filtered by sliprate > {threshold_sliprate:.4f} (quantile={filter_value}).")
    elif method == 'advanced_filtering':
        pdti_params = read_run_pdti_ta_params(config['working_directory'])
        _, filtered_df = stage_advanced_filtering(
            snap4_df.copy(),
            pdti_params['n_interval'],
            pdti_params['m_interval'],
            pdti_params['delta_tw'],
            filter_value
        )
    else:
        raise ValueError(f"Unknown filtering method '{method}'.")

    print(f"  Before filtering: {snap4_df.shape[0]} rows, After filtering: {filtered_df.shape[0]} rows")
    if filtered_df.empty:
        raise ValueError("Data became empty after filtering. Please check your filtering conditions.")

    return filtered_df


def stage_select_variables(df_snap_original, snap4_filtered, method, custom_columns_str=None):
    """Variable Selection"""
    print(f"\n--- 3. Variable Selection ---")

    print(f"Method: {method}")

    base_matrix_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    base_nv6_cols = ['n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)']
    common_cols_suffix = {
        "only": [], "nm": ['n', 'm'], "nmtw": ['n', 'm', 'tw'],
        "nmtw_sliprate": ['n', 'm', 'tw', 'sliprate']
    }
    data_for_clustering = pd.DataFrame()

    if method == 'custom_columns':
        if not custom_columns_str:
            raise ValueError("Custom columns selected but no columns were provided in config.")
        custom_columns = [col.strip() for col in custom_columns_str.split(',')]
        print(f"  Custom columns selected: {custom_columns}")

        snap4_temp_keys = snap4_filtered.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'})
        snap_y_cols_to_merge = [col for col in df_snap_original.columns if col not in snap4_temp_keys.columns or col in ['n','m','tw']]

        merged_df = pd.merge(snap4_temp_keys,
                             df_snap_original[snap_y_cols_to_merge],
                             on=['n', 'm', 'tw'], how='left')

        missing_cols = [col for col in custom_columns if col not in merged_df.columns]
        if missing_cols:
            raise ValueError(f"Some specified columns do not exist in the data: {missing_cols}")

        data_for_clustering = merged_df[custom_columns].copy()
        data_for_clustering.index = snap4_filtered.index

    elif method.startswith("matrix_components_") or method.startswith("nv6_components_"):
        if method.startswith("matrix_components_"):
            base_cols = base_matrix_cols
            prefix = "matrix_components_"
        else: # nv6_components
            base_cols = base_nv6_cols
            prefix = "nv6_components_"

        suffix_key = method.replace(prefix, "")
        if suffix_key not in common_cols_suffix:
            raise ValueError(f"Unknown suffix for {prefix}: {suffix_key}")

        additional_cols_keys = common_cols_suffix[suffix_key]
        
        if prefix == "nv6_components_":
            cols_to_use = base_cols[:]
            key_map = {'n': 'x_grid', 'm': 'y_grid', 'tw': 'tw_grid', 'sliprate': 'sliprate'}
            for key in additional_cols_keys:
                if key_map.get(key) in snap4_filtered.columns:
                    cols_to_use.append(key_map[key])
            data_for_clustering = snap4_filtered[cols_to_use].copy()
        else: # matrix components need a merge
            snap4_temp_keys = snap4_filtered.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'})
            cols_to_merge = base_cols + additional_cols_keys
            cols_to_merge = [c for c in cols_to_merge if c in df_snap_original.columns or c in ['n','m','tw']]
            
            merged_df = pd.merge(snap4_temp_keys[['n','m','tw']], df_snap_original[cols_to_merge], on=['n','m','tw'], how='left')
            final_cols = base_cols + additional_cols_keys
            data_for_clustering = merged_df[[c for c in final_cols if c in merged_df.columns]].copy()
            data_for_clustering.index = snap4_filtered.index
    else:
        raise ValueError(f"Unknown variable selection category: {method}")

    if data_for_clustering.isnull().values.any():
        print("Warning: Selected variables contain missing values. Dropping rows with NaNs.")
        data_for_clustering.dropna(inplace=True)
        if data_for_clustering.empty:
            raise ValueError("Data became empty after handling missing values in variable selection.")
    print(f"  Variables to be used: {data_for_clustering.columns.tolist()}")
    print(f"  Shape of data for clustering: {data_for_clustering.shape}")
    return data_for_clustering


def stage_standardize_data(df_for_clustering, apply_standardization):
    print(f"\n--- 4. Data Standardization ---")
    if not apply_standardization:
        print("  Method: Off")
        return df_for_clustering
    
    numeric_cols_for_std = df_for_clustering.select_dtypes(include=np.number).columns.tolist()
    if not numeric_cols_for_std:
        print("  Warning: No numerical columns to standardize. Skipping standardization.")
        return df_for_clustering
    
    print(f"  Method: On (Applying StandardScaler to columns {numeric_cols_for_std})")
    scaler = StandardScaler()
    df_for_clustering[numeric_cols_for_std] = scaler.fit_transform(df_for_clustering[numeric_cols_for_std])
    return df_for_clustering


def stage_pca(data_scaled, config, output_dir):
    """Principal Component Analysis (PCA)"""
    method = config['pca_method']
    print(f"\n--- 5. Principal Component Analysis (PCA) ---")
    print(f"Method: {method}")

    if method == 'pca_none':
        print("  PCA will not be performed.")
        return data_scaled

    data_for_pca_fit = data_scaled.select_dtypes(include=np.number)
    if data_for_pca_fit.empty:
        print("  Warning: No numerical data columns, PCA cannot be performed. Skipping PCA.")
        return data_scaled

    max_components = data_for_pca_fit.shape[1]
    n_components = 0

    if method in ['pca_80', 'pca_90']:
        pca_temp = PCA()
        pca_temp.fit(data_for_pca_fit)
        cumulative_variance_temp = np.cumsum(pca_temp.explained_variance_ratio_)
        current_pca_threshold = 0.80 if method == 'pca_80' else 0.90
        n_components = np.argmax(cumulative_variance_temp >= current_pca_threshold) + 1
        print(f"  Number of components for {current_pca_threshold*100:.0f}% cumulative variance: {n_components}")

    elif method == 'pca_n_components':
        n_components = config.get('pca_n_components_value')
        if n_components is None or n_components < 1:
            raise ValueError("Number of components is not specified or is invalid.")
        print(f"  User-specified number of components: {n_components}")
        if n_components > max_components:
            print(f"  Warning: Specified number of components ({n_components}) exceeds available features ({max_components}).")
            print(f"  Adjusting number of components to {max_components}.")
            n_components = max_components
    else:
        raise ValueError(f"Unknown PCA method: {method}")

    n_components = max(1, min(n_components, max_components))
    pca_model = PCA(n_components=n_components)
    data_after_pca_transformed = pca_model.fit_transform(data_for_pca_fit)

    data_after_pca = pd.DataFrame(data_after_pca_transformed, index=data_scaled.index,
                                  columns=[f'PC{i+1}' for i in range(n_components)])
    print(f"  Data shape after PCA: {data_after_pca.shape}")

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
        print(f"  PCA variance info saved to {pca_score_path}.")

    return data_after_pca


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
                print(f"Warning: Loop is too long in clset generation. ID {original_data_idx}, clsb {clsb_ids}")
                break
        clset_map[original_data_idx] = f"[{', '.join(map(str, clsb_ids))}]" if clsb_ids else "[0]"
    return clset_map

def stage_assign_and_save_results(
    df_snap_original_full,
    df_n_vector_clustered_data,
    linkage_result_matrix,
    final_labels,
    output_dir,
    run_config
):
    print(f"\n--- 8. Assigning Clusters and Saving Results ---")
    
    num_original_samples_in_clustering = df_n_vector_clustered_data.shape[0]

    unique_raw_labels = sorted(list(np.unique(final_labels)))
    label_mapping = {raw_label: i + 1 for i, raw_label in enumerate(unique_raw_labels)}
    final_mapped_labels = np.array([label_mapping[rl] for rl in final_labels])

    labels_df_temp = pd.DataFrame({'label': final_mapped_labels}, index=df_n_vector_clustered_data.index)
    
    df_for_output = df_snap_original_full.copy()
    
    if num_original_samples_in_clustering > 0 and linkage_result_matrix is not None and linkage_result_matrix.size > 0:
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
        df_for_output['clset'] = df_for_output.index.map(clset_map).fillna("[0]")
    else:
        df_for_output['clset'] = "[0]"

    df_for_output = df_for_output.merge(labels_df_temp, how='left', left_index=True, right_index=True)
    df_for_output['label'] = df_for_output['label'].fillna(0).astype(int)

    print(f"  Number of unique final cluster labels: {df_for_output['label'].nunique()}")
    print(f"  Number of members in each cluster (0 is unclassified):\n{df_for_output['label'].value_counts().sort_index()}")

    output_path_snap = os.path.join(output_dir, "clusteringSnap.dat")
    df_for_output.to_csv(output_path_snap, sep=' ', index=False, header=False)
    print(f"  Original format clusteringSnap.dat saved to {output_path_snap}.")
    output_path_snap2 = os.path.join(output_dir, "clusteringSnap2.dat")
    df_for_output.to_csv(output_path_snap2, sep='\t', index=False, header=False)
    print(f"  Original format clusteringSnap2.dat saved to {output_path_snap2}.")
    output_path_info = os.path.join(output_dir, "clusteringinfo.dat")
    df_for_output['label'].to_csv(output_path_info, sep='\t', index=False, header=False)
    print(f"  Original format clusteringinfo.dat saved to {output_path_info}.")

    if num_original_samples_in_clustering > 0 and linkage_result_matrix is not None and linkage_result_matrix.size > 0:
        fig_final_dendro, ax_final_dendro = plt.subplots(figsize=(18, 12))
        
        dendrogram(linkage_result_matrix, ax=ax_final_dendro, truncate_mode=None)
        
        title_str = f"Final Dendrogram (SHC Re-cluster)"
        ax_final_dendro.set_title(title_str)
        ax_final_dendro.set_xlabel("Sample index or (Cluster size)")
        ax_final_dendro.set_ylabel("Distance")
            
        info_text = "Re-Clustering Run Configuration:\n"
        info_text += f"  Source Dir: {run_config.get('source_dir')}\n"
        info_text += f"  New N Clusters: {run_config.get('n_clusters')}\n"

        fig_final_dendro.text(0.01, 0.98, info_text, transform=fig_final_dendro.transFigure,
                              fontsize=8, verticalalignment='top',
                              bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        dendrogram_final_path = os.path.join(output_dir, "dendrogram_recluster.png")
        plt.savefig(dendrogram_final_path)
        print(f"  Re-clustering dendrogram saved to {dendrogram_final_path}.")
        plt.close(fig_final_dendro)

    files_to_copy = ["snap_yr.dat", "snap2_y.dat", "fort.40",
                     "faultline.dat", "mrf.dat", "snap2.dat","rigid_amp.info"]
    print(f"\n  Starting to copy specified files...")
    for filename in files_to_copy:
        source_path = os.path.join(run_config['working_directory'], filename)
        destination_path = os.path.join(output_dir, filename)
        if os.path.exists(source_path):
            try:
                shutil.copy2(source_path, destination_path)
            except Exception as e:
                print(f"    Warning: An error occurred while copying {filename}: {e}")
        else:
            print(f"    Warning: Source file not found, skipping: {source_path}")


def read_config(filename='re_cluster_config.txt'):
    """Reads the configuration from a text file."""
    print(f"Reading configuration from '{filename}'...")
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Configuration file '{filename}' not found in the current directory.")
    
    config = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if '=' in line:
                separator = '='
            elif ':' in line:
                separator = ':'
            else:
                continue

            key, value = line.split(separator, 1)
            key = key.strip()
            value = value.strip()
            
            if value.lower() == 'true':
                config[key] = True
            elif value.lower() == 'false':
                config[key] = False
            else:
                try:
                    config[key] = int(value)
                except ValueError:
                    try:
                        config[key] = float(value)
                    except ValueError:
                        config[key] = value
    return config


def main():
    """Main execution function for re-clustering."""
    try:
        config = read_config('re_cluster_config.txt')
        
        source_dir = config.get('source_dir')
        new_n_clusters = config.get('n_clusters')
        
        if not source_dir or not new_n_clusters:
            raise ValueError("'source_dir' and 'n_clusters' must be specified in re_cluster_config.txt")
        
        linkage_matrix_path = os.path.join(source_dir, "shc_linkage_matrix.dat")
        metadata_path = os.path.join(source_dir, "metadata.txt")
        
        if not os.path.exists(linkage_matrix_path):
            raise FileNotFoundError(f"Linkage matrix not found at: {linkage_matrix_path}")
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"Metadata file not found at: {metadata_path}. Cannot reproduce pre-processing steps.")
            
        print(f"Loading original run configuration from {metadata_path}")
        original_config = read_config(metadata_path)
        
        new_config = original_config.copy()
        new_config['n_clusters'] = new_n_clusters
        
        # --- Reproduce pre-processing pipeline ---
        df_snap_orig, snap4_orig = stage_load_data(new_config['working_directory'])
        snap4_filtered = stage_filter_data(snap4_orig.copy(), new_config)
        data_for_clustering = stage_select_variables(df_snap_orig.copy(), snap4_filtered.copy(), new_config['variable_selection_method'], new_config.get('custom_columns'))
        data_std = stage_standardize_data(data_for_clustering.copy(), new_config['standardization'])
        data_pca = stage_pca(data_std.copy(), new_config, None) # Don't save PCA results again
        
        # --- Load linkage matrix and re-cluster ---
        print(f"Loading pre-computed linkage matrix from: {linkage_matrix_path}")
        linkage_mat = np.loadtxt(linkage_matrix_path)
        
        print(f"  Splitting into {new_n_clusters} clusters.")
        final_labels = fcluster(linkage_mat, new_n_clusters, criterion='maxclust')

        # --- Generate new output directory and save results ---
        output_dir_name = generate_output_dirname_from_config(new_config)
        output_dir = os.path.join(new_config['working_directory'], output_dir_name)
        ensure_dir(output_dir)
        print(f"Output for re-clustering will be saved to: {output_dir}")

        stage_assign_and_save_results(
            df_snap_original_full=df_snap_orig.copy(),
            df_n_vector_clustered_data=data_for_clustering.copy(),
            linkage_result_matrix=linkage_mat,
            final_labels=final_labels,
            output_dir=output_dir,
            run_config=new_config
        )
        
        save_metadata(output_dir, config=new_config, comments="This is a re-clustered result based on a previous SHC run.")
        print(f"=== Re-clustering completed successfully. Results saved in: {output_dir} ===")

    except Exception as e:
        import traceback
        print(f"\n--- RE-CLUSTERING FAILED ---")
        print(f"Error: {e}")
        print(traceback.format_exc())
        exit(1)


if __name__ == "__main__":
    main()

