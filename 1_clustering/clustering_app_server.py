#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# clustering_app_server.py
# サーバー環境での実行に最適化されたクラスタリングスクリプト。
# run_clustering.txt から設定を読み込み、SHCとウォード法を実行可能。
# SHCの計算結果を保存し、高速に再クラスタリングする機能を搭載。

print("Clustering App for Server (Complete) started.")

import pandas as pd
import numpy as np
import matplotlib
# サーバー環境でGUIバックエンドがなくても動作するように設定
matplotlib.use('Agg')
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

# 並列処理ライブラリをインポート
try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    print("Warning: joblib library not found. Parallel processing will be disabled.")
    print("You can install it by running: pip install joblib")


# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
# --- SHC (Stable Hierarchical Clustering) Module ---
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

class _SHCNode:
    """Internal class to represent the hierarchical tree structure for the SHC algorithm."""
    def __init__(self, node_id, left=None, right=None, count=0, distance=0):
        self.id = node_id
        self.left = left
        self.right = right
        self.count = count
        self.distance = distance

class StableHierarchicalClustering:
    """
    Implementation of the SHC algorithm with parallel processing and memory optimization.
    """
    def __init__(self, alpha: float, lambda_: float, max_depth: int, metric: str = 'euclidean', n_jobs: int = 1):
        if alpha <= 0:
            raise ValueError("alpha must be positive.")
        if lambda_ < 0:
            raise ValueError("lambda_ must be non-negative.")
        if max_depth <= 0:
            raise ValueError("max_depth must be positive.")
        if metric not in ['euclidean', 'cosine']:
            raise ValueError("metric must be either 'euclidean' or 'cosine'.")
            
        self.alpha = alpha
        self.lambda_ = lambda_
        self.max_depth = max_depth
        self.metric = metric
        self.n_jobs = n_jobs if JOBLIB_AVAILABLE else 1
        if self.n_jobs != 1 and JOBLIB_AVAILABLE:
            core_count_msg = 'all available' if self.n_jobs == -1 else str(self.n_jobs)
            print(f"  Parallel processing enabled with {core_count_msg} cores.")
        
        self._W = None
        self._n_samples = 0
        self._linkage_matrix_list = []

    def fit(self, X: np.ndarray):
        self._n_samples = X.shape[0]
        if self._n_samples < 2:
            return np.array([])
        
        print(f"  Building similarity graph using '{self.metric}' metric...")
        self._build_similarity_graph(X)
        
        initial_vertices = list(range(self._n_samples))
        
        print("  Building SHC tree (this may take a while)...")
        with tqdm(total=self._n_samples - 1, desc="SHC Total Splits") as pbar:
            root_node = self._recursive_build_tree(initial_vertices, 0, pbar)
        
        if root_node is None or (root_node.left is None and root_node.right is None):
            return np.array([])

        self._linkage_matrix_list = []
        self._build_linkage_recursively(root_node)
        
        return np.array(self._linkage_matrix_list)

    def _build_linkage_recursively(self, node):
        """
        Builds the linkage matrix by performing a post-order traversal of the SHC tree.
        This ensures the matrix is compliant with SciPy's requirements.
        """
        if node.left is None and node.right is None:
            return float(node.id)

        left_id = self._build_linkage_recursively(node.left)
        right_id = self._build_linkage_recursively(node.right)
        
        distance = float(self.max_depth - node.distance)
        count = float(node.count)
        
        self._linkage_matrix_list.append([left_id, right_id, distance, count])
        
        new_cluster_id = float(self._n_samples + len(self._linkage_matrix_list) - 1)
        
        return new_cluster_id

    def _build_similarity_graph(self, X: np.ndarray):
        if self.metric == 'euclidean':
            d_sq = squareform(pdist(X, 'sqeuclidean'))
        elif self.metric == 'cosine':
            d = squareform(pdist(X, 'cosine'))
            d_sq = d ** 2
        else:
            raise ValueError(f"Unsupported metric: {self.metric}")
        self._W = np.exp(-self.alpha * d_sq)

    def _recursive_build_tree(self, V_sub: list, current_depth: int, pbar: tqdm):
        if current_depth >= self.max_depth or len(V_sub) <= 1:
            if not V_sub: return None
            if len(V_sub) == 1:
                return _SHCNode(node_id=V_sub[0], count=1, distance=current_depth)
            else:
                return self._form_tree_from_leaves(V_sub, current_depth)

        S1_set = self._stable_sparsest_cut(V_sub)
        S1 = sorted(list(S1_set))
        S2 = sorted(list(set(V_sub) - S1_set))

        if not S1 or not S2:
            return self._form_tree_from_leaves(V_sub, current_depth)

        left_child = self._recursive_build_tree(S1, current_depth + 1, pbar)
        right_child = self._recursive_build_tree(S2, current_depth + 1, pbar)
        
        count = (left_child.count if left_child else 0) + (right_child.count if right_child else 0)
        
        pbar.update(1)
        
        return _SHCNode(node_id=-1, left=left_child, right=right_child, count=count, distance=current_depth)

    def _form_tree_from_leaves(self, V_sub: list, current_depth: int):
        if not V_sub: return None
        if len(V_sub) == 1: return _SHCNode(node_id=V_sub[0], count=1, distance=current_depth)
        
        left = _SHCNode(node_id=V_sub[0], count=1, distance=current_depth)
        right = _SHCNode(node_id=V_sub[1], count=1, distance=current_depth)
        parent = _SHCNode(node_id=-1, left=left, right=right, count=2, distance=current_depth)
        
        for i in range(2, len(V_sub)):
            new_node = _SHCNode(node_id=V_sub[i], count=1, distance=current_depth)
            parent = _SHCNode(node_id=-1, left=parent, right=new_node, count=parent.count + 1, distance=current_depth)
        return parent

    def _calculate_phi_only(self, pair, V_sub_set):
        """Helper function for parallel execution. Calculates ONLY the phi score for a pair."""
        i, j = pair
        S_ij = {k for k in V_sub_set if self._W[i, k] > self._W[j, k]}
        
        size_S = len(S_ij)
        size_S_complement = len(V_sub_set) - size_S
        
        if size_S == 0 or size_S_complement == 0:
            return np.inf
        
        S_complement_set = V_sub_set - S_ij
        cut_weight = self._W[np.ix_(list(S_ij), list(S_complement_set))].sum()
        return cut_weight / (size_S * size_S_complement)

    def _stable_sparsest_cut(self, V_sub: list):
        if len(V_sub) <= 1:
            return set(V_sub)

        V_sub_set = set(V_sub)
        pbar_desc = f"SSC on {len(V_sub)} pts"
        
        pairs_iterator = itertools.combinations(V_sub, 2)
        num_combinations = (len(V_sub) * (len(V_sub) - 1)) // 2

        quality_scores = Parallel(n_jobs=self.n_jobs)(
            delayed(self._calculate_phi_only)(pair, V_sub_set) 
            for pair in tqdm(pairs_iterator, total=num_combinations, desc=pbar_desc, leave=False)
        )
        quality_scores = np.array(quality_scores)
        
        if np.any(np.isfinite(quality_scores)):
            phi_min = np.min(quality_scores[np.isfinite(quality_scores)])
            shifted_scores = quality_scores - phi_min
            weights = np.exp(-self.lambda_ * shifted_scores)
            weights[np.isinf(quality_scores)] = 0
            total_weight = np.sum(weights)
        else:
            total_weight = 0

        if total_weight > 1e-9:
            probabilities = weights / total_weight
        else:
            valid_indices = np.where(np.isfinite(quality_scores))[0]
            if len(valid_indices) > 0:
                probabilities = np.zeros_like(quality_scores)
                probabilities[valid_indices] = 1.0 / len(valid_indices)
            else:
                probabilities = np.ones_like(quality_scores) / len(quality_scores)
                
        if np.sum(probabilities) > 0:
            probabilities /= np.sum(probabilities)
        else:
            probabilities = np.ones_like(quality_scores) / len(quality_scores)

        chosen_index = np.random.choice(len(quality_scores), p=probabilities)
        
        winning_pair_iterator = itertools.islice(itertools.combinations(V_sub, 2), chosen_index, chosen_index + 1)
        i, j = next(winning_pair_iterator)
        
        chosen_split = {k for k in V_sub_set if self._W[i, k] > self._W[j, k]}
        
        return chosen_split

# --- End of SHC Module ---


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

    clustering_method = config.get('clustering_method', 'ward')
    if clustering_method == 'ward':
        dist_metric = config.get('distance_metric', 'none')
        parts.append(f"ward-{dist_metric[:3]}")
        thresh_method = config.get('threshold_method', 'none')
        parts.append(f"thr-{thresh_method.replace('threshold_','')[:3]}")
    elif clustering_method == 'shc':
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


# --- 6. Clustering ---
def stage_linkage_clustering(data_for_linkage, distance_metric):
    """Performs traditional hierarchical clustering (Ward's method)."""
    print(f"\n--- 6. Hierarchical Clustering (Ward's Method) ---")
    print(f"  Distance metric: {distance_metric}")
    values_for_pdist = data_for_linkage
    if isinstance(data_for_linkage, pd.DataFrame):
        numeric_cols = data_for_linkage.select_dtypes(include=np.number).columns
        if numeric_cols.empty: raise ValueError("No numerical columns in data for clustering.")
        if len(numeric_cols) < data_for_linkage.shape[1]:
            print(f"  Warning: Linkage calculation will only use numerical columns: {list(numeric_cols)}")
        values_for_pdist = data_for_linkage[numeric_cols].values
    elif not isinstance(data_for_linkage, np.ndarray):
        raise TypeError("Data for clustering must be a DataFrame or ndarray.")
    if values_for_pdist.shape[0] < 2:
        raise ValueError(f"Too few data points ({values_for_pdist.shape[0]}) to perform hierarchical clustering.")
    if values_for_pdist.ndim == 1: values_for_pdist = values_for_pdist.reshape(-1, 1)
    try:
        distance_matrix_condensed = pdist(values_for_pdist, metric=distance_metric)
        linkage_matrix = linkage(distance_matrix_condensed, method='ward')
    except Exception as e:
        print(f"Error: An issue occurred during linkage calculation: {e}")
        print(f"Data shape: {values_for_pdist.shape}, Data type: {type(values_for_pdist)}")
        if isinstance(values_for_pdist, np.ndarray):
            print(f"Data contains NaN: {np.isnan(values_for_pdist).any()}, Data contains inf: {np.isinf(values_for_pdist).any()}")
        raise e
    print(f"  Shape of linkage matrix: {linkage_matrix.shape}")
    return linkage_matrix

def stage_determine_cluster_threshold(linkage_matrix, data_for_silhouette, config, output_dir):
    """Determines the cluster distance threshold for Ward's method."""
    method = config['threshold_method']
    print(f"\n--- 7. Cluster Distance Threshold (Ward's Method) ---")
    print(f"  Method: {method}")

    final_distance_threshold = None

    if method == 'threshold_arbitrary':
        target_k = config.get('target_cluster_count_for_arbitrary')
        if target_k is None:
            raise ValueError("For 'threshold_arbitrary', 'target_cluster_count_for_arbitrary' must be provided.")
        
        print(f"  Calculating distance threshold for target cluster count {target_k}.")
        if linkage_matrix.shape[0] + 1 < target_k or target_k < 2:
            raise ValueError(f"Target cluster count {target_k} is inappropriate. Number of samples: {linkage_matrix.shape[0]+1}")

        sorted_distances = np.sort(linkage_matrix[:, 2])
        if target_k > len(sorted_distances):
            final_distance_threshold = sorted_distances[-1] if len(sorted_distances) > 0 else 0.1
        else:
            final_distance_threshold = sorted_distances[-target_k]

        if final_distance_threshold == 0:
            non_zero_distances = sorted_distances[sorted_distances > 0]
            if len(non_zero_distances) > 0:
                if -target_k + 1 < 0 and sorted_distances[-target_k + 1] > 0:
                    final_distance_threshold = sorted_distances[-target_k + 1] / 2.0
                    if final_distance_threshold == 0: final_distance_threshold = 1e-6
                else:
                    final_distance_threshold = non_zero_distances[0] * 0.5 if non_zero_distances[0] * 0.5 > 0 else 1e-6
            else:
                final_distance_threshold = 1e-6
        print(f"  Calculated distance threshold (for target K={target_k}): {final_distance_threshold:.6f}")

    elif method == 'threshold_silhouette':
        print("  Determining optimal number of clusters based on silhouette score...")
        actual_data_for_sil_values = data_for_silhouette
        if isinstance(data_for_silhouette, pd.DataFrame):
            numeric_cols_sil = data_for_silhouette.select_dtypes(include=np.number).columns
            if numeric_cols_sil.empty: raise ValueError("No numerical columns in data for silhouette score calculation.")
            actual_data_for_sil_values = data_for_silhouette[numeric_cols_sil].values
        if actual_data_for_sil_values.ndim == 1:
            actual_data_for_sil_values = actual_data_for_sil_values.reshape(-1,1)
        if actual_data_for_sil_values.shape[0] < 2:
            raise ValueError("Not enough data points (< 2) for silhouette analysis.")

        upper_k_limit = min(21, actual_data_for_sil_values.shape[0]) # Test more K values on server
        if upper_k_limit <=2 :
            raise ValueError(f"Too few data points ({actual_data_for_sil_values.shape[0]}) to test multiple K values in silhouette analysis.")
        cluster_range_list = list(range(2, upper_k_limit))
        sil_scores = []
        for k_val in tqdm(cluster_range_list, desc="Silhouette Analysis", leave=False):
            try:
                labels = fcluster(linkage_matrix, k_val, criterion='maxclust')
                if len(np.unique(labels)) < 2 :
                    sil_scores.append(-1.0)
                    continue
                score = silhouette_score(actual_data_for_sil_values, labels, metric=config['distance_metric'])
                sil_scores.append(score)
            except ValueError as e:
                print(f"  Warning: Silhouette score calculation error for K={k_val}: {e}")
                sil_scores.append(-2.0)

        if not sil_scores or max(sil_scores) < -0.9 :
            raise ValueError("Could not determine optimal number of clusters from silhouette scores. No valid scores were obtained.")

        if output_dir:
            silhouette_scores_df = pd.DataFrame({
                'NumberOfClusters': cluster_range_list,
                'SilhouetteScore': sil_scores
            })
            slietscore_path = os.path.join(output_dir, "slietscore.csv")
            silhouette_scores_df.to_csv(slietscore_path, index=False)
            print(f"  Silhouette score data saved to {slietscore_path}.")

        fig_sil, ax_sil = plt.subplots(figsize=(8, 4))
        ax_sil.plot(cluster_range_list, sil_scores, marker='o')
        ax_sil.set_xlabel("Number of clusters")
        ax_sil.set_ylabel("Silhouette Score")
        ax_sil.set_title("Silhouette Analysis")
        ax_sil.grid(True)
        plt.tight_layout()
        silhouette_plot_path = os.path.join(output_dir, "silhouette_scores.png")
        plt.savefig(silhouette_plot_path)
        print(f"  Silhouette score plot saved to {silhouette_plot_path}.")
        plt.close(fig_sil)

        optimal_k_idx = np.argmax(sil_scores)
        optimal_k = cluster_range_list[optimal_k_idx]
        print(f"  Optimal number of clusters (max silhouette score): {optimal_k} (Score: {sil_scores[optimal_k_idx]:.4f})")

        sorted_distances = np.sort(linkage_matrix[:, 2])
        if optimal_k > len(sorted_distances):
            final_distance_threshold = sorted_distances[-1] if len(sorted_distances) > 0 else 0.1
        else:
            final_distance_threshold = sorted_distances[-optimal_k]
        print(f"  Recommended distance threshold (based on silhouette score): {final_distance_threshold:.6f}")

    else:
        raise ValueError(f"Unknown threshold determination method: {method}")

    return final_distance_threshold


def stage_shc_clustering_and_labeling(data_for_shc, config, output_dir):
    """Performs SHC clustering and labels data based on pre-defined cluster count."""
    print("\n--- 6. SHC (Stable Clustering) ---")
    
    alpha = config['shc_alpha']
    lambda_ = config['shc_lambda']
    metric = config['shc_metric']
    n_jobs = config.get('n_jobs', 1)
    n_clusters = config.get('n_clusters')
    max_depth = data_for_shc.shape[0]
    
    if n_clusters is None:
        raise ValueError("For SHC, 'n_clusters' must be provided in the config file.")
    
    print(f"  Parameters: alpha={alpha}, lambda={lambda_}, metric={metric}, n_jobs={n_jobs}, n_clusters={n_clusters}")
    
    shc = StableHierarchicalClustering(alpha=alpha, lambda_=lambda_, max_depth=max_depth, metric=metric, n_jobs=n_jobs)
    
    if isinstance(data_for_shc, pd.DataFrame):
        data_for_shc_np = data_for_shc.values
    else:
        data_for_shc_np = data_for_shc

    linkage_matrix = shc.fit(data_for_shc_np)
    
    if linkage_matrix.size == 0:
        raise ValueError("SHC clustering failed, resulting in an empty linkage matrix.")
        
    print(f"  Generated SHC linkage matrix. Shape: {linkage_matrix.shape}")

    # Save the linkage matrix for future re-clustering
    linkage_matrix_path = os.path.join(output_dir, "shc_linkage_matrix.dat")
    np.savetxt(linkage_matrix_path, linkage_matrix, fmt='%.8f')
    print(f"  SHC linkage matrix saved to {linkage_matrix_path} for re-clustering.")

    # Plot and save dendrogram without showing it
    fig_dendro, ax_dendro = plt.subplots(figsize=(16, 10))
    dendrogram(linkage_matrix, ax=ax_dendro)
    ax_dendro.set_title(f'SHC Dendrogram ({metric})')
    ax_dendro.set_xlabel("Sample index")
    ax_dendro.set_ylabel("Distance (Scipy format)")
    dendrogram_shc_path = os.path.join(output_dir, f"dendrogram_shc_full_{metric}.png")
    plt.savefig(dendrogram_shc_path)
    print(f"  SHC dendrogram saved to {dendrogram_shc_path}.")
    plt.close(fig_dendro)
            
    print(f"  Splitting into {n_clusters} clusters.")
    labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    if isinstance(data_for_shc, pd.DataFrame):
        final_labels_series = pd.Series(labels, index=data_for_shc.index)
    else:
        final_labels_series = pd.Series(labels)

    return linkage_matrix, final_labels_series


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
    
    if num_original_samples_in_clustering == 0:
        print("  Warning: No data to cluster, skipping cluster assignment and result saving.")
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
        
        df_for_output['clset'] = df_for_output.index.map(clset_map).fillna("[0]")
        
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

    if num_original_samples_in_clustering > 0:
        fig_final_dendro, ax_final_dendro = plt.subplots(figsize=(18, 12))
        
        color_threshold = -1
        if run_config.get('clustering_method') == 'ward' and 'final_cluster_threshold' in run_config:
            color_threshold = run_config['final_cluster_threshold']

        dendrogram(
            linkage_result_matrix, ax=ax_final_dendro,
            color_threshold=color_threshold, truncate_mode=None,
        )
        
        title_str = f"Final Dendrogram ({run_config.get('clustering_method')})"
        ax_final_dendro.set_title(title_str)
        ax_final_dendro.set_xlabel("Sample index or (Cluster size)")
        ax_final_dendro.set_ylabel("Distance")
        if color_threshold > 0:
            ax_final_dendro.axhline(y=color_threshold, color='r', linestyle='--')
            
        info_text = "Run Configuration:\n"
        info_text += f"  Mode: {run_config.get('mode', 'N/A')}\n"
        if run_config.get('mode') == 're-cluster':
            info_text += f"  Source Dir: {run_config.get('source_dir')}\n"
        info_text += f"  Filter: {run_config.get('filter_method')} (Value: {run_config.get('filter_value')})\n"
        info_text += f"  Variables: {run_config.get('variable_selection_method')}\n"
        info_text += f"  PCA: {run_config.get('pca_method')}\n"
        info_text += f"  Standardize: {run_config.get('standardization')}\n"
        info_text += f"  Clustering: {run_config.get('clustering_method')}\n"
        if run_config.get('clustering_method') == 'ward':
            info_text += f"    Distance: {run_config.get('distance_metric')}\n"
            info_text += f"    Threshold Method: {run_config.get('threshold_method')}\n"
        elif run_config.get('clustering_method') == 'shc':
            info_text += f"    Metric: {run_config.get('shc_metric')}\n"
            info_text += f"    alpha: {run_config.get('shc_alpha')}\n"
            info_text += f"    lambda: {run_config.get('shc_lambda')}\n"
            info_text += f"    n_clusters: {run_config.get('n_clusters')}\n"
            info_text += f"    n_jobs: {run_config.get('n_jobs', 1)}\n"

        fig_final_dendro.text(0.01, 0.98, info_text, transform=fig_final_dendro.transFigure,
                              fontsize=8, verticalalignment='top',
                              bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        dendrogram_final_path = os.path.join(output_dir, "dendrogram_final.png")
        plt.savefig(dendrogram_final_path)
        print(f"  Final dendrogram (with info) saved to {dendrogram_final_path}.")
        plt.close(fig_final_dendro)
    else:
        print("  No final dendrogram was saved because there was no data to cluster.")

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


# --- Main Processing Function ---
def run_single_pipeline(config, load_linkage_path=None):
    """Executes a single clustering pipeline with a given configuration."""
    print(f"\n=== Starting Pipeline Run ===")
    
    df_snap_orig, snap4_orig = stage_load_data(config['working_directory'])
    current_run_output_dir = ""

    try:
        snap4_filtered = stage_filter_data(snap4_orig.copy(), config)
        data_for_clustering = stage_select_variables(df_snap_orig.copy(), snap4_filtered.copy(), config['variable_selection_method'], config.get('custom_columns'))
        data_std = stage_standardize_data(data_for_clustering.copy(), config['standardization'])
        
        final_run_specific_dirname = generate_output_dirname_from_config(config)
        current_run_output_dir = os.path.join(config['working_directory'], final_run_specific_dirname)
        ensure_dir(current_run_output_dir)
        print(f"Output will be saved to: {current_run_output_dir}")

        data_pca = stage_pca(data_std.copy(), config, current_run_output_dir)

        linkage_mat, final_labels = None, None

        if config['clustering_method'] == 'ward':
            linkage_mat = stage_linkage_clustering(data_pca, config['distance_metric'])
            dist_thresh = stage_determine_cluster_threshold(linkage_mat, data_pca, config, current_run_output_dir)
            config['final_cluster_threshold'] = dist_thresh
            final_labels = fcluster(linkage_mat, dist_thresh, criterion='distance')

        elif config['clustering_method'] == 'shc':
            if load_linkage_path:
                print(f"Loading pre-computed linkage matrix from: {load_linkage_path}")
                linkage_mat = np.loadtxt(load_linkage_path)
                n_clusters = config['n_clusters']
                print(f"  Splitting into {n_clusters} clusters.")
                labels = fcluster(linkage_mat, n_clusters, criterion='maxclust')
                final_labels = pd.Series(labels, index=data_pca.index).values
            else:
                linkage_mat, final_labels_series = stage_shc_clustering_and_labeling(data_pca, config, current_run_output_dir)
                final_labels = final_labels_series.values
            
        if final_labels is not None:
            stage_assign_and_save_results(
                df_snap_original_full=df_snap_orig.copy(),
                df_n_vector_clustered_data=data_for_clustering.copy(),
                linkage_result_matrix=linkage_mat,
                final_labels=final_labels,
                output_dir=current_run_output_dir,
                run_config=config
            )
            save_metadata(current_run_output_dir, config=config)
            run_identifier = os.path.basename(current_run_output_dir)
            print(f"=== Pipeline completed successfully: {run_identifier} ===")
        else:
            raise ValueError("Cluster labels were not generated.")

    except Exception as e:
        import traceback
        print(f"\n--- PIPELINE FAILED ---")
        print(f"Error: {e}")
        print(traceback.format_exc())
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: {e}\n{traceback.format_exc()}")
        exit(1)


def recluster_from_saved_linkage(config):
    """Re-clusters data using a saved SHC linkage matrix and a new number of clusters."""
    print("\n=== Starting Re-Clustering from Saved Linkage Matrix ===")
    
    source_dir = config['source_dir']
    new_n_clusters = config['n_clusters']
    
    linkage_matrix_path = os.path.join(source_dir, "shc_linkage_matrix.dat")
    metadata_path = os.path.join(source_dir, "metadata.txt")
    
    if not os.path.exists(linkage_matrix_path):
        raise FileNotFoundError(f"Linkage matrix not found at: {linkage_matrix_path}")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found at: {metadata_path}. Cannot reproduce pre-processing steps.")
        
    print(f"Loading original run configuration from {metadata_path}")
    original_config = read_config(metadata_path)
    
    # Update config with new cluster count for the re-run
    new_config = original_config.copy()
    new_config['n_clusters'] = new_n_clusters
    
    run_single_pipeline(new_config, load_linkage_path=linkage_matrix_path)


def read_config(filename='run_clustering.txt'):
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
            
            key, value = line.split('=', 1)
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
                        config[key] = value # Treat as string
    return config


def main():
    """Main execution function driven by a configuration file."""
    config = read_config()
    
    mode = config.get('mode')
    if not mode:
        raise ValueError("'mode' not specified in run_clustering.txt")
        
    print(f"Execution mode: {mode}")

    if mode in ['run-shc', 'run-ward']:
        if mode == 'run-ward':
            config['clustering_method'] = 'ward'
        else:
            config['clustering_method'] = 'shc'
        run_single_pipeline(config)
        
    elif mode == 're-cluster':
        recluster_from_saved_linkage(config)
        
    else:
        raise ValueError(f"Unknown mode '{mode}' in configuration file.")


if __name__ == "__main__":
    main()

