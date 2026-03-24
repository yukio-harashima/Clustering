#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# clustering_app_lite.py
# 従来のクラスタリング手法に加え、SHC（安定クラスタリング）機能を統合。
# SHCの計算ボトルネックを並列化処理によって高速化し、メモリ効率を改善。

print("Clustering App Lite (Memory Optimized) script started.")

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

# ★ 並列処理ライブラリをインポート
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

        internal_nodes = []
        def find_internal_nodes(node):
            if node and node.left and node.right:
                internal_nodes.append(node)
                find_internal_nodes(node.left)
                find_internal_nodes(node.right)
        find_internal_nodes(root_node)

        internal_nodes.sort(key=lambda n: n.distance)

        node_to_new_id = {}
        for i, node in enumerate(internal_nodes):
            node_to_new_id[id(node)] = self._n_samples + i

        linkage_matrix = []
        for node in internal_nodes:
            if node.left.left is None and node.left.right is None:
                left_id = float(node.left.id)
            else:
                left_id = float(node_to_new_id[id(node.left)])

            if node.right.left is None and node.right.right is None:
                right_id = float(node.right.id)
            else:
                right_id = float(node_to_new_id[id(node.right)])
            
            linkage_matrix.append([left_id, right_id, node.distance, float(node.count)])

        return np.array(linkage_matrix)

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
                return _SHCNode(node_id=V_sub[0], count=1, distance=0)
            else:
                return self._form_tree_from_leaves(V_sub, current_depth)

        S1_set = self._stable_sparsest_cut(V_sub)
        S1 = sorted(list(S1_set))
        S2 = sorted(list(set(V_sub) - S1_set))

        if not S1 or not S2:
            return self._form_tree_from_leaves(V_sub, current_depth)

        left_child = self._recursive_build_tree(S1, current_depth + 1, pbar)
        right_child = self._recursive_build_tree(S2, current_depth + 1, pbar)
        
        distance = float(current_depth + 1)
        count = (left_child.count if left_child else 0) + (right_child.count if right_child else 0)
        
        pbar.update(1)
        
        return _SHCNode(node_id=-1, left=left_child, right=right_child, count=count, distance=distance)

    def _form_tree_from_leaves(self, V_sub: list, current_depth: int):
        if not V_sub: return None
        if len(V_sub) == 1: return _SHCNode(node_id=V_sub[0], count=1, distance=0)
        
        distance = float(current_depth + 0.5)
        left = _SHCNode(node_id=V_sub[0], count=1, distance=0)
        right = _SHCNode(node_id=V_sub[1], count=1, distance=0)
        parent = _SHCNode(node_id=-1, left=left, right=right, count=2, distance=distance)
        
        for i in range(2, len(V_sub)):
            new_node = _SHCNode(node_id=V_sub[i], count=1, distance=0)
            parent = _SHCNode(node_id=-1, left=parent, right=new_node, count=parent.count + 1, distance=distance)
        return parent
    
    # ★ 修正: メモリ効率改善のため、品質スコア(phi)のみを計算するヘルパー関数
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
        
        # ★ 修正: 巨大なリストを作らずイテレータを直接使う
        pairs_iterator = itertools.combinations(V_sub, 2)
        num_combinations = (len(V_sub) * (len(V_sub) - 1)) // 2

        # 1. 並列処理で品質スコア(phi)のみを計算・収集 (メモリ効率が良い)
        quality_scores = Parallel(n_jobs=self.n_jobs)(
            delayed(self._calculate_phi_only)(pair, V_sub_set) 
            for pair in tqdm(pairs_iterator, total=num_combinations, desc=pbar_desc, leave=False)
        )
        quality_scores = np.array(quality_scores)
        
        # 2. 収集したスコアから確率的に1つを選択
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
        
        # 3. 選択された1つのペアに対応する分割(S_ij)だけを再計算
        #    itertools.isliceで効率的に目的のペアを取得
        winning_pair_iterator = itertools.islice(itertools.combinations(V_sub, 2), chosen_index, chosen_index + 1)
        i, j = next(winning_pair_iterator)
        
        chosen_split = {k for k in V_sub_set if self._W[i, k] > self._W[j, k]}
        
        return chosen_split

# --- End of SHC Module ---


# --- Helper Functions ---
def ensure_dir(directory_path):
    """Creates a directory if it does not exist."""
    os.makedirs(directory_path, exist_ok=True)

def handle_overwrite_choice(target_path):
    """Handles user choice when a file or directory already exists."""
    if os.path.exists(target_path):
        while True:
            choice = input(f"'{target_path}' already exists. Overwrite, rename, or cancel? (o/r/c): ").lower()
            if choice == 'o':
                print(f"Overwriting '{target_path}'.")
                return 'overwrite'
            elif choice == 'r':
                print(f"Renaming and saving for '{target_path}'.")
                return 'rename'
            elif choice == 'c':
                print("Operation cancelled.")
                return 'cancel'
            else:
                print("Invalid choice. Please enter 'o', 'r', or 'c'.")
    return 'ok'

def generate_output_dirname_from_config(config, arbitrary_threshold_val=None):
    """Generates a concise, descriptive directory name from the configuration dictionary."""
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

    clustering_method = config.get('clustering_method', 'ward')
    if clustering_method == 'ward':
        dist_metric = config.get('distance_metric', 'none')
        if dist_metric == 'cosine': parts.append("distCos")
        elif dist_metric == 'euclidean': parts.append("distEuc")
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
    elif clustering_method == 'shc':
        metric = config.get('shc_metric', 'euc')
        alpha_val = config.get('shc_alpha', 0)
        lambda_val = config.get('shc_lambda_', 0)
        n_clusters = config.get('shc_n_clusters', 0)
        parts.append(f"SHC_{metric[:3]}_a{str(alpha_val).replace('.','p')}_l{str(lambda_val).replace('.','p')}_k{n_clusters}")

    return "_".join(parts)


def save_metadata(output_dir, filename="metadata.txt", config=None, comments=None):
    """Saves run settings and comments to a metadata file."""
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Execution Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("--- Configuration ---\n")
        if config:
            display_order = [
                'working_directory', 'filter_method', 'target_extraction_percentage', 'filter_quantile',
                'variable_selection_method', 'custom_columns', 'standardization_method',
                'pca_method', 'pca_n_components_value', 'clustering_method',
                'distance_metric', 'threshold_method', 'shc_metric',
                'shc_alpha', 'shc_lambda_', 'shc_n_clusters', 'shc_n_jobs'
            ]
            
            for key in display_order:
                if key in config:
                    f.write(f"{key}: {config[key]}\n")
            
            for key, value in config.items():
                if key not in display_order:
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
        print(f"Warning: {bash_file_path} not found.")
        while True:
            new_path = input(f"Please enter the correct path for run_PDTI_TA.bash (or press Enter to skip): ")
            if not new_path: raise FileNotFoundError("run_PDTI_TA.bash not found, skipping parameter reading.")
            if os.path.exists(new_path):
                bash_file_path = new_path; break
            else: print(f"The specified path '{new_path}' does not exist.")

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

def stage_filter_data(snap4_df, method, quantile_val=0.75, advanced_filter_params=None):
    """Data Filtering"""
    print(f"\n--- 2. Data Filtering ---")
    print(f"Method: {method}")
    df_processed_adv_filter_for_saving = None

    if method == 'arbitrary_percentage':
        threshold_sliprate = snap4_df['sliprate'].quantile(quantile_val)
        filtered_df = snap4_df[snap4_df['sliprate'] > threshold_sliprate].copy()
        print(f"  Filtered by sliprate > {threshold_sliprate:.4f} (quantile={quantile_val}).")
        print(f"  Before filtering: {snap4_df.shape[0]} rows, After filtering: {filtered_df.shape[0]} rows")
    elif method == 'skewness_kurtosis':
        print("  Warning: 'skewness_kurtosis' filtering is not implemented. No filtering will be performed.")
        filtered_df = snap4_df.copy()
    elif method == 'advanced_filtering':
        if advanced_filter_params is None:
            raise ValueError("advanced_filter_params are required for advanced filtering.")

        df_processed_adv_filter_for_saving, filtered_df = stage_advanced_filtering(
            snap4_df.copy(),
            advanced_filter_params['n_interval'],
            advanced_filter_params['m_interval'],
            advanced_filter_params['delta_tw'],
            advanced_filter_params['target_extraction_percentage']
        )
        print(f"  Advanced filtering applied. Before: {snap4_df.shape[0]} rows, After: {filtered_df.shape[0]} rows")
    else:
        print(f"  Warning: Unknown filtering method '{method}'. No filtering will be performed.")
        filtered_df = snap4_df.copy()

    if filtered_df.empty and method != 'advanced_filtering':
        raise ValueError("Data became empty after filtering. Please check your filtering conditions.")

    return filtered_df, df_processed_adv_filter_for_saving


def stage_select_variables(df_snap_original, snap4_filtered, method, custom_columns=None):
    """Variable Selection"""
    print(f"\n--- 3. Variable Selection ---")

    if isinstance(method, str):
        method = method.strip()
    else:
        raise TypeError(f"Variable selection method must be a string, but got {type(method)}.")

    print(f"Method: {method}")

    base_matrix_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    base_nv6_cols = ['n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)']
    common_cols_suffix = {
        "only": [], "nm": ['n', 'm'], "nmtw": ['n', 'm', 'tw'],
        "nmtw_sliprate": ['n', 'm', 'tw', 'sliprate']
    }
    data_for_clustering = pd.DataFrame()

    if method == 'custom_columns':
        if not custom_columns:
            raise ValueError("Custom columns selected but no columns were provided.")
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

    elif method.startswith("matrix_components_"):
        suffix_key = method.replace("matrix_components_", "")
        if suffix_key in common_cols_suffix:
            additional_cols_keys = common_cols_suffix[suffix_key]

            if suffix_key == "only":
                if not snap4_filtered.index.is_unique:
                    print("Warning: Index of snap4_filtered is not unique. Removing duplicates.")
                    snap4_filtered_unique_idx = snap4_filtered[~snap4_filtered.index.duplicated(keep='first')]
                else:
                    snap4_filtered_unique_idx = snap4_filtered

                common_indices = df_snap_original.index.intersection(snap4_filtered_unique_idx.index)
                df_snap_subset = df_snap_original.loc[common_indices]

                final_cluster_cols = [col for col in base_matrix_cols if col in df_snap_subset.columns]
                if not final_cluster_cols:
                        raise ValueError(f"For matrix_components_only, no usable columns found in df_snap_subset. Required: {base_matrix_cols}")
                data_for_clustering = df_snap_subset[final_cluster_cols].copy()
                data_for_clustering = data_for_clustering.reindex(snap4_filtered_unique_idx.index)
            else:
                columns_to_cluster_from_snap = base_matrix_cols + [col for col in additional_cols_keys if col in df_snap_original.columns]
                snap4_temp_keys = snap4_filtered[['x_grid', 'y_grid', 'tw_grid']].copy()
                snap4_temp_keys.rename(columns={'x_grid':'n', 'y_grid':'m', 'tw_grid':'tw'}, inplace=True)

                merge_on_keys = ['n', 'm', 'tw']

                if not all(key in df_snap_original.columns for key in merge_on_keys):
                    missing_keys = [key for key in merge_on_keys if key not in df_snap_original.columns]
                    raise ValueError(f"df_snap_original is missing merge keys: {missing_keys}")

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
                    raise ValueError("Merge failed or required columns are missing for matrix component variable selection.")
                data_for_clustering = merged_df[final_cluster_cols].copy()
                data_for_clustering.index = snap4_filtered.index
        else:
            raise ValueError(f"Unknown matrix component suffix: {suffix_key}")

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
                raise ValueError(f"No usable columns found in snap4_filtered for nv6 component selection.")
            data_for_clustering = snap4_filtered[columns_to_cluster].copy()
        else:
            raise ValueError(f"Unknown nv6 component suffix: {suffix_key}")
    else:
        raise ValueError(f"Unknown variable selection category: {method}")

    if data_for_clustering.isnull().values.any():
        print("Warning: Selected variables contain missing values. Dropping rows with NaNs.")
        data_for_clustering.dropna(inplace=True)
        if data_for_clustering.empty:
            raise ValueError("Data became empty after handling missing values in variable selection.")
    print(f"  Variables to be used: {data_for_clustering.columns.tolist()}")
    print(f"  Shape of data for clustering: {data_for_clustering.shape}")
    return data_for_clustering, data_for_clustering.columns.tolist()


def stage_standardize_data(df_for_clustering, columns_to_cluster, apply_standardization):
    print(f"\n--- 4. Data Standardization ---")
    scaler_model = None
    target_cols_for_std = [col for col in columns_to_cluster if col in df_for_clustering.columns]
    numeric_cols_for_std = df_for_clustering[target_cols_for_std].select_dtypes(include=np.number).columns.tolist()
    if not numeric_cols_for_std:
        print("  Warning: No numerical columns to standardize. Skipping standardization.")
        return df_for_clustering, None
    if apply_standardization:
        print(f"  Method: On (Applying StandardScaler to columns {numeric_cols_for_std})")
        scaler_model = StandardScaler()
        df_for_clustering[numeric_cols_for_std] = scaler_model.fit_transform(df_for_clustering[numeric_cols_for_std])
    else:
        print(f"  Method: Off")
    return df_for_clustering, scaler_model

def stage_pca(data_scaled, method, output_dir=None, config=None):
    """Principal Component Analysis (PCA)"""
    if config is None: config = {}
    print(f"\n--- 5. Principal Component Analysis (PCA) ---")
    print(f"Method: {method}")
    pca_model = None
    data_after_pca = data_scaled

    if method == 'pca_none':
        print("  PCA will not be performed.")
        return data_scaled, None

    data_for_pca_fit = data_scaled
    if isinstance(data_scaled, pd.DataFrame):
        numeric_cols = data_scaled.select_dtypes(include=np.number).columns
        if len(numeric_cols) == 0:
            print("  Warning: No numerical data columns, PCA cannot be performed. Skipping PCA.")
            return data_scaled, None
        if len(numeric_cols) < data_scaled.shape[1]:
            print(f"  Warning: PCA will only be performed on numerical columns: {list(numeric_cols)}")
        data_for_pca_fit = data_scaled[numeric_cols]

    if data_for_pca_fit.shape[1] < 1:
        print("  Warning: No data columns, PCA cannot be performed. Skipping PCA.")
        return data_scaled, None

    max_components = data_for_pca_fit.shape[1]

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

    if isinstance(data_scaled, pd.DataFrame):
        data_after_pca = pd.DataFrame(data_after_pca_transformed, index=data_scaled.index,
                                      columns=[f'PC{i+1}' for i in range(n_components)])
    else:
        data_after_pca = data_after_pca_transformed
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

    return data_after_pca, pca_model


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

def stage_determine_cluster_threshold(linkage_matrix, data_for_silhouette, method, distance_metric_for_silhouette, output_dir, config=None, columns_used=None):
    """Determines the cluster distance threshold for Ward's method."""
    if config is None: config = {}
    print(f"\n--- 7. Cluster Distance Threshold (Ward's Method) ---")
    print(f"  Method: {method}")

    fig_dendro_full, ax_dendro_full = plt.subplots(figsize=(16, 10))
    dendrogram(linkage_matrix, ax=ax_dendro_full, color_threshold=0)
    ax_dendro_full.set_title("Hierarchical Clustering Dendrogram (Full)")
    ax_dendro_full.set_xlabel("Sample index or (Cluster size)")
    ax_dendro_full.set_ylabel("Distance")
    dendrogram_full_path = os.path.join(output_dir, "dendrogram_full_for_thresholding.png")
    plt.savefig(dendrogram_full_path)
    print(f"  Full dendrogram saved to {dendrogram_full_path}.")

    if not config.get('run_all_combinations_mode', False):
        plt.show()
    plt.close(fig_dendro_full)

    final_distance_threshold = None
    user_input_threshold = None

    if method == 'threshold_arbitrary':
        target_k_for_arbitrary = config.get('target_cluster_count_for_arbitrary')
        if config.get('run_all_combinations_mode') and target_k_for_arbitrary is not None:
            print(f"  All-combinations mode: Calculating distance threshold for target cluster count {target_k_for_arbitrary}.")
            if linkage_matrix.shape[0] + 1 < target_k_for_arbitrary or target_k_for_arbitrary < 2:
                raise ValueError(f"Target cluster count {target_k_for_arbitrary} is inappropriate. Number of samples: {linkage_matrix.shape[0]+1}")

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
            print(f"  Calculated distance threshold (for target K={target_k_for_arbitrary}): {final_distance_threshold:.6f}")
        else:
            print(f"  Please refer to the dendrogram ({dendrogram_full_path}) and enter a distance threshold.")
            while True:
                try:
                    user_threshold_str = input("  Enter distance threshold: ")
                    final_distance_threshold = float(user_threshold_str)
                    if final_distance_threshold <= 0:
                        print("  Threshold must be a positive number.")
                        continue
                    user_input_threshold = final_distance_threshold
                    break
                except ValueError:
                    print("  Invalid number. Please try again.")
            print(f"  User-specified distance threshold: {final_distance_threshold}")

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

        upper_k_limit = min(11, actual_data_for_sil_values.shape[0])
        if upper_k_limit <=2 :
            raise ValueError(f"Too few data points ({actual_data_for_sil_values.shape[0]}) to test multiple K values in silhouette analysis.")
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

        if columns_used:
            display_text = "columns_used: " + ", ".join(columns_used)
            fig_sil.text(0.5, 0.01, display_text, ha='center', va='bottom', fontsize=8, wrap=True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

        silhouette_plot_path = os.path.join(output_dir, "silhouette_scores.png")
        plt.savefig(silhouette_plot_path)
        print(f"  Silhouette score plot saved to {silhouette_plot_path}.")
        if not config.get('run_all_combinations_mode', False): plt.show()
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

    elif method == 'threshold_elbow':
        raise NotImplementedError("Elbow method is not implemented.")
    elif method == 'threshold_hybrid':
        raise NotImplementedError("Hybrid method is not implemented.")
    else:
        raise ValueError(f"Unknown threshold determination method: {method}")

    return final_distance_threshold, user_input_threshold

def stage_shc_clustering_and_labeling(data_for_shc, config, output_dir):
    """Performs SHC clustering and labels data based on user input."""
    print("\n--- 6. SHC (Stable Clustering) ---")
    
    alpha = config['shc_alpha']
    lambda_ = config['shc_lambda_']
    metric = config['shc_metric']
    n_jobs = config.get('shc_n_jobs', 1)
    max_depth = data_for_shc.shape[0]
    
    print(f"  Parameters: alpha={alpha}, lambda={lambda_}, metric={metric}, n_jobs={n_jobs}")
    
    shc = StableHierarchicalClustering(alpha=alpha, lambda_=lambda_, max_depth=max_depth, metric=metric, n_jobs=n_jobs)
    
    if isinstance(data_for_shc, pd.DataFrame):
        data_for_shc_np = data_for_shc.values
    else:
        data_for_shc_np = data_for_shc

    linkage_matrix = shc.fit(data_for_shc_np)
    
    if linkage_matrix.size == 0:
        raise ValueError("SHC clustering failed, resulting in an empty linkage matrix.")
        
    print(f"  Generated SHC linkage matrix. Shape: {linkage_matrix.shape}")

    plt.figure(figsize=(16, 10))
    plt.title(f'SHC Dendrogram ({metric}) - Please determine the number of clusters')
    dendrogram(linkage_matrix)
    plt.xlabel("Sample index")
    plt.ylabel("Distance (depth-based)")
    dendrogram_shc_path = os.path.join(output_dir, f"dendrogram_shc_full_{metric}.png")
    plt.savefig(dendrogram_shc_path)
    print(f"  SHC dendrogram saved to {dendrogram_shc_path}.")
    print("  Please refer to the displayed dendrogram to decide on the number of clusters.")
    plt.show()

    while True:
        try:
            n_clusters_str = input("  Enter the desired number of clusters (e.g., 4): ")
            n_clusters = int(n_clusters_str)
            if 1 < n_clusters < data_for_shc.shape[0]:
                config['shc_n_clusters'] = n_clusters
                break
            else:
                print(f"  Number of clusters must be between 2 and the number of data points ({data_for_shc.shape[0]}).")
        except ValueError:
            print("  Invalid number. Please enter an integer.")
            
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
        info_text += f"  Filter: {run_config.get('filter_method')} ({run_config.get('filter_quantile', 'N/A') if run_config.get('filter_method') == 'arbitrary_percentage' else run_config.get('target_extraction_percentage','N/A') if run_config.get('filter_method') == 'advanced_filtering' else 'N/A'})\n"
        info_text += f"  Variables: {run_config.get('variable_selection_method')}\n"
        if run_config.get('variable_selection_method') == 'custom_columns':
            info_text += f"    Cols: {run_config.get('custom_columns')}\n"
        info_text += f"  PCA: {run_config.get('pca_method')}\n"
        info_text += f"  Standardize: {run_config.get('standardization_method')}\n"
        info_text += f"  Clustering: {run_config.get('clustering_method')}\n"
        if run_config.get('clustering_method') == 'ward':
            info_text += f"    Distance: {run_config.get('distance_metric')}\n"
            info_text += f"    Threshold Method: {run_config.get('threshold_method')}\n"
        elif run_config.get('clustering_method') == 'shc':
            info_text += f"    Metric: {run_config.get('shc_metric')}\n"
            info_text += f"    alpha: {run_config.get('shc_alpha')}\n"
            info_text += f"    lambda: {run_config.get('shc_lambda_')}\n"
            info_text += f"    n_clusters: {run_config.get('shc_n_clusters')}\n"
            info_text += f"    n_jobs: {config.get('shc_n_jobs', 1)}\n"

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
                print(f"    Successfully copied: {filename} -> {destination_path}")
            except Exception as e:
                print(f"    Warning: An error occurred while copying {filename}: {e}")
        else:
            print(f"    Warning: Source file not found, skipping: {source_path}")


# --- Main Processing Function ---
def run_single_pipeline(config, df_snap_orig, snap4_orig, pipeline_run_count=None):
    """Executes a single clustering pipeline with a given configuration."""
    print(f"\n=== Starting Pipeline Run (Count: {pipeline_run_count if pipeline_run_count is not None else 'N/A'}) ===")
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
            print("Warning: Data for clustering is empty after variable selection. Skipping this pipeline run.")
            if config.get('output_directory_path_override'):
                current_run_output_dir = config['output_directory_path_override']
                ensure_dir(current_run_output_dir)
                save_metadata(current_run_output_dir, config=config, comments="Skipped: No data after variable selection.")
            return

        data_std, scaler = stage_standardize_data(
            data_for_clustering.copy(), actual_cols_used,
            apply_standardization=(config['standardization_method'] == 'standardization_on')
        )

        data_pca, pca_model = stage_pca(
            data_std.copy(),
            config['pca_method'],
            output_dir=current_run_output_dir if current_run_output_dir and os.path.isdir(current_run_output_dir) else None,
            config=config
        )

        # --- Branching based on clustering method ---
        linkage_mat = None
        final_labels = None

        if config['clustering_method'] == 'ward':
            linkage_mat = stage_linkage_clustering(
                data_pca, config['distance_metric']
            )
            
            output_for_thresh_stage_final = os.path.join(config['working_directory'], f"_temp_ward_{datetime.now().strftime('%Y%m%d%H%M%S')}")
            ensure_dir(output_for_thresh_stage_final)

            sil_data_input = data_pca if config['pca_method'] != 'pca_none' else data_std

            dist_thresh, user_arb_thresh_val = stage_determine_cluster_threshold(
                linkage_mat, sil_data_input,
                method=config['threshold_method'],
                distance_metric_for_silhouette=config['distance_metric'],
                output_dir=output_for_thresh_stage_final,
                config=config,
                columns_used=actual_cols_used
            )
            if dist_thresh is None:
                print("Threshold determination was cancelled or failed. Terminating pipeline.")
                if os.path.exists(output_for_thresh_stage_final):
                    try: shutil.rmtree(output_for_thresh_stage_final)
                    except OSError as e: print(f"Error removing temp directory: {e}")
                return

            config['final_cluster_threshold'] = dist_thresh
            final_labels = fcluster(linkage_mat, dist_thresh, criterion='distance')

            if config['threshold_method'] == 'threshold_arbitrary' and user_arb_thresh_val is not None:
                arbitrary_threshold_value_for_filename = user_arb_thresh_val
            
            final_run_specific_dirname = generate_output_dirname_from_config(config, arbitrary_threshold_value_for_filename)
            current_run_output_dir = os.path.join(config['working_directory'], final_run_specific_dirname)
            
            overwrite_choice = handle_overwrite_choice(current_run_output_dir)
            if overwrite_choice == 'cancel': return
            elif overwrite_choice == 'rename':
                current_run_output_dir = f"{current_run_output_dir}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
            ensure_dir(current_run_output_dir)
            
            for item_name in os.listdir(output_for_thresh_stage_final):
                shutil.move(os.path.join(output_for_thresh_stage_final, item_name), os.path.join(current_run_output_dir, item_name))
            shutil.rmtree(output_for_thresh_stage_final)


        elif config['clustering_method'] == 'shc':
            temp_shc_output_dir = os.path.join(config['working_directory'], f"_temp_shc_{datetime.now().strftime('%Y%m%d%H%M%S')}")
            ensure_dir(temp_shc_output_dir)

            linkage_mat, final_labels_series = stage_shc_clustering_and_labeling(
                data_pca, config, temp_shc_output_dir
            )
            final_labels = final_labels_series.values
            
            final_run_specific_dirname = generate_output_dirname_from_config(config)
            current_run_output_dir = os.path.join(config['working_directory'], final_run_specific_dirname)

            overwrite_choice = handle_overwrite_choice(current_run_output_dir)
            if overwrite_choice == 'cancel': return
            elif overwrite_choice == 'rename':
                current_run_output_dir = f"{current_run_output_dir}_{datetime.now().strftime('%Y%m%d%H%M%S')}"

            ensure_dir(current_run_output_dir)
            for item_name in os.listdir(temp_shc_output_dir):
                shutil.move(os.path.join(temp_shc_output_dir, item_name), os.path.join(current_run_output_dir, item_name))
            shutil.rmtree(temp_shc_output_dir)

        # --- Save Results ---
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
            print("Error: Cluster labels were not generated.")


    except NotImplementedError as e:
        print(f"Error: An unimplemented feature was called: {e}")
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: NotImplementedError - {e}")
    except FileNotFoundError as e:
        print(f"Error: A file was not found: {e}")
    except ValueError as e:
        print(f"Error: A value-related issue occurred during data processing: {e}")
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: ValueError - {e}")
    except Exception as e:
        import traceback
        print(f"An unexpected error occurred: {e}")
        print(traceback.format_exc())
        if current_run_output_dir and os.path.isdir(current_run_output_dir):
            save_metadata(current_run_output_dir, config=config, comments=f"Error: Unexpected error - {e}\n{traceback.format_exc()}")


def get_user_config(default_wd, snap_cols, nvector_cols):
    """Interactively gets the configuration settings from the user."""
    config = {}
    print("\n--- Clustering Pipeline Configuration ---")
    config["working_directory"] = default_wd

    print("\n2. Data Filtering Method:")
    filter_options = {'1': 'arbitrary_percentage',
                      '2': 'skewness_kurtosis',
                      '3': 'advanced_filtering'}
    [print(f"  {k}: {v}{' (Not implemented)' if v == 'skewness_kurtosis' else ''}") for k, v in filter_options.items()]
    choice = input("    Select (1-3, default 3): ") or '3'
    config['filter_method'] = filter_options.get(choice, 'advanced_filtering')

    if config['filter_method'] == 'arbitrary_percentage':
        while True:
            q_val_str = input("    Enter sliprate filter percentile (0.0-1.0, e.g., 0.75 extracts top 25%, default 0.75): ") or "0.75"
            try:
                val = float(q_val_str); assert 0.0 <= val <= 1.0
                config['filter_quantile'] = val; break
            except (ValueError, AssertionError): print("    Value must be between 0.0 and 1.0.")
    elif config['filter_method'] == 'advanced_filtering':
        while True:
            perc_str = input("    Enter percentage of major data to extract (0.0-1.0, e.g., 0.25 for top 25%, default 0.25): ") or "0.25"
            try:
                val = float(perc_str); assert 0.0 < val <= 1.0
                config['target_extraction_percentage'] = val; break
            except (ValueError, AssertionError): print("    Value must be greater than 0.0 and at most 1.0.")


    print("\n3. Variable Selection Method:")
    var_methods_list = [
        "matrix_components_only", "matrix_components_nm", "matrix_components_nmtw", "matrix_components_nmtw_sliprate",
        "nv6_components_only", "nv6_components_nm", "nv6_components_nmtw", "nv6_components_nmtw_sliprate"
    ]
    var_options = {str(i+1): method for i, method in enumerate(var_methods_list)}

    custom_choice_num = len(var_options) + 1
    var_options[str(custom_choice_num)] = "Specify custom columns"

    default_var_choice = '8'
    for k, v in var_options.items(): print(f"  {k}: {v}")
    choice = input(f"    Select (1-{custom_choice_num}, default {default_var_choice}): ") or default_var_choice

    if choice == str(custom_choice_num):
        config['variable_selection_method'] = 'custom_columns'
        available_cols = sorted(list(set(snap_cols + nvector_cols)))
        print("\n    Available columns:")
        for i, col in enumerate(available_cols):
            print(f"      {i+1}: {col}")
        while True:
            try:
                col_indices_str = input("    Enter comma-separated column numbers (e.g., 1,2,8,9,10): ")
                selected_indices = [int(i.strip()) - 1 for i in col_indices_str.split(',')]
                if all(0 <= idx < len(available_cols) for idx in selected_indices):
                    config['custom_columns'] = [available_cols[i] for i in selected_indices]
                    break
                else:
                    print("      Invalid number included. Please try again.")
            except ValueError:
                print("      Invalid format. Please enter comma-separated numbers.")
    else:
        config['variable_selection_method'] = var_options.get(choice, var_methods_list[int(default_var_choice)-1])

    print("\n4. Data Standardization Method:")
    std_options = {'1': 'standardization_on', '2': 'standardization_off'}
    [print(f"  {k}: {'On' if v == 'standardization_on' else 'Off'}") for k,v in std_options.items()]
    choice = input("    Select (1-2, default 1): ") or '1'
    config['standardization_method'] = std_options.get(choice, 'standardization_on')

    print("\n5. Principal Component Analysis (PCA) Method:")
    pca_options = {
        '1': 'pca_80',
        '2': 'pca_90',
        '3': 'pca_none',
        '4': 'pca_n_components'
    }
    print("  1: 80% cumulative variance")
    print("  2: 90% cumulative variance")
    print("  3: None")
    print("  4: Specify number of components")
    choice = input("    Select (1-4, default 1): ") or '1'
    config['pca_method'] = pca_options.get(choice, 'pca_80')

    if config['pca_method'] == 'pca_n_components':
        while True:
            try:
                n_comp_str = input("    Enter number of principal components (e.g., 3): ")
                n_comp = int(n_comp_str)
                if n_comp > 0:
                    config['pca_n_components_value'] = n_comp
                    break
                else:
                    print("    Please enter an integer greater than 0.")
            except ValueError:
                print("    Invalid number. Please enter an integer.")

    print("\n6. Clustering Method:")
    clustering_options = {'1': 'ward', '2': 'shc'}
    print("  1: Hierarchical Clustering (Ward's Method)")
    print("  2: SHC (Stable Clustering)")
    
    choice = input("    Select (1-2, default 1): ") or '1'
    config['clustering_method'] = clustering_options.get(choice, 'ward')

    if config['clustering_method'] == 'ward':
        print("\n7. Distance Metric (for Ward's Method):")
        dist_options = {'1': 'cosine', '2': 'euclidean'}
        [print(f"  {k}: {v} distance") for k, v in dist_options.items()]
        choice = input("    Select (1-2, default 1): ") or '1'
        config['distance_metric'] = dist_options.get(choice, 'cosine')

        print("\n8. Cluster Threshold Method (for Ward's Method):")
        thresh_options = {
            '1': 'threshold_arbitrary', '2': 'threshold_silhouette',
            '3': 'threshold_elbow', '4': 'threshold_hybrid'
        }
        unimplemented_thresh = ['threshold_elbow', 'threshold_hybrid']
        for k, v_raw in thresh_options.items():
            v_display = v_raw.replace('threshold_','').replace('_',' ').capitalize()
            print(f"  {k}: {v_display}{' (Not implemented)' if v_raw in unimplemented_thresh else ''}")
        choice = input("    Select (1-4, default 2): ") or '2'
        config['threshold_method'] = thresh_options.get(choice, 'threshold_silhouette')
    
    elif config['clustering_method'] == 'shc':
        print("\n7. SHC Parameter Configuration:")
        
        if JOBLIB_AVAILABLE:
            print("  Parallel Processing:")
            while True:
                try:
                    import multiprocessing
                    cpu_count = multiprocessing.cpu_count()
                    n_jobs_str = input(f"    Enter number of CPU cores to use (Your system has {cpu_count} cores. -1 to use all, default: -1): ") or "-1"
                    n_jobs = int(n_jobs_str)
                    config['shc_n_jobs'] = n_jobs
                    break
                except ValueError:
                    print("    Invalid number. Please enter an integer.")
        
        print("  SHC Distance Metric:")
        shc_metric_options = {'1': 'euclidean', '2': 'cosine'}
        print("    1: Euclidean (default)")
        print("    2: Cosine")
        metric_choice = input("      Select (1-2, default 1): ") or '1'
        config['shc_metric'] = shc_metric_options.get(metric_choice, 'euclidean')

        while True:
            try:
                alpha_str = input("    Enter alpha (similarity sensitivity, positive, default 1.0): ") or "1.0"
                alpha = float(alpha_str)
                if alpha > 0:
                    config['shc_alpha'] = alpha
                    break
                else:
                    print("    alpha must be a positive number.")
            except ValueError:
                print("    Invalid number.")

        while True:
            try:
                lambda_str = input("    Enter lambda (stability control, non-negative, default 10.0): ") or "10.0"
                lambda_ = float(lambda_str)
                if lambda_ >= 0:
                    config['shc_lambda_'] = lambda_
                    break
                else:
                    print("    lambda must be a non-negative number.")
            except ValueError:
                print("    Invalid number.")

    return config

def run_all_combinations(base_working_dir, df_snap_orig, snap4_orig):
    """Executes the pipeline for all combinations of parameters."""
    print("\n--- All-Combinations Mode currently supports Ward's method only ---")
    pass

def main():
    """Main execution function."""
    print("Main function started.")
    default_working_directory = "./"

    working_directory = input(f"1. Enter working directory (press Enter for default: {default_working_directory}): ") or default_working_directory
    while not os.path.isdir(working_directory):
        print(f"Error: The specified directory '{working_directory}' does not exist or is not a directory.")
        working_directory = input(f"    Please enter a valid working directory (e.g., {default_working_directory}): ")
        if not working_directory:
            working_directory = default_working_directory
            if not os.path.isdir(working_directory):
                print("Default working directory is also invalid. Exiting program."); exit()

    print(f"Using working directory: {working_directory}")

    try:
        df_snap, df_nvector = stage_load_data(working_directory)
    except FileNotFoundError as e:
        print(f"Error: Failed to load data files: {e}")
        return
    except Exception as e:
        print(f"Error: Could not load data due to an unexpected issue: {e}")
        return

    print("\nSelect execution mode:")
    print("  1: Interactive mode for a single run")
    print("  2: Pre-defined settings for multiple runs (requires code modification)")
    print("  3: All-combinations mode (Ward's method only)")
    run_mode = input("Select (1-3, default 1): ") or '1'
    print(f"Selected execution mode: {run_mode}")

    if run_mode == '1':
        user_config = get_user_config(working_directory, df_snap.columns.tolist(), df_nvector.columns.tolist())
        if user_config:
            run_single_pipeline(user_config, df_snap_orig=df_snap, snap4_orig=df_nvector)
    elif run_mode == '2':
        print("Warning: Batch mode is experimental and does not support interactive SHC runs.")
        configs_to_run = [
            {
                "working_directory": working_directory, "filter_method": "advanced_filtering",
                "target_extraction_percentage": 0.25,
                "variable_selection_method": "nv6_components_nmtw_sliprate",
                "standardization_method": "standardization_on", "pca_method": "pca_90",
                "clustering_method": "ward", "distance_metric": "cosine", "threshold_method": "threshold_silhouette",
            }
        ]
        for i, cfg in enumerate(configs_to_run):
            print(f"\n--- Batch Run {i+1}/{len(configs_to_run)} ---")
            run_single_pipeline(cfg, df_snap_orig=df_snap, snap4_orig=df_nvector, pipeline_run_count=f"batch_{i+1:02d}")
    elif run_mode == '3':
        print("Starting all-combinations mode for Ward's method...")
        run_all_combinations(working_directory, df_snap_orig=df_snap, snap4_orig=df_nvector)
        print("All-combinations mode is extensive and currently disabled by default. Please check the code to run.")
    else:
        print("Invalid execution mode selected. Exiting.")

if __name__ == "__main__":
    main()

