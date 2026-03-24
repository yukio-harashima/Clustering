import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import itertools
import matplotlib.pyplot as plt

class StableHierarchicalClustering:
    """
    論文「階層クラスタリングの安定化」で提案されたSHCアルゴリズムを実装したクラス。

    Args:
        alpha (float): 類似度計算の感度パラメータ (α)。
        lambda_ (float): 安定性制御パラメータ (λ)。lambdaはPythonの予約語のためlambda_を使用。
        max_depth (int): クラスタリングの最大深度 (D)。
        metric (str): 距離計算に使用する尺度 ('euclidean' または 'cosine')。
    """
    def __init__(self, alpha: float, lambda_: float, max_depth: int, metric: str = 'euclidean'):
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
        
        # --- 内部変数 ---
        self._W = None  # 類似度行列
        self._n_samples = 0
        self._linkage_matrix = []
        self._next_cluster_id = 0

    def fit(self, X: np.ndarray):
        """
        データセットXに対してSHCアルゴリズムを実行し、リンケージ行列を構築する。

        Args:
            X (np.ndarray): データセット (n_samples, n_features)。

        Returns:
            np.ndarray: SciPy互換のリンケージ行列。
        """
        # --- 1. メイン関数: SHC ---
        self._n_samples = X.shape[0]
        self._next_cluster_id = self._n_samples
        self._linkage_matrix = []
        
        # 1.1. グラフ構築 (類似度計算)
        self._build_similarity_graph(X)
        
        # 1.2. 再帰的クラスタリングの呼び出し
        initial_vertices = list(range(self._n_samples))
        self._recursive_clustering(initial_vertices, 0)
        
        return np.array(self._linkage_matrix)

    def _build_similarity_graph(self, X: np.ndarray):
        """
        データ点から重み付きの隣接行列（類似度行列）Wを構築する。
        """
        # ユーザーが選択した距離尺度に基づいて距離を計算
        if self.metric == 'euclidean':
            # ユークリッド距離の2乗を計算
            d_sq = squareform(pdist(X, 'sqeuclidean'))
        elif self.metric == 'cosine':
            # コサイン距離を計算し、2乗する
            # コサイン距離は 0 (同一) から 2 (正反対) の範囲
            d = squareform(pdist(X, 'cosine'))
            d_sq = d ** 2
        else:
            # __init__でチェック済みだが、念のため
            raise ValueError(f"Unsupported metric: {self.metric}")

        # 類似度を計算 w = exp(-alpha * distance^2)
        self._W = np.exp(-self.alpha * d_sq)

    def _recursive_clustering(self, V_sub: list, current_depth: int):
        """
        グラフを受け取り、それを再帰的に分割して階層的な木構造（リンケージ行列）を構築する。
        """
        # --- 2. 再帰関数: REC ---
        
        # 2.1. 終了条件の判定 (ベースケース)
        if current_depth >= self.max_depth or len(V_sub) <= 1:
            # 葉ノードまたは単一要素のクラスタの場合、その頂点IDを返す
            return V_sub[0] if len(V_sub) == 1 else -1

        # 2.2. グラフの分割
        S1_set = self._stable_sparsest_cut(V_sub)
        S1 = sorted(list(S1_set))
        S2 = sorted(list(set(V_sub) - S1_set))

        # 分割が空になった場合のフォールバック
        if not S1 or not S2:
            # 非常に稀だが、分割が失敗した場合は現在のクラスタを葉として扱う
            return V_sub[0] if len(V_sub) == 1 else self._next_cluster_id -1 # 最後に作られたクラスタID

        # 2.3. 再帰呼び出し (部分グラフの作成は不要、頂点リストを渡すだけ)
        left_child_id = self._recursive_clustering(S1, current_depth + 1)
        right_child_id = self._recursive_clustering(S2, current_depth + 1)

        # 2.4. 親ノードの構築とリンケージ行列の作成
        parent_id = self._next_cluster_id
        self._next_cluster_id += 1
        
        # 距離は、階層の深さに基づいて定義する (論文では未定義)
        # 深いほど距離が小さくなるように (max_depth - current_depth) とする
        distance = float(self.max_depth - current_depth)
        
        # リンケージ行: [クラスタID1, クラスタID2, 距離, 新クラスタの点数]
        linkage_row = [
            float(left_child_id),
            float(right_child_id),
            distance,
            float(len(V_sub))
        ]
        self._linkage_matrix.append(linkage_row)
        
        return parent_id

    def _stable_sparsest_cut(self, V_sub: list):
        """
        与えられたグラフ(頂点部分集合)を、指数メカニズムに基づいて確率的に2つに分割する。
        """
        # --- 3. 安定化の核: SSC ---
        if len(V_sub) <= 1:
            return set(V_sub)

        # 3.1. 全分割候補の生成と評価
        candidate_splits = []
        quality_scores = []
        
        # V_sub内のインデックスに変換
        v_map = {original_idx: new_idx for new_idx, original_idx in enumerate(V_sub)}
        
        # 3.1.1. 全ての重心ペアの組み合わせを作成
        for i_idx, j_idx in itertools.combinations(range(len(V_sub)), 2):
            i = V_sub[i_idx]
            j = V_sub[j_idx]
            
            # 3.2.a. 分割候補 S_ij の生成
            S_ij = set()
            for k in V_sub:
                if self._W[i, k] > self._W[j, k]:
                    S_ij.add(k)
            
            # 3.2.b. 品質 phi の計算 (Sparsity)
            S_complement = set(V_sub) - S_ij
            size_S = len(S_ij)
            size_S_complement = len(S_complement)
            
            if size_S == 0 or size_S_complement == 0:
                phi = np.inf # エッジケース対応①: 無意味な分割
            else:
                # カットの重みを計算
                cut_weight = self._W[np.ix_(list(S_ij), list(S_complement))].sum()
                phi = cut_weight / (size_S * size_S_complement)
            
            # 3.2.c. 候補と品質の保存
            candidate_splits.append(S_ij)
            quality_scores.append(phi)

        # 3.3. 指数メカニズムによる確率的選択
        quality_scores = np.array(quality_scores)
        
        # 3.3.a. 数値的安定化 (log-sum-exp trick の考え方)
        phi_min = np.min(quality_scores[np.isfinite(quality_scores)]) if np.any(np.isfinite(quality_scores)) else 0
        shifted_scores = quality_scores - phi_min
        
        # 3.3.b. 重みの計算
        weights = np.exp(-self.lambda_ * shifted_scores)
        weights[np.isinf(quality_scores)] = 0 # phi=inf の場合は重み0
        
        # 3.3.c. 確率の計算
        total_weight = np.sum(weights)
        
        if total_weight > 0:
            probabilities = weights / total_weight
        else: # 全ての重みが0になった場合 (フォールバック)
            # 有限のphiを持つ候補から均等に選ぶ
            valid_indices = np.where(np.isfinite(quality_scores))[0]
            if len(valid_indices) > 0:
                probabilities = np.zeros_like(quality_scores)
                probabilities[valid_indices] = 1.0 / len(valid_indices)
            else: # 全てinfの場合（通常ありえない）
                probabilities = np.ones_like(quality_scores) / len(quality_scores)

        # 3.3.d. 最終的な分割の選択
        chosen_index = np.random.choice(len(candidate_splits), p=probabilities)
        chosen_split = candidate_splits[chosen_index]
        
        # 3.4. 結果の返却
        return chosen_split

# --- 使用例 ---
if __name__ == '__main__':
    # ダミーデータの生成
    from sklearn.datasets import make_blobs
    from scipy.cluster.hierarchy import dendrogram
    
    print("--- SHCアルゴリズムのテスト実行 ---")
    # X, y = make_blobs(n_samples=100, centers=4, n_features=4, random_state=42, cluster_std=0.8)
    wd = './'
    input_filename = 'testdata.dat'
    ps = wd + input_filename
    df = pd.read_csv(ps, sep='\t')
    dff = df[['col_1','col_2']]
    # クラスタリング用のNumPy配列に変換
    X = dff.values

    plt.figure(figsize=(8, 7))
    plt.scatter(X[:, 0], X[:, 1], marker='o', c=df['lab'], s=25, edgecolor='k')
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.show()


    # SHCのパラメータ
    alpha = 1.0
    lambda_ = 100 
    max_depth = 14
    metric_to_use = 'cosine' # 'euclidean' または 'cosine' を選択

    # SHCクラスタリングの実行
    print(f"パラメータ: alpha={alpha}, lambda={lambda_}, max_depth={max_depth}, metric={metric_to_use}")
    shc = StableHierarchicalClustering(alpha=alpha, lambda_=lambda_, max_depth=max_depth, metric=metric_to_use)
    linkage_matrix = shc.fit(X)



    # 結果のデンドログラム表示
    plt.figure(figsize=(10, 7))
    plt.title(f'Stable Hierarchical Clustering Dendrogram (metric={metric_to_use})')
    dendrogram(linkage_matrix)
    plt.show()

    # -------------------
    # 比較: 従来のウォード法
    from scipy.cluster.hierarchy import linkage
    # linkage_matrix_ward = linkage(X, method='ward', metric=metric_to_use)
    linkage_matrix_ward = linkage(pdist(X, metric=metric_to_use), method='ward')
    plt.figure(figsize=(10, 7))
    plt.title(f'Standard Ward Method Dendrogram (metric={metric_to_use})')
    dendrogram(linkage_matrix_ward)
    plt.show()
