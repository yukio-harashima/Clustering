#!/usr/bin/env python3
# クラスターごとの要素点の重心を求める (飛び地対応・出力形式維持版)
import numpy as np
import pandas as pd
from tqdm import tqdm
import os # エラー解消のためosモジュールをインポート

def find_sub_clusters(cluster_df):
    """
    クラスタ内のデータ点群から、空間的に連続したサブクラスタ（飛び地）を見つける。
    Args:
        cluster_df (pd.DataFrame): 単一のクラスタ番号(clno)に属するデータ。'n'と'm'列が必要。
    Returns:
        list[set]: サブクラスタのリスト。各サブクラスタは座標タプル(n, m)のセット。
    """
    if cluster_df.empty:
        return []
    
    # まだ訪問していない点のセットを作成
    points_to_visit = set(zip(cluster_df['n'], cluster_df['m']))
    sub_clusters = []

    # すべての点が訪問されるまでループ
    while points_to_visit:
        # 新しいサブクラスタを開始
        current_sub_cluster = set()
        # 訪問キュー（探索の開始点）
        queue = [points_to_visit.pop()]
        current_sub_cluster.add(queue[0])

        # 幅優先探索(BFS)で隣接するすべての点を見つける
        while queue:
            current_n, current_m = queue.pop(0)
            
            # 上下左右の隣接点をチェック
            for dn, dm in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (current_n + dn, current_m + dm)
                
                if neighbor in points_to_visit:
                    points_to_visit.remove(neighbor)
                    current_sub_cluster.add(neighbor)
                    queue.append(neighbor)
        
        sub_clusters.append(current_sub_cluster)
        
    return sub_clusters

# --- 1. データの読み込み ---
print("--- 1. データの読み込み ---")
od = "./" 
# od = "/Users/harashima-yukio/Desktop/2024_Noto/results_20240702-182819/" 
output_file = "clusteringSnap2.dat"
pss2 = os.path.join(od, output_file)

try:
    data = pd.read_table(pss2, sep='\t', header=None)
    data.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                   'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                        'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                            'trendp', 'trendt', 'trendb', 'plungp', 
                                'plungt', 'plungb', 'NDC', 'cllab', 'clno']

    datafile = "snap_yr.dat"
    ps = os.path.join(od, datafile)
    snp_yr = pd.read_table(ps, sep = r'\s+', header=None)
    snp_yr.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                   'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                        'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                            'trendp', 'trendt', 'trendb', 'plungp', 
                                'plungt', 'plungb', 'NDC']
    print("データの読み込みが完了しました。")
except FileNotFoundError as e:
    print(f"エラー: ファイルが見つかりません。パスを確認してください: {e}")
    exit()

# --- 2. データの前処理 ---
print("\n--- 2. データの前処理 ---")
data2 = data.copy()
if len(data2) == len(snp_yr):
    data2[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']] = snp_yr[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']]
    print("snap_yr.dat のモーメントテンソルデータで置き換えました。")
else:
    print("警告: clusteringSnap2.datとsnap_yr.datの行数が異なります。モーメントテンソルの置き換えをスキップします。")

# --- 3. 飛び地ごとの重心計算 ---
print("\n--- 3. 飛び地ごとの重心計算 ---")
# 各出力ファイルに対応する空のDataFrameリストを準備
df_list, centroids2_list, dfl_list, centroids2l_list, centroids4l_list = [], [], [], [], []

# clno > 0 のクラスタのみを対象
unique_clusters = sorted(data2[data2['clno'] > 0]['clno'].unique())

for clno in tqdm(unique_clusters, desc="Processing clusters"):
    cluster_df = data2[data2['clno'] == clno].copy()
    
    # クラスタ内の飛び地（サブクラスタ）を特定
    sub_clusters = find_sub_clusters(cluster_df)
    
    # 各サブクラスタについて重心を計算
    for sub_cluster_points in sub_clusters:
        sub_cluster_df = cluster_df[
            cluster_df[['n', 'm']].apply(tuple, axis=1).isin(sub_cluster_points)
        ]
        
        # --- ここから元コードの計算ロジックをサブグループDFに適用 ---
        # clusteringSnap3.dat用
        centroids = sub_cluster_df.groupby(['tw', 'clno'])[['dx', 'dy', 'sliprate']].mean().reset_index()
        moments = sub_cluster_df.groupby(['tw', 'clno'])[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].sum().reset_index()
        mo = moments.drop(moments.columns[[0, 1]], axis=1)
        df_list.append(pd.concat([centroids, mo], axis=1))

        # centroids.dat用
        centroids1 = sub_cluster_df.groupby(['tw', 'clno'])[['dx', 'dy', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].mean().reset_index()
        centroids2_list.append(centroids1)

        # clusteringSnap3l.dat用
        centroidsl = sub_cluster_df.groupby(['tw', 'clno'])[['lon','lat', 'sliprate']].mean().reset_index()
        momentsl = sub_cluster_df.groupby(['tw', 'clno'])[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].sum().reset_index()
        mol = momentsl.drop(momentsl.columns[[0, 1]], axis=1)
        dfl_list.append(pd.concat([centroidsl, mol], axis=1))
        
        # centroids2.dat用
        centroids1l = sub_cluster_df.groupby(['tw', 'clno'])[['lon','lat', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
                                                     'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 
                                                     'trendp', 'trendt', 'trendb', 'plungp','plungt', 'plungb', 'NDC']].mean().reset_index()
        centroids2l_list.append(centroids1l)
        
        # centroids4.dat用
        centroids3l = sub_cluster_df.groupby(['clno'])[['tw','lon','lat', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
                                                     'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 
                                                     'trendp', 'trendt', 'trendb', 'plungp','plungt', 'plungb', 'NDC']].mean().reset_index()
        centroids4l_list.append(centroids3l)

print("全サブクラスタの重心計算が完了しました。")

# --- 4. 結果の結合と出力 ---
print("\n--- 4. 結果の結合と出力 ---")

# 各リストが空でないことを確認してから結合
if df_list:
    df_final = pd.concat(df_list, ignore_index=True)
    output_file = "clusteringSnap3.dat"
    pss = os.path.join(od, output_file)
    df_final.to_csv(pss, sep='\t', index=False, header=False)
    print(f"  {output_file} を保存しました。")

if centroids2_list:
    centroids2_final = pd.concat(centroids2_list, ignore_index=True)
    output_file = "centroids.dat"
    pss = os.path.join(od, output_file)
    centroids2_final.to_csv(pss, sep='\t', index=False, header=False)
    print(f"  {output_file} を保存しました。")

if dfl_list:
    dfl_final = pd.concat(dfl_list, ignore_index=True)
    output_file = "clusteringSnap3l.dat"
    pss = os.path.join(od, output_file)
    dfl_final.to_csv(pss, sep='\t', index=False, header=False)
    print(f"  {output_file} を保存しました。")

if centroids2l_list:
    centroids2l_final = pd.concat(centroids2l_list, ignore_index=True)
    output_file = "centroids2.dat"
    pss = os.path.join(od, output_file)
    centroids2l_final.to_csv(pss, sep='\t', index=False, header=False)
    print(f"  {output_file} を保存しました。")

if centroids4l_list:
    centroids4l_final = pd.concat(centroids4l_list, ignore_index=True)
    output_file = "centroids4.dat"
    pss = os.path.join(od, output_file)
    centroids4l_final.to_csv(pss, sep='\t', index=False) # 元のコードの挙動に合わせてヘッダあり
    print(f"  {output_file} を保存しました。")

if not df_list: # いずれかのリストが空なら、処理対象がなかったことを示す
    print("警告: 計算対象となるクラスタが存在しなかったため、出力ファイルは生成されませんでした。")

print("\n処理が完了しました。")
