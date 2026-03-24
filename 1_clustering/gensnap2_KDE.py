#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# cluster_unifiedの出力結果から、最近傍法でsnap2.dat にクラスター番号を割り付ける

#/////////////////////////////////////////////////////////////////////////////////
import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
from tqdm import tqdm
import time

#///////////////////////////
# データ読み取り
# od = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241218-031730/const/" 
od = "./" 
output_file = "clusteringSnap2.dat"
pss2 = od + output_file
# dat.to_csv(pss, sep='\t', index=False, header=False)
data = pd.read_table(pss2, sep='\t', header=None)
data.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC', 'cllab', 'cluster']


datafile = "snap2_y.dat"
ps = od+datafile
snap2 = pd.read_table(ps, sep = r'\s+', header=None)
snap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

df_original = data
df_completion = snap2

#///////////////////////////////////////////
# 時間ごとにデータを補完する処理
cluster_results = []

# for tw_val in sorted(df_completion['tw'].unique()):
for tw_val in tqdm(sorted(df_completion['tw'].unique()), desc="gensnap2_KDE", ncols=100, ascii=True):
    # print(tw_val)
    # 補完データで現在のtw_valより小さい直近の元データの時間を選ぶ
    relevant_tw = df_original[df_original['tw'] <= tw_val]['tw'].max()

    # 現在の時間での元データをフィルタリングs
    original_subset = df_original[df_original['tw'] == relevant_tw]

    # KDEを実行するために座標のみを使用
    coordinates = original_subset[['dx', 'dy']].to_numpy()
    cluster_labels = original_subset['cluster'].to_numpy()

    # 座標に対してKDEを実行し、クラスター番号を補完
    kde = KernelDensity(bandwidth=5.0)  # 適宜バンド幅を調整
    kde.fit(coordinates)

    # 現在の補完データの座標を抽出
    completion_coords = df_completion[df_completion['tw'] == tw_val][['dx', 'dy']].to_numpy()

    # 各補完データ点に対する対数密度を計算
    log_densities = kde.score_samples(completion_coords)
    
    # 最も近い元データの点のクラスター番号を補完
    distances = np.linalg.norm(coordinates[:, np.newaxis] - completion_coords, axis=2)
    nearest_indices = np.argmin(distances, axis=0)
    completed_clusters = cluster_labels[nearest_indices]
    
    # 結果をリストに追加
    cluster_results.append(completed_clusters)
    time.sleep(0.01)

# 補完したクラスター番号を結合して新しい列として追加
df_completion['cluster'] = np.concatenate(cluster_results)

# 結果を表示
print(df_completion)

#///////////////////////////////////////////
# 結果のデータを出力
od = "./" 
output_file = "assignsnap2KDE.dat"
pss = od + output_file
df_completion.to_csv(pss, sep='\t', index=False, header=False)


#/////////////////////////////////////////////////////////////////////////////////
# # 元データの作成 (5間隔でx、3.8間隔でy、4秒間隔でtw)
# x_vals = np.arange(1, 101, 5)  # xは1から100まで5間隔
# y_vals = np.arange(1, 99, 3.8)  # yは1から98.8まで3.8間隔
# tw_vals = np.arange(1, 14, 4)  # twは1から13まで4秒間隔

# # 全組み合わせを生成
# xx, yy, tt = np.meshgrid(x_vals, y_vals, tw_vals, indexing='ij')
# xx = xx.flatten()
# yy = yy.flatten()
# tt = tt.flatten()

# # クラスター番号をランダムに生成
# cluster_vals = np.random.randint(0, 10, len(xx))

# # 元データフレーム
# df_original = pd.DataFrame({
#     'x': xx,
#     'y': yy,
#     'tw': tt,
#     'cluster': cluster_vals
# })

# # 補完データの作成 (1間隔でx、y、1秒間隔でtw)
# x_vals_comp = np.arange(1, 101, 1)  # xは1から100まで1間隔
# y_vals_comp = np.arange(1, 101, 1)  # yは1から100まで1間隔
# tw_vals_comp = np.arange(1, 14, 1)  # twは1から13まで1秒間隔

# # 全組み合わせを生成
# xx_comp, yy_comp, tt_comp = np.meshgrid(x_vals_comp, y_vals_comp, tw_vals_comp, indexing='ij')
# xx_comp = xx_comp.flatten()
# yy_comp = yy_comp.flatten()
# tt_comp = tt_comp.flatten()

# # 補完データフレーム
# df_completion = pd.DataFrame({
#     'x': xx_comp,
#     'y': yy_comp,
#     'tw': tt_comp
# })