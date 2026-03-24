#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# クラスターごとの要素点の重心を求める
import numpy as np
import pandas as pd
from scipy.spatial import distance

# データ読み取り
od = "./" 
# od = "/Users/harashima-yukio/Desktop/results_20240702-182819/" 
# od = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241207-002407_cp/" 
# od = "/Users/harashima-yukio/Desktop/2024_Noto/results_20240702-182819/" 
output_file = "clusteringSnap2.dat"
pss2 = od + output_file
# dat.to_csv(pss, sep='\t', index=False, header=False)
data = pd.read_table(pss2, sep='\t', header=None)
data.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC', 'cllab', 'clno']


datafile = "snap_yr.dat"
ps = od+datafile
snp_yr = pd.read_table(ps, sep = r'\s+', header=None)
snp_yr.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC']


#/////////////////////////////////////////////////////////////////

# dataの['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']をデータフレームsnap_yrの対応する列で置き換え
data2 = data.copy()
data2[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']] = snp_yr[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']]
                                                           
#/////////////////////////////////////////////////////////////////

# クラスターごとの重心を計算
# 各時刻 (tw) とクラスタ番号 (clno) ごとに、dx, dy, および sliprate の平均値を計算します。
centroids = data2.groupby(['tw', 'clno'])[['dx', 'dy', 'sliprate']].mean().reset_index()
# 各時刻 (tw) とクラスタ番号 (clno) ごとに、モーメントテンソルの要素 (Mrr, Mss, Mee, Mrs, Mre, Mse) を合計します。
moments = data2.groupby(['tw', 'clno'])[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].sum().reset_index()
# moments データフレームから、tw と clno 列を削除します。
mo=moments.drop(moments.columns[[0, 1]], axis=1)
# centroids データフレームと mo を列方向 (axis=1) に結合して、新しいデータフレームを作成します。
df=result = pd.concat([centroids, mo], axis=1)

# 各時刻 (tw) とクラスタ番号 (clno) ごとに、dx, dy, sliprate に加え、モーメントテンソルの各要素の平均値を計算します。
centroids1 = data2.groupby(['tw', 'clno'])[['dx', 'dy', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].mean().reset_index()
# moments = data.groupby(['tw', 'clno'])[[]].sum().reset_index()
# mo=moments.drop(moments.columns[[0, 1]], axis=1)
# df=result = pd.concat([centroids, mo], axis=1)
# クラスタ番号 (clno) が正の値のみを持つ行を抽出します。
centroids2 = centroids1[centroids1['clno'] > 0]

# 各時刻 (tw) とクラスタ番号 (clno) ごとに、経度 (lon), 緯度 (lat), および sliprate の平均値を計算します。
centroidsl = data.groupby(['tw', 'clno'])[['lon','lat', 'sliprate']].mean().reset_index()
# モーメントテンソルの要素の合計を計算。moments と同様の処理ですが、centroidsl 用。
momentsl = data.groupby(['tw', 'clno'])[['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']].sum().reset_index()
# モーメントテンソル以外の列 (tw と clno) を削除します。
mol=momentsl.drop(momentsl.columns[[0, 1]], axis=1)
# 経度・緯度重心 (centroidsl) とモーメントテンソル合計 (mol) を結合して結果を作成します。
dfl=result = pd.concat([centroidsl, mol], axis=1)

# 各時刻 (tw) とクラスタ番号 (clno) ごとに、経度・緯度情報 (lon, lat), sliprate およびモーメントテンソルの平均を計算します。
centroids1l = data.groupby(['tw', 'clno'])[['lon','lat', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
                                             'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 
                                             'trendp', 'trendt', 'trendb', 'plungp','plungt', 'plungb', 'NDC']].mean().reset_index()
# moments = data.groupby(['tw', 'clno'])[[]].sum().reset_index()
# mo=moments.drop(moments.columns[[0, 1]], axis=1)
# df=result = pd.concat([centroids, mo], axis=1)
# クラスタ番号 (clno) が正の値のみを持つ行を抽出します。
centroids2l = centroids1l[centroids1['clno'] > 0]

# クラスタ番号 (clno) ごとに、経度・緯度情報 (lon, lat), sliprate およびモーメントテンソルの平均を計算します。
centroids3l = data.groupby(['clno'])[['tw','lon','lat', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
                                             'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 
                                             'trendp', 'trendt', 'trendb', 'plungp','plungt', 'plungb', 'NDC']].mean().reset_index()
# moments = data.groupby(['tw', 'clno'])[[]].sum().reset_index()
# mo=moments.drop(moments.columns[[0, 1]], axis=1)
# df=result = pd.concat([centroids, mo], axis=1)
# クラスタ番号 (clno) が正の値のみを持つ行を抽出します。
centroids4l = centroids3l[centroids3l['clno'] > 0]

# 結果のデータを出力
od = "./" 
output_file = "clusteringSnap3.dat"
pss = od + output_file
df.to_csv(pss, sep='\t', index=False, header=False)

od = "./" 
output_file = "centroids.dat"
pss = od + output_file
centroids2.to_csv(pss, sep='\t', index=False, header=False)

od = "./" 
output_file = "clusteringSnap3l.dat"
pss = od + output_file
dfl.to_csv(pss, sep='\t', index=False, header=False)

od = "./" 
output_file = "centroids2.dat"
pss = od + output_file
# tw,    clno,      lon,      lat,  sliprate,    Mrr,     Mss,    Mee,    Mrs,     Mre,    Mse 
centroids2l.to_csv(pss, sep='\t', index=False, header=False)


od = "./" 
output_file = "centroids4.dat"
pss = od + output_file
# tw,    clno,      lon,      lat,  sliprate,    Mrr,     Mss,    Mee,    Mrs,     Mre,    Mse 
centroids4l.to_csv(pss, sep='\t', index=False)



#/////////////////////////////////////////////////////////////////
# # 各重心に最も近い格子点を求める
# def find_nearest_grid_point(centroid, points):
#     # 重心と全格子点の距離を計算
#     dist = distance.cdist([centroid], points, metric='euclidean')
#     # 最も距離が近い点のインデックスを取得
#     nearest_index = np.argmin(dist)
#     return points[nearest_index]

# nearest_points = []

# for _, centroid in centroids.iterrows():
#     # 重心の座標を取り出す
#     centroid_coords = [centroid['dx'], centroid['dy']]
#     # 重心に最も近い格子点を求める
#     nearest_point = find_nearest_grid_point(centroid_coords, points)
#     nearest_points.append(nearest_point)


# #/////////////////////////////


# # 結果を表示
# for cluster, point in zip(centroids['cluster'], nearest_points):
#     print(f"クラスター {cluster} の重心に最も近い格子点: {point}")