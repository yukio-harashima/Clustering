#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# ffi_MMcl4やPMcl4で使うassignsnap2KDE.datを形成する　  # 上記のプログラムの内部サブルーチンとして利用
import pandas as pd 
import numpy as np

#/////////////////////////////////////////////////////////////////////////////////
# データ読み取り
wd = "./" 
# wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241207-002407_cp/filt-adv25_var-nv6ONLY_pca80_stdOn_distCos_thrSil/" 
datafile = "mrf.dat"
ps = wd+datafile
mrf = pd.read_table(ps, sep = r'\s+', header=None)
                # 1     2           3    4    5    6    7    8
                # time, momentrate, Mrr, Mss, Mee, Mrs, Mre, Mse
mrf.columns = ['time', 'momentrate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']

datafile = "snap2.dat"
ps = wd+datafile
snap2 = pd.read_table(ps, sep = r'\s+', header=None)
                # 1  2  3   4   5   6         7    8    9
                # n, m, tw, dx, dy, sliprate, lat, lon, depth
snap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

datafile = "assignsnap2KDE.dat"
ps = wd+datafile
KDE = pd.read_table(ps, sep = r'\s+', header=None)
                # 1  2  3   4   5   6         7    8    9      10
                # n, m, tw, dx, dy, sliprate, lat, lon, depth, cluster
KDE.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth','cluster']

datafile = "centroids.dat"
ps = wd+datafile
centroids = pd.read_table(ps, sep='\t',header=None)
                    # 1     2      3     4    5         6    7    8    9    10   11
                    # tw,   clno,  dx,   dy, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse 
centroids.columns = ['tw', 'clno', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']

datafile = "centroids2.dat"
ps = wd+datafile
centroids2 = pd.read_table(ps, sep='\t',header=None)
                    # 1   2     3    4    5         6    7    8    9    10   11 12   13    14    15    16    17     18     19   20   21     22      23      24      25     
                    # tw, clno, lon, lat, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, trendp, trendt, trendb, plungp, plungt, plungb, NDC
centroids2.columns = ['tw', 'clno', 'lon','lat', 'sliprate','Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
                                             'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2', 
                                             'trendp', 'trendt', 'trendb', 'plungp','plungt', 'plungb', 'NDC']

# dataeditor.py実行後に生成されるclusteringSnap2c.datを使用
datafile = "clusteringSnap2c.dat"
ps = wd+datafile
clusteringSnap2 = pd.read_table(ps, sep='\t',header=None)
                            # 1   2     3    4    5         6    7    8    9    10   11 12   13    14    15    16    17     18     19   20   21     22      23      24      25     
                            # tw, clno, lon, lat, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, trendp, trendt, trendb, plungp, plungt, plungb, NDC
clusteringSnap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                            'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                        'trendp', 'trendt', 'trendb', 'plungp', 
                                            'plungt', 'plungb', 'NDC', 'clla', 'no']


#/////////////////////////////////////////////////////////////////////////////////

tmax = mrf['time'].max()
s2max = snap2['tw'].max()
s2ymax = KDE['tw'].max()
s2dtr = tmax/s2max
hs2 = s2dtr/2
s2ydtr = tmax/s2ymax
hs2y = s2ydtr/2

s2 = np.arange(hs2, tmax, s2dtr)
s2y = np.arange(hs2y, tmax, s2ydtr)

s2ar = np.arange(1, s2max+1, 1)
s2yar = np.arange(1, s2ymax+1, 1)

s2f = pd.DataFrame({'time':s2, 'tw':s2ar})
s2yf = pd.DataFrame({'time':s2y, 'tw':s2yar})

# ////////////////////////////////////////////////////////////////////////////////
# 最も近い時間を探す処理
def find_nearest_rows(s2f, s2yf):
    nearest_rows = []
    for t in s2f['time']:
        # 時間差の絶対値を計算
        diffs = np.abs(s2yf['time'] - t)
        # 最小差分のインデックスを取得
        nearest_idx = diffs.idxmin()
        # 対応する df2 の行を追加
        nearest_rows.append(s2yf.loc[nearest_idx])
    # 新しいデータフレームとして返す
    return pd.DataFrame(nearest_rows).reset_index(drop=True)

# df2 から最も近いデータを取得
nearest_s2 = find_nearest_rows(s2f, s2yf)



# ///////////////////////////////
# assignsnap2KDEの処理
# df3 から value2 に基づいてデータを取得
dfa = pd.merge(nearest_s2[['tw']], KDE, on='tw', how='left')
# 列の並び替え
dfa = dfa[KDE.columns]

# グループごとに連続する数値に置き換え
# 'value2' のユニークな値に対して連続する番号を割り振る
mapping = {v: i + 1 for i, v in enumerate(dfa['tw'].unique())}
dfa['tw'] = dfa['tw'].map(mapping)

# ///////////////////////////////
# centroidsの処理
dfa2 = pd.merge(nearest_s2[['tw']], centroids, on='tw', how='left')
# 列の並び替え
dfa2 = dfa2[centroids.columns]
# グループごとに連続する数値に置き換え
# 'value2' のユニークな値に対して連続する番号を割り振る
mapping = {v: i + 1 for i, v in enumerate(dfa2['tw'].unique())}
dfa2['tw'] = dfa2['tw'].map(mapping)

# ///////////////////////////////
# centroids2の処理
dfa3 = pd.merge(nearest_s2[['tw']], centroids2, on='tw', how='left')
# 列の並び替え
dfa3 = dfa3[centroids2.columns]
# グループごとに連続する数値に置き換え
# 'value2' のユニークな値に対して連続する番号を割り振る
mapping = {v: i + 1 for i, v in enumerate(dfa3['tw'].unique())}
dfa3['tw'] = dfa3['tw'].map(mapping)

# ///////////////////////////////
# clusteringSnap2の処理
dfa4 = pd.merge(nearest_s2[['tw']], clusteringSnap2, on='tw', how='left')
# 列の並び替え
dfa4 = dfa4[clusteringSnap2.columns]
# グループごとに連続する数値に置き換え
# 'value2' のユニークな値に対して連続する番号を割り振る
mapping = {v: i + 1 for i, v in enumerate(dfa4['tw'].unique())}
dfa4['tw'] = dfa4['tw'].map(mapping)



# ////////////////////////////////////////////////////////////////////////////////
# 結果のデータを出力
od = "./" 
# od = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241218-031730/" 
# od = wd

output_file = "assignsnap2KDE_v2.dat"
pss = od + output_file
dfa.to_csv(pss, sep='\t', index=False, header=False)

output_file = "centroids_v2.dat"
pss = od + output_file
dfa2.to_csv(pss, sep='\t', index=False, header=False)

output_file = "centroids2_v2.dat"
pss = od + output_file
dfa3.to_csv(pss, sep='\t', index=False, header=False)

output_file = "clusteringSnap2_v2.dat"
pss = od + output_file
dfa4.to_csv(pss, sep='\t', index=False, header=False)


print('data swapping done.')