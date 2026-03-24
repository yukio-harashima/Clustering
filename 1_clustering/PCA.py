#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# PDTIによって出力されたデータsnap.datを使ってクラスタリングを行う
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn
import time
import plotly.io as pio
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial.distance import pdist
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tqdm import tqdm



##/////////////////////////////////////////////////////////////////
        #データの読み込み
#/////////////////////////////////////////////////////////////////
wd = "./" 
datafile = "snap_y.dat"
ps = wd+datafile
data = pd.read_table(ps, sep = r'\s+', header=None)
dat = pd.DataFrame(data)
df = pd.DataFrame(data)
# 1  2  3   4   5   6         7    8    9    10   11   12   13    14    15    16    17     18     19   20   21     22      23      24      25      26      27      28
# n, m, tw, dx, dy, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, lat, lon, depth, trendp, trendt, trendb, plungp, plungt, plungb, NDC
df.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC']

datafile = "n_vector.dat"
ps = wd+datafile
snap4 = pd.read_table(ps, sep = r'\s+', header=None)
# 1  2  3  4  5  6  7       8       9       10      11      12      13
# n, m, tw,t, x, y, n1_(1), n1_(2), n1_(3), n2_(1), n2_(2), n2_(3), potency-rate
snap4.columns = ['n', 'm', 'tw','t', 'x', 'y',
                    'n1_(1)','n1_(2)','n1_(3)',
                         'n2_(1)','n2_(2)','n2_(3)',
                            'sliprate']

#/////////////////////////////////////////////////////////////////
        #データ成形
#/////////////////////////////////////////////////////////////////  スイッチポイント A
# #sanp.datを使う場合のデータ下処理(モーメントのスムージング)
# dfrow = df.shape[0]
# smooth = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
# row1 = ['Mrr', 'Mss', 'Mee']
# row2 = ['Mrs', 'Mre', 'Mse']
# r2 = 1 / np.sqrt(2)

# # M0を計算
# def rms(row):
#     return np.mean(np.square(row))

# M0 = pd.DataFrame()
# M0['Row1'] = df[row1].apply(rms, axis=1)   #3x3の対角成分の二乗和
# M0['Row2'] = 2*df[row2].apply(rms, axis=1) #3x3のi,j成分の二乗和
# M0['M0'] = r2*np.sqrt(M0['Row1']+M0['Row2']) # M0を計算

# # # デバッグ用計算
# # SMM = pd.DataFrame()
# # SMM['Mrr'] = df['Mrr'] / M0['M0']
# # print(df.iloc[10,6])
# # print(M0.iloc[10,4])
# # a = df.iloc[10,6] / M0.iloc[10,4]

# # 各モーメント成分に対してM0で割る
# for col in smooth:
#     df[col] = df[col] / M0['M0']
#/////////////////////////////////////////////////////////////////


# sliprate 上位n％以上(４分位数)
q = input('input sliprate quantile: ')
print('your input quantile is '+ q)
qq = float(q)




# ////////////////////    snap.datのsliprateの絞り込み    ///////////////////////  スイッチポイントA
# quant= df['sliprate'].quantile(qq) 
# qu = float(quant)
# ddd = df[df['sliprate'] > qu]

# ////////////////////    n_vector.datのsliprateの絞り込み    ///////////////////////  スイッチポイントB
qu = snap4['sliprate'].quantile(qq)   # 四分位数
# qu = snap4['sliprate'].max() * qq       # sliprate の最大値からの閾値

qu = float(qu)
ddd = snap4[snap4['sliprate'] > qu]

#/////////////////////////////////////////////////////////////////
        #クラスタリング
#/////////////////////////////////////////////////////////////////

# クラスタリングに使用する列を指定
# columns_to_cluster = [ 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']    #スイッチポイント A
columns_to_cluster = ['tw', 'x', 'y','n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)']     #スイッチポイントB
columns_to_cluster = [ 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']    #スイッチポイント c


# クラスタリング 
# # //////////////////    #コサイン類似度から距離を求める #ウォード法（最小分散法）   ///////////////////////  スイッチポイント
# print('The distance calculation used is "cosine".')
# cldf = linkage(pdist(ddd[columns_to_cluster], metric='cosine'), method='ward')
# # cldf = linkage(pdist(snap4[columns_to_cluster], metric='cosine'), method='ward')

# 指定した列のデータを標準化
scaler = StandardScaler()
# ddd_scaled = pd.DataFrame(scaler.fit_transform(ddd[columns_to_cluster]), columns=columns_to_cluster)
ddd[columns_to_cluster] = scaler.fit_transform(ddd[columns_to_cluster])
ddd_scaled = ddd[columns_to_cluster]
# 平均と分散の計算
nrm_dat = pd.DataFrame({
    'Mean': ddd_scaled.mean(),
    'Variance': ddd_scaled.var()
})

# インデックス名を列名として追加
nrm_dat.reset_index(inplace=True)
nrm_dat.rename(columns={'index': 'Column'}, inplace=True)

print(nrm_dat)

# コサイン類似度から距離を求め、ウォード法でクラスタリング
print('The distance calculation used is "cosine".')
cldf = linkage(pdist(ddd[columns_to_cluster], metric='cosine'), method='ward')


# # ////////////////////   #ユークリッド距離を求める #ウォード法（最小分散法）   ///////////////////////   スイッチポイント
# print('The distance calculation used is "ward".')
# cldf = linkage(df[columns_to_cluster], metric='euclidean', method='ward')

# 樹形図作成
#図の作成
fig = plt.figure(figsize=(16, 10))
# 樹形図作成
dendrogram(cldf,color_threshold=0.84)
# fig.update_layout(width=1400, height=900)
#図の保存
# plt.savefig("dendrogram.pdf") 
# 樹形図表示
plt.show()

#散布図を描画するクラスター番号の決定
n = input('input distance: ')
print('your input distance is '+ n)


#///////////////////////////////////////////////////////////////////////
# ct = input('input distance: ')
# ct = int(ct)
# print('your input distance is '+ ct)

# dendrogram(cldf,color_threshold=ct)

# #図の保存
# plt.savefig("dendrogram.pdf") 
# # 樹形図表示
# plt.show()

#////////////////////////////////////

#データ出力
cl = pd.DataFrame(cldf,columns=['ind1', 'ind2', 'dist', 'NoD'])
row=ddd.shape[0]
cl['clNo'] = range(row, len(cl.index) + row)
cl.to_csv("index.csv", index = False)
print(cl)
