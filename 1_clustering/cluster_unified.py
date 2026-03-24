#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# PDTIによって出力されたデータsnap.datを使ってクラスタリングを行う

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage,dendrogram
from matplotlib import pyplot as plt
import plotly.io as pio
import plotly.express as px
from scipy.spatial.distance import pdist
import plotly.figure_factory as ff
import plotly.graph_objects as go
from tqdm import tqdm
import time
from sklearn.preprocessing import StandardScaler


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
columns_to_cluster = ['n', 'm', 'tw','n1_(1)','n1_(2)','n1_(3)', 'n2_(1)','n2_(2)','n2_(3)', 'sliprate']     #スイッチポイントB

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

# # コサイン類似度から距離を求め、ウォード法でクラスタリング
print('The distance calculation used is "cosine".')
cldf = linkage(pdist(ddd_scaled, metric='cosine'), method='ward')


# # ////////////////////   #ユークリッド距離を求める #ウォード法（最小分散法）   ///////////////////////   スイッチポイント
# print('The distance calculation used is "ward".')
# cldf = linkage(df[columns_to_cluster], metric='euclidean', method='ward')
# cldf = linkage(pdist(ddd[columns_to_cluster], metric='euclidean'), method='ward')

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


#dummy/assign
#/////////////////////////////////////////////////////////////////
# # 繰り返し処理
dfrow = ddd.shape[0] #dfの行数
clmax = cl['clNo'].max() 
clset = [] #全てのデータ対応のクラスター番号記録

# dfrow = 2
# i = 2
# for i in range(dfrow):
for i in tqdm(range(dfrow), desc="assign cluster NO", ncols=100, ascii=True):
    clsb = [] #i行のクラスター番号の記録
    # print(i)
    abcd = cl[(cl['ind1'] == i )| (cl['ind2'] == i)]  #iの含まれる行を検索
    d = abcd.iat[0,4] #クラスター番号の抽出
    clsb.append(d) #クラスター番号の記録
    
    while d != clmax: #dがクラスター番号の最大値でない時
        # print(d)
        abcd = cl[(cl['ind1'] == d )| (cl['ind2'] == d)]  #dの含まれる行を検索
        d = abcd.iat[0,4] #クラスター番号の抽出
        clsb.append(d) #クラスター番号の記録
    clset.append(clsb)
    # print(clsb)
    time.sleep(0.01)

ddd['clset'] = clset #元データにクラスター番号のセットを追加
ddd['clset'] = ddd['clset'].astype('str')
df['clset'] = ddd['clset']
df['clset'] = df['clset'].fillna("[0]")

dat['clset'] = ddd['clset']
dat['clset'] = dat['clset'].astype('str')
#/////////////////////////////////////////////////////////////////

#データ出力
ass = pd.DataFrame (df)
ass.to_csv("assign.csv", index = False)
print(ass)


#dummy/scattergram
#/////////////////////////////////////////////////////////////////
# #散布図を描画するクラスター番号の決定
# n = input('input distance: ')
# print('your input distance is '+ n)
nn = float(n)
abcd = cl[(cl['dist'] >= nn )] #クラスター距離から任意の距離以上のクラスターを検索
dm = []
dm.extend(abcd['clNo']) #上のクラスター番号の取得

dn = []
dn.extend(abcd['ind1']) #上のクラスターを構成するクラスター番号を取得
dn.extend(abcd['ind2']) #同上

lab = list(set(dn) - set(dm)) #dnとdmクラスター番号の差集合を取得
lal = len(lab)
la = list(range(lal)) #整数番号割り付け用の連番生成
lab.sort() #クラスター番号の小さい順にソート

#/////////////////////////////////////////////////////////////////
#散布図用のデータにクラスター番号のラベル割り付け
n = len(lab) #差集合の要素数を取得
collab = []
ass['label'] = ''  #からの列の生成

for i in range(n):
    # print (lab[i])
    setNo = int(lab[i]) #整数型に変換
    sen = int(la[i]+1)
    ll = str(setNo) #文字列に変換
    asset = ass[ass['clset'].str.contains(ll)] #i(文字)の含まれる行を検索
    ass.loc[asset.index, 'label'] = sen #整数番号割り付け
    # print (len(asset))

dat['label'] = ass['label'] #整数番号割り付け
dat['label'].replace('', 0, inplace=True) #sliprate下位n%にタグ付け（空欄項目を文字列置換）
# dat['label'] = dat['label'].astype('int')
ass['label'] = dat['label']

# i = 0
# setNo = int(lab[i]) #整数型に変換
# ll = str(setNo) #文字列に変換

# asset = ass[ass['clset'].str.contains(ll)]#i(文字)の含まれる行を検索
# ass.loc[asset.index, 'label'] = setNo
# print (len(asset))


#/////////////////////////////////////////////////////////////////
# 散布図を描画　2次元
x = ass.iloc[:,18]  # x軸のデータ
y = ass.iloc[:,19]  # y軸のデータ
z = ass['label']   # 点の色をデータから設定
cm = plt.get_cmap('rainbow')  # カラーマップを設定
# plt.scatter( x ,y, c=z)
# plt.scatter( x ,y, cmap=cm)
# plt.legend(loc='upper right')
# plt.show()

#/////////////////////////////////
# # 散布図を描画　３次元
# x = ass.iloc[:,0]  # x軸のデータ
# y = ass.iloc[:,1]  # y軸のデータ
# z = ass.iloc[:,2]  # z軸のデータ
# value = ass['label'] # 点(x, y, z)がもつ量(クラスター番号)
# cm = plt.get_cmap('rainbow')# カラーマップを生成

# fig = plt.figure()# figureを生成する
# ax = fig.add_subplot(1, 1, 1, projection='3d')# axをfigureに設定する
 
# # axに散布図を描画、戻り値にPathCollectionを得る
# mappable = ax.scatter(x, y, z, c=value, cmap=cm)
# fig.colorbar(mappable, ax=ax)
 
# # 表示する
# plt.show()

#/////////////////////////////////
# 散布図を描画　３次元(インタラクティブ)
x = ass.iloc[:, 0]
y = ass.iloc[:, 1]
z = ass.iloc[:, 2]
value = ass['label']

# fig = px.scatter_3d(ass, x='n', y='m', z='tw', color='label')
# fig.update_traces(marker_size=3)
# pio.write_html(fig, file='scatter_3d.html', auto_open=True)

#/////////////////////////////////
# インタラクティブなデンドログラムを作成
dendro = ff.create_dendrogram(data, linkagefun=lambda x: cldf)

# 分岐点の座標を取得する
def get_dendrogram_branch_points(dendro_fig):
    branch_points = []
    for trace in dendro_fig.data:
        if 'line' in trace.mode:
            for i in range(len(trace.x)):
                branch_points.append((trace.x[i], trace.y[i]))
    return branch_points

branch_points = get_dendrogram_branch_points(dendro)
# デンドログラムをプロット
fig = dendro

# 分岐点を追加
x_points, y_points = zip(*branch_points)
points = pd.DataFrame()
points['x_points'] = x_points
points['y_points'] = y_points

# データのフィルタリングと整形
#/////////////////////////////////
# 1. y座標の値が0より上であること。
pt_filtered = points[points['y_points'] > 0]

# 2. y座標の値が同じ点では、x座標の値が中央値を求める
pt_filtered = pt_filtered.groupby('y_points')['x_points'].median().reset_index()

# クラスター番号の生成
xlen = row
pt_filtered = pt_filtered.sort_values(['y_points'])
pxlen = int(len(pt_filtered['x_points']))
pt_filtered['label'] = list(range(xlen,pxlen+xlen))

# planeMap等に表示されるクラスター番号のみ置き換える。
dictio = pd.DataFrame()
new = list(range(1,lal+1))
dictio['old'] = lab
dictio['new'] = new
replace_dict = dictio.set_index('old')['new'].to_dict()

# pt_filteredの列['label']の値をdictの['new']に置き換え
pt_filtered['label'] = pt_filtered['label'].replace(replace_dict)

#/////////////////////////////////
def determine_color(label):
    if label in new:
        return 'red'
    else:
        return 'black'

# pt_filteredのラベルに基づいて色を設定
pt_filtered['color'] = pt_filtered['label'].apply(determine_color)

# 分岐点のプロット
fig.add_trace(go.Scatter(x=pt_filtered['x_points'], y=pt_filtered['y_points'], mode = 'markers+text', 
                         marker=dict(color=pt_filtered['color']),
                         text =  pt_filtered['label'], textposition="top center",
                         textfont=dict(color=pt_filtered['color'])))

# fig.add_trace(go.Scatter(x=x_points, y=y_points, mode = 'markers', 
#                          marker=dict(color='black')))

fig.update_layout(width=1400, height=900)
# fig.update_layout(autosize=True) 

# グラフを表示
# fig.show()
# グラフを保存
pio.write_html(fig, file='dendrogram.html')


# 樹形図作成
#図の作成
fig = plt.figure()
# 樹形図作成
dendrogram(cldf,color_threshold=nn)
# fig.update_layout(width=1400, height=900)
#図の保存
plt.savefig("dendrogram.pdf") 
# 樹形図表示
# plt.show()

#////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# 結果のデータを再出力
od = "./" 
output_file = "clusteringSnap.dat"
pss = od + output_file
df.to_csv(pss, sep=' ', index=False, header=False)

# 結果のデータを再出力 (非スムーシング)
od = "./" 
output_file = "clusteringSnap2.dat"
pss = od + output_file
dat.to_csv(pss, sep='\t', index=False, header=False)


# 結果の情報
od = "./" 
output_file = "clusteringinfo.dat"
pss = od + output_file
dat['label'].to_csv(pss, sep='\t', index=False, header=False)