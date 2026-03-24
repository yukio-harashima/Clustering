#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# 出力データの相違を調べるプログラム　スタンドアローンかつ絶対パスを指定すること

import pandas as pd


#///////////////////////////
# データ読み取り
wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241218-031730/" 
datafile = "snap2.dat"
ps = wd+datafile
df1 = pd.read_table(ps, sep = r'\s+', header=None)
df1.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

datafile = "snap2_y.dat"
ps = wd+datafile
df2 = pd.read_table(ps, sep = r'\s+', header=None)
df2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

# datafile = "mrf.dat"
# ps = wd+datafile
# mrf = pd.read_table(ps, sep = r'\s+', header=None)
#                 # 1     2           3    4    5    6    7    8
#                 # time, momentrate, Mrr, Mss, Mee, Mrs, Mre, Mse
# mrf.columns = ['time', 'momentrate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']

# datafile = "mrf_1.dat"
# ps = wd+datafile
# mrf_1 = pd.read_table(ps, sep = r'\s+', header=None)
#                 # 1     2           3    4    5    6    7    8
#                 # time, momentrate, Mrr, Mss, Mee, Mrs, Mre, Mse
# mrf_1.columns = ['time', 'momentrate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']

datafile = "snap_y.dat"
ps = wd+datafile
snp_y = pd.read_table(ps, sep = r'\s+', header=None)
snp_y.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC']

datafile = "snap_yr.dat"
ps = wd+datafile
snp_yr = pd.read_table(ps, sep = r'\s+', header=None)
snp_yr.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
               'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                        'trendp', 'trendt', 'trendb', 'plungp', 
                            'plungt', 'plungb', 'NDC']



#////////////////////    snap2.dat用     ///////////////////////
# 各行をタプルに変換
df1_tuples = df1.apply(tuple, axis=1)
df2_tuples = df2.apply(tuple, axis=1)

# インデックスを効率よく検索するための辞書を作成
df2_index_map = {row: idx for idx, row in enumerate(df2_tuples)}

# 一致したインデックスを保存するリスト
matches = [(idx1, df2_index_map[row]) for idx1, row in enumerate(df1_tuples) if row in df2_index_map]

# 結果をデータフレームに整形
result = pd.DataFrame(matches, columns=['df1_index', 'df2_index'])
    
# 出力
print("一致する行のインデックス:")
print(result)
#/////////////////////////////////////////////////////////////////


#////////////////////    mrf.dat用     ///////////////////////
# mrf_diff = mrf_1 - mrf
# mrfres = mrf_diff.sum()
# print(mrfres)

#/////////////////////////////////////////////////////////////////


#////////////////////    snap_y.dat用     ///////////////////////
snp_diff =snp_y  - snp_yr

print(snp_diff)

#/////////////////////////////////////////////////////////////////
