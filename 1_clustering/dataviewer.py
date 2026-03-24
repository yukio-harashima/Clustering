#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# vscode上のインタラクティブウィンドウで各種データを確認する用

import pandas as pd

#/////////////////////////////////////////////////////////////////////////////////
# データ読み取り
wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241207-002407_cp/" 

# datafile = "snap.dat"
# ps = wd+datafile
# snap1 = pd.read_table(ps, sep = r'\s+', header=None)

# datafile = "snap2.dat"
# ps = wd+datafile
# snap2 = pd.read_table(ps, sep = r'\s+', header=None)

datafile = "snap_y.dat"
ps = wd+datafile
snap3 = pd.read_table(ps, sep = r'\s+', header=None)

# datafile = "snap2_y.dat"
# ps = wd+datafile
# snap4 = pd.read_table(ps, sep = r'\s+', header=None)

datafile = "n_vector.dat"
ps = wd+datafile
snap5 = pd.read_table(ps, sep = r'\s+', header=None)

# datafile = "centroids2_v2.dat"
# ps = wd+datafile
# KDE = pd.read_table(ps, sep='\t',header=None)



# # dat = pd.DataFrame(data)
# # df = pd.DataFrame(data)
# # # 1  2  3   4   5   6         7    8    9    10   11   12   13    14    15    16    17     18     19   20   21     22      23      24      25      26      27      28
# # # n, m, tw, dx, dy, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, lat, lon, depth, trendp, trendt, trendb, plungp, plungt, plungb, NDC
# snap1.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
#                'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
#                     'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
#                         'trendp', 'trendt', 'trendb', 'plungp', 
#                             'plungt', 'plungb', 'NDC']



# snap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

snap3.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                   'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                        'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                            'trendp', 'trendt', 'trendb', 'plungp', 
                                'plungt', 'plungb', 'NDC']

# snap4.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate','lat', 'lon', 'depth']

snap5.columns = ['n', 'm', 'tw','t', 'x', 'y',
                    'n1_(1)','n1_(2)','n1_(3)',
                         'n2_(1)','n2_(2)','n2_(3)',
                            'sliprate']

# # s = ['n', 'm', 'tw']
# # diff = pd.DataFrame(snap3[s]-snap4[s])
# diff = pd.DataFrame(snap2-snap4)

#/////////////////////////////////////////////////////////////////////////////////

# # 結果のデータを再出力
# od = "./" 
# output_file = "clusteringSnap.dat"
# pss1 = wd + output_file
# # df.to_csv(pss, sep=' ', index=False, header=False)

# data1 = pd.read_table(pss1,  sep=' ',  header=None)

# # 結果のデータを再出力 (非スムーシング)
# od = "./" 
# output_file = "clusteringSnap2c.dat"
# pss2 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data2 = pd.read_table(pss2, sep='\t', header=None)

# output_file = "clusteringSnap2_v2.dat"
# pss3 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data3 = pd.read_table(pss3, sep='\t', header=None)




# output_file = "centroids4.dat"
# pss4 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data4 = pd.read_table(pss4, sep='\t', header=None)
# # datafile = "mrf.dat"
# # ps = wd+datafile
# # mrf = pd.read_table(ps, sep = r'\s+', header=None)
# # # 1     2           3    4    5    6    7    8
# # # time, momentrate, Mrr, Mss, Mee, Mrs, Mre, Mse

# output_file = "clusteringSnap2c.dat"
# pss5 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data5 = pd.read_table(pss5, sep='\t', header=None)
# data5.columns = [ # 出力ファイルの列構造
#         'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
#         'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
#         'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
#         'trendp', 'trendt', 'trendb', 'plungp', 
#         'plungt', 'plungb', 'NDC', 'clla', 'no'
#     ]


# output_file = "clmeca_fault.dat"
# pss6 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data6 = pd.read_table(pss6, sep='\t', header=0)




# output_file = "clarea.dat"
# pss7 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data7 = pd.read_table(pss7, sep='\t', header=None)
# data7.columns = [ # 出力ファイルの列構造
#         'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
#         'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
#         'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
#         'trendp', 'trendt', 'trendb', 'plungp', 
#         'plungt', 'plungb', 'NDC', 'clla', 'no'
#     ]




# output_file = "clave_mod.dat"
# pss8 = wd + output_file
# # dat.to_csv(pss, sep='\t', index=False, header=False)
# data8 = pd.read_table(pss8, sep='\t', header=0)