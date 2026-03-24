#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# vscode上のインタラクティブウィンドウで各種データを確認する用

import pandas as pd

#/////////////////////////////////////////////////////////////////////////////////
# データ読み取り
# wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241207-002407_cp/filt-adv25_var-nv6ONLY_pca80_stdOn_distCos_thrSil/"
# wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20250730-100041/filt-adv25_var-nv6NMTW_pcaN3_stdOn_distCos_thrSil/"
wd = "./"

datafile = "clusteringSnap2.dat"
ps = wd+datafile
clusteringSnap2 = pd.read_table(ps, sep='\t',header=None)
clusteringSnap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                            'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                        'trendp', 'trendt', 'trendb', 'plungp', 
                                            'plungt', 'plungb', 'NDC', 'clla', 'no']


cli = clusteringSnap2['no']

# データ出力
od = wd
output_file = "clusteringSnap2.dat"
pss = od + output_file
clusteringSnap2.to_csv(pss, sep='\t', index=False, header=False)

output_file = "clusteringinfo.dat"
pss = od + output_file
cli.to_csv(pss, sep='\t', index=False, header=False)
