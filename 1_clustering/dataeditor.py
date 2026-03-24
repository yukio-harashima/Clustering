#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# vscode上のインタラクティブウィンドウで各種データを確認する用

import pandas as pd

#/////////////////////////////////////////////////////////////////////////////////
# データ読み取り
# wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20241207-002407_cp/filt-adv25_var-nv6ONLY_pca80_stdOn_distCos_thrSil/"\t
# wd = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20250730-100041/filt-arbP75_distCUSTw10-10-10_nClust8_20250927-220248/relabeled_nClust8/"
wd = "./"

datafile = "clusteringSnap2.dat"
ps = wd+datafile
clusteringSnap2 = pd.read_table(ps, sep='\t',header=None)
clusteringSnap2.columns = ['n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
                            'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
                                    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
                                        'trendp', 'trendt', 'trendb', 'plungp', 
                                            'plungt', 'plungb', 'NDC', 'clla', 'no']

# 数字と書き換え先の対応辞書を作成
mapping = {
0: 0,
1: 2,
2: 3,
3: 4,
4: 1,
5: 5,

}


# 'no'列の数字を対応する数字に書き換える
clusteringSnap2['no'] = clusteringSnap2['no'].map(mapping)

# データ出力
od = wd
output_file = "clusteringSnap2c.dat"
pss = od + output_file
clusteringSnap2.to_csv(pss, sep='\t', index=False, header=False)
