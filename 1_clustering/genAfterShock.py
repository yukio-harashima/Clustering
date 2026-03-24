#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
import numpy as np
import pandas as pd

# データ読み取り
# id="/Users/harashima-yukio/Desktop/"

id = "./" 
input_file = "地震リスト (2).csv"
pss1 = id + input_file
# dat.to_csv(pss, sep='\t', index=False, header=False)
data = pd.read_csv(pss1)

depth = pd.DataFrame(data['深さ'])
depth [['depth', 'km']]= depth['深さ'].str.split(' ', expand=True)
data = pd.concat(
    [data.iloc[:, :6], depth['depth'], data.iloc[:, 6:]],
    axis=1
)


lat= pd.DataFrame(data['緯度'])
# 正規表現を使って度、分、方位を抽出
lat[['latdeg', 'latmin', 'latdir']]= lat['緯度'].str.extract(r'(\d+)°(\d+\.\d+)′([NS])')
# 抽出結果を数値型に変換
lat['latdeg'] = lat['latdeg'].astype(float)  
lat['latmin'] = lat['latmin'].astype(float)
# 方位に基づいて度を正負に変換
lat['latdeg'] = lat.apply(lambda row: -row['latdeg'] if row['latdir'] == 'S' else row['latdeg'], axis=1)
# 度と分を統合して小数点形式に変換
lat['decimal_latdeg'] = lat['latdeg'] + lat['latmin'] / 60 
lat['decimal_latdeg'] = lat['decimal_latdeg'].astype(float)
# 緯度データのインデックスをリセット
lat = lat.reset_index(drop=True)
# データフレームに挿入
data = pd.concat(
    [data.loc[:, :"緯度"], lat[['latdeg', 'latmin', 'decimal_latdeg']], data.loc[:, "緯度":]],
    axis=1
)



lon= pd.DataFrame(data['経度'])
# 正規表現を使って度、分、方位を抽出
lon[['londeg', 'lonmin', 'londir']]= lon['経度'].str.extract(r'(\d+)°(\d+\.\d+)′([EW])')
# 抽出結果を数値型に変換
lon['londeg'] = lon['londeg'].astype(float)  
lon['lonmin'] = lon['lonmin'].astype(float)
# 方位に基づいて度を正負に変換
lon['londeg'] = lon.apply(lambda row: -row['londeg'] if row['londir'] == 'W' else row['londeg'], axis=1)
# 度と分を統合して小数点形式に変換
lon['decimal_londeg'] = lon['londeg'] + lon['lonmin'] / 60 
lon['decimal_londeg'] = lon['decimal_londeg'].astype(float)
# 経度データのインデックスをリセット
lon = lon.reset_index(drop=True)
# データフレームに挿入
data = pd.concat(
    [data.loc[:, :"経度"], lon[['londeg', 'lonmin', 'decimal_londeg']], data.loc[:,"経度":]],
    axis=1
)


data=data.drop(['震央地名', '緯度','経度', '深さ', '最大震度'], axis=1)

# data.columns = ['ymd','time', 'lontime', 'londeg', 'lon10', 'lattime', 'latdeg', 'lat10', 'depth', 'M']


while True:  # 無限ループで入力待機を繰り返す
    # 入力の待機
    n = input('Input version or location (leave empty for default): ').strip()

    # 入力がある場合の処理
    if n:
        print(f'Filename is aftershock_{n}.txt')
        output_file = f"aftershock_{n}.txt"
        break  # 入力がある場合はループを終了
    else:
        # 入力がない場合の処理
        confirm = input('WARNING!! You did not define a version. Overwrite the old file? (y/n): ').strip().lower()
        if confirm == 'y':
            output_file = "aftershock.txt"
            break  # ユーザーがファイル作成を許可した場合はループを終了
        elif confirm == 'n':
            print("Please provide a version or location.")
            continue  # 再度バージョン入力の待機に戻る
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
            continue  # 無効な入力の場合は再試行



# ファイルパスの生成
od = "/Users/harashima-yukio/modelplane/"
pss = od + output_file

# データを保存
try:
    data.to_csv(pss, sep='\t', index=False, header=False)
    print(f"File successfully saved as: {pss}")
except Exception as e:
    print(f"An error occurred while saving the file: {e}")


# n = input('input version or location: ')
# print('filename is aftershock_'+ n)

# od = "./" 
# output_file = "aftershock_"
# version = n
# file = ".txt"
# pss = od + output_file + version + file
# data.to_csv(pss, sep='\t', index=False, header=False)



# print('WARNING!! you not difine version. overwrite old file.')

# od = "./" 
# output_file = "aftershock.txt"
# pss = od + output_file
# data.to_csv(pss, sep='\t', index=False, header=False)