#!/usr/bin/env python3
# ----------------------------------------------------------------
# 1. ライブラリのインポート
# ----------------------------------------------------------------
import pandas as pd
import numpy as np
from geographiclib.geodesic import Geodesic
import os

# ----------------------------------------------------------------
# 2. グローバル設定
# ----------------------------------------------------------------

# 測地線計算のためのインスタンスを生成 (WGS84楕円体モデル)
geod = Geodesic.WGS84

# ----------------------------------------------------------------
# 3. 関数の定義
# ----------------------------------------------------------------

def get_coords_tsunami_fault_model(lat, lon, strike, dip, length, width, upper_depth, lower_depth):
    """
    断層パラメータから地表投影された断層の矩形座標と深さを計算する。
    
    Args:
        lat (float): 基準点の緯度
        lon (float): 基準点の経度
        strike (float): 走向 (度)
        dip (float): 傾斜角 (度)
        length (float): 長さ (km)
        width (float): 幅 (km)
        upper_depth (float): 上端深度 (km)
        lower_depth (float): 下端深度 (km)
        
    Returns:
        (np.array, np.array, np.array): 矩形の頂点のx座標(経度), y座標(緯度), z座標(深さ)
    """
    # 傾斜を考慮した見かけの幅を計算
    width_on_map = width * np.cos(np.deg2rad(dip))
    
    xs = np.zeros(5)
    ys = np.zeros(5)

    # 基準点(左上:LT)から右上(RT)へ
    tmp0 = geod.Direct(lat, lon, strike, length * 1e3) # 距離はメートルに変換
    xs[0], ys[0] = tmp0['lon1'], tmp0['lat1'] # LT
    xs[1], ys[1] = tmp0['lon2'], tmp0['lat2'] # RT
    xs[4], ys[4] = tmp0['lon1'], tmp0['lat1'] # 閉じた矩形にするため始点を終点にも設定

    # 右上(RT)から右下(RB)へ (走向から+90度の方向)
    tmp = geod.Direct(tmp0['lat2'], tmp0['lon2'], strike + 90, width_on_map * 1e3)
    xs[2], ys[2] = tmp['lon2'], tmp['lat2'] # RB

    # 基準点(左上:LT)から左下(LB)へ
    tmp1 = geod.Direct(lat, lon, strike + 90, width_on_map * 1e3)
    xs[3], ys[3] = tmp1['lon2'], tmp1['lat2'] # LB
    
    # 各頂点に対応する深さを設定
    # [左上, 右上, 右下, 左下, 左上] の順
    zs = np.array([upper_depth, upper_depth, lower_depth, lower_depth, upper_depth])
    
    return xs, ys, zs

# ----------------------------------------------------------------
# 4. メイン処理
# ----------------------------------------------------------------

def main():
    """
    メインの実行関数
    """
    # --- 4.1 データの読み込み ---
    try:
        # 断層モデルの名前とセグメント数の読み込み
        df_seg = pd.read_excel('/Users/harashima-yukio/Desktop/日本海における大規模地震に関する調査検討会_津波断層モデルのパラメータ_断層モデルの名前とセグメントの数.xlsx')

        # 各断層セグメントのパラメータの読み込み
        # ★注意: 'upper_depth', 'lower_depth' カラムがExcelファイルに含まれている必要があります。
        df = pd.read_excel('/Users/harashima-yukio/Desktop/日本海における大規模地震に関する調査検討会_津波断層モデルのパラメータ_緯度経度など (1).xlsx')
    except FileNotFoundError as e:
        print(f"エラー: 入力ファイルが見つかりません。パスを確認してください。")
        print(e)
        return

    # --- 4.2 出力ディレクトリの作成 ---
    output_dir = 'output_fault_data'
    os.makedirs(output_dir, exist_ok=True)
    print(f"座標データは '{output_dir}' ディレクトリに出力されます。")

    # --- 4.3 断層セグメントの座標計算とファイル出力 ---
    seg_before = 0
    file_count = 0
    # モデルごとにループ
    for model, num_seg in zip(df_seg['model'], df_seg['num_seg']):
        index_b = seg_before
        index_e = seg_before + num_seg - 1
        
        tmpdf = df.loc[index_b:index_e].copy()
        
        # モデルを構成するセグメントごとにループ
        cols_to_zip = ['lat', 'lon', 'strike', 'dip', 'length', 'width', 'upper_depth', 'lower_depth', 'model_sub_segment']
        for row_data in zip(*[tmpdf[col] for col in cols_to_zip]):
            lat, lon, strike, dip, length, width, u_depth, l_depth, sub_segment = row_data

            # 座標と深さを計算
            xs, ys, zs = get_coords_tsunami_fault_model(lat, lon, strike, dip, length, width, u_depth, l_depth)
            
            # 出力ファイル名の生成
            output_filename = f"{sub_segment}.dat"
            output_path = os.path.join(output_dir, output_filename)

            # ファイルにタブ区切りで書き込み
            with open(output_path, 'w') as f:
                # ヘッダー情報の書き込み (任意)
                # f.write(f"# Model: {model}, Segment: {sub_segment}\n")
                # f.write("# Lat\tLon\tDepth(km)\n")
                # 座標データの書き込み
                for i in range(5):
                    f.write(f"{ys[i]:.3f}\t{xs[i]:.3f}\t{zs[i]:.3f}\n")
            
            file_count += 1

        seg_before += num_seg
        
    print(f"\n処理が完了しました。合計 {file_count} 個のファイルが生成されました。")


if __name__ == '__main__':
    main()
