#!/usr/bin/env python3
# calc_fault_v4_detail.py
# 目的: 
# 1. モデル断層面に近い面を判定し、stk0側に整列させる（v3と同様）
# 2. その際の法線ベクトル成分、内積の絶対値などの詳細情報を vec_inner.dat に出力する
# 3. モデル断層面の法線ベクトルを model_vec.dat に出力する
#
# 法線ベクトルの座標系: [North, East, Down] (Fortranコード準拠)

import numpy as np
import pandas as pd
import sys

# --- モデル断層面の設定 ---
MODEL_STRIKE_DEG = 179.0  # 走向
MODEL_DIP_DEG = 55.0    # 傾斜

# --- 入力データの列名定義 ---
INPUT_COLUMNS_30 = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'stk0', 'stk1', 'dip0', 'dip1', 
    'rake0', 'rake1', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]

# --- 出力する列の定義 (vec_inner.dat) ---
OUTPUT_COLUMNS = [
    'n', 'm', 'tw', 
    'stk0', 'stk1', 'dip0', 'dip1', 'rake0', 'rake1', 
    'no',
    'near_N', 'near_E', 'near_D',  # モデルに近い面の法線ベクトル
    'far_N', 'far_E', 'far_D',     # モデルに遠い面の法線ベクトル
    'inner_abs_near',              # 近い方の内積絶対値
    'inner_abs_far'                # 遠い方の内積絶対値
]

def faultnormalvec(stk, dip):
    """
    断層面の法線ベクトルを計算する。
    戻り値: [North, East, Down]
    """
    rad_stk = np.deg2rad(stk)
    rad_dip = np.deg2rad(dip)
    
    # 北成分 (North)
    nn = -np.sin(rad_stk) * np.sin(rad_dip)
    # 東成分 (East)
    ne =  np.cos(rad_stk) * np.sin(rad_dip)
    # 下成分 (Down)
    nd = -np.cos(rad_dip)
    
    return np.array([nn, ne, nd])

def main():
    input_filename = "clusteringSnap2.dat"
    output_vec_filename = "vec_inner.dat"
    output_model_filename = "model_vec.dat"

    # --- データ読み込み ---
    print(f"読み込み開始: {input_filename}")
    try:
        df = pd.read_csv(
            input_filename,
            delimiter='\t',
            header=None,
            names=INPUT_COLUMNS_30,
            skipinitialspace=True
        )
    except FileNotFoundError:
        print(f"エラー: '{input_filename}' が見つかりません。")
        sys.exit(1)

    # --- 1. モデル断層面のベクトル計算と出力 ---
    vec_model = faultnormalvec(MODEL_STRIKE_DEG, MODEL_DIP_DEG)
    print(f"モデル断層面 (Strike={MODEL_STRIKE_DEG}, Dip={MODEL_DIP_DEG})")
    print(f"モデル法線ベクトル [North, East, Down]: {vec_model}")

    # model_vec.dat の作成
    try:
        df_model = pd.DataFrame([vec_model], columns=['model_N', 'model_E', 'model_D'])
        df_model.to_csv(output_model_filename, sep='\t', index=False, header=True)
        print(f"モデルベクトルを '{output_model_filename}' に保存しました。")
    except Exception as e:
        print(f"エラー: {output_model_filename} の保存に失敗しました: {e}")

    # --- 2. 各行の計算処理 ---
    print("各データの断層面判定とベクトル計算を実行中...")
    
    # 入れ替えカウント用
    swap_count = 0

    def process_row(row):
        nonlocal swap_count
        
        # 既存データの取得
        n, m, tw, no = row['n'], row['m'], row['tw'], row['no']
        s0, d0, r0 = row['stk0'], row['dip0'], row['rake0']
        s1, d1, r1 = row['stk1'], row['dip1'], row['rake1']

        # ベクトル計算 [North, East, Down]
        v0 = faultnormalvec(s0, d0)
        v1 = faultnormalvec(s1, d1)

        # モデルとの内積 (絶対値)
        score0 = abs(np.inner(vec_model, v0))
        score1 = abs(np.inner(vec_model, v1))

        # 判定とデータの整列
        # score1 (面1) の方が大きければ入れ替える
        if score1 > score0:
            swap_count += 1
            # 近い方 = 面1
            stk_near, dip_near, rake_near = s1, d1, r1
            stk_far,  dip_far,  rake_far  = s0, d0, r0
            
            vec_near = v1
            vec_far  = v0
            
            score_near = score1
            score_far  = score0
        else:
            # 近い方 = 面0 (そのまま)
            stk_near, dip_near, rake_near = s0, d0, r0
            stk_far,  dip_far,  rake_far  = s1, d1, r1
            
            vec_near = v0
            vec_far  = v1
            
            score_near = score0
            score_far  = score1

        # 指定された列順序でデータを返す
        return pd.Series([
            n, m, tw,
            stk_near, stk_far, dip_near, dip_far, rake_near, rake_far,
            no,
            vec_near[0], vec_near[1], vec_near[2], # near_N, near_E, near_D
            vec_far[0],  vec_far[1],  vec_far[2],  # far_N, far_E, far_D
            score_near,
            score_far
        ], index=OUTPUT_COLUMNS)

    # 適用して新しいDataFrameを作成
    df_result = df.apply(process_row, axis=1)

    print(f"処理完了: 入れ替えられた行数 {swap_count} / {len(df)}")

    # --- 3. 結果の出力 ---
    try:
        # floatのフォーマットを指定して綺麗に出力することも可能ですが、
        # ここでは標準のまま出力します。
        df_result.to_csv(output_vec_filename, sep='\t', index=False, header=True)
        print(f"詳細データを '{output_vec_filename}' に保存しました。")
        
        print("\n[出力列の構成]")
        print(", ".join(OUTPUT_COLUMNS))
        
    except Exception as e:
        print(f"エラー: {output_vec_filename} の保存に失敗しました: {e}")

if __name__ == '__main__':
    main()