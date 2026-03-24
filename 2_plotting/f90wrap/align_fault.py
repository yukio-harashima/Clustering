#!/usr/bin/env python3
# align_fault.py
# calc_fault_v2.py/v3.py をベースにした複数ファイル対応版
# 
# 目的:
# 1. snap_y.dat (角度データ) を読み込み、モデル断層面に近い方を面1 (str1等) に揃えて上書き保存
# 2. n_vector.dat (ベクトルデータ) を読み込み、モデル断層面に近い方を面1 (n1等) に揃えて上書き保存
#
# モデル断層面設定: 走向 179.0, 傾斜 55.0
# ベクトル座標系: [北, 東, 下]

import numpy as np
import pandas as pd

# --- モデル断層面の設定 ---
MODEL_STRIKE = 179.0
MODEL_DIP = 55.0

# --- ファイル1: snap_y.dat の設定 ---
FILE1_NAME = 'snap_y.dat'
FILE1_COLS = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse',
    'str1', 'str2', 'dip1', 'dip2', 'rake1', 'rake2',
    'lat', 'lon', 'depth', 'trendp', 'trendt', 'trendb',
    'plungp', 'plungt', 'plungb', 'NDC'
]

# --- ファイル2: n_vector.dat の設定 ---
FILE2_NAME = 'n_vector.dat'
FILE2_COLS = [
    'n', 'm', 'tw', 't', 'x_coord', 'y_coord',
    'n1_(1)', 'n1_(2)', 'n1_(3)',  # 面0のベクトル (北, 東, 下)
    'n2_(1)', 'n2_(2)', 'n2_(3)',  # 面1のベクトル (北, 東, 下)
    'sliprate'
]


def faultnormalvec(stk, dip):
    """
    走向(stk)と傾斜(dip)から法線ベクトルを計算する関数
    戻り値の順序: [北(North), 東(East), 下(Down)]
    """
    rad_stk = np.deg2rad(stk)
    rad_dip = np.deg2rad(dip)
    
    nn = -np.sin(rad_stk) * np.sin(rad_dip) # 北
    ne =  np.cos(rad_stk) * np.sin(rad_dip) # 東
    nd = -np.cos(rad_dip)                   # 下
    
    return np.array([nn, ne, nd])


def process_snap_y():
    print(f"--- {FILE1_NAME} の処理を開始します ---")
    try:
        # sep=r'\s+' は空白区切り(スペース複数も可)に対応
        df = pd.read_csv(FILE1_NAME, sep=r'\s+', header=None, names=FILE1_COLS)
    except FileNotFoundError:
        print(f"エラー: {FILE1_NAME} が見つかりません。")
        print(f"fgensnap_y(gensnap_y.f90)または、fgensnap_y_mf(gen_snap_ymf.f90)を実行してください")
        return

    # モデル断層面の法線ベクトル (北, 東, 下)
    vec_model = faultnormalvec(MODEL_STRIKE, MODEL_DIP)
    
    # 入れ替えカウント用変数
    swap_count = 0

    def apply_sort_planes(row):
        nonlocal swap_count # 外側の変数を更新するために必要

        # 面0 (str1) の法線ベクトル計算
        v1 = faultnormalvec(row['str1'], row['dip1'])
        # 面1 (str2) の法線ベクトル計算
        v2 = faultnormalvec(row['str2'], row['dip2'])

        # モデルとの内積 (絶対値が大きい方が向きが近い)
        score1 = abs(np.inner(vec_model, v1))
        score2 = abs(np.inner(vec_model, v2))

        # 面1の方がモデルに近い場合、値を入れ替える
        # 近い方を (str1, dip1, rake1) にする
        if score2 > score1:
            swap_count += 1
            return pd.Series([
                row['str2'], row['dip2'], row['rake2'],
                row['str1'], row['dip1'], row['rake1']
            ])
        else:
            # そのまま
            return pd.Series([
                row['str1'], row['dip1'], row['rake1'],
                row['str2'], row['dip2'], row['rake2']
            ])

    # 適用して値を更新
    target_cols = ['str1', 'dip1', 'rake1', 'str2', 'dip2', 'rake2']
    df[target_cols] = df.apply(apply_sort_planes, axis=1)

    # 上書き保存 (区切り文字はタブとして出力します。読み込み時は\s+なので互換性あり)
    df.to_csv(FILE1_NAME, sep='\t', index=False, header=False)
    print(f"{FILE1_NAME} の処理が完了し、上書き保存しました。")
    print(f"入れ替えられた行数: {swap_count} / {len(df)}")


def process_n_vector():
    print(f"--- {FILE2_NAME} の処理を開始します ---")
    try:
        df = pd.read_csv(FILE2_NAME, sep=r'\s+', header=None, names=FILE2_COLS)
    except FileNotFoundError:
        print(f"エラー: {FILE2_NAME} が見つかりません。")
        print(f"fgensnap_y(gensnap_y.f90)または、fgensnap_y_mf(gen_snap_ymf.f90)を実行してください")
        return

    # モデル断層面の法線ベクトル (北, 東, 下)
    vec_model = faultnormalvec(MODEL_STRIKE, MODEL_DIP)

    # 入れ替えカウント用変数
    swap_count = 0

    def apply_sort_vectors(row):
        nonlocal swap_count # 外側の変数を更新するために必要
        
        # データ内のベクトルを取得 (既に [北, 東, 下] の順であることを前提)
        v1 = np.array([row['n1_(1)'], row['n1_(2)'], row['n1_(3)']])
        v2 = np.array([row['n2_(1)'], row['n2_(2)'], row['n2_(3)']])

        # モデルとの内積
        score1 = abs(np.inner(vec_model, v1))
        score2 = abs(np.inner(vec_model, v2))

        # 面1の方がモデルに近い場合、ベクトル成分を入れ替える
        # 近い方を n1系 にする
        if score2 > score1:
            swap_count += 1
            return pd.Series([
                v2[0], v2[1], v2[2],
                v1[0], v1[1], v1[2]
            ])
        else:
            return pd.Series([
                v1[0], v1[1], v1[2],
                v2[0], v2[1], v2[2]
            ])

    # 適用して値を更新
    target_cols = ['n1_(1)', 'n1_(2)', 'n1_(3)', 'n2_(1)', 'n2_(2)', 'n2_(3)']
    df[target_cols] = df.apply(apply_sort_vectors, axis=1)

    # 上書き保存
    df.to_csv(FILE2_NAME, sep='\t', index=False, header=False)
    print(f"{FILE2_NAME} の処理が完了し、上書き保存しました。")
    print(f"入れ替えられた行数: {swap_count} / {len(df)}")


def main():
    print("処理を開始します...")
    print(f"参照面 strike:{MODEL_STRIKE},dip:{MODEL_DIP}")
    process_snap_y()
    print("-" * 30)
    process_n_vector()
    print("全ての処理が終了しました。")

if __name__ == '__main__':
    main()