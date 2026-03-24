#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# schimdt_cal.py
# n_vector.dat を読み込み、シュミットネット（下半球等積投影）座標を計算して出力する
# ffi_schimdt_plot_nv.bash用データ生成コード
# 基コード　gen_schimdt_pr.py

import pandas as pd
import numpy as np
import sys
import os

# 入力ファイルの設定
INPUT_FILENAME = "n_vector.dat"
OUTPUT_FILENAME = "schmidt_vector.dat"

# n_vector.dat の列構成 (align_fault_data.py の仕様に基づく)
# n1_(1): North, n1_(2): East, n1_(3): Down
INPUT_COLS = [
    'n', 'm', 'tw', 't', 'x_coord', 'y_coord',
    'n1_N', 'n1_E', 'n1_D',  # n1_(1)..(3)
    'n2_N', 'n2_E', 'n2_D',  # n2_(1)..(3)
    'sliprate'
]

def calculate_schmidt(nx, ny, nz):
    """
    シュミットネット（下半球等積投影）座標の計算
    
    引数:
      nx: East 成分
      ny: North 成分
      nz: Down 成分
    
    戻り値:
      X, Y: 投影面上の座標
    """
    # 単位ベクトル化
    norm = np.sqrt(nx**2 + ny**2 + nz**2)
    # 0除算対策
    norm = np.where(norm == 0, 1e-9, norm)
    nx, ny, nz = nx/norm, ny/norm, nz/norm
    
    # 下半球（nz >= 0）に統一。nzが負ならベクトルを反転させる
    mask = nz < 0
    nx = np.where(mask, -nx, nx)
    ny = np.where(mask, -ny, ny)
    nz = np.where(mask, -nz, nz)
    
    # シュミット投影計算 (Lambert Equal-Area Projection)
    # 投影式: R = 2 * sin(theta/2) などの変形等価式
    # 係数 k = sqrt( 2 / (1 + nz) ) * (1/sqrt(2)) = sqrt( 1 / (1 + nz) )
    k = np.sqrt(1.0 / (1.0 + nz))
    
    # X軸をEast, Y軸をNorthとする投影座標
    X = nx * k
    Y = ny * k
    
    return X, Y

def main():
    print("先にalign_fault.py(節面選択用コード)を実行することを推奨します") 
    print(f"処理を開始します: {INPUT_FILENAME} -> {OUTPUT_FILENAME}")

    # 1. 入力ファイルの存在確認 (実行ディレクトリ)
    if not os.path.exists(INPUT_FILENAME):
        print(f"エラー: 実行ディレクトリに入力ファイル '{INPUT_FILENAME}' が見つかりません。処理を中止します。")
        sys.exit(1)

    try:
        # データ読み込み (空白区切り)
        df = pd.read_csv(INPUT_FILENAME, sep=r'\s+', header=None, names=INPUT_COLS)
        print(f"データ読み込み完了: {len(df)} 行")

        # --- 計算処理 ---
        print("投影座標を計算中...")

        # n_vector.dat の仕様: [North, East, Down]
        # calculate_schmidt の仕様: (East, North, Down)
        # そのため、引数の順序を入れ替えて渡す必要があります。
        
        # 面1 (n1) の計算
        # nx=East(n1_E), ny=North(n1_N), nz=Down(n1_D)
        df['n1_X'], df['n1_Y'] = calculate_schmidt(df['n1_E'], df['n1_N'], df['n1_D'])

        # 面2 (n2) の計算
        df['n2_X'], df['n2_Y'] = calculate_schmidt(df['n2_E'], df['n2_N'], df['n2_D'])

        # --- データ出力 ---
        # 必要な列を選択して保存
        # 出力構成: n, m, tw, n1_X, n1_Y, n2_X, n2_Y, sliprate
        output_cols = ['n', 'm', 'tw', 'n1_X', 'n1_Y', 'n2_X', 'n2_Y', 'sliprate']
        
        df_output = df[output_cols]
        
        df_output.to_csv(OUTPUT_FILENAME, sep='\t', index=False, header=False)
        
        print(f"--- 成功 ---")
        print(f"投影データを '{OUTPUT_FILENAME}' に保存しました。")
        print("出力フォーマット: n, m, tw, n1_X, n1_Y, n2_X, n2_Y, sliprate")

    except Exception as e:
        print(f"実行中に予期せぬエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()