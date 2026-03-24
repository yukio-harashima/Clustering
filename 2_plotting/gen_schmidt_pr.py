#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# gen_schimdt_pr.py
# 下半球等積投影のための座標変換コード
# 修正版: 入力データのベクトル順序を [北, 東, 下] として正しく処理する
# 2025/12/30

import pandas as pd
import numpy as np
import os
import sys
import subprocess

def find_file(filename):
    """カレントディレクトリまたは1つ上のディレクトリからファイルを探す"""
    paths = [os.path.join(".", filename), os.path.join("..", filename)]
    for path in paths:
        if os.path.exists(path):
            return path
    return None

def calculate_schmidt(nx, ny, nz):
    """
    シュミットネット（下半球等積投影）座標の計算
    
    引数:
      nx: East 成分 (X軸方向)
      ny: North 成分 (Y軸方向)
      nz: Down 成分
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
    # 係数 k = sqrt( 1 / (1 + nz) )
    k = np.sqrt(1.0 / (1.0 + nz))
    X = nx * k
    Y = ny * k
    
    return X, Y

def main():
    input_file = "clustering_vector.dat"
    generator_script = "gen_clustering_vector.py"
    output_filename = "schmidt_projection.dat"
    
    # 1. 入力ファイルの存在確認
    path = find_file(input_file)
    
    if not path:
        print(f"'{input_file}' が見つかりません。")
        # 2. 生成スクリプトの存在確認と実行
        gen_path = find_file(generator_script)
        if gen_path:
            print(f"'{generator_script}' を実行してデータを生成します...")
            try:
                subprocess.run([sys.executable, gen_path], check=True)
                path = find_file(input_file)
                if not path:
                    print(f"エラー: データ生成を試みましたが '{input_file}' が見つかりません。")
                    sys.exit(1)
            except subprocess.CalledProcessError as e:
                print(f"エラー: '{generator_script}' の実行に失敗しました: {e}")
                sys.exit(1)
        else:
            print(f"エラー: '{input_file}' も '{generator_script}' も見つかりません。処理を中断します。")
            sys.exit(1)
            
    print(f"データ読み込み開始: {path}")
    
    try:
        # clustering_vector.dat はタブ区切り、ヘッダーなし
        df = pd.read_table(path, sep='\t', header=None)
        
        # --- 修正箇所: データの読み込み順序 ---
        # カラム構成の前提（n, m, tw, t, x, y, n1(N,E,D), n2(N,E,D), sliprate, clset, label）
        # 列インデックス:
        # 6: n1_North, 7: n1_East, 8: n1_Down
        # 9: n2_North, 10: n2_East, 11: n2_Down
        
        n1_n, n1_e, n1_d = df[6], df[7], df[8]
        n2_n, n2_e, n2_d = df[9], df[10], df[11]
        
        print("投影座標を計算中...")
        
        # --- 修正箇所: 引数の渡し順序 ---
        # calculate_schmidt(nx, ny, nz) は (East, North, Down) を期待しているため、
        # 引数を (n1_e, n1_n, n1_d) の順に入れ替えて渡す。
        
        df['n1_X'], df['n1_Y'] = calculate_schmidt(n1_e, n1_n, n1_d)
        df['n2_X'], df['n2_Y'] = calculate_schmidt(n2_e, n2_n, n2_d)
        
        # 出力用データの抽出
        # 列構成: n(0), m(1), tw(2), n1_X, n1_Y, n2_X, n2_Y, sliprate(12), label(14)
        # ※ sliprateが12列目、labelが14列目にある前提（入力ファイルの構造による）
        output_cols = [0, 1, 2, 'n1_X', 'n1_Y', 'n2_X', 'n2_Y', 12, 14]
        df_output = df[output_cols]
        
        df_output.to_csv(output_filename, sep='\t', index=False, header=False)
        
        print(f"--- 成功 ---")
        print(f"投影データを '{output_filename}' に保存しました。")
        print("\n[出力フォーマット (9列)]")
        print("n, m, tw, n1_X, n1_Y, n2_X, n2_Y, sliprate, label")
            
    except Exception as e:
        print(f"実行中にエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()