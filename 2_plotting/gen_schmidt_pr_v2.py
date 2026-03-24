#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# gen_schimdt_pr_v2.py
# 下半球等積投影のための座標変換コード
# 修正版: 入力データのベクトル順序を [北, 東, 下] として正しく処理する
# 追加機能: vec_inner.dat, model_vec.dat の処理を追加
# 2026/1/1

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

def process_original_clustering_data():
    """既存のclustering_vector.datに対する処理"""
    input_file = "clustering_vector.dat"
    generator_script = "gen_clustering_vector.py"
    output_filename = "schmidt_projection.dat"
    
    print(f"\n--- {input_file} の処理 ---")
    
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
                    return
            except subprocess.CalledProcessError as e:
                print(f"エラー: '{generator_script}' の実行に失敗しました: {e}")
                return
        else:
            print(f"エラー: '{input_file}' も '{generator_script}' も見つかりません。この処理をスキップします。")
            return
            
    try:
        # clustering_vector.dat はタブ区切り、ヘッダーなし
        df = pd.read_table(path, sep='\t', header=None)
        
        # カラム構成の前提（n, m, tw, t, x, y, n1(N,E,D), n2(N,E,D), sliprate, clset, label）
        n1_n, n1_e, n1_d = df[6], df[7], df[8]
        n2_n, n2_e, n2_d = df[9], df[10], df[11]
        
        # 計算: calculate_schmidt(East, North, Down)
        df['n1_X'], df['n1_Y'] = calculate_schmidt(n1_e, n1_n, n1_d)
        df['n2_X'], df['n2_Y'] = calculate_schmidt(n2_e, n2_n, n2_d)
        
        # 出力用データの抽出
        output_cols = [0, 1, 2, 'n1_X', 'n1_Y', 'n2_X', 'n2_Y', 12, 14]
        df_output = df[output_cols]
        
        df_output.to_csv(output_filename, sep='\t', index=False, header=False)
        print(f"完了: '{output_filename}' を保存しました。")
            
    except Exception as e:
        print(f"エラー: {e}")

def process_vec_inner():
    """vec_inner.datに対する処理"""
    input_file = "vec_inner.dat"
    output_file = "vec_inner_sch.dat"
    
    print(f"\n--- {input_file} の処理 ---")
    
    path = find_file(input_file)
    if not path:
        print("vec_inner_sch.datが生成されません")
        return

    try:
        # vec_inner.dat はヘッダーありと想定 (calc_fault_v4_detail.pyの出力)
        df = pd.read_csv(path, sep='\t')
        
        # 法線ベクトル成分の取得 (calc_fault_v4_detail.pyの出力列名に基づく)
        # near_N, near_E, near_D / far_N, far_E, far_D
        
        # 計算: calculate_schmidt(East, North, Down)
        # 引数の順序に注意: East, North, Down の順で渡す
        df['n1_X'], df['n1_Y'] = calculate_schmidt(df['near_E'], df['near_N'], df['near_D'])
        df['n2_X'], df['n2_Y'] = calculate_schmidt(df['far_E'], df['far_N'], df['far_D'])
        
        # 出力データの作成
        # 要求列: n, m, tw, n1_X, n1_Y, n2_X, n2_Y, no, inner_abs_near, inner_abs_far
        output_cols = [
            'n', 'm', 'tw', 
            'n1_X', 'n1_Y', 
            'n2_X', 'n2_Y', 
            'no', 
            'inner_abs_near', 'inner_abs_far'
        ]
        
        # 指定された列が存在するか確認しつつ抽出
        df_output = df[output_cols]
        
        # ヘッダーなしで保存
        df_output.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"完了: '{output_file}' を保存しました。")
        
    except Exception as e:
        print(f"エラー: {e}")

def process_model_vec():
    """model_vec.datに対する処理"""
    input_file = "model_vec.dat"
    output_file = "model_vec_sch.dat"
    
    print(f"\n--- {input_file} の処理 ---")
    
    path = find_file(input_file)
    if not path:
        print("model_vec_sch.datが生成されません")
        return

    try:
        # model_vec.dat はヘッダーありと想定
        # 列名: model_N, model_E, model_D
        df = pd.read_csv(path, sep='\t')
        
        # 計算: calculate_schmidt(East, North, Down)
        df['n1_X'], df['n1_Y'] = calculate_schmidt(df['model_E'], df['model_N'], df['model_D'])
        
        # 出力データの作成
        # 要求列: モデル断層面の法線ベクトル成分（North, East, Down）, n1_X, n1_Y
        output_cols = ['model_N', 'model_E', 'model_D', 'n1_X', 'n1_Y']
        
        df_output = df[output_cols]
        
        # ヘッダーなしで保存
        df_output.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"完了: '{output_file}' を保存しました。")
        
    except Exception as e:
        print(f"エラー: {e}")

def main():
    # 既存の処理
    process_original_clustering_data()
    
    # 追加された処理
    process_vec_inner()
    process_model_vec()
    
    print("\n全ての処理が終了しました。")

if __name__ == "__main__":
    main()