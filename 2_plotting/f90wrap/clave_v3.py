#!/usr/bin/env python3
# clave_v3.py

'''
clusteringSnap2_mod.dat(calc_fault_v2.py) ファイルを読み込み、30列目の no の値でデータをグループ化します。
各グループ内で、モーメントテンソル6成分(Mrr, Mss, Mee, Mrs, Mre, Mse)をそれぞれ合計します。
合計したモーメントテンソルを外部のFortranプログラム `reverse_mrf_v2` に渡し、
断層パラメータ(sliprate, strike, dip, rake)を再計算します。
最終的に、クラスタ番号、再計算された断層パラメータ、合計したモーメントテンソルを
指定されたフォーマットで clave_mod2.dat に出力します。
'''

import numpy as np
import pandas as pd
import sys
import os
import subprocess
import re
from collections import Counter

def run_reverse_mrf(m_rr, m_ss, m_ee, m_rs, m_re, m_se, executable_path):
    """
    コンパイル済みのFortranプログラムを実行し、モーメントテンソルから物理量を計算する。
    """
    if not os.path.exists(executable_path):
        sys.stderr.write(f"エラー: Fortran実行可能ファイルが見つかりません: {executable_path}\n")
        return None

    args = [
        executable_path,
        f"{m_rr:.10e}", f"{m_ss:.10e}", f"{m_ee:.10e}",
        f"{m_rs:.10e}", f"{m_re:.10e}", f"{m_se:.10e}"
    ]

    try:
        result = subprocess.run(args, capture_output=True, text=True, check=True, encoding='utf-8')
        output_params = {}
        pattern = re.compile(r"([A-Z0-9]+):\s*(-?\d+\.\d+)")
        for line in result.stdout.strip().split('\n'):
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = float(match.group(2))
                output_params[key] = value
        
        # 必要なキーが全て揃っているか確認
        required_keys = ['SLIPRATE', 'STRIKE1', 'DIP1', 'RAKE1', 'STRIKE2', 'DIP2', 'RAKE2']
        if not all(key in output_params for key in required_keys):
            sys.stderr.write(f"警告: Fortranからの出力に必要なキーが不足しています。出力: {result.stdout}\n")
            return None
            
        return output_params

    except FileNotFoundError:
        sys.stderr.write(f"エラー: コマンドが見つかりません: {executable_path}\n")
        return None
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"エラー: Fortranプログラムの実行に失敗しました (リターンコード: {e.returncode})\n")
        sys.stderr.write(f"  コマンド: {' '.join(args)}\n")
        sys.stderr.write(f"  標準エラー出力:\n{e.stderr}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"エラー: Fortranプログラムの呼び出し中に予期せぬエラーが発生しました: {e}\n")
        return None

def main():
    input_filename = "clusteringSnap2_mod.dat" # calc_fault_v3.py で出力
    output_filename = "clave_mod2.dat"
    fortran_executable = "/Users/harashima-yukio/program/SnippetsGMT6_m/bin/reverse_mrf_v2" # Fortran実行可能ファイルのパス

    # --- データ読み込み ---
    try:
        df = pd.read_csv(
            input_filename,
            delimiter='\t',
            header=0, # ヘッダーありと指定
            skipinitialspace=True
        )
        print(f"'{input_filename}' から {len(df)} 行のデータを読み込みました。")
    except FileNotFoundError:
        print(f"エラー: 入力ファイル '{input_filename}' が見つかりません。")
        return
    except Exception as e:
        print(f"エラー: '{input_filename}' の読み込み中にエラーが発生しました: {e}")
        return

    # --- noの値毎にモーメントを合計し、Fortranで再計算 ---
    print("グループごとにモーメントを合計し、Fortranで再計算を開始します...")
    
    # 'no' 列とモーメント列が存在するか確認
    moment_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    required_cols = ['no'] + moment_cols
    for col in required_cols:
        if col not in df.columns:
            print(f"エラー: 入力ファイルに必須列 '{col}' が見つかりません。")
            return
            
    # 関連列を数値型に変換 (エラー時はNaN)
    for col in required_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # 'no'がNaNの行は除外
    df.dropna(subset=['no'], inplace=True)
    df['no'] = df['no'].astype(int)

    final_results_list = []
    
    # 'no' でグループ化してループ処理
    for cluster_no, group in df.groupby('no'):
        # 1. モーメントテンソル6成分を合計
        summed_moments = group[moment_cols].sum()
        
        # 2. Fortranプログラム `reverse_mrf_v2` を呼び出す
        recalculated_params = run_reverse_mrf(
            summed_moments['Mrr'], summed_moments['Mss'], summed_moments['Mee'],
            summed_moments['Mrs'], summed_moments['Mre'], summed_moments['Mse'],
            fortran_executable
        )
        
        # 3. Fortran実行が失敗した場合はこのクラスタの処理をスキップ
        if recalculated_params is None:
            print(f"警告: クラスタ no={cluster_no} の計算に失敗しました。スキップします。")
            continue
            
        # 4. 出力用の行データを作成
        output_row = {
            'no': cluster_no,
            'sliprate': recalculated_params.get('SLIPRATE'),
            'mrr': summed_moments['Mrr'],
            'mss': summed_moments['Mss'],
            'mee': summed_moments['Mee'],
            'mrs': summed_moments['Mrs'],
            'mre': summed_moments['Mre'],
            'mse': summed_moments['Mse'],
            'stk0': recalculated_params.get('STRIKE1'),
            'dip0': recalculated_params.get('DIP1'),
            'rake0': recalculated_params.get('RAKE1'),
            'stk1': recalculated_params.get('STRIKE2'),
            'dip1': recalculated_params.get('DIP2'),
            'rake1': recalculated_params.get('RAKE2')
        }
        final_results_list.append(output_row)
        
    df_final = pd.DataFrame(final_results_list)
    print("全クラスタの再計算処理が完了しました。")

    # --- データ出力 ---
    # 指定された列順序
    output_columns_order = [
        'no', 'sliprate', 'mrr', 'mss', 'mee', 'mrs', 'mre', 'mse', 
        'stk0', 'dip0', 'rake0', 'stk1', 'dip1', 'rake1'
    ]
    
    if df_final.empty:
        print("警告: 処理結果が空のため、ファイルは作成されませんでした。")
        return

    # DataFrameの列をこの順序に並べ替える
    df_final = df_final[output_columns_order]

    try:
        df_final.to_csv(output_filename, sep='\t', index=False, header=True, float_format='%.6f')
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()
