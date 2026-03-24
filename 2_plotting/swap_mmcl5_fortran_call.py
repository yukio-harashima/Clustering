#!/usr/bin/env python3

'''
swap_mmcl5_fortran_call.py (v5):
環境変数 CLNO の読み取りとそれに基づくフィルタリング処理を削除しました。
新しい時間窓ごとに、その時間窓内のデータを空間グリッド (n, m) でグループ化します。
各グリッドグループ内で、no != 0 のクラスタ番号の最頻値を求め、代表クラスタを決定します。
改良点:
- clusteringSnap2c.dat が存在しない場合、代わりに clusteringSnap2.dat を使用します。
- (v5) ルールBの変更:
  - 代表クラスタのデータをsliprateで選ぶ代わりに、モーメント6成分(Mrr, Mss, Mee, Mrs, Mre, Mse)を合計します。
  - 合計したモーメントを外部のFortranプログラム `reverse_mrf_v2` に渡し、各種物理量(sliprate, strike, dip, rake等)を再計算します。
  - 再計算された値と、元のデータから引き継いだ位置情報、clla, noを結合して最終的な代表データとします。
'''

import pandas as pd
import numpy as np
import sys
import os
import subprocess # Fortranプログラム呼び出しのためインポート
import re         # 結果の解析のためインポート
from collections import Counter

# clusteringSnap2c.dat および出力ファイルの列定義
CLUSTER_COLUMNS = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]

def run_reverse_mrf(m_rr, m_ss, m_ee, m_rs, m_re, m_se, executable_path):
    """
    コンパイル済みのFortranプログラムを実行し、モーメントテンソルから物理量を計算する。

    Args:
        m_rr, m_ss, m_ee, m_rs, m_re, m_se (float): 6つのモーメントテンソル成分。
        executable_path (str): Fortran実行可能ファイルのパス。

    Returns:
        dict: 計算された物理量の辞書。エラー時はNoneを返す。
    """
    if not os.path.exists(executable_path):
        sys.stderr.write(f"エラー: Fortran実行可能ファイルが見つかりません: {executable_path}\n")
        return None

    # コマンドライン引数を文字列のリストとして作成
    args = [
        executable_path,
        f"{m_rr:.10e}", f"{m_ss:.10e}", f"{m_ee:.10e}",
        f"{m_rs:.10e}", f"{m_re:.10e}", f"{m_se:.10e}"
    ]

    try:
        # Fortranプログラムをサブプロセスとして実行
        result = subprocess.run(args, capture_output=True, text=True, check=True, encoding='utf-8')

        # 標準出力を解析して結果を辞書に格納
        output_params = {}
        # 正規表現パターン: キー(英字):値(浮動小数点数)
        pattern = re.compile(r"([A-Z0-9]+):\s*(-?\d+\.\d+)")
        for line in result.stdout.strip().split('\n'):
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = float(match.group(2))
                output_params[key] = value
        
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


def revised_swap_process():
    wd = "./"
    ford = "/Users/harashima-yukio/program/SnippetsGMT6_m/bin/"
    fortran_executable = os.path.join(ford, "reverse_mrf_v2") # Fortran実行可能ファイルのパス

    try:
        # (データ読み込み部分は元のコードと同じ)
        mrf_df = pd.read_csv(
            wd + "mrf.dat", sep=r'\s+', header=None,
            names=['time', 'momentrate', 'Mrr_mrf', 'Mss_mrf', 'Mee_mrf', 'Mrs_mrf', 'Mre_mrf', 'Mse_mrf'],
            engine='python', skipinitialspace=True, dtype={'time': float}
        )
        if mrf_df.empty or 'time' not in mrf_df.columns or mrf_df['time'].isna().all():
            sys.stderr.write("エラー: mrf.dat が空か、有効な 'time' 列がありません。\n")
            return 1
        t_max = mrf_df['time'].max()
        if pd.isna(t_max) or t_max <= 0:
            sys.stderr.write(f"エラー: mrf.dat から有効な t_max を取得できませんでした (t_max={t_max})。\n")
            return 1

        cluster_file_c = wd + "clusteringSnap2c.dat"
        cluster_file_plain = wd + "clusteringSnap2.dat"
        
        cluster_data_file = ""
        if os.path.exists(cluster_file_c):
            cluster_data_file = cluster_file_c
            print(f"情報: {cluster_data_file} を使用します。")
        elif os.path.exists(cluster_file_plain):
            cluster_data_file = cluster_file_plain
            print(f"情報: {cluster_file_c} が見つかりませんでした。代わりに {cluster_data_file} を使用します。")
        else:
            sys.stderr.write(f"エラー: {cluster_file_c} も {cluster_file_plain} も見つかりません。\n")
            return 1

        cluster_df_raw = pd.read_csv(
            cluster_data_file, sep='\t', header=None,
            names=CLUSTER_COLUMNS, dtype=str, engine='python', skipinitialspace=True
        )

        if cluster_df_raw.empty:
            sys.stderr.write(f"情報: {cluster_data_file} が空です。\n")
            pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
                 wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
            )
            return 0
        
        cluster_df = pd.DataFrame()
        for col in CLUSTER_COLUMNS:
            if col in cluster_df_raw.columns:
                if col in ['clla']: 
                    cluster_df[col] = cluster_df_raw[col].replace({'NaN': np.nan, 'nan': np.nan, '':np.nan}).astype(object)
                else:
                    cluster_df[col] = pd.to_numeric(cluster_df_raw[col], errors='coerce')
            else:
                cluster_df[col] = np.nan
        
        essential_cols = ['tw', 'no', 'n', 'm'] + ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
        cluster_df.dropna(subset=essential_cols, inplace=True)
        if cluster_df.empty:
            sys.stderr.write(f"情報: 必須列のNaN除去後、{cluster_data_file} に有効なデータがありません。\n")
            pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
                 wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
            )
            return 0

        cluster_max_original_tw = cluster_df['tw'].max()
        if pd.isna(cluster_max_original_tw) or cluster_max_original_tw <= 0:
            sys.stderr.write(f"エラー: {cluster_data_file} から有効な cluster_max_original_tw を取得できませんでした。\n")
            return 1

        new_time_windows_df = pd.read_csv(
            wd + "new_time_windows.dat", sep=r'\s+', header=None,
            names=['new_tw_id', 'start_time', 'end_time'], engine='python', skipinitialspace=True,
            dtype={'new_tw_id': int, 'start_time': float, 'end_time': float}
        )
        if new_time_windows_df.empty:
            sys.stderr.write("情報: new_time_windows.dat が空です。\n")
            pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
                 wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
            )
            return 0

    except FileNotFoundError as e:
        sys.stderr.write(f"エラー: 入力ファイルが見つかりません - {e}\n")
        return 1
    except Exception as e:
        sys.stderr.write(f"エラー: データ読み込みまたは初期処理中にエラー - {e}\n")
        return 1

    time_per_original_tw = t_max / cluster_max_original_tw
    cluster_df['calc_real_time'] = (cluster_df['tw'].astype(float) - 0.5) * time_per_original_tw
    
    overall_start_time = new_time_windows_df['start_time'].min()
    overall_end_time = new_time_windows_df['end_time'].max()
    
    cluster_df.dropna(subset=['calc_real_time'], inplace=True)

    initial_filtered_df = cluster_df[
        (cluster_df['calc_real_time'] >= overall_start_time) &
        (cluster_df['calc_real_time'] < overall_end_time) 
    ].copy() 

    if initial_filtered_df.empty:
        sys.stderr.write(f"情報: 指定された全体時間範囲に一致するデータが {cluster_data_file} にありません。\n")
        pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
             wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
        )
        return 0
        
    final_output_rows_list = []

    for _, new_window_row in new_time_windows_df.iterrows():
        current_new_tw_id = int(new_window_row['new_tw_id'])
        new_window_start_real_time = new_window_row['start_time']
        new_window_end_real_time = new_window_row['end_time']

        data_in_current_new_window = initial_filtered_df[
            (initial_filtered_df['calc_real_time'] >= new_window_start_real_time) &
            (initial_filtered_df['calc_real_time'] < new_window_end_real_time)
        ]

        if data_in_current_new_window.empty:
            continue

        for grid_keys, group_at_grid in data_in_current_new_window.groupby(['n', 'm']):
            valid_no_series = group_at_grid[group_at_grid['no'] != 0]['no']
            if valid_no_series.empty:
                continue

            no_counts = Counter(valid_no_series.dropna().astype(int))
            if not no_counts: continue
            
            most_common_nos_with_freq = no_counts.most_common()
            if not most_common_nos_with_freq: continue

            max_freq = most_common_nos_with_freq[0][1]
            representative_no_candidates = [item[0] for item in most_common_nos_with_freq if item[1] == max_freq]
            if not representative_no_candidates: continue
            
            representative_no = min(representative_no_candidates)
            data_with_representative_no = group_at_grid[group_at_grid['no'] == representative_no]
            
            if data_with_representative_no.empty: continue

            # ルールB: モーメント成分を合計し、Fortranで再計算
            moment_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
            summed_moments = data_with_representative_no[moment_cols].sum()

            # Fortranプログラムを実行して新しい物理量を取得
            recalculated_params = run_reverse_mrf(
                summed_moments['Mrr'], summed_moments['Mss'], summed_moments['Mee'],
                summed_moments['Mrs'], summed_moments['Mre'], summed_moments['Mse'],
                fortran_executable
            )
            
            if recalculated_params is None:
                sys.stderr.write(f"警告: グリッド(n={grid_keys[0]}, m={grid_keys[1]})の計算に失敗しました。スキップします。\n")
                continue

            # --- ▼▼▼ ここからがご要望の反映箇所です ▼▼▼ ---
            #
            # 結果を格納するための新しい行を作成
            # 代表クラスタの最初の行をテンプレートとしてコピーする。
            # これにより、no, clla, lat, lon, depth等が自動的に引き継がれる。
            selected_row = data_with_representative_no.iloc[0].copy()

            # 1. テンプレートの値を、合計したモーメントと再計算した物理量で上書き
            selected_row.update(summed_moments) # 合計したモーメントで更新
            
            selected_row['sliprate'] = recalculated_params.get('SLIPRATE', np.nan)
            selected_row['str1']     = recalculated_params.get('STRIKE1', np.nan)
            selected_row['dip1']     = recalculated_params.get('DIP1', np.nan)
            selected_row['rake1']    = recalculated_params.get('RAKE1', np.nan)
            selected_row['str2']     = recalculated_params.get('STRIKE2', np.nan)
            selected_row['dip2']     = recalculated_params.get('DIP2', np.nan)
            selected_row['rake2']    = recalculated_params.get('RAKE2', np.nan)
            selected_row['trendp']   = recalculated_params.get('TRENDP', np.nan)
            selected_row['plungp']   = recalculated_params.get('PLUNGP', np.nan)
            selected_row['trendt']   = recalculated_params.get('TRENDT', np.nan)
            selected_row['plungt']   = recalculated_params.get('PLUNGT', np.nan)
            selected_row['trendb']   = recalculated_params.get('TRENDB', np.nan)
            selected_row['plungb']   = recalculated_params.get('PLUNGB', np.nan)
            selected_row['NDC']      = recalculated_params.get('NDC', np.nan)
            
            # 2. 最後に新しい時間窓IDを設定
            selected_row['tw'] = current_new_tw_id

            # 3. 最終結果リストに追加
            final_output_rows_list.append(selected_row)
            #
            # --- ▲▲▲ ここまでがご要望の反映箇所です ▲▲▲ ---

    if not final_output_rows_list:
        sys.stderr.write("情報: 全ての処理後、有効なデータが残りませんでした。\n")
        final_result_df = pd.DataFrame(columns=CLUSTER_COLUMNS)
    else:
        # リストからDataFrameを作成
        final_result_df = pd.DataFrame(final_output_rows_list)
        if 'calc_real_time' in final_result_df.columns:
            final_result_df.drop(columns=['calc_real_time'], inplace=True)
        
        # 列の並び順を整える
        for col in CLUSTER_COLUMNS:
            if col not in final_result_df.columns:
                final_result_df[col] = np.nan
        final_result_df = final_result_df[CLUSTER_COLUMNS]

    output_filename = wd + "clusteringSnap2_v2.dat"
    try:
        final_result_df.to_csv(output_filename, sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN')
        print(f"情報: {output_filename} が正常に生成されました。行数: {len(final_result_df)}")
    except Exception as e:
        sys.stderr.write(f"エラー: {output_filename} の書き込み中にエラー - {e}\n")
        return 1
        
    return 0

if __name__ == '__main__':
    print("swap_mmcl5_fortran_call.py (v5) を実行中...")
    exit_code = revised_swap_process()
    if exit_code == 0:
        print('swap_mmcl5_fortran_call.py: データ処理完了。')
    else:
        print('swap_mmcl5_fortran_call.py: エラーが発生しました。')
    sys.exit(exit_code)

