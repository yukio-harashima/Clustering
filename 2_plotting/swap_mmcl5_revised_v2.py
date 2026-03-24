#!/usr/bin/env python3

'''
swap_mmcl5_revised_v2.py:
環境変数 CLNO の読み取りとそれに基づくフィルタリング処理を削除しました。
新しい時間窓ごとに、その時間窓内のデータを空間グリッド (n, m) でグループ化します。
各グリッドグループ内で、no != 0 のクラスタ番号の最頻値を求め、そのクラスタ番号を持つデータを抽出します。最頻値が複数ある場合は、no の値が最も小さいものを代表とします。
**改良点:** sliprateによる最終フィルタリングを削除し、代表クラスタ番号を持つ全てのデータを残すように変更しました。
**改良点:** clusteringSnap2c.dat が存在しない場合、代わりに clusteringSnap2.dat を使用するように変更しました。
'''


import pandas as pd
import numpy as np
import sys
import os  # ファイルの存在確認のためにosライブラリをインポート
from collections import Counter

# clusteringSnap2c.dat および出力ファイルの列定義
CLUSTER_COLUMNS = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]

def revised_swap_process():
    wd = "./"

    try:
        # 1. 入力ファイルの読み込み
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

        # 使用するクラスタデータファイルを決定する
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
            # どちらのファイルも見つからない場合はエラー終了
            sys.stderr.write(f"エラー: {cluster_file_c} も {cluster_file_plain} も見つかりません。\n")
            return 1

        # 決定したファイルパスを使ってデータを読み込む
        cluster_df_raw = pd.read_csv(
            cluster_data_file, sep='\t', header=None,
            names=CLUSTER_COLUMNS, dtype=str, engine='python', skipinitialspace=True
        )

        if cluster_df_raw.empty:
            # メッセージ内のファイル名を、実際に使用したファイル名に動的に変更
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
        
        essential_cols = ['tw', 'no', 'n', 'm', 'sliprate']
        cluster_df.dropna(subset=essential_cols, inplace=True)
        if cluster_df.empty:
            # メッセージ内のファイル名を、実際に使用したファイル名に動的に変更
            sys.stderr.write(f"情報: 必須列のNaN除去後、{cluster_data_file} に有効なデータがありません。\n")
            pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
                 wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
            )
            return 0

        cluster_max_original_tw = cluster_df['tw'].max()
        if pd.isna(cluster_max_original_tw) or cluster_max_original_tw <= 0:
            # メッセージ内のファイル名を、実際に使用したファイル名に動的に変更
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

    # 2. 実時間への変換
    time_per_original_tw = t_max / cluster_max_original_tw
    cluster_df['calc_real_time'] = (cluster_df['tw'].astype(float) - 0.5) * time_per_original_tw
    
    # 3. 初期時間フィルタリング
    overall_start_time = new_time_windows_df['start_time'].min()
    overall_end_time = new_time_windows_df['end_time'].max()
    
    cluster_df.dropna(subset=['calc_real_time'], inplace=True)

    initial_filtered_df = cluster_df[
        (cluster_df['calc_real_time'] >= overall_start_time) &
        (cluster_df['calc_real_time'] < overall_end_time) 
    ].copy() 

    if initial_filtered_df.empty:
        # メッセージ内のファイル名を、実際に使用したファイル名に動的に変更
        sys.stderr.write(f"情報: 指定された全体時間範囲に一致するデータが {cluster_data_file} にありません。\n")
        pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(
             wd + "clusteringSnap2_v2.dat", sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN'
        )
        return 0
        
    final_output_rows_list = []

    # 4. 新しい時間窓ごとの処理
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

        selected_rows_for_this_new_tw = []

        for grid_keys, group_at_grid in data_in_current_new_window.groupby(['n', 'm']):
            valid_no_series = group_at_grid[group_at_grid['no'] != 0]['no']
            if valid_no_series.empty:
                continue

            no_counts = Counter(valid_no_series.dropna().astype(int))
            if not no_counts:
                continue
            
            most_common_nos_with_freq = no_counts.most_common()
            
            if not most_common_nos_with_freq:
                continue

            max_freq = most_common_nos_with_freq[0][1]
            
            representative_no_candidates = [
                item[0] for item in most_common_nos_with_freq if item[1] == max_freq
            ]
            if not representative_no_candidates:
                continue
            representative_no = min(representative_no_candidates)

            data_with_representative_no = group_at_grid[group_at_grid['no'] == representative_no]
            
            # --- ▼▼▼ ここからが変更箇所です ▼▼▼ ---
            # sliprateでのフィルタリングを行わず、代表noを持つ全てのデータを追加します
            if not data_with_representative_no.empty:
                for _, row in data_with_representative_no.iterrows():
                    selected_rows_for_this_new_tw.append(row.copy())
            # --- ▲▲▲ ここまでが変更箇所です ▲▲▲ ---

        if selected_rows_for_this_new_tw:
            temp_df = pd.DataFrame(selected_rows_for_this_new_tw)
            temp_df['tw'] = current_new_tw_id
            final_output_rows_list.append(temp_df)

    if not final_output_rows_list:
        sys.stderr.write("情報: 全ての処理後、有効なデータが残りませんでした。\n")
        final_result_df = pd.DataFrame(columns=CLUSTER_COLUMNS)
    else:
        final_result_df = pd.concat(final_output_rows_list, ignore_index=True)
        if 'calc_real_time' in final_result_df.columns:
            final_result_df.drop(columns=['calc_real_time'], inplace=True)
        
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
    print("swap_mmcl5_revised.py (v4 - sliprateフィルタなし) を実行中...")
    exit_code = revised_swap_process()
    if exit_code == 0:
        print('swap_mmcl5_revised.py: データ処理完了。')
    else:
        print('swap_mmcl5_revised.py: エラーが発生しました。')
    sys.exit(exit_code)
