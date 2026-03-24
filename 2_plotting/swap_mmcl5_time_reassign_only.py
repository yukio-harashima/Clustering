#!/usr/bin/env python3

'''
swap_mmcl5_time_reassign_only.py:
入力されたクラスタデータに対して、物理的な実時間を計算し、
'new_time_windows.dat' に基づいて新しい時間窓IDを再割り当てします。
代表データの抽出（フィルタリング）処理は行いません。
**改良点:** clusteringSnap2c.dat が存在しない場合、代わりに clusteringSnap2.dat を使用するように変更しました。
'''


import pandas as pd
import numpy as np
import sys
import os

# clusteringSnap2c.dat および出力ファイルの列定義
CLUSTER_COLUMNS = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]

def time_reassign_process():
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
            sys.stderr.write(f"エラー: {cluster_file_c} も {cluster_file_plain} も見つかりません。\n")
            return 1

        # 決定したファイルパスを使ってデータを読み込む
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
        
        essential_cols = ['tw', 'no', 'n', 'm']
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
            sys.stderr.write("情報: new_time_windows.dat が空です。時間再割り当てができません。\n")
            # 元のデータをそのまま出力することも考えられるが、ここでは空ファイルを出力
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
    
    # --- ▼▼▼ ここからが変更箇所です ▼▼▼ ---

    # 3. 新しい時間窓IDの再割り当て
    
    # 各実時間がどの新しい時間窓に属するかを判定するための関数
    def get_new_tw_id(real_time):
        for _, window in new_time_windows_df.iterrows():
            if window['start_time'] <= real_time < window['end_time']:
                return window['new_tw_id']
        return np.nan

    # 全てのデータ行に対して新しい時間窓IDを適用
    cluster_df['new_tw'] = cluster_df['calc_real_time'].apply(get_new_tw_id)

    # どの新しい時間窓にも属さないデータを削除
    final_result_df = cluster_df.dropna(subset=['new_tw']).copy()
    
    if final_result_df.empty:
        sys.stderr.write(f"情報: 指定された全体時間範囲に一致するデータが {cluster_data_file} にありません。\n")
        final_result_df = pd.DataFrame(columns=CLUSTER_COLUMNS)
    else:
        # tw列を新しいIDで上書きし、不要な列を削除
        final_result_df['tw'] = final_result_df['new_tw'].astype(int)
        final_result_df.drop(columns=['calc_real_time', 'new_tw'], inplace=True)
    
        # 出力列の整備
        for col in CLUSTER_COLUMNS:
            if col not in final_result_df.columns:
                final_result_df[col] = np.nan
        final_result_df = final_result_df[CLUSTER_COLUMNS]

    # --- ▲▲▲ ここまでが変更箇所です ▲▲▲ ---

    # 4. 結果の出力
    output_filename = wd + "clusteringSnap2_v2.dat"
    try:
        final_result_df.to_csv(output_filename, sep='\t', index=False, header=False, float_format='%.6f', na_rep='NaN')
        print(f"情報: {output_filename} が正常に生成されました。行数: {len(final_result_df)}")
    except Exception as e:
        sys.stderr.write(f"エラー: {output_filename} の書き込み中にエラー - {e}\n")
        return 1
        
    return 0

if __name__ == '__main__':
    print("swap_mmcl5_time_reassign_only.py を実行中...")
    exit_code = time_reassign_process()
    if exit_code == 0:
        print('swap_mmcl5_time_reassign_only.py: データ処理完了。')
    else:
        print('swap_mmcl5_time_reassign_only.py: エラーが発生しました。')
    sys.exit(exit_code)
