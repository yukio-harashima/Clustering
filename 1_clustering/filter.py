#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import Counter
import sys # sysモジュールをインポート

# --- 設定値 ---
# clusteringSnap2c.dat の列名 (30列を想定)
INPUT_COLUMNS = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no' # 30列目が 'no' (クラスター番号)
]

# tw の1ステップあたりの時間間隔 (秒)
TIME_INTERVAL_SEC = 0.5

# クラスター番号ごとに抽出したい開始秒数と終了秒数を指定
# 例: {クラスター番号: (開始秒数, 終了秒数)}
CLUSTER_EXTRACTION_PARAMS = {
    1: (0.0, 10.0),  # クラスター1: 9秒 から 11秒 まで
    2: (6.0, 10.5), # クラスター2: 14.5秒 から 15.5秒 まで
    3: (10.0, 18.0), 
    4: (16.0, 30.0), 
    # 1: (0.0, 10.0),  # クラスター1: 9秒 から 11秒 まで
    # 2: (0.0, 10.0), # クラスター2: 14.5秒 から 15.5秒 まで
    # 3: (10.0, 30.0), 
    # 4: (10.0, 30.0), 


    # 必要に応じて他のクラスターの設定を追加
}

INPUT_FILE = "clusteringSnap2c.dat"
OUTPUT_FILE = "clarea.dat"
# --- 設定値ここまで ---

def load_and_preprocess_data(filepath, columns, time_interval):
    """データを読み込み、実時間を計算する"""
    try:
        df = pd.read_csv(filepath, sep='\t', header=None, names=columns, dtype=str, engine='python', skipinitialspace=True)
    except FileNotFoundError:
        sys.stderr.write(f"エラー: 入力ファイル '{filepath}' が見つかりません。\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"エラー: ファイル '{filepath}' の読み込み中にエラーが発生しました: {e}\n")
        sys.exit(1)

    for col in df.columns:
        if col not in ['clla']: 
            df[col] = pd.to_numeric(df[col], errors='coerce')

    if 'tw' not in df.columns or 'no' not in df.columns:
        sys.stderr.write("エラー: 'tw' または 'no' 列が入力ファイルに見つかりません。\n")
        sys.exit(1)

    df.dropna(subset=['tw', 'no', 'n', 'm'], inplace=True) 
    df['actual_time'] = df['tw'] * time_interval
    df['no'] = df['no'].astype(int) 
    df['n'] = df['n'].astype(int)
    df['m'] = df['m'].astype(int)
    return df

def extract_data_by_cluster(df, cluster_params):
    """クラスターごとに指定された時間範囲のデータを抽出する"""
    all_extracted_data = []
    for cluster_no, (time_start, time_end) in cluster_params.items(): # パラメータのアンパックを変更
        cluster_specific_data = df[df['no'] == cluster_no].copy()
        
        extracted = cluster_specific_data[
            (cluster_specific_data['actual_time'] >= time_start) &
            (cluster_specific_data['actual_time'] <= time_end)
        ]
        if not extracted.empty:
            # 抽出されたデータに、抽出条件となったクラスター番号と時間範囲を記録しておく (重複処理で使用)
            extracted['extraction_cluster_no'] = cluster_no
            extracted['extraction_time_start'] = time_start
            extracted['extraction_time_end'] = time_end
            all_extracted_data.append(extracted)

    if not all_extracted_data:
        return pd.DataFrame(columns=list(df.columns) + ['extraction_cluster_no', 'extraction_time_start', 'extraction_time_end'])
    
    return pd.concat(all_extracted_data, ignore_index=True)

def filter_by_highest_frequency_per_second_per_knot(df):
    """
    各ノット(n, m)において、最も秒間あたりの出現頻度が高い 'no' を持つデータのみを残す。
    秒間頻度が同じ場合は、no の値が小さい方を優先する。
    """
    if df.empty:
        return df

    final_df_list = []
    for knot_coords, group_data_for_knot in df.groupby(['n', 'm']):
        if group_data_for_knot.empty:
            continue
        
        cluster_info = [] # (秒間頻度, extraction_cluster_no, データフレーム) を格納

        # このノットに含まれる可能性のある各抽出元クラスターについて処理
        for extraction_no, data_from_extraction_no in group_data_for_knot.groupby('extraction_cluster_no'):
            if data_from_extraction_no.empty:
                continue

            # この抽出元クラスターの出現頻度
            occurrence_count = len(data_from_extraction_no[data_from_extraction_no['no'] == extraction_no])
            
            # この抽出元クラスターの抽出時間幅
            time_start = data_from_extraction_no['extraction_time_start'].iloc[0]
            time_end = data_from_extraction_no['extraction_time_end'].iloc[0]
            duration = time_end - time_start
            if duration <= 0: # 期間が0または負の場合は頻度を0とするか、エラー処理
                freq_per_second = 0
            else:
                freq_per_second = occurrence_count / duration
            
            # このノット内で、現在の extraction_no に実際に属するデータのみを対象とする
            # (filter_by_most_frequent_no_per_knot からの変更点として、
            #  ここでは extraction_cluster_no のデータのみを考える)
            actual_data_for_this_extraction_no = data_from_extraction_no[data_from_extraction_no['no'] == extraction_no]
            
            if not actual_data_for_this_extraction_no.empty:
                 cluster_info.append({'freq_per_sec': freq_per_second, 
                                      'cluster_no': extraction_no, 
                                      'data': actual_data_for_this_extraction_no})
        
        if not cluster_info:
            continue
            
        # 秒間頻度が最も高く、それが同じ場合はクラスター番号が小さいものを選択
        cluster_info_sorted = sorted(cluster_info, key=lambda x: (-x['freq_per_sec'], x['cluster_no']))
        
        best_cluster_data = cluster_info_sorted[0]['data']
        final_df_list.append(best_cluster_data)
        
    if not final_df_list:
        return pd.DataFrame(columns=df.columns) # 元の列構造で空のDataFrameを返す
        
    return pd.concat(final_df_list, ignore_index=True)

def save_data(df, filepath, columns_order):
    """結果を指定された列順でタブ区切りで保存する"""
    if df.empty:
        open(filepath, 'w').close()
        sys.stderr.write(f"情報: 出力データが空のため、空のファイル '{filepath}' を作成しました。\n")
        return

    # 不要になった一時列を削除
    cols_to_drop = ['extraction_cluster_no', 'extraction_time_start', 'extraction_time_end', 'actual_time']
    for col in cols_to_drop:
        if col in df.columns:
            df = df.drop(columns=[col])
            
    # 出力する列が存在することを確認し、なければNaNで埋める
    # また、columns_order は INPUT_COLUMNS を使う (元の30列に戻す)
    final_output_columns = INPUT_COLUMNS 
    for col in final_output_columns:
        if col not in df.columns:
            df[col] = np.nan
            
    df_to_save = df[final_output_columns] 
    
    try:
        df_to_save.to_csv(filepath, sep='\t', header=False, index=False, float_format='%.6f', na_rep='NaN')
    except Exception as e:
        sys.stderr.write(f"エラー: ファイル '{filepath}' の書き込み中にエラーが発生しました: {e}\n")
        sys.exit(1)

def main():
    """メイン処理"""
    df_raw = load_and_preprocess_data(INPUT_FILE, INPUT_COLUMNS, TIME_INTERVAL_SEC)
    
    if df_raw.empty:
        sys.stderr.write("情報: 入力データが空か、前処理後にデータが残りませんでした。\n")
        save_data(pd.DataFrame(columns=INPUT_COLUMNS), OUTPUT_FILE, INPUT_COLUMNS) 
        return

    extracted_df = extract_data_by_cluster(df_raw, CLUSTER_EXTRACTION_PARAMS)

    if extracted_df.empty:
        sys.stderr.write("情報: 指定されたパラメータで抽出されたデータはありませんでした。\n")
        save_data(pd.DataFrame(columns=INPUT_COLUMNS), OUTPUT_FILE, INPUT_COLUMNS) 
        return
        
    final_result_df = filter_by_highest_frequency_per_second_per_knot(extracted_df)

    save_data(final_result_df, OUTPUT_FILE, INPUT_COLUMNS) # 保存時の列順は INPUT_COLUMNS を使用
    sys.stderr.write(f"処理が完了しました。結果は '{OUTPUT_FILE}' に保存されました。\n")

if __name__ == '__main__':
    main()
