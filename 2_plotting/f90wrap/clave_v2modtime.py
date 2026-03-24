#!/usr/bin/env python3
# clave_v2.py

'''
clusteringSnap2_mod.dat(calc_faul0t_v2.py) ファイルを読み込み、指定された時間区間内のデータを抽出します。
その後、30列目の no の値でデータをグループ化し、各グループ内で指定された列の平均値を計算します。
特に角度データ（走向、傾斜、すべり角）については、ベクトル平均法を用いて適切に平均値を算出します。
最後に、各グループの平均断層パラメータ2組に対して selectplane 関数を適用し、
モデル断層面に近い方を主断層面として整理し、結果を指定されたフォーマットで動的なファイル名で出力します。
'''

import numpy as np
import pandas as pd
import math

# --- 使用するデータ直接書き込み (モデル断層面のパラメータ) ---
# これらは selectplane 関数で参照面として使用されます。
MODEL_STRIKE_DEG = 179.0  # モデル(参照)断層面の走向 (度)
MODEL_DIP_DEG = 55.0    # モデル(参照)断層面の傾斜 (度)
DT_INTERVAL_SEC = 0.5     # 時間窓(tw)あたりの時間間隔 (秒)
START_TIME_SEC = 31     # 切り出す開始時間 (秒)。Noneの場合はデータの最初から
END_TIME_SEC = 46       # 切り出す終了時間 (秒)。Noneの場合はデータの最後まで

# --- 平均値計算用の関数定義 ---

def simple_ave(series):
    """Pandas Series の単純平均を計算します。NaNは無視されます。"""
    return series.mean()

def angle_average_degrees(angles_deg_list, result_range_0_360=True, normalize_to_plus_minus_180=False):
    """
    角度のリスト（度単位）からベクトル平均を計算します。
    """
    # NaNやNoneをフィルタリング
    valid_angles_deg = [a for a in angles_deg_list if pd.notna(a)]
    
    if not valid_angles_deg: # 有効な角度がない場合はNaNを返す
        return np.nan
        
    angles_rad = np.deg2rad(valid_angles_deg) # 度をラジアンに変換
    
    sum_cos = np.sum(np.cos(angles_rad))
    sum_sin = np.sum(np.sin(angles_rad))
    
    num_angles = len(valid_angles_deg)
    
    mean_cos = sum_cos / num_angles
    mean_sin = sum_sin / num_angles
    
    # 平均ベクトルが非常に小さい場合（角度が均等に散らばっているなど）
    if np.isclose(mean_cos, 0.0) and np.isclose(mean_sin, 0.0):
        return np.nan 

    mean_angle_rad = math.atan2(mean_sin, mean_cos)
    mean_angle_deg = math.degrees(mean_angle_rad)
    
    if normalize_to_plus_minus_180:
        # 結果を-180～180度の範囲に正規化
        mean_angle_deg = (mean_angle_deg + 180) % 360 - 180
        if np.isclose(mean_angle_deg, -180.0): # -180.0 を 180.0 にする (慣習による)
            mean_angle_deg = 180.0
    elif result_range_0_360:
        # 結果を0～360度の範囲に正規化
        if mean_angle_deg < 0:
            mean_angle_deg += 360.0
            
    return mean_angle_deg

def stk_ave(series):
    """走向 (stk0, stk1) 列用の角度平均計算 (0-360度範囲)"""
    return angle_average_degrees(series.tolist(), result_range_0_360=True)

def dip_ave(series):
    """傾斜 (dip0, dip1) 列用の角度平均計算"""
    avg_angle = angle_average_degrees(series.tolist(), result_range_0_360=True) 
    if pd.notna(avg_angle):
        if avg_angle > 90 and avg_angle <= 180:
            avg_angle = 180 - avg_angle
        elif avg_angle > 180 and avg_angle <= 270:
             avg_angle = avg_angle - 180
        elif avg_angle > 270 and avg_angle < 360:
            avg_angle = 360 - avg_angle
    return avg_angle


def rake_ave(series):
    """すべり角 (rake0, rake1) 列用の角度平均計算 (-180～180度範囲)"""
    return angle_average_degrees(series.tolist(), result_range_0_360=False, normalize_to_plus_minus_180=True)

# --- faultnormalvec 関数の定義 ---
def faultnormalvec(stk, dip):
    """
    断層面の走向 (stk) と傾斜 (dip) から、その断層面の法線ベクトルを計算します。
    """
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip)) # North component
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip)) # East component
    nd = -np.cos(np.deg2rad(dip))                         # Down component
    return np.array([ne, nn, nd]) # East, North, Down

# --- selectplane 関数の定義 ---
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    """
    2つの候補となる断層面のうち、モデル断層面により「近い」方を優先的な断層面として選択します。
    """
    if any(pd.isna(x) for x in [stk0, dip0, stk1, dip1]):
        return stk0, dip0, rake0, stk1, dip1, rake1
        
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    
    tmp0 = np.inner(vecmodelplane, vecplane0)
    tmp1 = np.inner(vecmodelplane, vecplane1)
    
    if abs(tmp0) >= abs(tmp1): 
        return stk0, dip0, rake0, stk1, dip1, rake1
    else: 
        return stk1, dip1, rake1, stk0, dip0, rake0

# --- メイン処理 ---
def main():
    input_filename = "clusteringSnap2_mod.dat"
    base_output_filename = "clave_mod"

    # --- データ読み込み ---
    try:
        df = pd.read_csv(
            input_filename,
            delimiter='\t',
            header=0, 
            skipinitialspace=True
        )
        print(f"'{input_filename}' から {len(df)} 行のデータを読み込みました。")
    except FileNotFoundError:
        print(f"エラー: 入力ファイル '{input_filename}' が見つかりません。")
        return
    except Exception as e:
        print(f"エラー: '{input_filename}' の読み込み中にエラーが発生しました: {e}")
        return

    # --- 時間区間の指定 (定数から読み込み) ---
    if 'tw' not in df.columns:
        print("エラー: 入力ファイルに 'tw' (時間窓) 列が見つかりません。")
        return
        
    tw_min_data, tw_max_data = df['tw'].min(), df['tw'].max()
    
    # 定数から開始・終了時間窓を計算
    if START_TIME_SEC is not None:
        # 秒数をtwインデックスに変換 (切り上げ)
        tw_start = math.ceil(START_TIME_SEC / DT_INTERVAL_SEC)
    else:
        tw_start = tw_min_data

    if END_TIME_SEC is not None:
        # 秒数をtwインデックスに変換 (切り捨て)
        tw_end = math.floor(END_TIME_SEC / DT_INTERVAL_SEC)
    else:
        tw_end = tw_max_data

    if tw_start > tw_end:
        print(f"エラー: 設定された開始時間({START_TIME_SEC}秒)が終了時間({END_TIME_SEC}秒)より後です。処理を中断します。")
        return

    print("\n--- 時間区間の指定 ---")
    print(f"設定された時間間隔(dt): {DT_INTERVAL_SEC} 秒")
    print(f"設定された処理時間: {START_TIME_SEC if START_TIME_SEC is not None else '最初'}秒 から {END_TIME_SEC if END_TIME_SEC is not None else '最後'}秒 まで")
    print(f"対応する時間窓(tw)の範囲: {tw_start} から {tw_end}")
    
    output_filename = f"{base_output_filename}_tw{tw_start}-{tw_end}.dat"

    # --- 時間によるデータフィルタリング ---
    df_time_filtered = df[(df['tw'] >= tw_start) & (df['tw'] <= tw_end)].copy()
    if df_time_filtered.empty:
        print(f"指定された時間区間 ({tw_start} - {tw_end}) にデータが存在しません。")
        return
    print(f"時間フィルタリング後のデータ数: {len(df_time_filtered)} 行")


    # --- noの値毎にデータのグループ化と平均値計算 ---
    print("\nグループごとの平均値計算を開始します...")
    
    agg_funcs = {
        'sliprate': simple_ave, 'Mrr': simple_ave, 'Mss': simple_ave,
        'Mee': simple_ave, 'Mrs': simple_ave, 'Mre': simple_ave, 'Mse': simple_ave,
        'stk0': stk_ave, 'stk1': stk_ave, 'dip0': dip_ave, 'dip1': dip_ave,
        'rake0': rake_ave, 'rake1': rake_ave
    }
    
    if 'no' not in df_time_filtered.columns:
        print(f"エラー: 入力ファイルに 'no' 列が見つかりません。")
        return

    df_time_filtered['no'] = pd.to_numeric(df_time_filtered['no'], errors='coerce')
    df_time_filtered.dropna(subset=['no'], inplace=True)
    df_time_filtered['no'] = df_time_filtered['no'].astype(int)

    cols_to_convert = list(agg_funcs.keys())
    for col in cols_to_convert:
        if col in df_time_filtered.columns:
            df_time_filtered[col] = pd.to_numeric(df_time_filtered[col], errors='coerce')
        else:
            print(f"警告: 列 '{col}' が入力ファイルに存在しません。スキップします。")
            del agg_funcs[col]
            
    df_averaged = df_time_filtered.groupby('no').agg(agg_funcs).reset_index()
    
    print("グループごとの平均値計算が完了しました。")

    # --- selectplane の処理を平均化されたデータに適用 ---
    print("平均化された断層パラメータに対する selectplane 処理を開始します...")
    
    results_after_selectplane = []
    for index, row in df_averaged.iterrows():
        pref_s, pref_d, pref_r, other_s, other_d, other_r = selectplane(
            MODEL_STRIKE_DEG, MODEL_DIP_DEG,
            row['stk0'], row['dip0'], row['rake0'],
            row['stk1'], row['dip1'], row['rake1']
        )
        
        output_row = {
            'no': row['no'], 'sliprate': row['sliprate'],
            'Mrr': row['Mrr'], 'Mss': row['Mss'], 'Mee': row['Mee'],
            'Mrs': row['Mrs'], 'Mre': row['Mre'], 'Mse': row['Mse'],
            'stk0': pref_s, 'dip0': pref_d, 'rake0': pref_r,
            'stk1': other_s, 'dip1': other_d, 'rake1': other_r
        }
        results_after_selectplane.append(output_row)
        
    df_final = pd.DataFrame(results_after_selectplane)
    print("selectplane 処理が完了しました。")

    # --- データ出力 ---
    output_columns_order = [
        'no', 'sliprate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse', 
        'stk0', 'dip0', 'rake0', 'stk1', 'dip1', 'rake1'
    ]
    final_output_columns = [col for col in output_columns_order if col in df_final.columns]

    try:
        df_final[final_output_columns].to_csv(output_filename, sep='\t', index=False, header=True, float_format='%.6f')
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()
