#!/usr/bin/env python3
# clave_v2.py

'''
clusteringSnap2_mod.dat (calc_fault_v2.py) ファイルを読み込み、30列目の no の値でデータをグループ化します。
その後、各グループ内で指定された列の平均値を計算します。特に角度データ（走向、傾斜、すべり角）については、ベクトル平均法を用いて適切に平均値を算出します。
最後に、各グループの平均断層パラメータ2組に対して selectplane 関数を適用し、モデル断層面に近い方を主断層面として整理し、結果を指定されたフォーマットで clave_mod.dat に出力します。
'''

import numpy as np
import pandas as pd
import math

# --- 使用するデータ直接書き込み (モデル断層面のパラメータ) ---
# これらは selectplane 関数で参照面として使用されます。
MODEL_STRIKE_DEG = 179.0  # モデル(参照)断層面の走向 (度)
MODEL_DIP_DEG = 55.0    # モデル(参照)断層面の傾斜 (度)

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
        # print(f"警告: グループ内で角度が均等に分散しているため、平均角度は未定義です。NaNを返します。データ: {angles_deg_list}")
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
    # 傾斜の平均は0-90度の範囲が物理的に期待される。
    # ベクトル平均の結果がこの範囲に収まらない場合、解釈に注意が必要。
    # ここではベクトル平均を計算し、結果をそのまま返す。
    # 必要であれば、結果を0-90度にクリップするなどの後処理を検討。
    avg_angle = angle_average_degrees(series.tolist(), result_range_0_360=True) # まず0-360で計算
    if pd.notna(avg_angle):
        # 0-90度への変換ロジック (例: 350度なら10度、170度なら10度など)
        # これは単純なベクトル平均の解釈とは異なるため、今回は行わない。
        # 物理的な傾斜として、0-90の範囲に収まるように値を調整するか、
        # あるいは単純平均の方が適切かを検討する余地がある。
        # 今回は指示通りベクトル平均の結果を返す。
        if avg_angle > 90 and avg_angle <= 180:
            avg_angle = 180 - avg_angle
        elif avg_angle > 180 and avg_angle <= 270:
             avg_angle = avg_angle - 180
        elif avg_angle > 270 and avg_angle < 360:
            avg_angle = 360 - avg_angle
        # 0度または90度に近い場合はそのまま
    return avg_angle


def rake_ave(series):
    """すべり角 (rake0, rake1) 列用の角度平均計算 (-180～180度範囲)"""
    return angle_average_degrees(series.tolist(), result_range_0_360=False, normalize_to_plus_minus_180=True)

# --- faultnormalvec 関数の定義 (pyutils_2024Noto.py より) ---
def faultnormalvec(stk, dip):
    """
    断層面の走向 (stk) と傾斜 (dip) から、その断層面の法線ベクトルを計算します。
    走向と傾斜は度単位で入力します。
    法線ベクトルは [東成分, 北成分, 下成分] の順で返されます。
    """
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip)) # North component
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip)) # East component
    nd = -np.cos(np.deg2rad(dip))                         # Down component
    return np.array([ne, nn, nd]) # East, North, Down

# --- selectplane 関数の定義 (pyutils_2024Noto.py より) ---
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    """
    2つの候補となる断層面のうち、モデル断層面により「近い」方を優先的な断層面として選択します。
    「近さ」は法線ベクトルの内積の絶対値で判断します。
    選択された断層面の走向、傾斜、すべり角を返します。
    """
    # 入力値がNaNの場合のフォールバック
    if any(pd.isna(x) for x in [stk0, dip0, stk1, dip1]):
        # print(f"警告: selectplaneへの入力にNaNが含まれるため、デフォルトで候補0を返します。")
        # NaNが含まれる場合、比較ができないため、ここでは仮に候補0のパラメータを優先面とする
        # (またはエラーを出すか、NaNを返すなどの処理も考えられる)
        return stk0, dip0, rake0 # 優先面
        
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    
    tmp0 = np.inner(vecmodelplane, vecplane0)
    tmp1 = np.inner(vecmodelplane, vecplane1)
    
    selected_stk, selected_dip, selected_rake = np.nan, np.nan, np.nan
    other_stk, other_dip, other_rake = np.nan, np.nan, np.nan
    
    if abs(tmp0) >= abs(tmp1): # tmp0が大きいか等しい場合は候補0を選択
        selected_stk, selected_dip, selected_rake = stk0, dip0, rake0
        other_stk, other_dip, other_rake = stk1, dip1, rake1
    else: # abs(tmp0) < abs(tmp1) の場合 (tmp1が大きい場合) は候補1を選択
        selected_stk, selected_dip, selected_rake = stk1, dip1, rake1
        other_stk, other_dip, other_rake = stk0, dip0, rake0
        
    # この関数は元々選択された1組だけを返していたが、
    # 今回のclave_v2.pyの要件では、選択された方とされなかった方の両方が必要。
    # なので、(優先面stk, dip, rake), (非優先面stk, dip, rake) の形で返す。
    return selected_stk, selected_dip, selected_rake, other_stk, other_dip, other_rake


# --- メイン処理 ---
def main():
    input_filename = "clusteringSnap2_mod.dat"
    output_filename = "clave_mod.dat"

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

    # --- noの値毎にデータのグループ化と平均値計算 ---
    print("グループごとの平均値計算を開始します...")
    
    # 平均化する列と適用する関数を定義
    agg_funcs = {
        'sliprate': simple_ave,
        'Mrr': simple_ave,
        'Mss': simple_ave,
        'Mee': simple_ave,
        'Mrs': simple_ave,
        'Mre': simple_ave,
        'Mse': simple_ave,
        'stk0': stk_ave,
        'stk1': stk_ave,
        'dip0': dip_ave,
        'dip1': dip_ave,
        'rake0': rake_ave,
        'rake1': rake_ave
    }
    
    # 'no' 列が存在するか確認
    if 'no' not in df.columns:
        print(f"エラー: 入力ファイルに 'no' 列が見つかりません。")
        return

    # 'no' 列を数値型に変換 (エラー時はNaNにし、その後整数型に変換)
    df['no'] = pd.to_numeric(df['no'], errors='coerce')
    df.dropna(subset=['no'], inplace=True) # 'no'がNaNの行は除外
    df['no'] = df['no'].astype(int)

    # 各列を数値型に変換 (エラー時はNaN)
    cols_to_convert = list(agg_funcs.keys())
    for col in cols_to_convert:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        else:
            print(f"警告: 列 '{col}' が入力ファイルに存在しません。スキップします。")
            # agg_funcs から存在しない列を削除
            del agg_funcs[col]
            
    # グループ化して集計
    # reset_index() を使って 'no' を列に戻す
    df_averaged = df.groupby('no').agg(agg_funcs).reset_index()
    
    print("グループごとの平均値計算が完了しました。")

    # --- selectplane の処理を平均化されたデータに適用 ---
    print("平均化された断層パラメータに対する selectplane 処理を開始します...")
    
    results_after_selectplane = []
    for index, row in df_averaged.iterrows():
        # selectplane は (優先stk,優先dip,優先rake, 非優先stk,非優先dip,非優先rake) を返す
        pref_s, pref_d, pref_r, other_s, other_d, other_r = selectplane(
            MODEL_STRIKE_DEG, MODEL_DIP_DEG,
            row['stk0'], row['dip0'], row['rake0'],
            row['stk1'], row['dip1'], row['rake1']
        )
        
        # 結果を辞書として格納
        output_row = {
            'no': row['no'],
            'sliprate': row['sliprate'],
            'Mrr': row['Mrr'],
            'Mss': row['Mss'],
            'Mee': row['Mee'],
            'Mrs': row['Mrs'],
            'Mre': row['Mre'],
            'Mse': row['Mse'],
            'stk0': pref_s,    # モデルに近い方
            'dip0': pref_d,
            'rake0': pref_r,
            'stk1': other_s,   # そうでない方
            'dip1': other_d,
            'rake1': other_r
        }
        results_after_selectplane.append(output_row)
        
    df_final = pd.DataFrame(results_after_selectplane)
    print("selectplane 処理が完了しました。")

    # --- データ出力 ---
    # 指定された列順序
    output_columns_order = [
        'no', 'sliprate', 'Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse', 
        'stk0', 'dip0', 'rake0', 'stk1', 'dip1', 'rake1'
    ]
    # DataFrameの列をこの順序に並べ替える (存在しない列があればエラーになるので注意)
    # 念のため、実際に存在する列のみで順序を構成
    final_output_columns = [col for col in output_columns_order if col in df_final.columns]

    try:
        df_final[final_output_columns].to_csv(output_filename, sep='\t', index=False, header=True, float_format='%.6f')
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()
