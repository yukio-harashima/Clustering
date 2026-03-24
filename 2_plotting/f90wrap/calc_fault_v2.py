#!/usr/bin/env python3
# calc_fault_v2.py
# clusteringSnap2c.dat(dataeditor.py)を読み込み設定した節面に近い面に数値を入れ替える
import numpy as np
import pandas as pd
import math # faultnormalvec内でnp.deg2radを使うため、mathは直接不要かも

# --- 使用するデータ直接書き込み (モデル断層面のパラメータ) ---
# これらは selectplane 関数で参照面として使用されます。
# 必要に応じてこれらの値を変更してください。
MODEL_STRIKE_DEG = 179.0  # モデル(参照)断層面の走向 (度)
MODEL_DIP_DEG = 55.0    # モデル(参照)断層面の傾斜 (度)



# --- 入力データの列名定義 ---
INPUT_COLUMNS_30 = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'stk0', 'stk1', 'dip0', 'dip1', 
    'rake0', 'rake1', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]



# ------ for calc_fault.py ------
'''
# 目的: 
    2つの候補となる断層面（(stk0, dip0, rake0) と (stk1, dip1, rake1)）のうち、モデルとして与えられた断層面（modelstk, modeldip）により「近い」方を優先的な断層面として選択します。
# 引数:
    modelstk, modeldip: モデル（参照）断層面の走向と傾斜。
    stk0, dip0, rake0: 候補1の断層面の走向、傾斜、すべり角。
    stk1, dip1, rake1: 候補2の断層面の走向、傾斜、すべり角。
# 処理ロジック:
    faultnormalvec 関数（後述）を使って、モデル断層面、候補1、候補2のそれぞれの法線ベクトルを計算します。
    モデル断層面の法線ベクトルと、候補1および候補2の法線ベクトルとの内積 (np.inner) をそれぞれ計算します (tmp0, tmp1)。
    内積の絶対値が大きい方が、法線ベクトルの向きがより近い（またはより平行）と判断し、その候補断層面のパラメータ（走向、傾斜、すべり角）を選択します。
    内積の絶対値が等しい場合は、候補0を優先します。
# 返り値: 
    選択された断層面の走向 (stk_s)、傾斜 (dip_s)、すべり角 (rake_s)。
'''

# select a preferred fault plane
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    tmp0 = np.inner(vecmodelplane, vecplane0)
    tmp1 = np.inner(vecmodelplane, vecplane1)
    if abs(tmp0) > abs(tmp1):
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    elif abs(tmp0) < abs(tmp1):
        stk_s = stk1
        dip_s = dip1
        rake_s = rake1
    else:
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    return stk_s, dip_s, rake_s



# ------ for calc_fault.py ------
'''
# 目的:
    断層面の走向 (stk) と傾斜 (dip) から、その断層面の法線ベクトルを計算します。
# 引数:
    stk: 走向 (度)。
    dip: 傾斜 (度)。
# 処理ロジック:
    走向と傾斜をラジアンに変換します (np.deg2rad)。
    地震学で一般的に用いられる規約（例: Aki & Richards）に基づいて、法線ベクトルの北(N)、東(E)、下(D)成分を計算します。
        nn = -sin(strike) * sin(dip) (北成分)
        ne =  cos(strike) * sin(dip) (東成分)
        nd = -cos(dip) (下成分、法線が下向きを指す場合)
    計算された3成分をNumPy配列として返します。注意: コード内では np.array([ne, nn, nd]) となっており、東成分が先に来ています。これは使用する座標系の規約に依存します。
# 返り値: 
    法線ベクトルの3成分を格納したNumPy配列。
'''

# define a fault-normal vector
def faultnormalvec(stk, dip):
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    nd = -np.cos(np.deg2rad(dip))
    return np.array([ne, nn, nd])




# ------ for calc_fault.py ------
'''
# 目的: 6成分のモーメントテンソル（USE規約、例: GCMTのpsmecaフォーマット）を入力として、P軸、N(B)軸、T軸のプランジ（傾伏角）とアジマス（方位角）を計算
# 引数:
    focmec: 6成分のモーメントテンソル (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp の順で、USE座標系: Up-South-East)。
# 処理ロジック:
    入力された6成分とスカラーモーメント（ここでは26を仮定、実際には単位やスケールに依存）から obspy.imaging.beachball.MomentTensor オブジェクトを作成します。
    mt.mt (3x3のモーメントテンソル行列) に対して np.linalg.eigh を使い、固有値 (d) と固有ベクトル (v) を計算します。eigh は対称（エルミート）行列用で、固有値は昇順にソートされて返されます。
    固有ベクトル v の各列がP, N, T軸の方向ベクトルに対応します（固有値の順序に依存）。ObsPyの規約では、固有ベクトルは (Up, South, East) 座標系で得られます。
        pl = np.arcsin(-v[0]): プランジを計算。v[0] はUp成分なので、-v[0] はDown成分に相当し、arcsin で水平面からの角度を求めます。
        az = np.arctan2(v[2], -v[1]): アジマスを計算。v[2] はEast成分、v[1] はSouth成分なので、-v[1] はNorth成分に相当し、arctan2(East, North) で北からの角度を求めます。
    プランジが下向き正 (0-90度) になるように、またアジマスが (0-360度) になるように調整します。
    計算されたプランジとアジマスを度に変換 (np.rad2deg) します。
    固有値の順序に基づき、P軸、N軸、T軸のプランジとアジマスをリストに格納して返します。ObsPyの MomentTensor クラスの内部処理や eigh の返り値の規約に依存するため、どの固有ベクトルがP,N,Tに対応するかは注意が必要です。一般的に、固有値が最大、中間、最小のものがそれぞれT軸、B(N)軸、P軸に対応します（圧縮を正とする場合）。
# 返り値: 
    [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth] のリスト（角度は度単位）。
'''

def get_plunge_azimuth(focmec):
    '''
    input: 6 components of  mooment tensor in USE convention (e.g., default GCMT, psmeca format)
    output: A list of [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth] angles in degree
    '''
    m1,m2,m3,m4,m5,m6 = focmec[0],focmec[1],focmec[2],focmec[3],focmec[4],focmec[5]
    mt = obspy.imaging.beachball.MomentTensor(m1,m2,m3,m4,m5,m6, 26)

    (d, v) = np.linalg.eigh(mt.mt)
    pl = np.arcsin(-v[0])
    az = np.arctan2(v[2], -v[1])
    for i in range(0, 3):
        if pl[i] <= 0:
            pl[i] = -pl[i]
            az[i] += np.pi
        if az[i] < 0:
            az[i] += 2 * np.pi
        if az[i] > 2 * np.pi:
            az[i] -= 2 * np.pi

    pl = np.rad2deg(pl)
    az = np.rad2deg(az)

    p_plunge = pl[0]
    p_azimuth = az[0]
    n_plunge = pl[1]
    n_azimuth = az[1]
    t_plunge = pl[2]
    t_azimuth = az[2]

    return [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth]



# --- メイン処理 ---
def main():
    input_filename = "clusteringSnap2c.dat"
    output_filename = "clusteringSnap2_mod.dat"

    # --- データ読み込み ---
    try:
        df = pd.read_csv(
            input_filename,
            delimiter='\t',
            header=None,
            names=INPUT_COLUMNS_30,
            skipinitialspace=True
        )
        print(f"'{input_filename}' から {len(df)} 行のデータを読み込みました。")
    except FileNotFoundError:
        print(f"エラー: 入力ファイル '{input_filename}' が見つかりません。")
        return
    except Exception as e:
        print(f"エラー: '{input_filename}' の読み込み中にエラーが発生しました: {e}")
        return

    # --- selectplane の処理を各行に適用 ---
    print("断層面の選択処理を開始します...")
    
    # 結果を格納するための新しい列を準備 (または既存列を上書きするための準備)
    # DataFrameが大きくなる場合、新しい列を作ってから最後にリネームする方が効率的な場合もある
    
    # apply関数を使って各行に関数を適用
    def apply_select_plane_to_row(row):
        s0, d0, r0 = row['stk0'], row['dip0'], row['rake0']
        s1, d1, r1 = row['stk1'], row['dip1'], row['rake1'] # 指示通り 'stk1' を使用

        # モデルに近い方の断層面パラメータを取得
        selected_stk, selected_dip, selected_rake = selectplane(
            MODEL_STRIKE_DEG, MODEL_DIP_DEG,
            s0, d0, r0,
            s1, d1, r1
        )

        # 選択されたのが元々の (s0,d0,r0) か (s1,d1,r1) かを判断
        # 浮動小数点数の比較なので、ある程度の許容誤差を設けるのが安全だが、
        # selectplane の返り値は入力のどちらかなので、直接比較で問題ないはず
        if (selected_stk == s0 and selected_dip == d0 and selected_rake == r0):
            # 候補0が選択された場合
            new_stk0, new_dip0, new_rake0 = s0, d0, r0
            new_stk1, new_dip1, new_rake1 = s1, d1, r1
        elif (selected_stk == s1 and selected_dip == d1 and selected_rake == r1):
            # 候補1が選択された場合
            new_stk0, new_dip0, new_rake0 = s1, d1, r1
            new_stk1, new_dip1, new_rake1 = s0, d0, r0
        else:
            # selectplaneが入力のどちらでもない値を返した場合 (通常は起こらないはず)
            print(f"警告: selectplaneの返り値が元の候補と一致しませんでした。行: {row.name}")
            # この場合のデフォルト処理 (元の値を維持など) を決める必要がある
            new_stk0, new_dip0, new_rake0 = s0, d0, r0
            new_stk1, new_dip1, new_rake1 = s1, d1, r1
            
        return pd.Series([new_stk0, new_dip0, new_rake0, new_stk1, new_dip1, new_rake1])

    # 新しい列名 (または上書きする列名)
    # 指示: 「モデル面に近い断層面をstk0, dip0, rake0に、そうでない断層面をstk1, dip1, rake1に置き換える」
    # なので、既存の列を上書きする
    result_columns = ['stk0', 'dip0', 'rake0', 'stk1', 'dip1', 'rake1']
    
    df[result_columns] = df.apply(apply_select_plane_to_row, axis=1)
    
    print("断層面の選択処理が完了しました。")

    # --- データ出力 ---
    try:
        df.to_csv(output_filename, sep='\t', index=False, header=True)
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()