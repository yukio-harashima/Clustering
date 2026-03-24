#!/usr/bin/env python3
# calc_fault_v3_count.py
# 修正点: 
# 1. faultnormalvecの戻り値をFortranコードに合わせて [North, East, Down] の順に変更
# 2. 断層面の入れ替えが発生した行数をカウントして表示する機能を追加
# clusteringSnap2c.dat(dataeditor.py)を読み込み設定した節面に近い面に数値を入れ替える

import numpy as np
import pandas as pd
import math

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


# select a preferred fault plane
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    """
    モデル断層面に「より近い」断層面を選択して返す。
    内部で faultnormalvec を使用して法線ベクトルを計算し、内積で比較する。
    """
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    
    # 内積計算 (成分順序が変わっても、全てのベクトルで統一されていれば結果は同じ)
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


# define a fault-normal vector
def faultnormalvec(stk, dip):
    """
    断層面の法線ベクトルを計算する。
    
    修正: Fortranコード (get_n_vector) に合わせ、戻り値の順序を
          [North, East, Down] に統一しました。
          
    stk: 走向 (度)
    dip: 傾斜 (度)
    """
    # 北成分 (North) = -sin(dip) * sin(stk)
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    
    # 東成分 (East) = sin(dip) * cos(stk)
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    
    # 下成分 (Down) = -cos(dip)
    nd = -np.cos(np.deg2rad(dip))
    
    # 以前のバージョン: return np.array([ne, nn, nd]) # [East, North, Down]
    # 今回の修正: Fortranに合わせて [North, East, Down] の順で返す
    return np.array([nn, ne, nd])


def get_plunge_azimuth(focmec):
    '''
    input: 6 components of  mooment tensor in USE convention (e.g., default GCMT, psmeca format)
    output: A list of [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth] angles in degree
    
    Note: この関数を使用するには 'import obspy' が必要です。現在はmainから呼ばれていません。
    '''
    # import obspy # 必要であればここでインポート、または冒頭に追加
    m1,m2,m3,m4,m5,m6 = focmec[0],focmec[1],focmec[2],focmec[3],focmec[4],focmec[5]
    # 実行には obspy が必要ですが、インストールされていない環境を考慮し、
    # 呼び出されない限りエラーにならないようそのままにしてあります。
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
    
    # 入れ替えカウント用変数
    swap_count = 0

    def apply_select_plane_to_row(row):
        nonlocal swap_count # 外側の変数を更新するために必要
        
        s0, d0, r0 = row['stk0'], row['dip0'], row['rake0']
        s1, d1, r1 = row['stk1'], row['dip1'], row['rake1']

        # モデルに近い方の断層面パラメータを取得
        selected_stk, selected_dip, selected_rake = selectplane(
            MODEL_STRIKE_DEG, MODEL_DIP_DEG,
            s0, d0, r0,
            s1, d1, r1
        )

        # 選択されたのが元々の (s0,d0,r0) か (s1,d1,r1) かを判断
        if (selected_stk == s0 and selected_dip == d0 and selected_rake == r0):
            # 候補0が選択された場合 (入れ替えなし)
            new_stk0, new_dip0, new_rake0 = s0, d0, r0
            new_stk1, new_dip1, new_rake1 = s1, d1, r1
        elif (selected_stk == s1 and selected_dip == d1 and selected_rake == r1):
            # 候補1が選択された場合 (入れ替えあり)
            swap_count += 1
            new_stk0, new_dip0, new_rake0 = s1, d1, r1
            new_stk1, new_dip1, new_rake1 = s0, d0, r0
        else:
            # 予期しないケース (とりあえず入れ替えなしとして扱う)
            print(f"警告: selectplaneの返り値が元の候補と一致しませんでした。行: {row.name}")
            new_stk0, new_dip0, new_rake0 = s0, d0, r0
            new_stk1, new_dip1, new_rake1 = s1, d1, r1
            
        return pd.Series([new_stk0, new_dip0, new_rake0, new_stk1, new_dip1, new_rake1])

    # 結果列の定義
    result_columns = ['stk0', 'dip0', 'rake0', 'stk1', 'dip1', 'rake1']
    
    # 適用
    df[result_columns] = df.apply(apply_select_plane_to_row, axis=1)
    
    print("断層面の選択処理が完了しました。")
    print(f"入れ替えられた行数: {swap_count} / {len(df)}")

    # --- データ出力 ---
    try:
        df.to_csv(output_filename, sep='\t', index=False, header=True)
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()