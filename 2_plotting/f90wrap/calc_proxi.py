#!/usr/bin/env python3
# calc_proxi.py
'''
RMSD (Root Mean Square Distance) 指標の解釈
意味: RMSDは、グループ内の各データ点から、そのグループの時空間的な重心（平均位置）までの距離の二乗平均平方根です。
解釈:
値が小さい: グループ内のデータ点が重心の近くに密集していることを意味します。つまり、そのグループのイベントは時空間的に非常に**近接している（まとまっている）**と言えます。
値が大きい: グループ内のデータ点が重心から広範囲に散らばっていることを意味します。つまり、そのグループのイベントは時空間的に**あまり近接していない（散らばっている）**と言えます。
特徴:
直感的に「平均的なばらつきの半径」のようなイメージで捉えることができます。
計算が比較的シンプルです。
データの散らばりの異方性（特定の方向に伸びているなど）は直接的には表現しません。あくまで重心からの平均的な距離です。
'''
'''
PCA (Principal Component Analysis) 指標の解釈
calc_proxi.py で計算しているPCA指標は、具体的には「標準化された座標における主成分の固有値の幾何平均」です。

PCAの基本: 主成分分析は、データのばらつきが最も大きい方向（第1主成分）、次に大きい方向（第2主成分）…というように、データの構造を捉える手法です。各主成分の「広がり具合」は固有値で表されます。
固有値の意味:
3次元のデータ（標準化された 'n', 'm', 'tw'）に対してPCAを行うと、3つの固有値 (λ1,λ2,λ3) が得られます。
これらの固有値は、それぞれ対応する主成分方向へのデータの分散の大きさを表します。つまり、データがその方向にどれだけ広がっているかを示します。
幾何平均の意味:
固有値の幾何平均 (λ1⋅λ2⋅λ3)1/3
  は、データの散らばりが作る3次元の楕円体の「代表的な半径」あるいは「体積的なスケール」のようなものと解釈できます。
解釈:
値が小さい: 3つの主成分方向への広がりがいずれも小さい、あるいは少なくともいくつかの方向への広がりが非常に小さいことを意味します。これは、データが時空間的に小さな「ボリューム」内に密集していることを示唆し、近接性が高いと言えます。
値が大きい: 少なくとも1つ以上の主成分方向へデータが大きく広がっていることを意味します。これは、データが時空間的に大きな「ボリューム」を占めていることを示唆し、**近接性が低い（散らばっている）**と言えます。
特徴:
データの散らばりの異方性（特定の方向に細長く伸びている、あるいはある平面上で広がっているなど）を考慮した指標となります。
RMSDが同じような値を示すグループでも、PCA指標（や個々の固有値の分布）を見ることで、散らばりの形状の違いを区別できる可能性があります。
'''
'''
2つの指標の比較と利用
RMSD は、グループの「平均的なコンパクトさ」を測るのに適しています。
PCA指標 は、グループの「時空間的な占有ボリュームのスケール」や「散らばりの異方性を含めた総合的な広がり」を評価するのに適しています。
'''

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# --- 使用するデータ直接書き込み ---
# これらの値は、'n', 'm', 'tw' の各単位が実際の物理量としてどれだけのスケールを持つかを示します。
# 必要に応じてこれらの値を変更してください。
KM_PER_GRID_N_UNIT = 5.0  # 'n'の1単位あたりの距離 (km)
KM_PER_GRID_M_UNIT = 5.0  # 'm'の1単位あたりの距離 (km)
SEC_PER_TW_UNIT = 0.5     # 'tw'の1単位あたりの時間 (秒)

# --- 入力データの列名定義 (clusteringSnap2_mod.dat のヘッダーを想定) ---
# 実際のファイルヘッダーと一致しているか確認してください。
# もしファイルにヘッダーがない場合は、names=... で指定する必要があります。

# --- 指標計算用の関数定義 ---

def calculate_rmsd(group_df_std, columns=['n_std', 'm_std', 'tw_std']):
    """
    グループ内の標準化されたデータ点群の重心からのRMS距離を計算します。
    group_df_std: 標準化された座標を持つグループのDataFrame
    columns: RMSD計算に使用する列名のリスト
    """
    if group_df_std.empty or group_df_std[columns].isnull().all().all():
        return np.nan
    
    # NaNを含む行を削除してから計算 (またはNaNを0で埋めるなどの処理も検討可)
    valid_data = group_df_std[columns].dropna()
    if valid_data.empty:
        return np.nan

    centroid = valid_data.mean(axis=0).values # 各列の平均 (重心)
    
    # 各点から重心までの距離の二乗
    squared_distances = np.sum((valid_data.values - centroid)**2, axis=1)
    rmsd = np.sqrt(np.mean(squared_distances))
    return rmsd

def calculate_pca_metric(group_df_std, columns=['n_std', 'm_std', 'tw_std']):
    """
    グループ内の標準化されたデータ点群に対してPCAを行い、
    固有値の幾何平均を散らばりの指標として計算します。
    group_df_std: 標準化された座標を持つグループのDataFrame
    columns: PCAに使用する列名のリスト
    """
    if group_df_std.empty or group_df_std[columns].isnull().all().all():
        return np.nan

    valid_data = group_df_std[columns].dropna()
    if len(valid_data) < len(columns): # PCAを適用するにはサンプル数が特徴量数以上必要
        # print(f"警告: グループ内の有効なデータ点数 ({len(valid_data)}) が特徴量数 ({len(columns)}) より少ないためPCA指標はNaNになります。")
        return np.nan

    try:
        pca = PCA(n_components=len(columns))
        pca.fit(valid_data.values)
        eigenvalues = pca.explained_variance_
        
        # 固有値が非常に小さい（またはゼロ）場合の処理
        # 幾何平均は正の数に対して定義されるため、非負の値に調整
        eigenvalues_processed = np.maximum(eigenvalues, 1e-12) # 非常に小さい正の値に置き換え
        
        # 幾何平均: (λ1 * λ2 * λ3)^(1/k)
        pca_metric = np.exp(np.mean(np.log(eigenvalues_processed)))
        
    except Exception as e:
        # print(f"警告: PCA計算中にエラーが発生しました ({e})。PCA指標はNaNになります。")
        return np.nan
        
    return pca_metric

# --- メイン処理 ---
def main():
    input_filename = "clusteringSnap2_mod.dat"
    output_filename = "cl_proxi.dat"
    plot_filename = "proximity_metrics_plot.png"

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

    # --- 必要な列の確認と型変換 ---
    required_coord_cols = ['n', 'm', 'tw', 'no']
    for col in required_coord_cols:
        if col not in df.columns:
            print(f"エラー: 入力ファイルに必要な列 '{col}' が見つかりません。")
            return
        try:
            df[col] = pd.to_numeric(df[col], errors='raise') # 数値でない場合はエラー
        except ValueError as e:
            print(f"エラー: 列 '{col}' の数値への変換に失敗しました: {e}")
            print(f"       '{col}' 列に数値以外のデータが含まれている可能性があります。確認してください。")
            return
            
    df.dropna(subset=required_coord_cols, inplace=True) # 必須列にNaNがある行は除外
    df['no'] = df['no'].astype(int)

    # --- 各グリッドと時間窓のズレを物理単位に変換 ---
    df['n_phys'] = df['n'] * KM_PER_GRID_N_UNIT
    df['m_phys'] = df['m'] * KM_PER_GRID_M_UNIT
    df['tw_phys'] = df['tw'] * SEC_PER_TW_UNIT
    
    phys_cols_to_std = ['n_phys', 'm_phys', 'tw_phys']

    # --- データ全体の標準化 ---
    print("物理単位に変換されたデータの標準化を開始します...")
    scaler = StandardScaler()
    try:
        df_std_values = scaler.fit_transform(df[phys_cols_to_std])
    except ValueError as e:
        print(f"エラー: データの標準化中にエラーが発生しました。データが少なすぎるか、分散がゼロの可能性があります。: {e}")
        # データが1点しかない場合や、全ての値が同じ場合に発生しうる
        if len(df) <=1 :
            print("データが1行以下のため、標準化と指標計算はスキップします。")
            return
        # 3列とも同じ値の場合など
        is_constant = True
        for col in phys_cols_to_std:
            if df[col].nunique() > 1:
                is_constant = False
                break
        if is_constant:
            print("n_phys, m_phys, tw_phys の全ての値が定数のため、標準化できません。指標計算はスキップします。")
            return
        raise e # それ以外のValueErrorは再スロー

    df_std = pd.DataFrame(df_std_values, columns=['n_std', 'm_std', 'tw_std'], index=df.index)
    df = pd.concat([df, df_std], axis=1)
    print("データの標準化が完了しました。")

    # --- noの値毎にデータのグループ化と指標計算 ---
    print("グループごとの指標計算を開始します...")
    
    results = []
    # 'no' 列でグループ化する前に、標準化された列にNaNがないか確認
    # (StandardScalerは通常NaNを生成しないが、入力にNaNが多かった場合の影響を考慮)
    df.dropna(subset=['n_std', 'm_std', 'tw_std'], inplace=True)

    if df.empty:
        print("有効なデータがないため、グループ化と指標計算はスキップします。")
    else:
        for group_no, group_data in df.groupby('no'):
            print(f"  グループ no={group_no} の指標を計算中...")
            
            rmsd_val = calculate_rmsd(group_data, columns=['n_std', 'm_std', 'tw_std'])
            pca_val = calculate_pca_metric(group_data, columns=['n_std', 'm_std', 'tw_std'])
            
            results.append({'no': group_no, 'rmsd': rmsd_val, 'pca': pca_val})
    
    print("グループごとの指標計算が完了しました。")

    if not results:
        print("計算結果がありません。")
        return

    df_results = pd.DataFrame(results)

    # --- グラフ化 ---
    if not df_results.empty:
        print(f"指標のグラフをプロットして '{plot_filename}' に保存します...")
        try:
            fig, ax1 = plt.subplots(figsize=(10, 6))

            color1 = 'tab:red'
            ax1.set_xlabel('no (グループ番号)')
            ax1.set_ylabel('RMSD 指標 (標準化後)', color=color1)
            ax1.plot(df_results['no'], df_results['rmsd'], color=color1, marker='o', linestyle='-', label='RMSD')
            ax1.tick_params(axis='y', labelcolor=color1)
            # NaN値をプロットでどう扱うか（現状は線が途切れる）

            ax2 = ax1.twinx()  # 共通のx軸を持つ2番目のy軸
            color2 = 'tab:blue'
            ax2.set_ylabel('PCA 指標 (標準化後固有値の幾何平均)', color=color2)
            ax2.plot(df_results['no'], df_results['pca'], color=color2, marker='s', linestyle='--', label='PCA')
            ax2.tick_params(axis='y', labelcolor=color2)

            fig.tight_layout() # ラベルが重ならないように調整
            plt.title('各グループの時空間的近接性指標')
            # 凡例を一つにまとめる (少し工夫が必要)
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2, loc='best')
            
            plt.savefig(plot_filename)
            print(f"グラフが '{plot_filename}' に保存されました。")
            # plt.show() # 対話的に表示する場合はコメントを外す
            plt.close(fig) # メモリ解放
        except Exception as e:
            print(f"エラー: グラフのプロット中にエラーが発生しました: {e}")
    else:
        print("プロットするデータがありません。")


    # --- データ出力 ---
    # 指定された列順序
    output_columns_order = ['no', 'rmsd', 'pca']
    
    try:
        df_results[output_columns_order].to_csv(output_filename, sep='\t', index=False, header=True, float_format='%.6g')
        print(f"処理結果が '{output_filename}' に保存されました。")
    except Exception as e:
        print(f"エラー: '{output_filename}' の書き込み中にエラーが発生しました: {e}")

if __name__ == '__main__':
    main()
