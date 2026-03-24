#!/usr/bin/env python3
'''

'''


import sys
import os
import argparse
import pandas as pd
import numpy as np
# math モジュールは ar2plp.py 内でインポート・使用されるため、
# calc_fault.py で直接使用していなければ、ここでのインポートは必須ではありません。

# --- calc_fault.py があるディレクトリをPythonのモジュール検索パスに追加 ---
# これにより、calc_fault.py と同じディレクトリにある ar2plp.py を
# カレントディレクトリに関わらずインポートできるようになる。
try:
    # スクリプトとして実行された場合、__file__ はそのスクリプトの絶対パスを指す
    # (シンボリックリンク経由で実行された場合でも、リンク先の実際のパス)
    script_path = os.path.abspath(__file__)
    script_directory = os.path.dirname(script_path)
    if script_directory not in sys.path:
        sys.path.insert(0, script_directory)
except NameError:
    # 対話的インタプリタなど、__file__ が定義されていない環境で実行された場合
    # (このスクリプトは通常コマンドラインから実行されるので、あまり考慮不要かもしれないが念のため)
    print("警告: スクリプトのディレクトリを特定できませんでした (__file__ が未定義)。")
    print("      ar2plp.py がPythonの検索パス上にあることを確認してください。")
    script_directory = os.getcwd() # フォールバックとしてカレントディレクトリを仮定

try:
    import ar2plp as ps # ar2plp.py を ps という名前でインポート
except ImportError:
    print(f"エラー: ar2plp.py モジュールが見つかりません。")
    print(f"calc_fault.py と同じディレクトリ ({script_directory}) に ar2plp.py を配置してください。")
    print(f"現在のPython検索パスの先頭: {sys.path[0] if sys.path else '空'}")
    exit(1)
except Exception as e:
    print(f"ar2plp.py のインポート中に予期せぬエラーが発生しました: {e}")
    exit(1)


# 入力データの期待される列名 (ユーザー提供のもの)
INPUT_COLUMNS_30 = [
    'n', 'm', 'tw', 'dx', 'dy', 'sliprate', 'Mrr', 'Mss',
    'Mee', 'Mrs', 'Mre', 'Mse', 'str1', 'str2', 'dip1',
    'dip2', 'rake1', 'rake2', 'lat', 'lon', 'depth',
    'trendp', 'trendt', 'trendb', 'plungp',
    'plungt', 'plungb', 'NDC', 'clla', 'no'
]

def convert_harvard_to_aki_richards_tensor(mrr, mss, mee, mrs, mre, mse):
    """
    ハーバードCMT規約 (Up, South, East) のモーメントテンソル成分を
    Aki & Richards規約 (North, East, Down) の3x3テンソルに変換します。
    """
    am = np.zeros((3, 3), dtype=np.float64)
    am[0, 0] = mss  # Mxx (North-North) = M_South-South
    am[1, 1] = mee  # Myy (East-East)   = M_East-East
    am[2, 2] = mrr  # Mzz (Down-Down)   = M_Up-Up
    am[0, 1] = -mse # Mxy (North-East)  = -M_South-East
    am[1, 0] = am[0, 1]
    am[0, 2] = -mrs # Mxz (North-Down)  = -M_Up-South
    am[2, 0] = am[0, 2]
    am[1, 2] = mre  # Myz (East-Down)   = M_Up-East
    am[2, 1] = am[1, 2]
    return am

def calculate_parameters_from_mt_components(mrr, mss, mee, mrs, mre, mse, ref_strike_deg, ref_dip_deg):
    """
    モーメントテンソル成分から断層パラメータ一式を計算するヘルパー関数。
    ar2plp.py (psエイリアス) 内のPython関数を使用。
    """
    am_tensor = convert_harvard_to_aki_richards_tensor(mrr, mss, mee, mrs, mre, mse)
    
    try:
        outputs_ar2plp = ps.ar2plp_py(am_tensor)
    except Exception as e:
        # エラー発生時にどのテンソルで問題が起きたか分かるように情報を追加
        # print(f"デバッグ情報: am_tensor = \n{am_tensor}")
        raise RuntimeError(f"ps.ar2plp_py の呼び出し中にエラー: {e}")

    if len(outputs_ar2plp) != 20:
        raise ValueError(f"ar2plp_pyが予期しない出力数 ({len(outputs_ar2plp)}) を返しました。期待値20。")
    
    am0,am1,e_iso,am0b,phia,deltaa,alama,slipa,phib,deltab,alamb,slipb, \
    trendp,plungp,trendt,plungt,trendb,plungb,eta_clvd,ierr_fps = outputs_ar2plp

    try:
        nodal_p = ps.selection_nodal_py(
            phia, deltaa, phib, deltab,
            ref_strike_deg, ref_dip_deg 
            # ps.RAD_CONST は引数として渡さない (selection_nodal_py内部で参照すると仮定)
        )
        s1, d1, r1, s2, d2, r2 = ps.swap_fault_py(
            phia, deltaa, alama,
            phib, deltab, alamb,
            nodal_p
        )
    except Exception as e:
        raise RuntimeError(f"節面選択/交換処理中にエラー (Python版): {e}")

    return {
        'str1_calc': s1, 'dip1_calc': d1, 'rake1_calc': r1,
        'str2_calc': s2, 'dip2_calc': d2, 'rake2_calc': r2,
        'trendp_calc': trendp, 'plungp_calc': plungp,
        'trendt_calc': trendt, 'plungt_calc': plungt,
        'trendb_calc': trendb, 'plungb_calc': plungb,
        'fps_M0': am0, 
        'fps_M1': am1,
        'fps_e_iso': e_iso,
        'fps_M0b': am0b, 
        'fps_slipa': slipa, 
        'fps_slipb': slipb,
        'fps_eta_clvd': eta_clvd, 
        'fps_ierr': ierr_fps
    }

def process_fault_data_row(row_series, ref_strike_deg, ref_dip_deg):
    """
    DataFrameの1行を処理して、新しい断層パラメータを計算し、元の行にマージします。
    """
    try:
        mrr = float(row_series['Mrr'])
        mss = float(row_series['Mss'])
        mee = float(row_series['Mee'])
        mrs = float(row_series['Mrs'])
        mre = float(row_series['Mre'])
        mse = float(row_series['Mse'])
    except KeyError as e:
        raise ValueError(f"必要なモーメントテンソル列が見つかりません: {e}. 利用可能な列: {list(row_series.index)}")
    except ValueError as e:
        raise ValueError(f"モーメントテンソル値 '{row_series.get('Mrr', 'N/A')}', etc. の数値変換エラー: {e}")

    calculated_params = calculate_parameters_from_mt_components(
        mrr, mss, mee, mrs, mre, mse, ref_strike_deg, ref_dip_deg
    )
    
    updated_row = row_series.copy()
    for key, value in calculated_params.items():
        updated_row[key] = value
        
    return updated_row

def process_specific_clustering_files(ref_strike_deg, ref_dip_deg):
    """
    clusteringSnap2c.dat と clmeca_fault.dat に関連する追加処理を実行します。
    """
    input_filename = "clusteringSnap2c.dat"
    output_r_filename = "clusteringSnap2r.dat"
    output_clmeca_filename = "clmeca_fault.dat"
    
    print(f"\n### 追加処理開始: {input_filename} ###")

    try:
        # ファイルが存在するかまず確認 (よりPythonicな方法)
        if not os.path.exists(input_filename):
            raise FileNotFoundError(f"{input_filename} が見つかりません。")
            
        df_c = pd.read_csv(input_filename, delimiter='\t', skipinitialspace=True, names=INPUT_COLUMNS_30, header=None)
        print(f"{input_filename} から {len(df_c)} 行のデータを読み込みました (機能1)。")
    except FileNotFoundError:
        print(f"エラー: {input_filename} が見つかりません。追加処理をスキップします。")
        return
    except Exception as e:
        print(f"エラー: {input_filename} の読み込み中にエラー (機能1): {e}")
        return

    required_mt_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
    actual_columns_c = [str(col).strip() for col in df_c.columns] 
    df_c.columns = actual_columns_c
    missing_cols_c = [col for col in required_mt_cols if col not in df_c.columns]
    if missing_cols_c:
        print(f"エラー: {input_filename} に必要なモーメントテンソル列がありません: {missing_cols_c}")
        print(f"利用可能な列: {list(df_c.columns)}")
        return

    processed_rows_for_r = []
    for index, row in df_c.iterrows():
        print(f"機能1処理中: 行 {index + 1}/{len(df_c)}", end="\r")
        try:
            params = calculate_parameters_from_mt_components(
                float(row['Mrr']), float(row['Mss']), float(row['Mee']),
                float(row['Mrs']), float(row['Mre']), float(row['Mse']),
                ref_strike_deg, ref_dip_deg
            )
            updated_row = row.copy()
            updated_row['str1'] = params['str1_calc']
            updated_row['dip1'] = params['dip1_calc']
            updated_row['rake1'] = params['rake1_calc']
            updated_row['str2'] = params['str2_calc']
            updated_row['dip2'] = params['dip2_calc']
            updated_row['rake2'] = params['rake2_calc']
            processed_rows_for_r.append(updated_row)
        except Exception as e:
            print(f"\nエラー: 機能1 行 {index + 1} 処理中: {e}")
            error_row = row.copy()
            error_row['processing_error_func1'] = str(e)
            processed_rows_for_r.append(error_row)
    print("\n機能1の全行処理完了。")

    if processed_rows_for_r:
        df_r_output = pd.DataFrame(processed_rows_for_r)
        try:
            df_r_output.to_csv(output_r_filename, sep='\t', index=False, header=True)
            print(f"機能1結果が {output_r_filename} に保存されました。")
        except Exception as e:
            print(f"エラー: {output_r_filename} の書き込み中にエラー: {e}")
    
    # 'no', 'sliprate', 'clla' 列の存在確認と型変換 (機能2用)
    required_cols_func2 = ['no', 'sliprate', 'clla']
    missing_cols_func2 = [col for col in required_cols_func2 if col not in df_c.columns]
    if missing_cols_func2:
        print(f"エラー: {input_filename} に機能2で必要な列がありません: {missing_cols_func2}。機能2をスキップします。")
        return
        
    try:
        df_c['no'] = pd.to_numeric(df_c['no'], errors='coerce').fillna(-1).astype(int)
        df_c['sliprate'] = pd.to_numeric(df_c['sliprate'], errors='coerce')
        # cllaは文字列のまま扱うか、必要に応じて型変換
        df_c['clla'] = df_c['clla'].astype(str) 
    except Exception as e:
        print(f"エラー: 'no', 'sliprate', 'clla' 列の型変換中にエラー: {e}")
        return

    print(f"\n機能2処理開始: グループ化と平均モーメントテンソルからの計算...")
    grouped_clmeca_data = []
    
    # 'no' 列でグループ化する前に、NaNや変換失敗で -1 になったものを除外するかどうか検討
    # ここでは -1 グループも処理されるが、通常は有効な 'no' 値のみが対象となるはず
    for no_val, group in df_c.groupby('no'):
        if no_val == -1 and group['no'].iloc[0] == -1 : # 'no'が無効だった行のグループ
            print(f"機能2処理中: 'no'列が無効だったデータのグループ (no={no_val}) をスキップ")
            continue
        print(f"機能2処理中: グループ no={no_val}")
        try:
            # モーメントテンソル成分を数値に変換してから平均を計算
            mt_group_numeric = group[required_mt_cols].apply(pd.to_numeric, errors='coerce')
            avg_mt_series = mt_group_numeric.mean()
            
            if avg_mt_series.isnull().any():
                print(f"警告: グループ no={no_val} で平均モーメントテンソル成分にNaNが含まれます。スキップします。")
                print(f"  NaNを含む平均MT成分: {avg_mt_series[avg_mt_series.isnull()]}")
                grouped_clmeca_data.append({'no': no_val, 'processing_error_func2': 'NaN in avg MT components'})
                continue

            representative_params = calculate_parameters_from_mt_components(
                avg_mt_series['Mrr'], avg_mt_series['Mss'], avg_mt_series['Mee'],
                avg_mt_series['Mrs'], avg_mt_series['Mre'], avg_mt_series['Mse'],
                ref_strike_deg, ref_dip_deg
            )
            
            clla_representative = np.nan # デフォルト値
            valid_sliprates = group['sliprate'].dropna() # NaNを除外
            if not valid_sliprates.empty:
                # idxmax() はNaNを除外した上で最大値のインデックスを返す
                clla_representative = group.loc[valid_sliprates.idxmax()]['clla']
            
            record = {'no': no_val, 'clla_representative': clla_representative}
            record.update(representative_params)
            for mt_col in required_mt_cols:
                record[f'avg_{mt_col}'] = avg_mt_series[mt_col]
            
            grouped_clmeca_data.append(record)
        except Exception as e:
            print(f"\nエラー: 機能2 グループ no={no_val} 処理中: {e}")
            # エラー発生時のデバッグ情報
            # print(f"  エラー発生時のグループデータ (先頭5行):\n{group.head()}")
            # print(f"  エラー発生時の平均MTシリーズ:\n{avg_mt_series if 'avg_mt_series' in locals() else '未計算'}")
            grouped_clmeca_data.append({'no': no_val, 'processing_error_func2': str(e)})
            
    print("機能2のグループ処理完了。")

    if grouped_clmeca_data:
        df_clmeca_output = pd.DataFrame(grouped_clmeca_data)
        # 出力列の順序を定義 (存在しない列は無視される)
        desired_cols_order = ['no', 'clla_representative'] + \
                             [f'avg_{col}' for col in required_mt_cols] + \
                             ['str1_calc', 'dip1_calc', 'rake1_calc', 
                              'str2_calc', 'dip2_calc', 'rake2_calc',
                              'trendp_calc', 'plungp_calc', 'trendt_calc', 'plungt_calc', 'trendb_calc', 'plungb_calc',
                              'fps_M0', 'fps_M1', 'fps_e_iso', 'fps_M0b', 
                              'fps_slipa', 'fps_slipb', 'fps_eta_clvd', 'fps_ierr',
                              'processing_error_func2'] # エラー列も考慮
        
        # 実際に存在する列のみで順序を再構成
        final_cols_order = [col for col in desired_cols_order if col in df_clmeca_output.columns]
        # desired_cols_order にないがDataFrameには存在する列を末尾に追加
        for col in df_clmeca_output.columns:
            if col not in final_cols_order:
                final_cols_order.append(col)

        df_clmeca_output = df_clmeca_output[final_cols_order]

        try:
            df_clmeca_output.to_csv(output_clmeca_filename, sep='\t', index=False, header=True)
            print(f"機能2結果が {output_clmeca_filename} に保存されました。")
        except Exception as e:
            print(f"エラー: {output_clmeca_filename} の書き込み中にエラー: {e}")

    print("### 追加処理終了 ###")

def main():
    parser = argparse.ArgumentParser(
        description="指定されたデータファイルから断層パラメータを再計算します。また、特定のクラスタリングデータファイルに対する追加処理も実行します。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", nargs='?', default=None, 
                        help="[汎用処理用] 入力データファイルへのパス (CSV形式を想定)。省略すると汎用処理はスキップ。")
    parser.add_argument("output_file", nargs='?', default=None,
                        help="[汎用処理用] 処理結果を保存する出力ファイルへのパス (CSV形式)。input_file指定時必須。")
    parser.add_argument(
        "--ref_strike", 
        type=float, 
        default=0.0, 
        help="節面選択のための参照走向 (度)。全処理で共通。"
    )
    parser.add_argument(
        "--ref_dip", 
        type=float, 
        default=90.0, 
        help="節面選択のための参照傾斜 (度)。全処理で共通。"
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default=",",
        help="[汎用処理用] 入力CSVファイルの区切り文字。"
    )
    
    args = parser.parse_args()

    try:
        ps.fpsset_py() 
        print("Python版の定数初期化関数 (ps.fpsset_py) が呼び出されました。")
    except AttributeError: # ar2plp.py に fpsset_py がない場合など
        print("エラー: ar2plp.py に fpsset_py 関数が見つかりません。ar2plp.pyの実装を確認してください。")
        return
    except Exception as e:
        print(f"ps.fpsset_py の呼び出し中にエラーが発生しました: {e}")
        return

    # --- 汎用ファイル処理 ---
    if args.input_file and args.output_file:
        print(f"\n### 汎用ファイル処理開始: {args.input_file} -> {args.output_file} ###")
        try:
            if not os.path.exists(args.input_file):
                raise FileNotFoundError(f"{args.input_file} が見つかりません。")
            df = pd.read_csv(args.input_file, delimiter=args.delimiter, skipinitialspace=True)
            print(f"{args.input_file} から {len(df)} 行のデータを読み込みました。")
        except FileNotFoundError:
            print(f"エラー: 入力ファイルが見つかりません: {args.input_file}")
            # 汎用処理はスキップし、追加処理に進む
        except Exception as e:
            print(f"入力ファイルの読み込み中にエラーが発生しました: {e}")
            return # 致命的なエラーの場合は全体を終了
        else: # ファイル読み込み成功時のみ実行
            required_mt_cols = ['Mrr', 'Mss', 'Mee', 'Mrs', 'Mre', 'Mse']
            actual_columns = [str(col).strip() for col in df.columns]
            df.columns = actual_columns
            missing_cols = [col for col in required_mt_cols if col not in df.columns]
            if missing_cols:
                print(f"エラー: 入力ファイルに必要なモーメントテンソル列がありません: {missing_cols}")
            else:
                processed_rows_list = []
                for index, row in df.iterrows():
                    print(f"汎用処理中: 行 {index + 1}/{len(df)}", end="\r")
                    try:
                        updated_row_series = process_fault_data_row(row, args.ref_strike, args.ref_dip)
                        processed_rows_list.append(updated_row_series)
                    except Exception as e:
                        print(f"\nエラー: 汎用処理 行 {index + 1} 処理中: {e}")
                        error_info_row = row.copy()
                        error_info_row['processing_error_main'] = str(e)
                        processed_rows_list.append(error_info_row)
                print("\n汎用ファイルの全行処理が完了しました。")

                if processed_rows_list:
                    output_df = pd.DataFrame(processed_rows_list)
                    try:
                        output_df.to_csv(args.output_file, index=False) 
                        print(f"汎用処理結果が {args.output_file} に保存されました。")
                    except Exception as e:
                        print(f"出力ファイルの書き込み中にエラーが発生しました: {e}")
        print("### 汎用ファイル処理終了 ###")
    elif args.input_file or args.output_file: # どちらか一方だけ指定された場合
        print("警告: 汎用ファイル処理のためには input_file と output_file の両方を指定する必要があります。汎用処理はスキップされます。")

    # --- clusteringSnap2c.dat に対する追加処理 ---
    process_specific_clustering_files(args.ref_strike, args.ref_dip)

if __name__ == "__main__":
    main()
