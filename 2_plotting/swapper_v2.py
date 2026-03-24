#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def get_env_var_as_float(var_name, script_name="swapper_v2.py"):
    """環境変数を読み込み、浮動小数点数に変換する。失敗時はエラー終了。"""
    val_str = os.environ.get(var_name)
    if val_str is None:
        sys.stderr.write(f"エラー ({script_name}): 環境変数 {var_name} が設定されていません。\n")
        sys.exit(1)
    try:
        return float(val_str)
    except ValueError:
        sys.stderr.write(f"エラー ({script_name}): 環境変数 {var_name} ('{val_str}') を数値に変換できません。\n")
        sys.exit(1)

def process_data_file(input_filename, output_filename, snap_interval, plot_starttime, plot_endtime, tw_map):
    """
    データファイルをフィルタリングし、twを再割り当てして出力する。
    tw_map: フィルタリング後の original_tw を新しい sequential_tw にマッピングする辞書。
            キーが original_tw、値が new_sequential_tw。
    """
    processed_lines = []
    try:
        with open(input_filename, 'r') as f_in:
            for line_num, line_content in enumerate(f_in, 1):
                fields = line_content.strip().split()
                if not fields or len(fields) < 3:
                    # sys.stderr.write(f"警告 ({input_filename}:{line_num}): 列数が3未満です。スキップします。\n")
                    continue
                try:
                    original_tw = int(float(fields[2])) # 3列目 (0-indexed)
                except ValueError:
                    # sys.stderr.write(f"警告 ({input_filename}:{line_num}): 3列目のTW ('{fields[2]}') が数値に変換できません。スキップします。\n")
                    continue

                # この original_tw がフィルタリング対象として tw_map に含まれているか確認
                if original_tw in tw_map:
                    new_tw = tw_map[original_tw]
                    fields[2] = str(new_tw) # 新しいTWで置き換え
                    processed_lines.append(" ".join(fields)) # 元のスペース区切りを維持

    except FileNotFoundError:
        sys.stderr.write(f"情報 ({input_filename}): ファイルが見つかりません。空のファイルとして処理します。\n")
        # 出力ファイルは空になる
    except Exception as e:
        sys.stderr.write(f"エラー ({input_filename}): ファイル処理中にエラー: {e}\n")
        # エラーが発生した場合、出力ファイルは不完全か空になる可能性がある

    try:
        with open(output_filename, 'w') as f_out:
            for line in processed_lines:
                f_out.write(line + "\n")
        print(f"情報: {output_filename} を更新しました。処理後の行数: {len(processed_lines)}")
    except Exception as e:
        sys.stderr.write(f"エラー: {output_filename} の書き込み中にエラー: {e}\n")
        sys.exit(1)


def main():
    print("swapper_v2.py を実行中...")
    snap_interval = get_env_var_as_float('SWAPPER_SNAP_INTERVAL')
    plot_starttime = get_env_var_as_float('SWAPPER_PLOT_STARTTIME')
    plot_endtime = get_env_var_as_float('SWAPPER_PLOT_ENDTIME')

    if snap_interval <= 0:
        sys.stderr.write("エラー: SWAPPER_SNAP_INTERVAL は正である必要があります。\n")
        sys.exit(1)
    if plot_starttime > plot_endtime: # start == end は許容 (データなしになる)
        sys.stderr.write("エラー: SWAPPER_PLOT_STARTTIME は SWAPPER_PLOT_ENDTIME より小さいか等しくする必要があります。\n")
        sys.exit(1)

    original_tws_to_keep = set()
    # まず snap2.dat (より多くの行を含む可能性が高い) をスキャンして対象となる元のTWを特定
    # snap.dat のTWも snap2.dat に含まれるはず (fgensnapの構造上)
    files_to_determine_tws = ["snap2.dat", "snap.dat"] # snap2.dat を優先的にスキャン

    for filename_to_scan in files_to_determine_tws:
        if not os.path.exists(filename_to_scan):
            sys.stderr.write(f"情報: TW特定のためのファイル {filename_to_scan} が見つかりません。スキップします。\n")
            continue
        try:
            with open(filename_to_scan, 'r') as f:
                for line in f:
                    fields = line.strip().split()
                    if not fields or len(fields) < 3: continue
                    try:
                        original_tw = int(float(fields[2]))
                    except ValueError: continue
                    
                    actual_time_start = (original_tw - 1) * snap_interval
                    actual_time_end = original_tw * snap_interval

                    if (actual_time_start < plot_endtime) and (actual_time_end > plot_starttime):
                        original_tws_to_keep.add(original_tw)
        except Exception as e:
            sys.stderr.write(f"エラー: {filename_to_scan} のTW特定スキャン中にエラー: {e}\n")
            # 続行するが、結果が不正確になる可能性

    if not original_tws_to_keep:
        sys.stderr.write("情報: 指定された時間範囲に該当する元の時間窓が見つかりませんでした。\n")
        # 元のファイルを空にする
        for f_to_empty in ["snap.dat", "snap2.dat"]:
            if os.path.exists(f_to_empty):
                open(f_to_empty, 'w').close()
                print(f"情報: {f_to_empty} は空として出力されました（対象データなし）。")
        sys.exit(0)

    sorted_original_tws = sorted(list(original_tws_to_keep))
    # original_tw を新しい連続したtw (1-indexed) にマッピング
    original_tw_to_new_sequential_tw_map = {
        orig_tw: new_seq_tw + 1 for new_seq_tw, orig_tw in enumerate(sorted_original_tws)
    }
    print(f"情報: フィルタリング対象の元のTW: {sorted_original_tws}")
    print(f"情報: TWマッピング (元TW -> 新TW): {original_tw_to_new_sequential_tw_map}")


    process_data_file("snap.dat", "snap.dat", snap_interval, plot_starttime, plot_endtime, original_tw_to_new_sequential_tw_map)
    process_data_file("snap2.dat", "snap2.dat", snap_interval, plot_starttime, plot_endtime, original_tw_to_new_sequential_tw_map)

    print("swapper_v2.py: 処理完了。")

if __name__ == '__main__':
    main()
