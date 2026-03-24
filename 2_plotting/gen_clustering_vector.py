#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# n_vector.dat + クラスタリング結果
# -*- coding: utf-8 -*-

import pandas as pd
import os
import sys

def find_file(filename):
    """
    カレントディレクトリまたは1つ上のディレクトリからファイルを探す
    """
    paths = [
        os.path.join(".", filename),
        os.path.join("..", filename)
    ]
    for path in paths:
        if os.path.exists(path):
            return path
    return None

def main():
    print("--- クラスタリング結果結合プロセスを開始します ---")

    # 1. ファイルの探索
    n_vector_name = "n_vector.dat"
    snap2_name = "clusteringSnap2.dat"

    path_n_vector = find_file(n_vector_name)
    path_snap2 = find_file(snap2_name)

    # 存在チェック
    error_occurred = False
    if not path_n_vector:
        print(f"エラー: {n_vector_name} が見つかりません。")
        error_occurred = True
    if not path_snap2:
        print(f"エラー: {snap2_name} が見つかりません。")
        error_occurred = True

    if error_occurred:
        print("必要なファイルが見つからないため、処理を終了します。")
        sys.exit(1)

    print(f"読み込み中: {path_n_vector}")
    print(f"読み込み中: {path_snap2}")

    try:
        # 2. データの読み込み
        # n_vector.dat は元のコードに合わせ、空白区切り（1つ以上のスペース）で読み込み
        df_vector = pd.read_table(path_n_vector, sep=r'\s+', header=None)
        
        # clusteringSnap2.dat はタブ区切り
        df_snap2 = pd.read_table(path_snap2, sep='\t', header=None)

        # 3. 行数の整合性チェック
        if len(df_vector) != len(df_snap2):
            print(f"警告: 行数が一致しません。")
            print(f"  {n_vector_name}: {len(df_vector)}行")
            print(f"  {snap2_name}: {len(df_snap2)}行")
            print("データの整合性が取れない可能性があるため、処理を中断します。")
            sys.exit(1)

        # 4. clusteringSnap2.dat の末尾2列（clset と label）を抽出
        # 列のインデックスで指定 (最後から2つ目と1つ目)
        results_to_append = df_snap2.iloc[:, -2:]

        # 5. n_vector のデータに結合
        df_combined = pd.concat([df_vector, results_to_append], axis=1)

        # 6. 保存 (clusteringSnap2.dat と同じ設定: タブ区切り、ヘッダーなし、インデックスなし)
        output_filename = "clustering_vector.dat"
        df_combined.to_csv(output_filename, sep='\t', index=False, header=False)

        print(f"--- 成功 ---")
        print(f"結果を {output_filename} に保存しました。")
        print(f"出力列数: {df_combined.shape[1]} 列 (元の13列 + 結果2列)")

        # 7. 変数名の一覧表示
        variable_names = [
            "n", "m", "tw", "t", "x_coord", "y_coord",
            "n1_(1)", "n1_(2)", "n1_(3)",
            "n2_(1)", "n2_(2)", "n2_(3)",
            "sliprate", "clset", "label"
        ]
        
        print("\n--- clustering_vector.dat 変数名一覧 ---")
        for i, name in enumerate(variable_names, 1):
            print(f"列 {i:02d}: {name}")
        print("---------------------------------------")

    except Exception as e:
        print(f"予期せぬエラーが発生しました: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()