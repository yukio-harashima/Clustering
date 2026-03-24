#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# vscode上のインタラクティブウィンドウで各種データを確認する用
import os
import shutil

def copy_single_file(source_dir, dest_dir, filename):
    """
    一つのファイルをsource_dirからdest_dirへコピーします。

    Args:
        source_dir (str): コピー元のディレクトリパス
        dest_dir (str): コピー先のディレクトリパス
        filename (str): コピーするファイル名
    """
    # ファイルパスを組み立てる
    source_path = os.path.join(source_dir, filename)
    dest_path = os.path.join(dest_dir, filename)

    print(f"処理中のファイル: {filename}")

    # コピー元ファイルが存在するか確認
    if not os.path.exists(source_path):
        print(f"  -> スキップ: ファイルが見つかりません: {source_path}\n")
        return

    # ファイルをコピー
    try:
        # shutil.copy2 はタイムスタンプなどのメタデータも一緒にコピーします
        shutil.copy2(source_path, dest_path)
        print(f"  -> 成功: コピーが完了しました。\n")
    except Exception as e:
        print(f"  -> エラー: ファイルのコピーに失敗しました。\n     理由: {e}\n")

# --- ここから設定 ---

# コピー元のディレクトリ
SOURCE_DIRECTORY = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20250730-100041/filt-arbP75_distCUSTw10-10-10_nClust8_20250930-072648"

# コピー先のディレクトリ
DESTINATION_DIRECTORY = "/Users/harashima-yukio/Desktop/1993_Hokkaido/results_20250730-100041/filt-arbP75_distCUSTw10-10-10_nClust8_20250930-072648/relabeled_nClust6"

# コピーしたいファイル名のリスト（ここにファイル名を追加・削除してください）
FILES_TO_COPY = [
    "snap_yr.dat",
    "snap2_y.dat",
    "fort.40",
    "faultline.dat",
    "mrf.dat",
    "snap2.dat",
    "rigid_amp.info" # このファイルはスキップされます
]

# --- 設定ここまで ---

# --- メイン処理 ---
print(f"コピー元: '{SOURCE_DIRECTORY}'")
print(f"コピー先: '{DESTINATION_DIRECTORY}'")
print("-" * 40)

# コピー先ディレクトリが存在しない場合は作成
try:
    os.makedirs(DESTINATION_DIRECTORY, exist_ok=True)
except OSError as e:
    print(f"致命的なエラー: コピー先ディレクトリの作成に失敗しました。\n理由: {e}")
    # ディレクトリが作成できない場合は処理を中断
    exit()

# リスト内のファイルを一つずつコピー
for file_name in FILES_TO_COPY:
    copy_single_file(SOURCE_DIRECTORY, DESTINATION_DIRECTORY, file_name)

print("すべてのコピー処理が完了しました。")