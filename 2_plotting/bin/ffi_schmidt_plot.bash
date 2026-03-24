#!/bin/bash

# --- 設定 ---
# 入力データファイル名
INFILE="schmidt_projection.dat"
#前提データファイル名
INFILE2="clustering_vector.dat"
# データ生成用Pythonスクリプト
GEN_SCRIPT="gen_schmidt_pr.py"
# 出力画像ファイル名
OUTFILE="ffi_schmidt_vector_plot"
# 【追加設定】図の縮尺 (1単位あたりのcm)
SCALE=5

# --- 1. データファイルの存在確認と自動生成 ---
if [ ! -f "$INFILE" ]; then
  echo "警告: データファイル $INFILE が見つかりません。"
    if [ ! -f "$INFILE2" ]; then
      echo "警告: データファイル $INFILE2 が見つかりません。"
      gen_clustering_vector.py
      gen_schmidt_pr.py
    else
      echo "スクリプト $GEN_SCRIPT を実行してデータを生成します..."
      # Pythonスクリプトを実行 (環境に合わせて python3 に変更してください)
      gen_schmidt_pr.py
    fi

    # 生成後の再確認
    if [ ! -f "$INFILE" ]; then
      echo "エラー: $GEN_SCRIPT を実行しましたが $INFILE が生成されませんでした。"
      exit 1
    fi
    echo "データ生成完了。"

fi

# --- 2. labelmax の自動取得 ---
# labelは9列目なので -i8 (0始まりのインデックス) を指定
info=($(gmt info "$INFILE" -i8 -C))
labelmax=${info[1]}

# もし最大値が1未満、あるいは取得失敗した場合はデフォルト値を設定
if [ -z "$labelmax" ] || [ $(echo "$labelmax < 1" | bc -l) -eq 1 ]; then
    labelmax=10
    echo "Info: labelmaxをデフォルト値(10)に設定しました。"
else
    echo "Label Max: ${labelmax}"
fi

# --- 3. カラーパレット (CPT) の作成 ---
# 【修正点】最大値(4)を含む区間を作るため、上限を labelmax + 0.5 に設定します
# labelmax=4 の場合、0.5 ～ 4.5 の範囲となり、1, 2, 3, 4 が各区間の中心に来ます。
cpt_max=$(echo "$labelmax + 0.5" | bc)

gmt makecpt -Cvik -T0.5/${cpt_max}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > vector_color.cpt

# --- 4. プロット開始 ---
# 【修正】png形式で出力し、上下左右に0.5cmの余白(A+m0.5c)を設定
gmt begin $OUTFILE
  gmt figure $OUTFILE png A+m0.5c

  # ベースマップの設定 (線形投影 -Jx)
  # 範囲は半径1の円が収まる -1.2 ～ 1.2 に設定
  # 【修正】縮尺を変数 SCALE で指定
  # 【修正】外枠線を描画しないため -B0 を -B+n に変更
  gmt basemap -R-1.2/1.2/-1.2/1.2 -Jx${SCALE}c -B+n

  # --- 外枠と方位の描画 ---
  # 半径1（直径2）の円を描画
  # 【修正】縮尺(5cm) x 直径(2) = 10cm のサイズを指定
  echo 0 0 | gmt plot -Sc$((${SCALE}*2))c -W1p,black

  # 方位記号 (N, E, S, W)
  echo 0 1.1 N | gmt text -F+f12p,Helvetica-Bold,black -N
  echo 1.1 0 E | gmt text -F+f12p,Helvetica-Bold,black -N
  echo 0 -1.1 S | gmt text -F+f12p,Helvetica-Bold,black -N
  echo -1.1 0 W | gmt text -F+f12p,Helvetica-Bold,black -N

  # --- nベクトルのプロット (丸) ---
  # 変数: 4:n1_X, 5:n1_Y, 9:label, 8:sliprate
  # sort -k8 -n : 8列目(sliprate)で昇順ソート
  # awk : label($9) >= 1 のデータを抽出
  sort -k8 -n "$INFILE" | \
  awk -F'\t' '{if ($9 >= 1) print $4, $5, $9}' | \
  gmt plot -Sc0.2c -W0.2p,black -Cvector_color.cpt

  # --- vベクトルのプロット (四角) ---
  # 変数: 6:n2_X, 7:n2_Y, 9:label, 8:sliprate
  sort -k8 -n "$INFILE" | \
  awk -F'\t' '{if ($9 >= 1) print $6, $7, $9}' | \
  gmt plot -Ss0.2c -W0.2p,black -Cvector_color.cpt

  # --- 仕上げ ---
  # 凡例の追加 (手動設定)
  gmt legend -DjBR+o0.1c/0.5c -F+gwhite+p0.5p << EOF
S 0.2c c 0.2c white 0.2p,black 0.4c n vector
S 0.2c s 0.2c white 0.2p,black 0.4c v vector
EOF

  # カラーバーの追加
  # 【修正】位置を下にずらし、補助目盛りを削除
  # -DJBC...: 図の外側(J)の下中央(BC)基準。+o0c/1.5c で下に1.5cmオフセット。
  # -Ba1f1: アノテーション(a)と目盛り(f)を1間隔に設定（補助目盛りなし）。
  # -B+l...: タイトルを設定（数値の下に表示）。
  gmt colorbar -Cvector_color.cpt -DJBC+w8c/0.3c+o0c/0c+h -Ba1f1 -B+l"Cluster Number"

gmt end show

# 一時ファイルの削除
rm -f vector_color.cpt gmt.conf gmt.history