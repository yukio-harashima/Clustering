#!/bin/bash

# --- 設定 ---
# 入力データファイル名
INFILE="schmidt_projection.dat"
#前提データファイル名
INFILE2="clustering_vector.dat"
# データ生成用Pythonスクリプト
GEN_SCRIPT="gen_schmidt_pr.py"
# 出力画像ファイル名
OUTFILE="ffi_schmidt_plot_v2"
# 【設定】図の縮尺 (1単位あたりのcm)
SCALE=5

# --- 1. データファイルの存在確認と自動生成 ---
if [ ! -f "$INFILE" ]; then
  echo "警告: データファイル $INFILE が見つかりません。"
    if [ ! -f "$INFILE2" ]; then
      echo "警告: データファイル $INFILE2 が見つかりません。"
      gen_clustering_vector.py
      gen_schmidt_pr_v2.py
    else
      echo "スクリプト $GEN_SCRIPT を実行してデータを生成します..."
      # Pythonスクリプトを実行 (環境に合わせて python3 に変更してください)
      gen_schmidt_pr_v2.py
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

if [ -z "$labelmax" ] || [ $(echo "$labelmax < 1" | bc -l) -eq 1 ]; then
    labelmax=10
    echo "Info: labelmaxをデフォルト値(10)に設定しました。"
else
    echo "Label Max: ${labelmax}"
fi

# --- 3. カラーパレット (CPT) の作成 ---
cpt_max=$(echo "$labelmax + 0.5" | bc)
gmt makecpt -Cvik -T0.5/${cpt_max}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > vector_color.cpt

# --- 4. プロット開始 ---
gmt begin $OUTFILE
  gmt figure $OUTFILE png A+m0.5c

  # サブプロット開始 (1行2列)
  gmt subplot begin 1x2 -Fs12c/12c -M0.1c -R-1.2/1.2/-1.2/1.2 -B+n

    # === 左側 (0): n vector ===
    gmt subplot set 0
    
    # タイトルを表示
    echo 0 1.25 "n vector" | gmt text -F+f18p,Helvetica-Bold,black -N

    # 外枠の円
    echo 0 0 | gmt plot -Sc$((${SCALE}*2))c -W1p,black

    # 方位記号
    echo 0 1.1 N | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 1.1 0 E | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 0 -1.1 S | gmt text -F+f12p,Helvetica-Bold,black -N
    echo -1.1 0 W | gmt text -F+f12p,Helvetica-Bold,black -N

    # nベクトル (丸) のプロット
    # 1. sort -k8 -n : 8列目(sliprate)で昇順(小さい順)にソート
    # 2. awk : n1_X($4), n1_Y($5), label($9) を使用 (label >= 1 のみ)
    sort -k8 -n "$INFILE" | \
    awk -F'\t' '{if ($9 >= 1) print $4, $5, $9}' | \
    gmt plot -Sc0.2c -W0.2p,black -Cvector_color.cpt

    # === 右側 (1): v vector ===
    gmt subplot set 1
    
    # タイトルを表示
    echo 0 1.25 "v vector" | gmt text -F+f18p,Helvetica-Bold,black -N

    # 外枠の円
    echo 0 0 | gmt plot -Sc$((${SCALE}*2))c -W1p,black

    # 方位記号
    echo 0 1.1 N | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 1.1 0 E | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 0 -1.1 S | gmt text -F+f12p,Helvetica-Bold,black -N
    echo -1.1 0 W | gmt text -F+f12p,Helvetica-Bold,black -N

    # vベクトル (四角) のプロット
    # 1. sort -k8 -n : 8列目(sliprate)で昇順(小さい順)にソート
    # 2. awk : n2_X($6), n2_Y($7), label($9) を使用
    sort -k8 -n "$INFILE" | \
    awk -F'\t' '{if ($9 >= 1) print $6, $7, $9}' | \
    gmt plot -Ss0.2c -W0.2p,black -Cvector_color.cpt

  gmt subplot end

  # カラーバー (共通)
  gmt colorbar -Cvector_color.cpt -DJBC+w15c/0.3c+o0c/0c+h -Ba1f1 -B+l"Cluster Number"

gmt end show

# 一時ファイルの削除
rm -f vector_color.cpt gmt.conf gmt.history