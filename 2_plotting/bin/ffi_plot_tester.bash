#!/bin/bash

# --- 変数の設定 ---
# 地図の幅 (例: 15cm)
mapwidth=15


# --- データ範囲の自動取得 ---
info=($(gmt gmtinfo Multiple_fault.dat -C))
xmin=$(echo "${info[4]}-10" | bc); xmax=$(echo "${info[5]}+10" | bc)
ymin=$(echo "${info[6]}-10" | bc); ymax=$(echo "${info[7]}+10" | bc)
# -I0.1/0.1 は、データの最小/最大値から少しマージン（余白）をとる設定
r=${xmin}/${xmax}/${ymin}/${ymax}
echo ${r}
# r=-100/100/-100/100

# --- GMTスクリプト開始 ---
gmt begin ffi_plot_Mf png

  # ベースマップの描画
  # -R${r}: gmt infoで取得したデータ範囲
  # -JM${mapwidth}: メルカトル図法、幅を指定
  # -BWSne: 外枠と目盛りを描画 (W:西, S:南, N:北, E:東)
  # -Bxaf+l"X軸 (3列目)": X軸の目盛り(a)とグリッド(f)、ラベル(+l)を設定
  # -Byaf+l"Y軸 (4列目)": Y軸の目盛り(a)とグリッド(f)、ラベル(+l)を設定
  # -B+t"Multiple_fault.dat プロット": プロット全体のタイトルを設定
  gmt basemap -R${r} -JX${mapwidth} -BWSne -Ba0.8
  
  # データのプロット
  # awkで3列目(x)と4列目(y)を抽出してgmt plotに渡す
  awk '{print $3,$4}' Multiple_fault.dat | gmt plot -Sc0.15 -W0.4

gmt end show

# --- 一時ファイルの削除 ---
rm -rf gmt.history gmt.conf