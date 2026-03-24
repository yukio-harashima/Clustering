#!/bin/bash

#--- 出力ファイル名
outfile="stereo_vector_plot"

#--- データファイル
datafile="n_vector.dat"

#--- スタイルの設定
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 18p
gmt set FONT_ANNOT_PRIMARY 18p
gmt set MAP_FRAME_PEN 1p
gmt set MAP_TICK_PEN 1p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
gmt set PLOT_DEGREE_FORMAT D

#--- GMTスクリプト開始
gmt begin ${outfile} pdf

  #--- ステレオネットの設定
  PROJ="-Jp1c"
  # 修正点：仰角の範囲を-90から90に広げる
  # REGION="-R0/360/-90/90"
  REGION="-R0/360/-90/90"
  FRAME="-Bpa3g3 -Bnwse"

  #--- ステレオネットの基盤（グリッド）を描画
  gmt basemap ${PROJ} ${REGION} ${FRAME} --MAP_GRID_PEN_PRIMARY=0.5p,gray

  #--- 方角のラベルを追加 ---
  gmt text ${PROJ} ${REGION} -F+f16p,Helvetica-Bold+jMC << EOF
0   -8  N
90  -8  E
180 -8  S
270 -8  W
EOF

  #--- ベクトルデータをプロット ---

  # 1. 下半球のベクトル (z >= 0) を塗りつぶし丸でプロット
  # awkで行を選択し、gmt mathで計算する方式に変更
  awk '$9 >= 0 {print $7, $8, $9}' ${datafile} | \
  gmt math -o 0,1 \
      T 0 COL 1 COL ATAN2 R2D DUP 0 LT { 360 ADD } IFT = \
      T 2 COL ASIN R2D = \
  | gmt plot ${PROJ} ${REGION} -Sc0.2c -Gblack

  # 2. 上半球のベクトル (z < 0) を白抜き丸でプロット
  # awkで行を選択し、gmt mathで計算する方式に変更
  awk '$9 < 0 {print $7, $8, $9}' ${datafile} | \
  gmt math -o 0,1 \
      T 0 COL 1 COL ATAN2 R2D DUP 0 LT { 360 ADD } IFT = \
      T 2 COL ASIN R2D NEG = \
  | gmt plot ${PROJ} ${REGION} -Sc0.2c -W0.5p,black

gmt end show

#--- 生成された一時ファイルを削除
rm -rf gmt.history gmt.conf