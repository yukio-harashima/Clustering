#!/bin/bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=0.8      # To change the size of the mechanism solution, change the value of mecascale.
#---
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
source FileChecker.bash
Fort40Checker
RigidChecker
fgenpdtdis
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 18p
gmt set FONT_ANNOT_PRIMARY 18p
gmt set MAP_FRAME_PEN 1p
gmt set MAP_TICK_PEN 1p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
gmt set PLOT_DEGREE_FORMAT D
#---
info=($(gmt gmtinfo pddis.dat -C))
#+++++  2025/3/30 Yuji Yagi  ++++++++
DX=$(echo "${info[3]}-(${info[2]})+0.01"|bc)
DY=$(echo "${info[1]}-(${info[0]})+0.01"|bc)
echo ${DY} ${DX} 
result=`echo "${DY}/${DX} > 1.0 " | bc`
if [ $result -eq 1 ]; then
  wk_x=$(echo "(${DY}-(${DX}))*0.5+${DY}*0.025"|bc)
  wk_y=$(echo "${DY}*0.05"|bc)
else
  wk_x=$(echo "${DX}*0.05"|bc)
  wk_y=$(echo "(${DX}-(${DY}))*0.5+${DX}*0.025"|bc)
fi
echo ${wk_x} ${wk_y}
xmin=$(echo "scale=2;((${info[2]}-(${wk_x}))/1.)" | bc); xmax=$(echo "scale=2;((${info[3]}+${wk_x})/1.)" | bc)
ymin=$(echo "scale=2;((${info[0]}-(${wk_y}))/1.)" | bc); ymax=$(echo "scale=2;((${info[1]}+${wk_y})/1.)" | bc)
#++++++++++++++++++++++++++++++++++++
r=${xmin}/${xmax}/${ymin}/${ymax}
slipmax=${info[9]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
half=$(echo "scale=${len};${slipmax}/2" | bc)
mapwidth=14; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.6*'${mecascale}'}') 
interval=$(echo "scale=${len};$slipmax/5" | bc)
dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
echo max PD: ${slipmax} m,  interval: ${interval}  m
slipmax=$(echo "${slipmax}*1.1" |bc )


#++++++++++++++++++++++++++++++++++++ 7/17 modernized harahima
# --- 1. 初期設定 ---
# 元となるデータのおおよその描画範囲
r_initial=${r}

# 目標とする図の寸法 (この比率に合わせる)
result=`echo "${DY}/${DX} > 1.0 " | bc`
if [ $result -eq 1 ]; then
  target_width=10  # cm
  # target_height=16.18 # gold
  target_height=14.14213 # silver
  wk_x=$(echo "(${DY}-(${DX}))*0.5+${DY}*0.025"|bc)
  wk_y=$(echo "${DY}*0.05"|bc)
else
  # target_width=16.18  # gold
  target_width=14.14213  # silver

  target_height=10 # cm
  wk_x=$(echo "${DX}*0.05"|bc)
  wk_y=$(echo "(${DX}-(${DY}))*0.5+${DX}*0.025"|bc)
fi

echo "元の描画範囲 (r_initial): ${r_initial}"
echo "目標の図の寸法: ${target_width}cm x ${target_height}cm"


# --- 2. 現状の分析：現在の成り行きの高さを計cd 算 ---
map_corners_initial=$(echo ${r_initial} | awk -F/ '{print $2, $4}')
current_height=$(echo ${map_corners_initial} | gmt mapproject -R${r_initial} -JM${target_width} | awk '{print $2}')
echo "現在の成り行きの高さ: ${current_height} cm"


# --- 3. 比較と判断：緯度(Y)と経度(X)のどちらを拡大するか決定 ---
# 目標の高さと現在の高さを比較
height_diff=`echo "${target_height}/${current_height} > 1.0" | bc`
# if (( $(echo "${height_diff} > 0" | bc -l) )); then
if [ $height_diff -eq 1 ]; then
    # --- ケースA: 現在の高さが足りない -> 南北(緯度)方向を拡大 ---
    echo "判断: 地図が目標より横長です。南北(緯度)方向の範囲を拡大します。"
    height_diff=$(echo "scale=5; ${target_height} - ${current_height}" | bc)
    # 4. 拡大量の計算
    # 緯度1度あたりのcmを計算
    lat_range=$(echo ${r_initial} | awk -F/ '{print $4 - $3}')
    cm_per_degree_lat=$(echo "scale=5; ${current_height} / ${lat_range}" | bc)
    # 不足している高さを緯度量に変換
    lat_padding_degrees=$(echo "scale=5; ${height_diff} / ${cm_per_degree_lat}" | bc)

    # 5. 新しい描画範囲の生成
    # 元の範囲を読み込み
    xmin=$(echo ${r_initial} | awk -F/ '{print $1}')
    xmax=$(echo ${r_initial} | awk -F/ '{print $2}')
    ymin=$(echo ${r_initial} | awk -F/ '{print $3}')
    ymax=$(echo ${r_initial} | awk -F/ '{print $4}')
    # 上下に均等に拡大量を加える
    new_ymin=$(echo "scale=5; ${ymin} - ${lat_padding_degrees} / 2" | bc)
    new_ymax=$(echo "scale=5; ${ymax} + ${lat_padding_degrees} / 2" | bc)
    # 最終的な描画範囲を定義
    r_final="${xmin}/${xmax}/${new_ymin}/${new_ymax}"

else
    # --- ケースB: 現在の高さが過剰 -> 東西(経度)方向を拡大 ---
    echo "判断: 地図が目標より縦長です。東西(経度)方向の範囲を拡大します。"

    # 4. 拡大量の計算
    # 目標の幅と現在の幅から、必要な経度の拡大量を計算
    width_diff=$(echo "scale=5; (${current_height} * ${target_width} / ${target_height}) - ${target_width}" | bc)
    # 経度1度あたりのcmを計算
    lon_range=$(echo ${r_initial} | awk -F/ '{print $2 - $1}')
    cm_per_degree_lon=$(echo "scale=5; ${target_width} / ${lon_range}" | bc)
    # 不足している幅を経度量に変換
    lon_padding_degrees=$(echo "scale=5; ${width_diff} / ${cm_per_degree_lon}" | bc)

    # 5. 新しい描画範囲の生成
    xmin=$(echo ${r_initial} | awk -F/ '{print $1}')
    xmax=$(echo ${r_initial} | awk -F/ '{print $2}')
    ymin=$(echo ${r_initial} | awk -F/ '{print $3}')
    ymax=$(echo ${r_initial} | awk -F/ '{print $4}')
    # 左右に均等に拡大量を加える
    new_xmin=$(echo "scale=5; ${xmin} - ${lon_padding_degrees} / 2" | bc)
    new_xmax=$(echo "scale=5; ${xmax} + ${lon_padding_degrees} / 2" | bc)
    # 最終的な描画範囲を定義
    r_final="${new_xmin}/${new_xmax}/${ymin}/${ymax}"
fi

echo "--- 計算結果 ---"
echo "最終的な描画範囲 (r_final): ${r_final}"
echo "----------------"

r=${r_final}
# --- 検証と描画 ---
# 新しい描画範囲で高さを再計算して、目標値と一致するか確認
map_corners_final=$(echo ${r} | awk -F/ '{print $2, $4}')
final_height=$(echo ${map_corners_final} | gmt mapproject -R${r} -JM${target_width} | awk '{print $2}')
echo "検証: 新しい描画範囲での高さは ${final_height} cm です (目標: ${target_height} cm)"
mapwidth=${target_width}
mapheight=${final_height}
#++++++++++++++++++++++++++++++++++++


#---
gmt begin ffi_MapMeca pdf
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JM${mapwidth} -BWSne -Bx -By  
  gmt coast -Df  -Glightgray
  #--
  #--  Remove "#" if you want to write the distribution.
  awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${interval}
  gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  path_bin="$(dirname "$0")"
  path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  awk '{print $7,$8}' ../aftershock.dat | gmt plot -Sc0.1c -W1p
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W0.5p
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency Density (m)"
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf