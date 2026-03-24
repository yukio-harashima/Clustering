#!/bin/bash
# 複数断面対応のffi_MapMeca.bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=1.8      # To change the size of the mechanism solution, change the value of mecascale.
#---
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
source FileChecker.bash
Fort40Checker
RigidChecker
fgenpdtdis_mf

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
echo ${r}
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

#--- 1. カラーパレットの準備 (gmt beginの前) ---
# 地形用のCPTファイルを作成
gmt grdcut @earth_relief_15s -R120/150/20/50 -Gdem.nc
gmt makecpt -Cgray -T-10000/0/100 -Z > topo.cpt
# 断層すべり用のCPTファイルを作成
gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > slip.cpt

#---
gmt begin ffi_MapMeca_Mf png
  # gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > slip.cpt
  gmt basemap -R${r} -JM${mapwidth} -BWSne -Bx -By -Ba0.5 
  # gmt coast -Df  -Glightgray
  # 地形用プロット
  gmt grdgradient dem.nc -Ggrad.grd -A45 -Ne0.8
  gmt grdimage dem.nc -Igrad.grd -Ctopo.cpt -t20
  gmt coast -Df -W0.25 -Glightgray

  #--
  #--  Remove "#" if you want to write the distribution.
  # awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  # gmt grdimage tmp.grd -C  -t35
  # gmt grdcontour tmp.grd -C${interval}
  # gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -Cslip.cpt  # -t20

  
  # 1. 全体を細い線で描画
  # awk '{if($0 ~ />/) print ">"; else if(NF>=2) print $2,$1}' faultline.dat | gmt plot -W0.6

  # # 2. 断層の上端（各セグメントの浅い2点）を太線で上書き描画
  # ---------------------------------------------------------
  # 工程1: upper.dat の作成
  # ---------------------------------------------------------

  # # 1. 【重要】データのクリーニング
  # #    元データの区切り行（>を含む行）を正規化して、一時ファイル clean_fault.dat に保存
  # awk '{if($0 ~ />/) print ">"; else print $0}' faultline.dat > clean_fault.dat

  # # 2. 分割
  # #    きれいになった clean_fault.dat を使って分割を実行
  # gmt convert clean_fault.dat -D"tmp_seg_%d.txt"

  # # 3. upper.dat を初期化
  # rm -f upper.dat

  # # 4. 分割したファイルをループ処理
  # #    (tmp_seg_*.txt ではなく tmp_seg_* と指定することで生成された全ファイルを対象にする)
  # for file in tmp_seg_*; do
  #   # セグメント開始記号を書き込み
  #   echo ">" >> upper.dat
    
  #   # 並べ替えと抽出
  #   # 3列目でソート(昇順) -> 上位2行(浅い点) -> 座標入れ替え
  #   sort -k3 -n "$file" | head -n 2 | awk '{print $2, $1}' >> upper.dat
  # done

  # # 5. お掃除 (一時ファイルの削除)
  # rm -f clean_fault.dat tmp_seg_*

  # ---------------------------------------------------------
  # 工程2: upper.dat の描画
  # ---------------------------------------------------------

  # 作成した upper.dat を太線で描画
  gmt plot upper.dat -W1.0


  path_bin="$(dirname "$0")"
  path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  # gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  # gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  # awk '{print $7,$8}' ../aftershock.dat | gmt plot -Sc0.1c -W1p
  # awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W0.5p
  gmt colorbar -Cslip.cpt -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency Density (m)"
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf


