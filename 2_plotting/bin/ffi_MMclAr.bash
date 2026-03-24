#!/bin/bash
# 地図上にクラスターごとの分布を点で示す。

#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=8     # To change the size of the mechanism solution, change the value of mecascale.
#---
source FileChecker.bash
Fort40Checker
RigidChecker
fgenpdtdis
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 20p
gmt set FONT_ANNOT_PRIMARY 20p
gmt set MAP_FRAME_PEN 1p
gmt set MAP_TICK_PEN 1p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
gmt set PLOT_DEGREE_FORMAT D
#---
info=($(gmt gmtinfo pddis.dat -C))

xmin=$(echo "${info[2]}-0.7" | bc); xmax=$(echo "${info[3]}+0.7" | bc)
ymin=$(echo "${info[0]}-0.1" | bc); ymax=$(echo "${info[1]}+0.1" | bc)
xmin=$(echo "scale=2;${xmin}/1." | bc); xmax=$(echo "scale=2;${xmax}/1." | bc)
ymin=$(echo "scale=2;${ymin}/1." | bc); ymax=$(echo "scale=2;${ymax}/1." | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
slipmax=${info[9]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
half=$(echo "scale=${len};${slipmax}/2" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.6*'${mecascale}'}') 
interval=$(echo "scale=${len};$slipmax/5" | bc)
dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
echo max PD: ${slipmax} m,  interval: ${interval}  m
slipmax=$(echo "${slipmax}*1.1" |bc )

echo ${r}
#---
gmt begin ffi_MapMecaArea pdf
  gmt makecpt -Cvik -T0.5/6.5/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JM${mapwidth} -BWSne -Ba1
  gmt coast -Df -G200
  #--
  #--  Remove "#" if you want to write the distribution.
  # awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  #gmt grdimage tmp.grd -C  -t35
  # gmt grdcontour tmp.grd -C${interval}
  # gmt coast -Df -W0.25
  #--
  # sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
  #awk -F'\t' '{OFS="\t"; print $0}' clarea.dat | sort -k 6 -n -r > tmp.dat
  awk -F'\t' '{OFS="\t"; print $0}' clarea.dat | sort -k 30 -n -r > tmp.dat
  #awk -F'\t' '{ print $20,$19,$30,$7,$8,$9,$10,$11,$12,23}' tmp.dat | gmt meca -Sm0.2 -T0 -C
  #awk -F'\t' '{ print $20,$19,$30,$7,$8,$9,$10,$11,$12,23}' tmp.dat | gmt meca -Sm${mecasize}c+l -T0 -C

  #awk -F'\t' '{OFS="\t"; print $20,$19,$30,$6}' clarea.dat | sort -k 30 -n -r | sort -k 3 -n -r > _tmp.dat
  awk -F'\t' '{OFS="\t"; print $20,$19,$30,$6}' clarea.dat | sort -k 6 -n -r | sort -k 4 -n > _tmp.dat
  #awk -F'\t' '{ print $1,$2,$3,$4*0.5}' _tmp.dat | gmt plot -Sc -W -C -t10 
  awk -F'\t' '{ print $1,$2,$3,$4*0.4}' _tmp.dat | gmt plot -Sc -W -C -t10 

  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  path_bin="$(dirname "$0")"
  path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  # gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  # gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W1p
  # gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba1+l"Cluster Number"
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf #tmp.dat
