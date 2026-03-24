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
#---
gmt begin ffi_MapMeca pdf
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JM${mapwidth} -BWSne -Bx -By -Ba0.5 
  gmt coast -Df  -Glightgray
  #--
  #--  Remove "#" if you want to write the distribution.
  awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${interval}
  gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.6
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  path_bin="$(dirname "$0")"
  path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  # gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  # gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  # awk '{print $7,$8}' ../aftershock.dat | gmt plot -Sc0.1c -W1p
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W0.5p
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency Density (m)"
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf