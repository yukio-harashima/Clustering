#!/bin/bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=1.4      # To change the size of the mechanism solution, change the value of mecascale.
#---
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
source FileChecker.bash
Fort40Checker
RigidChecker
if [ $# -ne 3 ]; then
  echo -e \\n " MapMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1]" 
  fgenmrf
  info=($(gmt gmtinfo mrf.dat -C))
  tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
  interval=$(echo "scale=0;${tmax}/10." | bc)
  if [ $interval -le 0 ]; then
      interval=1
  fi
  echo -e  " Default fun: ffi_MapMecaSnap.bash ${interval}  ${tmax}  0 " \\n
  echo ${interval}  ${tmax}  0  | fgensnap  > .tmpout  ; rm -rf .tmpout
  i1=${interval}; i2=${tmax}; i3=0
else
  echo ${1} ${2} ${3} | fgensnap
  i1=${1}; i2=${2}; i3=${3}
fi
#---
gmt set MAP_TICK_LENGTH 1p
gmt set FONT_LABEL 15p
gmt set FONT_ANNOT_PRIMARY 15p
gmt set MAP_FRAME_PEN 0.4p
gmt set MAP_TICK_PEN 0.4p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 1p
# gmt set PLOT_DEGREE_FORMAT D
gmt set FORMAT_GEO_MAP D
#---
info=($(gmt gmtinfo snap2.dat -C))
#+++++  2025/3/30 Yuji Yagi  ++++++++
DX=$(echo "${info[15]}-(${info[14]})+0.01"|bc)
DY=$(echo "${info[13]}-(${info[12]})+0.01"|bc)
echo ${DY} ${DX} 
result=`echo "${DY}/${DX} > 1.0 " | bc`
if [ $result -eq 1 ]; then
  wk_x=$(echo "(${DY}-${DX})*0.5+${DY}*0.025"|bc)
  wk_y=$(echo "${DY}*0.05"|bc)
else
  wk_x=$(echo "${DX}*0.05"|bc)
  wk_y=$(echo "(${DX}-${DY})*0.5+${DX}*0.025"|bc)
fi
xmin=$(echo "scale=2;((${info[14]}-${wk_x})/1.)" | bc); xmax=$(echo "scale=2;((${info[15]}+${wk_x})/1.)" | bc)
ymin=$(echo "scale=2;((${info[12]}-${wk_y})/1.)" | bc); ymax=$(echo "scale=2;((${info[13]}+${wk_y})/1.)" | bc)
#++++++++++++++++++++++++++++++++++++
# r=${xmin}/${xmax}/${ymin}/${ymax}
r=138.34/140.19/41.54/43.48
num_tw=${info[5]}
slipmax=${info[11]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
half=$(echo "scale=${len};${slipmax}/2" | bc)
interval=$(echo "scale=2;$slipmax/5" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
column=10; row=$(echo "scale=2;${num_tw}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')
mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.6*'${mecascale}'}' )
dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
interval=$(echo "scale=2;$slipmax/5" | bc)
echo max PRD: ${slipmax} m/s,  interval: ${interval}  m/s
slipmax=$(echo "${slipmax}*1.1" |bc )
#---
path_bin="$(dirname "$0")"
path_data=${path_bin/bin/data}
#---
for tw in $(seq 1 ${num_tw}); do echo -n "+"; done; echo -e " " ${num_tw}
gmt begin ffi_MapMecaSnap png
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt subplot begin ${row}x${column} -Fs${mapwidth} -M0.10/0.3 -SCb -SRl -Bxa0.5 -Bya0.5 -BWSne -R${r} -JM${mapwidth}
    for tw in $(seq 1 ${num_tw}); do
      echo -n "+"
      gmt subplot set
      #---
      awk '{if ($3=='${tw}') print $8,$7,$6}' snap2.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r}
      gmt grdimage tmp.grd -C  -t35.
      gmt grdcontour tmp.grd -C${interval}
      gmt coast -Df -W0.25 -Glightgray
      # awk '{print $2,$1}' faultline.dat | gmt plot
      if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
      awk '{if ($3=='${tw}' &&  $6 > 0.00001) print $20,$19,$6,$7,$8,$9,$10,$11,$12,23}' snap.dat |sort -k 3 -n  | gmt meca -Sm${mecasize}c+l -T0 -C
      #---
      
#      awk '{print $2,$1}' faultline.dat | gmt plot
#      if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
      # gmt set IO_SEGMENT_MARKER "***"
      # gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
      # gmt set IO_SEGMENT_MARKER ">"
      # gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
      awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.2c  
      #---
      t1=$(echo "scale=1;(${tw}-1)*${i1}" | bc)
      t2=$(echo "scale=1;${tw}*${i1}" | bc)
      if [ ${i3} -eq 1 ]; then sec=${t2}" s"; else sec=${t1}"-"${t2}" s"; fi
      gmt text -Dj0.0/-0.35 -F+cLT+jTL+t"${sec}" -N
      #---
      if [ ${tw} -eq ${num_tw} ]; then gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency-rate Density (m/s)"; fi
    done
  gmt subplot end
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
echo -e " done !"
rm -rf gmt.history gmt.conf