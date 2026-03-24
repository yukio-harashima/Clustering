#!/bin/bash
#---
# ffi_RupStrDip.bashの改良版
mecascale=1.0 

source FileChecker.bash
Fort40Checker
RigidChecker
#---
#----
fgenmrf
fgentpdt
fgenpdtdis
#---
if [ $# -ne 1 ] || [ $1 -le 0 ]; then
    echo -e \\n " RupStrDip.bash  [duration (duration > 0.)]" 
    echo -e     "   default run" \\n
    vr=$(awk 'NR==2 {print $7}' fort.40 )
    duration=$(awk 'NR==6 {print ($1*$3)/'${vr}'+$7*$8}' fort.40 )
    echo ${duration} | fgenstrdiprup
else
    echo $1 | fgenstrdiprup
fi

#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 15p
gmt set FONT_ANNOT_PRIMARY 15p
gmt set MAP_FRAME_PEN 1p
gmt set MAP_TICK_PEN 1p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
gmt set PLOT_DEGREE_FORMAT D
gmt set PAGE_COLOR 254/254/254
#---
info=($(gmt gmtinfo pddis.dat -C))

xmin=$(echo "${info[2]}-0.5" | bc); xmax=$(echo "${info[3]}+0.5" | bc)
ymin=$(echo "${info[0]}-0.1" | bc); ymax=$(echo "${info[1]}+0.1" | bc)
xmin=$(echo "scale=2;${xmin}/1." | bc); xmax=$(echo "scale=2;${xmax}/1." | bc)
ymin=$(echo "scale=2;${ymin}/1." | bc); ymax=$(echo "scale=2;${ymax}/1." | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
mslipmax=${info[9]}
wk=$(echo " scale=0;${mslipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
mhalf=$(echo "scale=${len};${mslipmax}/2" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
mecasize=$(echo ${mslipmax} | awk '{print 5/$1*0.6*'${mecascale}'}') 
minterval=$(echo "scale=${len};$mslipmax/5" | bc)
dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
echo max PD: ${mslipmax} m,  interval: ${minterval}  m
mslipmax=$(echo "${mslipmax}*1.1" |bc )

#---
#---
aspect=2.00 # larger value means higher
strinfo=($(gmt gmtinfo vr_str.dat -C))
dipinfo=($(gmt gmtinfo vr_dip.dat -C))
strmin=${strinfo[0]}; strmax=${strinfo[1]}
dipmin=${dipinfo[0]}; dipmax=${dipinfo[1]}
tmax=${strinfo[3]}
strr=0/${tmax}/${strmin}/${strmax}; dipr=0/${tmax}/${dipmin}/${dipmax}

ratio=$(echo " scale=4;${tmax}/(${strmax}-(${strmin}))" | bc)
strwidth=10; mapheight=$(echo "${strwidth}*${ratio}*${aspect}" | bc)
ratio=$(echo " scale=4;${tmax}/(${dipmax}-(${dipmin}))" | bc)
dipwidth=$(echo "(${mapheight}/${ratio})/${aspect}" | bc)

echo ${strinfo[5]} > tmp; echo ${dipinfo[5]} >> tmp; info=($(gmt gmtinfo tmp -C))
slipmax=${info[1]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
interval=$(echo "scale=2;$slipmax/5" | bc)
half=$(echo "scale=2;$slipmax/2" | bc)
xmv=$(echo "${strwidth}" | bc)
cororheight=$(echo "${mapheight}+${dipwidth}*1.2" | bc)

#---
# yshift=$(echo "-(${mapheight}+3)" | bc)
# yshift=$(echo "${mapheight}+${dipwidth}*1.2" | bc)
yshift=$(echo "${dipwidth}*1.2" | bc)

#---
#---
gmt begin ffi_MRSD pdf
  gmt makecpt -Chot -T0/${mslipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  # # map
  gmt basemap -R${r} -JX${mapwidth}/${cororheight} -BWSne -Ba1
  gmt coast -Df  -Glightgray
  #--
  #--  Remove "#" if you want to write the distribution.
  awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${minterval}
  gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  # path_bin="$(dirname "$0")"
  # path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  # gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"s
  # gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W1p
  gmt basemap -Bnswe
  gmt colorbar -DjBR+w${cororheight}/0.2+jBL+o0.5/0.0 -Ba${mhalf}+l"Potency-rate Density (m/s)"

  # strike
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${strr} -JX${strwidth}/-${mapheight} -BWSne -By+l"Strike (km)" -X${xmv} -Y${yshift}
  awk '{print $1,$2,$3}' vr_str.dat | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${strr}
  gmt grdimage tmp.grd -C
  gmt grdcontour tmp.grd -C${interval}
  gmt basemap -Bnswe

  # dip
  # gmt basemap -R${dipr} -JX${strwidth}/${dipwidth} -BWSne -Bx+l"Time (s)" -By+l"Dip (km)" -X${xmv}
  gmt basemap -R${dipr} -JX${strwidth}/${dipwidth} -BWSne -Bx+l"Time (s)" -By+l"Dip (km)"   -Y-${yshift}
  awk '{print $1,$2,$3}' vr_dip.dat | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${dipr}
  gmt grdimage tmp.grd -C
  gmt grdcontour tmp.grd -C${interval}
  gmt colorbar -DjBR+w${cororheight}/0.2+jBL+o0.5/0.0 -Ba${half}+l"Potency-rate Density (m/s)"
  gmt basemap -Bnswe
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.5
gmt end show
rm -rf gmt.history gmt.conf tmp tmp.grd
