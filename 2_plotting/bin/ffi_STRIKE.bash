#!/bin/bash
#---
source FileChecker.bash
Fort40Checker
RigidChecker
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
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 11p
gmt set FONT_ANNOT_PRIMARY 11p
gmt set MAP_FRAME_PEN 0.5p
gmt set MAP_TICK_PEN 0.5p
gmt set PAGE_COLOR 254/254/254

#---
aspect=2.0 # larger value means higher
strinfo=($(gmt gmtinfo vr_str.dat -C))
dipinfo=($(gmt gmtinfo vr_dip.dat -C))
strmin=${strinfo[0]}; strmax=${strinfo[1]}
dipmin=${dipinfo[0]}; dipmax=${dipinfo[1]}
tmax=${strinfo[3]}
# strr=0/${tmax}/${strmin}/${strmax}; dipr=0/${tmax}/${dipmin}/${dipmax}
strr=0/${tmax}/${strmin}/${strmax}; dipr=0/${tmax}/${dipmin}/${dipmax}

ratio=$(echo " scale=4;${tmax}/(${strmax}-(${strmin}))" | bc)
strwidth=10; mapheight=$(echo "${strwidth}*${ratio}*${aspect}" | bc)
ratio=$(echo " scale=4;${tmax}/(${dipmax}-(${dipmin}))" | bc)
dipheight=$(echo "(${mapheight}/${ratio})/${aspect}" | bc)
echo ${strinfo[5]} > tmp; echo ${dipinfo[5]} >> tmp; info=($(gmt gmtinfo tmp -C))
slipmax=${info[1]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
interval=$(echo "scale=2;$slipmax/5" | bc)
half=$(echo "scale=2;$slipmax/2" | bc)
xmv=$(echo "${strwidth}+0.05" | bc)

#---
cororheight=$(echo "${mapheight}+${dipheight}*1.2" | bc)
echo cororheight${cororheight}, mapheight${mapheight}, dipheight ${dipheight}

# gold ratio
# widthratio=1.618
# shilver ratio
widthratio=1.414
rat=$(echo "scale=5;${strwidth}/(${widthratio}*${cororheight})"|bc)
cororheight=$(echo "${rat}*${cororheight}" | bc)
mapheight=$(echo "${strwidth}/${widthratio}" | bc)

dipheight=$(echo "${rat}*${dipheight}" | bc)
echo cororheight${cororheight}, mapheight${mapheight}, dipheight ${dipheight}

yshift=$(echo "${dipheight}*1.2" | bc)

#---
gmt begin ffi_RSD png
  # strike
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${strr} -JX${strwidth}/-${mapheight} -BWSne -Bx+l"Time (s)" -By+l"Strike (km)"
  awk '{print $1,$2,$3}' vr_str.dat | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${strr}
  gmt grdimage tmp.grd -C
  gmt grdcontour tmp.grd -C${interval}
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency-rate Density (m/s)"
  gmt basemap -Bnswe
  
  # dip
  # gmt basemap -R${dipr} -JX${strwidth}/${dipheight} -BWSne -Bx+l"Time (s)" -By+l"Dip (km)" -Y-${yshift}
  # awk '{print $1,$2,$3}' vr_dip.dat | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${dipr}
  # gmt grdimage tmp.grd -C
  # gmt grdcontour tmp.grd -C${interval}

  # gmt basemap -Bnswe
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf tmp tmp.grd
