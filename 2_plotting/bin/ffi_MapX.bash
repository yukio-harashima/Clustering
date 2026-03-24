#!/bin/bash
if [ $# = 0 ]; then
    echo -e \\n " ffi_MapX.bash [ s p t b ]" \\n
    exit 1
fi
p_s=false; p_p=false; p_t=false; p_b=false
while (( $# > 0 ))
do
  case $1 in
    "s" | "S" ) p_s=true ;;
    "p" | "P" ) p_p=true ;;
    "t" | "T" ) p_t=true ;;
    "b" | "B" ) p_b=true ;;
    * ) echo -e \\n " Invalid input" \\n; exit 1 ;;
  esac
  shift
done
#---
source FileChecker.bash
Fort40Checker
RigidChecker
fgenpdtdis
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 11p
gmt set FONT_ANNOT_PRIMARY 11p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
#---
info=($(gmt gmtinfo pddis.dat -C))
xmin=$(echo "${info[2]}-0.15" | bc); xmax=$(echo "${info[3]}+0.15" | bc)
ymin=$(echo "${info[0]}-0.15" | bc); ymax=$(echo "${info[1]}+0.15" | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
slipmax=${info[9]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
half=$(echo "scale=2;$slipmax/2" | bc)
mapwidth=12; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
#---
gmt begin ffi_MapX pdf
  gmt makecpt -Cplasma -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JM${mapwidth} -BWSne  -Baf
  if "${p_s}"; then # plot strike
    awk '{print $6,$5,$8,$16,$8/'${slipmax}'*0.75}' pdtdis.dat | gmt plot -SV0c+jc -W1.5+cl -C
    awk '{print $6,$5,$8,$17,$8/'${slipmax}'*0.75}' pdtdis.dat | gmt plot -SV0c+jc -W1.5+cl -C
  fi
  if "${p_p}"; then # plot P axis
    sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$22,$8/'${slipmax}'*0.75}' | gmt plot -SV0c+jc -W1.5+cl -C
  fi
  if "${p_t}"; then # plot T axis
    sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$23,$8/'${slipmax}'*0.75}' | gmt plot -SV0c+jc -W1.5+cl -C
  fi
  if "${p_b}"; then # plot B axis
    sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$24,$8/'${slipmax}'*0.75}' | gmt plot -SV0c+jc -W1.5+cl -C
  fi
  gmt coast -Df -W0.25
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.3c -W1p
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.15/0.0 -Ba${half}+l"Potency Density (m)"
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
gmt end
rm -rf gmt.history gmt.conf
