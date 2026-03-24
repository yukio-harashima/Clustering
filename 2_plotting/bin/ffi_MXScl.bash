#!/bin/bash
#---
source FileChecker.bash
Fort40Checker
RigidChecker
#---
if [ $# -lt 4 ]; then
  echo -e \\n " ffi_MapXSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1] [ s p t b ]" \\n
  fgenmrf
  info=($(gmt gmtinfo mrf.dat -C))
  tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
  dur=$(echo "scale=0;${tmax}/10." | bc)
  if [ $dur -le 0 ]; then
      dur=1
  fi
  echo -e  " Default fun: ffi_MapXSnap.bash ${dur}  ${tmax}  0  p " \\n
  typ=0; axi=p
else
  
  dur=${1}; tmax=${2}; typ=${3}; axi=${4}
fi
#---
#---
p_s=false; p_p=false; p_t=false; p_b=false
case ${axi} in
  "s" | "S" ) p_s=true ;;
  "p" | "P" ) p_p=true ;;
  "t" | "T" ) p_t=true ;;
  "b" | "B" ) p_b=true ;;
esac
#---
echo $dur $tmax $typ | fgensnap
swapper.py
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 9p
gmt set FONT_ANNOT_PRIMARY 10p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
#---
info=($(gmt gmtinfo snap2.dat -C))
xmin=$(echo "${info[14]}-0.1" | bc); xmax=$(echo "${info[15]}+0.1" | bc)
ymin=$(echo "${info[12]}-0.1" | bc); ymax=$(echo "${info[13]}+0.1" | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
num_tw=${info[5]}
slipmax=${info[11]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
interval=$(echo "scale=2;$slipmax/5" | bc)
half=$(echo "scale=2;$slipmax/2" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
column=5; row=$(echo "scale=2;${num_tw}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')
#---


info2=($(gmt gmtinfo clusteringinfo.dat -C))
labmax=$(echo "${info2[1]}+0.5" | bc)
echo ${labmax}


for tw in $(seq 1 ${num_tw}); do echo -n "+"; done; echo -e " " ${num_tw}
gmt begin ffi_MapXSnap_${4} pdf
  # gmt makecpt -Chawaii -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt makecpt -Cvik -T-0.5/${labmax}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt subplot begin ${row}x${column} -Fs${mapwidth} -M0.15/0.3 -SCb -SRl -Bx -By -Bwsne -R${r} -JM${mapwidth}
    for tw in $(seq 1 ${num_tw}); do
      echo -n "+"
      gmt subplot set

      # awk  -F'\t' '{if ($3=='${tw}' && $6 > 0.00001) print $8,$7,$10}' assignsnap2KDE_v2.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r}
      # gmt grdimage tmp.grd -C  -t35.
      #---
      if "${p_s}"; then # plot strike
        awk '{if ($1=='${tw}') print $3,$4,$2,$12,$6/'${slipmax}'*0.4}' centroids2_v2.dat | gmt plot -SV0c+jc -W1.5+cl -C
        awk '{if ($1=='${tw}') print $3,$4,$2,$13,$6/'${slipmax}'*0.4}' centroids2_v2.dat | gmt plot -SV0c+jc -W1.5+cl -C
      fi
      if "${p_p}"; then # plot P axis
        awk '{if ($1=='${tw}') print $3,$4,$2,$18,$6/'${slipmax}'*0.4}' centroids2_v2.dat | gmt plot -SV0c+jc -W1.5+cl -C
      fi
      if "${p_t}"; then # plot T axis
        awk '{if ($1=='${tw}') print $3,$4,$2,$19,$6/'${slipmax}'*0.4}' centroids2_v2.dat | gmt plot -SV0c+jc -W1.5+cl -C
      fi
      if "${p_b}"; then # plot B axis
        awk '{if ($1=='${tw}') print $3,$4,$2,$20,$6/'${slipmax}'*0.4}' centroids2_v2.dat | gmt plot -SV0c+jc -W1.5+cl -C
      fi
      #---
      gmt coast -Df -W0.25
      awk '{print $2,$1}' faultline.dat | gmt plot
      awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.2c
      #---
      t1=$(echo "scale=1;(${tw}-1)*$dur" | bc)
      t2=$(echo "scale=1;${tw}*$dur" | bc)
      if [ $typ -eq 1 ]; then sec=${t2}" s"; else sec=${t1}"-"${t2}" s"; fi
      gmt text -Dj0.0/-0.35 -F+cLT+jTL+t"${sec}" -N
      #---
      # if [ ${tw} -eq ${num_tw} ]; then gmt colorbar -DjBR+w${mapheight}/0.1+jBL+o0.15/0.0 -Ba${half}+l"Potency-rate density (m/s)"; fi
      if [ ${tw} -eq ${num_tw} ]; then gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba1+l"Cluster Number"; fi
    done
  gmt subplot end
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
echo -e " done !"
rm -rf gmt.history gmt.conf
