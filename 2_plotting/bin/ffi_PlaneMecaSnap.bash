#!/bin/bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=0.8      # To change the size of the mechanism solution, change the value of mecascale.
#---
source FileChecker.bash
Fort40Checker
RigidChecker
#-----
if [ $# -ne 4 ]; then
  echo -e \\n " ffi_PlaneMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1] [Rotate-Focal-Mechanism : 0, No-Rotation : 1]" 
  fgenmrf
  info=($(gmt gmtinfo mrf.dat -C))
  tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
  interval=$(echo "scale=0;${tmax}/10." | bc)
  if [ $interval -le 0 ]; then
      interval=1
  fi
  i1=${interval}
  i2=${tmax}
  i3=0
  i4=0
  echo -e  " Default fun: ffi_PlaneMecaSnap.bash ${i1} ${i2} ${i3} ${i4} " \\n
else
  i1=${1}; i2=${2}; i3=${3}; i4=${4}
fi
#----  Rotation 
if [ ${i4}  -eq  0 ]; then
  strike=$(awk 'NR==4 {print $1}' fort.40 )
  dip=$(awk 'NR==4 {print $2}' fort.40 )
  str=$(echo "scale=1; 90-$strike" |bc)
  echo  ${str}   -${dip}  > meca_rotation.info
  #  see "RotFiveTensors" in Git-hub (SELT)
else
  rm -rf meca_rotation.info
fi
#-----
echo ${i1} ${i2} ${i3} | fgensnap
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 11p
gmt set FONT_ANNOT_PRIMARY 11p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
#---
xx=$(awk '{if(NR==6) print $1}' fort.40); yy=$(awk '{if(NR==6) print $2}' fort.40)
info=($(gmt gmtinfo snap.dat -C))
xmin=$(echo "${info[6]}-${xx}" | bc); xmax=$(echo "${info[7]}+${xx}" | bc)
ymin=$(echo "${info[8]}-${yy}" | bc); ymax=$(echo "${info[9]}+${yy}" | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
num_tw=${info[5]}
info=($(gmt gmtinfo snap2.dat -C))
slipmax=${info[11]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
half=$(echo "scale=${len};${slipmax}/2" | bc)

# interval=$(echo "scale=2;$slipmax/5" | bc)
interval=0.1

ratio=$(echo " scale=4;(${ymax}-(${ymin}))/(${xmax}-(${xmin}))" | bc)
mapwidth=8; mapheight=$(echo "${mapwidth}*${ratio}" | bc)
wk=$(echo " scale=0;10*(${ymax}-(${ymin}))/(${xmax}-(${xmin}))" | bc)
column=5; if [ ${wk} -le  3 ]; then column=1 ;fi
row=$(echo "scale=2;${num_tw}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')
mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.5*'${mecascale}'}' )
echo max PRD: ${slipmax} m/s,  interval: ${interval}  m/s
slipmax=$(echo "${slipmax}*1.1" |bc )
#---
for tw in $(seq 1 ${num_tw}); do echo -n "+"; done; echo -e " " ${num_tw}
gmt begin ffi_PlaneMecaSnap pdf
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt subplot begin ${row}x${column} -Fs${mapwidth}/${mapheight} -M0.10/0.2 -SCb -SRl -Bx+l"Strike (km)" -By+l"Dip (km)" -Bwsne -R${r} -JX${mapwidth}/${mapheight}
    for tw in $(seq 1 ${num_tw}); do
      echo -n "+"
      gmt subplot set
      #---
      #---
      awk '{if ($3=='${tw}') print $4,$5,$6}' snap2.dat | gmt xyz2grd -Gtmp.grd -I0.5 -R${r}
      gmt grdimage tmp.grd -C  -t0
      #gmt grdimage tmp.grd -C  -t35
      gmt grdcontour tmp.grd -C${interval}
      awk '{if ($3=='${tw}' && $6 > 0.00001) print $4,$5,$6,$7,$8,$9,$10,$11,$12,23}' snap.dat |sort -k 3 -n | gmt meca -Sm${mecasize}c+l -T0 -C -N
      #---
      echo 0 0 | gmt plot -Sa0.2c
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
