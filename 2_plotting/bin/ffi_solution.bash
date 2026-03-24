#!/bin/bash
#-----
mecascale=0.8      # To change the size of the mechanism solution, change the value of mecascale.
#-----
if [ $# -ne 1 ]; then
  echo -e \\n "ffi_solution.bash [Rotate-Focal-Mechanism : 0, No Rotation : 1]" 
  echo -e     "  default run [ Rotate-Focal-Mechanism : 0]" \\n
  meca=0
else
  meca=${1}
fi
#----
source FileChecker.bash
Fort40Checker
RigidChecker
#----  Rotation 
if [ ${meca}  -eq  0 ]; then
  strike=$(awk 'NR==4 {print $1}' fort.40 )
  dip=$(awk 'NR==4 {print $2}' fort.40 )
  str=$(echo "scale=1; 90-$strike" |bc)
  echo  ${str}   -${dip}  > meca_rotation.info
  #  see "RotFiveTensors" in Git-hub (SELT)
else
  rm -rf meca_rotation.info
fi
#----
fgenmrf
fgentpdt
fgenpdtdis
#---
gmt set PROJ_LENGTH_UNIT cm
gmt set MAP_TICK_LENGTH 1p
gmt set FONT_LABEL 10p
gmt set FONT_ANNOT_PRIMARY 10p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
#---
#---
gmt begin ffi_solution png
#------ infomation
  Smoment=$(awk 'NR==2 {print $1}' fort.40 )
  Mw=$(awk 'NR==2 {print $2}' fort.40 )
  alat=$(awk 'NR==2 {print $4}' fort.40 )
  alon=$(awk 'NR==2 {print $5}' fort.40 )
  depth=$(awk 'NR==2 {print $6}' fort.40 )
  variance=$(awk 'NR==8 {print $1}' fort.40 )
  echo 0 2  Seismic Moment = ${Smoment} Nm, Mw = ${Mw} | gmt text -F+f10p,0,black+jLB  -N -R0/10/0/4  -JX10/1.8  -Y25
  echo 0 1  '(Strike, Dip)=('${strike}, ${dip} ')', Hypo.='('${alat}, ${alon}, ${depth} km')' | gmt text -F+f10p,0,black+jLB  -N
  echo 0 0  Variance = ${variance} | gmt text -F+f10p,0,black+jLB  -N

#------ Moment tensor solution
  awk '{if (NR==1) print 0.5,0.5,1,$0,20}' tpdt.dat | gmt meca -R0/1/0/1 -JX4 -Sm4.5 -M -T0 -G195/51/44  -Y-5 -N

#------ Moment rate function
  info=($(gmt gmtinfo mrf.dat -C))
  xmin=${info[0]}; xmax=${info[1]}
  ymin=${info[2]}; ymax=$(echo "${info[3]}*1.05" | bc)
  gmt basemap -JX5/4 -R${xmin}/${xmax}/${ymin}/${ymax} -Bx+l"Time (s)" -By+l"Moment rate (10@+18@+ Nm/s)" -BWS -X6 -Y0.3
  # gmt basemap -JX5/4 -R${xmin}/${xmax}/${ymin}/${ymax} -Bx+l"Time (s)" -By+l"Moment rate (10@+18@+ Nm/s)" -BWS -Y1
  awk '{print $1, $2}' mrf.dat | gmt plot -G195/51/44 -L -N
  #gmt basemap -Bws

#------- Potency density tensor plane map
  xx=$(awk '{if(NR==6) print $1}' fort.40); yy=$(awk '{if(NR==6) print $2}' fort.40)
  info=($(gmt gmtinfo pdtdis.dat -C))
  xmin=$(echo "${info[6]}-${xx}" | bc); xmax=$(echo "${info[7]}+${xx}" | bc)
  ymin=$(echo "${info[4]}-${yy}" | bc); ymax=$(echo "${info[5]}+${yy}" | bc)
  r=${xmin}/${xmax}/${ymin}/${ymax}
  info=($(gmt gmtinfo pddis.dat -C))
  slipmax=${info[9]}
  wk=$(echo " scale=0;${slipmax}*10/1" | bc)
  len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
  half=$(echo "scale=${len};${slipmax}/2" | bc)
  ratio=$(echo " scale=4;(${ymax}-(${ymin}))/(${xmax}-(${xmin}))" | bc)
  yshift=$(echo "-1.5-${ratio}*12" | bc)
  mapwidth=12; mapheight=$(echo "${mapwidth}*${ratio}" | bc)
  mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.55*'${mecascale}'}' )
  interval=$(echo "scale=${len};$slipmax/5" | bc)
  echo max PD: ${slipmax} m,  interval: ${interval}  m
  slipmax=$(echo "${slipmax}*1.1" |bc )
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JX${mapwidth}/${mapheight} -BWSne -Bx+l"Strike (km)" -By+l"Dip (km)"  -Y${yshift}  -X-6
  #---
  awk '{print $3,$4,$5}' pddis.dat | gmt xyz2grd -Gtmp.grd -I0.5 -R${r}
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${interval}
  #---
  sort -k 8 -n pdtdis.dat | awk '{print $4,$3,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -C -T0 -W0.6
  echo 0 0 | gmt plot -Sa0.2c
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency Density (m)"

#------- Potency density tensor distribution


  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf  meca_rotation.info
