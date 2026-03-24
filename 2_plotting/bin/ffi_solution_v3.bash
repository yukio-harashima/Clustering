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
# #----  Rotation 
# if [ ${meca}  -eq  0 ]; then
#   strike=$(awk 'NR==4 {print $1}' fort.40 )
#   dip=$(awk 'NR==4 {print $2}' fort.40 )
#   str=$(echo "scale=1; 90-$strike" |bc)
#   echo  ${str}   -${dip}  > meca_rotation.info
#   #  see "RotFiveTensors" in Git-hub (SELT)
# else
#   rm -rf meca_rotation.info
# fi
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
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
#---






#--- GMT6
gmt begin ffi_solution_v2 png
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
  awk '{if (NR==1) print 0.5,0.5,1,$0,20}' tpdt.dat | gmt meca -R0/1/0/1 -JX4 -Sm3.5 -M -T0 -G195/51/44  -Y-5 -N

#------ Moment rate function
  #----  / Moment rate function variable /
  info=($(gmt gmtinfo mrf.dat -C))
  xmin=${info[0]}; xmax=${info[1]}
  ymin=${info[2]}; ymax=$(echo "${info[3]}*1.05" | bc)

  #----  drow
  gmt basemap -JX5/4 -R${xmin}/${xmax}/${ymin}/${ymax} -Bx+l"Time (s)" -By+l"Moment rate (10@+18@+ Nm/s)" -BWS -X6
  awk '{print $1, $2}' mrf.dat | gmt plot -G195/51/44 -L -N
  #gmt basemap -Bws

#------- Potency density tensor distribution / map meca/
  #----  / map meca variable /
  info=($(gmt gmtinfo pddis.dat -C))
  #+++++  2025/3/30 Yuji Yagi  ++++++++
  DX=$(echo "${info[3]}-(${info[2]})+0.01"|bc)
  DY=$(echo "${info[1]}-(${info[0]})+0.01"|bc)
  echo ${DY} ${DX} 
  result=`echo "${DY}/${DX} > 1.0 " | bc`
  if [ $result -eq 1 ]; then
    wk_x=$(echo "(${DY}-${DX})*0.5+${DY}*0.025"|bc)
    wk_y=$(echo "${DY}*0.05"|bc)
  else
    wk_x=$(echo "${DX}*0.05"|bc)
    wk_y=$(echo "(${DX}-${DY})*0.5+${DX}*0.025"|bc)
  fi
  xmin=$(echo "scale=2;((${info[2]}-${wk_x})/1.)" | bc); xmax=$(echo "scale=2;((${info[3]}+${wk_x})/1.)" | bc)
  ymin=$(echo "scale=2;((${info[0]}-${wk_y})/1.)" | bc); ymax=$(echo "scale=2;((${info[1]}+${wk_y})/1.)" | bc)
  #++++++++++++++++++++++++++++++++++++
  r=${xmin}/${xmax}/${ymin}/${ymax}
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

  #----  drow
  gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${r} -JM${mapwidth} -BWSne -Bx -By  -X6
  gmt coast -Df  -Glightgray
  #--
  #--  Remove "#" if you want to write the distribution.
  awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r} 
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${interval}
  gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  path_bin="$(dirname "$0")"
  path_data=${path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  gmt psxy -W1p $path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  gmt psxy -W1p $path_data/gem_active_faults.gmt.txt -t60
  awk '{print $7,$8}' ../aftershock.dat | gmt plot -Sc0.1c -W1p
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W0.5p
  gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency Density (m)"


  #------- Potency density tensor distribution / plane map meca/
  #----  Rotation 
  strike=$(awk 'NR==4 {print $1}' fort.40 )
  dip=$(awk 'NR==4 {print $2}' fort.40 )
  str=$(echo "scale=1; 90-$strike" |bc)
  echo  ${str}   -${dip}  > meca_rotation.info
  #  see "RotFiveTensors" in Git-hub (SELT)
  fgenmrf
  fgentpdt
  fgenpdtdis

  #----  / plane map meca variable /
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

  #----  drow
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

  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf  meca_rotation.info
