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
gmt set PLOT_DEGREE_FORMAT D
gmt set PAGE_COLOR 254/254/254
# gmt set DEGREE FORMAT 1
#---


#----     /     variable     /      ----#
#----  / Moment rate function variable /a
a_info=($(gmt gmtinfo mrf.dat -C))
a_xmin=${a_info[0]}; a_xmax=${a_info[1]}
a_ymin=${a_info[2]}; a_ymax=$(echo "${a_info[3]}*1.05" | bc)
a_r=${a_xmin}/${a_xmax}/${a_ymin}/${a_ymax}


#----  / map meca variable /b
b_info=($(gmt gmtinfo pddis.dat -C))
#+++++  2025/3/30 Yuji Yagi  ++++++++
b_DX=$(echo "${b_info[3]}-(${b_info[2]})+0.01"|bc)
b_DY=$(echo "${b_info[1]}-(${b_info[0]})+0.01"|bc)
echo ${b_DY} ${b_DX} 
b_result=`echo "${b_DY}/${b_DX} > 1.0 " | bc`
if [ $b_result -eq 1 ]; then
  b_wk_x=$(echo "(${b_DY}-${b_DX})*0.5+${b_DY}*0.025"|bc)
  b_wk_y=$(echo "${b_DY}*0.05"|bc)
else
  b_wk_x=$(echo "${b_DX}*0.05"|bc)
  b_wk_y=$(echo "(${b_DX}-${b_DY})*0.5+${b_DX}*0.025"|bc)
fi
b_xmin=$(echo "scale=2;((${b_info[2]}-${b_wk_x})/1.)" | bc); b_xmax=$(echo "scale=2;((${b_info[3]}+${b_wk_x})/1.)" | bc)
b_ymin=$(echo "scale=2;((${b_info[0]}-${b_wk_y})/1.)" | bc); b_ymax=$(echo "scale=2;((${b_info[1]}+${b_wk_y})/1.)" | bc)
#++++++++++++++++++++++++++++++++++++
b_r=${b_xmin}/${b_xmax}/${b_ymin}/${b_ymax}
b_slipmax=${b_info[9]}
b_wk=$(echo " scale=0;${b_slipmax}*10/1" | bc)
b_len=1;  if [ ${b_wk} -le 1  ]; then b_len=3; fi;  if [ ${b_wk} -ge 100 ]; then b_len=0; fi
b_half=$(echo "scale=${b_len};${b_slipmax}/2" | bc)
# b_mapwidth=14; b_mapheight=$(echo ${b_xmax} ${b_ymax} | gmt mapproject -JM${b_mapwidth} -R${b_r} | awk '{print $2}')
b_mecascale=0.8 
b_mecasize=$(echo ${b_slipmax} | awk '{print 5/$1*0.6*'${b_mecascale}'}') 
b_interval=$(echo "scale=${b_len};$b_slipmax/5" | bc)
b_dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
echo max PD: ${b_slipmax} m,  interval: ${b_interval}  m
b_slipmax=$(echo "${b_slipmax}*1.1" |bc )


#----  / plane map meca variable /
c_xx=$(awk '{if(NR==6) print $1}' fort.40); c_yy=$(awk '{if(NR==6) print $2}' fort.40)
c_info=($(gmt gmtinfo pdtdis.dat -C))
c_xmin=$(echo "${c_info[6]}-${c_xx}" | bc); c_xmax=$(echo "${c_info[7]}+${c_xx}" | bc)
c_ymin=$(echo "${c_info[4]}-${c_yy}" | bc); c_ymax=$(echo "${c_info[5]}+${c_yy}" | bc)
c_r=${c_xmin}/${c_xmax}/${c_ymin}/${c_ymax}
c_info=($(gmt gmtinfo pddis.dat -C))
c_slipmax=${c_info[9]}
c_wk=$(echo " scale=0;${c_slipmax}*10/1" | bc)
c_len=1;  if [ ${c_wk} -le 1  ]; then c_len=3; fi;  if [ ${c_wk} -ge 100 ]; then c_len=0; fi
c_half=$(echo "scale=${c_len};${c_slipmax}/2" | bc)
c_ratio=$(echo " scale=4;(${c_ymax}-(${c_ymin}))/(${c_xmax}-(${c_xmin}))" | bc)
c_yshift=$(echo "-1.5-${c_ratio}*12" | bc)
c_mapwidth=12; c_mapheight=$(echo "${c_mapwidth}*${c_ratio}" | bc)
c_mecascale=0.8 
c_mecasize=$(echo ${c_slipmax} | awk '{print 5/$1*0.55*'${c_mecascale}'}' )
c_interval=$(echo "scale=${c_len};$c_slipmax/5" | bc)
echo max PD: ${c_slipmax} m,  interval: ${c_interval}  m
c_slipmax=$(echo "${c_slipmax}*1.1" |bc )

#----  /   shift   /   ----#
b_single_width=10; b_mapheight_unit=$(echo ${b_xmax} ${b_ymax} | gmt mapproject -JM${b_single_width} -R${b_r} | awk '{print $2}')
c_ydelta=$(echo "${c_yshift}" | sed 's/^-//')
b_mapheight=$(echo "${c_ydelta}+5" | bc)
b_mapwidth=$(echo "scale=4;(${b_mapheight})/((${b_mapheight_unit})/10)" | bc)
echo c_yshift ${c_yshift} c_ydelta ${c_ydelta} b_ratio ${b_mapheight_unit}, b_mapheight ${b_mapheight}, b_mapwidth ${b_mapwidth}





#--- GMT6
gmt begin ffi_solution_v2 png
#------ infomation
  Smoment=$(awk 'NR==2 {print $1}' fort.40 )
  Mw=$(awk 'NR==2 {print $2}' fort.40 )
  strike=$(awk 'NR==4 {print $1}' fort.40 )
  dip=$(awk 'NR==4 {print $2}' fort.40 )
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
  # #----  / Moment rate function /
  gmt basemap -JX5/4 -R${a_r} -Bx+l"Time (s)" -By+l"Moment rate (10@+18@+ Nm/s)" -BWS -X6
  awk '{print $1, $2}' mrf.dat | gmt plot -G195/51/44 -L -N
  #gmt basemap -Bws

#------- Potency density tensor distribution / map meca/
  #----  / map meca /
  gmt makecpt -Chot -T0/${b_slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${b_r} -JM${b_mapwidth} -BWSne -Bx -By -X7.75 -Y${c_yshift}
  gmt coast -Df  -Glightgray
  #--
  #--  Remove "#" if you want to write the distribution.
  awk '{print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${b_r} 
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${b_interval}
  gmt coast -Df -W0.25
  #--
  sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${b_mecasize}c+l -T0 -C  # -t20
  awk '{print $2,$1}' faultline.dat | gmt plot -W0.2
  if [ ${b_dip} -gt 0 ]; then sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1 ; fi
  b_path_bin="$(dirname "$0")"
  b_path_data=${b_path_bin/bin/data}
  gmt set IO_SEGMENT_MARKER "***"
  gmt psxy -W1p $b_path_data/plateboundary_Bird_2003.dat -t20
  gmt set IO_SEGMENT_MARKER ">"
  gmt psxy -W1p $b_path_data/gem_active_faults.gmt.txt -t60
  awk '{print $7,$8}' ../aftershock.dat | gmt plot -Sc0.1c -W1p
  awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W0.5p
  gmt colorbar -DjBR+w${b_mapheight}/0.2+jBL+o0.05/0.0 -Ba${b_half}+l"Potency Density (m)"


  #------- Potency density tensor distribution / plane map meca/c
  #----  Rotation 
  c_strike=$(awk 'NR==4 {print $1}' fort.40 )
  c_dip=$(awk 'NR==4 {print $2}' fort.40 )
  c_str=$(echo "scale=1; 90-$c_strike" |bc)
  echo  ${c_str}   -${c_dip}  > meca_rotation.info
  #  see "RotFiveTensors" in Git-hub (SELT)
  fgenmrf
  fgentpdt
  fgenpdtdis

  gmt makecpt -Chot -T0/${c_slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${c_r} -JX${c_mapwidth}/${c_mapheight} -BWSne -Bx+l"Strike (km)" -By+l"Dip (km)" -X-13.75
  #---
  awk '{print $3,$4,$5}' pddis.dat | gmt xyz2grd -Gtmp.grd -I0.5 -R${c_r}
  gmt grdimage tmp.grd -C  -t35
  gmt grdcontour tmp.grd -C${c_interval}
  #---
  sort -k 8 -n pdtdis.dat | awk '{print $4,$3,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${c_mecasize}c+l -C -T0 -W0.6
  echo 0 0 | gmt plot -Sa0.2c
  gmt colorbar -DjBR+w${c_mapheight}/0.2+jBL+o0.05/0.0 -Ba${c_half}+l"Potency Density (m)"

  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
gmt end show
rm -rf gmt.history gmt.conf  meca_rotation.info
