#!/bin/bash
#-----
mecascale=1.2     # To change the size of the mechanism solution, change the value of mecascale.
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
#----
fgenmrf
fgentpdt
fgenpdtdis
# i1=1; i2=280; i3=0
# echo $i1 $i2 $i3 | fgensnap
#---Waveform
if [ ! -e ./Wave/ ] ; then
  for dir in wave.obs_f wave.syn .station.abic ; do
    if [ ! -e ${dir} ]; then
      if [ -e ./const/${dir} ]; then
        cp -r ./const/${dir} .
      else
        echo -e \\n " Stop. "${dir}" is needed." \\n
        exit 1
      fi
    fi
  done
  GetASCIIwave.bash
fi
#---
if [ ! -e plot_cwave_f.list ] && [ -e ./const/plot_cwave_f.list ]; then
  cp ./const/plot_cwave_f.list  .
fi
if [ ! -e plot_cwave_f.list ] && [ -e ./const/plot_cwave_f.csh ]; then
  cp ./const/plot_cwave_f.csh  plot_cwave_f.list
fi
#---
if [ ! -e wave.obs_f ]; then
  if [ -e ./const/wave.obs_f ]; then
    cp -r ./const/wave.obs_f .
  else
    echo -e \\n " Stop. wave.obs_f is needed." \\n
    exit 1
  fi
fi
#---
gmt set PROJ_LENGTH_UNIT cm
gmt set MAP_TICK_LENGTH 1p
gmt set FONT_LABEL 20p
gmt set FONT_ANNOT_PRIMARY 20p
gmt set MAP_FRAME_TYPE plain
# gmt set PLOT_DEGREE_FORMAT F
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
gmt set IO_SEGMENT_MARKER "***"
# gmt set MAP_ANNOT_ORTHO 45
# gmt set GMT_THEME cookbook
# #---
# #---
# info=($(gmt gmtinfo pddis.dat -C))
# xmin=$(echo "${info[2]}+0.1" | bc); xmax=$(echo "${info[3]}-0.5" | bc)
# ymin=$(echo "${info[0]}+0.1" | bc); ymax=$(echo "${info[1]}" | bc)
# r=${xmin}/${xmax}/${ymin}/${ymax}
# xleft=$(echo "${xmin}" | bc)
# yleft=$(echo "${ymin}+0.2" | bc)
# xright=$(echo "${xmax}+1.0" | bc)
# yright=$(echo "${ymax}-0.4" | bc)
# r_oblique=${xleft}/${yleft}/${xright}/${yright}r
# xmin=$(echo "scale=2;${xmin}/1." | bc); xmax=$(echo "scale=2;${xmax}/1." | bc)
# ymin=$(echo "scale=2;${ymin}/1." | bc); ymax=$(echo "scale=2;${ymax}/1." | bc)
# r_neighbor=${xmin}/${xmax}/${ymin}/${ymax}
# slipmax=${info[9]}
# wk=$(echo " scale=0;${slipmax}*10/1" | bc)
# len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
# half=$(echo "scale=2;${slipmax}/2" | bc)
# mapwidth=7; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
# mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.6*'${mecascale}'}') 
# interval=$(echo "scale=${len};$slipmax/5" | bc)
# dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
# # echo max PD: ${slipmax} m,  interval: ${interval}  m
# # slipmax=$(echo "${slipmax}*1.1" |bc )
# #---
# #---
# gmt begin FIG2_TOTAL_v3 pdf
# #------ Map meca distribution
#   gmt makecpt -Cbilbao -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
#   gmt basemap -R${r_oblique} -JOc-25/-58/15/75/${mapwidth}c -Bwsne -Bx2 -By1
#   # gmt basemap -R${r_oblique} -JM${mapwidth} -BWSne -Bxa1 -Bya1 
#   gmt coast -R${r_oblique} -Df -Glightgray
#   #--
#   awk '{if($5 > 0) print $2,$1,$5}' pddis.dat | gmt nearneighbor -Gtmp.grd -I0.01 -S0.05 -R${r_neighbor} 
#   gmt grdimage tmp.grd -C  -t35
#   gmt coast -Df -W0.25 -Ggray
#   #--
#   gmt plot -W0.5 plateboundary_Bird_2003.dat
#   gmt plot -W1.0 -Sf0.75/0.05+l+t sandwich_subduction_boundary.dat
#   gmt basemap -R${r_oblique} -JOc-25/-58/15/75/${mapwidth}c -BWSne -Bxa2g2 -Bya1g1
#   sort -k 8 -n pdtdis.dat | awk '{print $6,$5,$8,$9,$10,$11,$12,$13,$14,23}' | gmt meca -Sm${mecasize}c+l -T0 -C  # -t20
#   awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.4c -W1p
#   cb_height=$(echo "${mapheight}*0.15" | bc)
#   gmt colorbar -DjBR+w${cb_height}/0.2+jBL+o-2.15/1.35 -Ba${half}+l"Potency Density (m)" --FONT_LABEL=42p --FONT_ANNOT_PRIMARY=42p
#   # gmt basemap -R${r_oblique} -JOc-25/-58/15/75/${mapwidth}c -BWSne -Bxa2 -Bya1
#   gmt text -Dj0.0/0.2 -F+f20p+cTL+jBL+t"(c)" -N
# #------ Total Moment Tensor
#   mapwidth_meca=2.5
#   # yoff_meca=$(echo "${mapheight}-${mapwidth_meca}" | bc)
#   yoff_meca=13
#   # awk '{if (NR==1) print 0.5,0.5,1,$0,20}' tpdt.dat | gmt meca -R0/1/0/1 -JX${mapwidth_meca} -Sm${mapwidth_meca} -M -T0 -Ggray -N -X4.2 -Y${yoff_meca}
#   # gmt text -Dj0.0/0.2 -F+f20p+cTL+jBL+t"(d)" -N
# #------ Stations
#   st_list=($(grep BHZ plot_cwave_f.list | awk '{print $3}'))
#   amp_list=($(grep @~m@~m/s plot_cwave_f.list | awk '{print $3}'))
#   azi_list=($(grep Az plot_cwave_f.list | sed -e 's/=/ /g' | awk '{print $4}'))
#   del_list=($(grep Del plot_cwave_f.list | sed -e 's/=/ /g' | awk '{print $4}'))
#   ExStationList=("PAB" "MSEY" "SNZO" "OTAV")
#   mapwidth_azi=7
#   #mapwidth_wave=5; mapheight_wave=1.8
#   #xoff_wave=$(echo "${mapwidth_mrf}-${mapwidth_wave}" | bc) #mapwidth_mrf=10
#   #yoff_wave=$(echo "${mapheight}-${mapheight_wave}-(${mapheight_mrf}+1.0)" | bc)
#   LAT=$(awk '{{if(NR==2) print $4}}' fort.40)
#   LON=$(awk '{{if(NR==2) print $5}}' fort.40)
#   #yoff_azi=$(echo "${mapheight_mrf}+1.0-${yoff_meca}" | bc)
#   gmt coast -Rg -JE$LON/$LAT/110/$mapwidth_azi -B110 -Ggray -Y17.5
#   gmt grdmath -Rg -I1 $LON $LAT SDIST = tmp.grd
#   gmt grdcontour tmp.grd -C3335.85 -L3000/4000 -W0.3,black,-
#   gmt grdcontour tmp.grd -C11111.11 -L10000/12000 -W0.3,black,-
#   for station in "${st_list[@]}"; do
#     awk '{if (NR==1) print $5,$4}' "wave.obs_f/"$station"BHZ" | \
#     gmt psxy -St0.5c -W0.3,white -Gblack
#   done
#   # for station in "${ExStationList[@]}"; do
#   #   awk '{if (NR==1) print $5,$4}' "wave.obs_f/"$station"BHZ" | \
#   #   gmt psxy -St0.4c -W0.6,black -Gwhite
#   #   awk '{if (NR==1) print $5,$4,"'"$station"'"}' "wave.obs_f/"$station"BHZ" | \
#   #   gmt text -N -F+jTC+f10p,1,black=~1.2,white -D0/-0.2
#   # done
#   echo $LON $LAT | gmt plot -Sa0.4c -Gblack -W0.3,white
# #   gmt text -F+f << END
# #    $(echo $(echo "scale=4;0.5*$mapwidth_azi" | bc) $(echo "scale=4;0.37*$mapwidth_azi" | bc) | gmt mapproject -Rg -JE$LON/$LAT/110/$mapwidth_azi -I ) 10p,0,black=~1.0,white 30@.
# #    $(echo $(echo "scale=4;0.5*$mapwidth_azi" | bc) $(echo "scale=4;0.05*$mapwidth_azi" | bc) | gmt mapproject -Rg -JE$LON/$LAT/110/$mapwidth_azi -I ) 10p,0,black=~1.0,white 100@.
# # END
#   gmt text -Dj0.0/0.2 -F+f20p+cTL+jBL+t"(a)" -N
# #------ Moment rate function
#   info=($(gmt gmtinfo mrf.dat -C))
#   xmin=${info[0]}; xmax=280
#   ymin=${info[2]}; ymax=$(echo "${info[3]}*1.1" | bc)
#   mapheight_mrf=4.5; mapwidth_mrf=$(echo "${mapwidth_meca}+${mapwidth_azi}+8" | bc)
#   gmt basemap -JX${mapwidth_mrf}/7 -R${xmin}/${xmax}/${ymin}/${ymax} -Bxa50f10 -Bya40+l"Moment rate (10@+18@+ Nm/s)" -BwsnE -X8
# #---EP1
#   gmt plot -G247/247/247 -L << END
# 0 0
# 0 ${ymax}
# 30 ${ymax}
# 30 0
# END
# #---EP2
#   gmt plot -G211/211/211 -L << END
# 30 0
# 30 ${ymax}
# 100 ${ymax}
# 100 0
# END
# #---EP3
#   gmt plot -G166/166/166 -L << END
# 100 0
# 100 ${ymax}
# 145 ${ymax}
# 145 0
# END
# #---EP4
#   gmt plot -G122/122/122 -L << END
# 145 0
# 145 ${ymax}
# 280 ${ymax}
# 280 0
# END
# #---Text
#   gmt text -N -F+jTL+f20p,black=~1.0,white -D0/-0.1 << END
# 0,${ymax},Ep1
# 30,${ymax},Ep2
# 100,${ymax},Ep3
# 145,${ymax},Ep4
# END
#   awk '{print $1, $2}' mrf.dat | gmt plot -W2,black -L
#   gmt text -Dj0.0/0.2 -F+f20p+cTL+jBL+t"(b)" -N
#------
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
# aspect=2.00 # larger value means higher
strinfo=($(gmt gmtinfo vr_str.dat -C))
dipinfo=($(gmt gmtinfo vr_dip.dat -C))
strmin=${strinfo[0]}; strmax=${strinfo[1]}
tmax=${strinfo[3]}
strr=0/${tmax}/${strmin}/${strmax}
# ratio=$(echo " scale=4;${tmax}/(${strmax}-(${strmin}))" | bc)
# mapheight=$(echo "${mapwidth}*${ratio}*${aspect}" | bc)
echo ${strinfo[5]} > tmp; echo ${dipinfo[5]} >> tmp; info=($(gmt gmtinfo tmp -C))
slipmax=${info[1]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
interval=$(echo "scale=2;$slipmax/5" | bc)
cont_min=$(echo "scale=2;(${slipmax}/5)+0.0001" | bc)
half=$(echo "scale=2;${slipmax}/2" | bc)
yoffset=$(echo "5.8-${yoff_meca}" | bc)
#---
  gmt makecpt -Cbilbao -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt basemap -R${strr} -JX${mapwidth_mrf}/-16.3 -BwSnE -Bxa50f10+l"Time (s)" -Bya100+l"Strike (km)" 
  awk '{print $1,$2,$3}' vr_str.dat | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${strr}
  # gmt xyz2grd -Gtmp.grd -I0.1 -R${strr} vr_str_pro.dat
  gmt grdimage tmp.grd -C
  gmt grdcontour tmp.grd -L${cont_min}/${slipmax} -C${interval}
  echo ${cont_min} ${interval}
  # awk '{print $1,$2}' tw_mec_xy.dat | gmt plot -W2,darkblue
  cb_height=2.8 #7*0.2
  gmt colorbar -DjTR+w${cb_height}/0.2+jTR+o2.5/3.0 -Ba${half}+l"Potency-rate Density (m/s)" --FONT_LABEL=42p --FONT_ANNOT_PRIMARY=42p
  gmt plot -W1.5,3_3:0 << END
30,${strmin}
30,80
***
100,${strmin}
100,220
***
145,${strmin}
145,320
END
#---USGS-E2
  gmt plot -Sa1.0c -W1.0,black << END
145,90
END
# echo 145 90 "USGS @%3%M@%%@-w@- 8.1" |  gmt text -N -F+jTR+f20p,1,black=~1.0,white -D-0.3/0
#---USGS-E3
  gmt plot -Sa0.6c -W1,0,black << END
247,285
END
# echo 247 285 "USGS mb 6.7" |  gmt text -N -F+jTR+f20p,1,black=~1.0,white -D-0.3/0
# echo 247 283 "USGS mb 6.7" | gmt text -N -F+jTR+f8p,black=~1.0,white -D-0.3/0
#---velocity
  gmt plot -W1.5 << END
20,150
50,210
***
20,150
50,240
***
20,150
50,270
END
#---Text
  gmt text -N -F+jML+f20p,black=~1.0,white -D0.1/0 << END
50,215,2.0 km/s
50,245,3.0 km/s
50,275,4.0 km/s
END
#---
  gmt text -N -F+jTL+f20p,black=~1.0,white -D0/-0.1 << END
0,${strmin},Ep1
30,${strmin},Ep2
100,${strmin},Ep3
145,${strmin},Ep4
END
  gmt basemap -R${strr} -JX${mapwidth_mrf}/-16.3 -BwSnE -Bxa50f10+l"Time (s)" -Bya100+l"Strike (km)"
  gmt text -Dj0.0/0.2 -F+f20p+cTL+jBL+t"(d)" -N
  gmt text -Dj0.1/-0.2 -F+f20p+cTR+jBL+t"NNE" -N
  gmt text -Dj0.1/-0.2 -F+f20p+cBR+jTL+t"SSW" -N
#-----
gmt end show
rm -rf gmt.conf gmt.history tmp.grd