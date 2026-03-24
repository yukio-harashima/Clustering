#!/bin/bash
source FileChecker.bash
Fort40Checker
#---
if [ ! -e plot_cwave_f.list ]; then
  if [ -e ./const/plot_cwave_f.list ]; then
    cp ./const/plot_cwave_f.list .
  else
    echo -e \\n " Stop. plot_cwave_f.list is needed." \\n
    exit 1
  fi
fi
st_list=($(grep BHZ plot_cwave_f.list | awk '{{print $3}}'))
# station_count=${#st_list[@]} # ステーションの合計数を取得
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
gmt set MAP_FRAME_PEN thinner
gmt set FONT_ANNOT_PRIMARY 11p
gmt set FONT_LABEL 11p
gmt set MAP_TICK_LENGTH 2p
#---
mapwidth=12
#---
gmt begin ffi_AziEquiStaMap_v2 png
  LAT=$(awk '{{if(NR==2) print $4}}' fort.40)
  LON=$(awk '{{if(NR==2) print $5}}' fort.40)
  gmt coast -Rg -JE$LON/$LAT/110/$mapwidth -B110 -Ggray
  gmt grdmath -Rg -I1 $LON $LAT SDIST = tmp.grd
  gmt grdcontour tmp.grd -C3335.85 -L3000/4000 -W0.3,black,-
  gmt grdcontour tmp.grd -C10007.55 -L9000/11000 -W0.3,black,-
  for station in "${st_list[@]}"; do
    awk '{if (NR==1) print $5,$4}' "wave.obs_f/"$station"BHZ" | \
    gmt psxy -St0.8c -W0.3,white -Gblack
    # awk '{if (NR==1) print $5,$4,"'"$station"'"}' "wave.obs_f/"$station"BHZ" | \
    # gmt text -N -F+jML+f0p,0,black=~1,white -D0.1/0
  done
  echo $LON $LAT | gmt plot -Sa0.7c -Gblack -W0.3,white
  gmt text -F+f << END
  $(echo $(echo "scale=4;0.5*$mapwidth" | bc) $(echo "scale=4;0.37*$mapwidth" | bc) | gmt mapproject -Rg -JE$LON/$LAT/110/$mapwidth -I ) 11p,0,black=~0.5,white 30@.
  $(echo $(echo "scale=4;0.5*$mapwidth" | bc) $(echo "scale=4;0.09*$mapwidth" | bc) | gmt mapproject -Rg -JE$LON/$LAT/110/$mapwidth -I ) 11p,0,black=~0.5,white 90@.
END
  # ステーションの合計数を表示
  # gmt text -F+f << END
  #  $(echo $(echo "scale=4;0.9*$mapwidth" | bc) $(echo "scale=4;0.8*$mapwidth" | bc) | gmt mapproject -Rg -JE$LON/$LAT/110/$mapwidth -I ) 11p,0,black=~0.5,white "${station_count}"
END
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.0
gmt end show
rm -rf gmt.history gmt.conf tmp.grd
