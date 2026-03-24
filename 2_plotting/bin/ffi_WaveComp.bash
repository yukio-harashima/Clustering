#!/bin/bash
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
gmt set PROJ_LENGTH_UNIT cm
gmt set MAP_FRAME_PEN 0.6p
gmt set FONT_ANNOT_PRIMARY 9p
gmt set FONT_LABEL 9p
gmt set MAP_TICK_LENGTH 1p
#---
st_list=($(grep BHZ plot_cwave_f.list | awk '{print $3}'))
amp_list=($(grep @~m@~m/s plot_cwave_f.list | awk '{print $3}'))
azi_list=($(grep Az plot_cwave_f.list | sed -e 's/=/ /g' | awk '{print $4}'))
del_list=($(grep Del plot_cwave_f.list | sed -e 's/=/ /g' | awk '{print $4}'))
#---
tmax=0.00; ampmax=0.00; num_station=0
for station in "${st_list[@]}"; do
  info=($(gmt gmtinfo ./Wave/${station}"BHZ_obs" -C))
  tmax=$(echo ${tmax} ${info[1]} | awk '{if ($1>=$2) print $1; else print $2}')
  ampmax=$(echo ${ampmax} ${info[3]} | awk '{if ($1>=$2) print $1; else print $2}')
  num_station=$(echo "$num_station+1" | bc)
done
echo " # of station: "${num_station}" / Max time: "${tmax}" (s)"
#---
tmax=$(echo "scale=2;${tmax}+0.99" | bc | awk '{printf("%d\n", $1)}')
ampmax=$(echo "scale=2;${ampmax}+0.99" | bc | awk '{printf("%d\n", $1)}')
ampmin=$(echo "scale=2;${ampmax}*(-1.00)" | bc | awk '{printf("%d\n", $1)}')
#r=0/${tmax}/${ampmin}/${ampmax}
r=0/${tmax}/-1.2/1.2
mapwidth=4; mapheight=2
column=4; row=$(echo "scale=2;${num_station}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')
half=$(echo "scale=4;${mapwidth}/2.2" | bc)
quarter=$(echo "scale=4;${mapwidth}/4" | bc)
#---
for i in $(seq 1 ${num_station}); do echo -n "+"; done; echo -e " " ${num_station}
i=0
gmt begin ffi_WaveComp png
  gmt subplot begin ${row}x${column} -Fs${mapwidth}/${mapheight} -M${quarter}/0.0 -B -R${r} -JX${mapwidth}/${mapheight}
    for station in "${st_list[@]}"; do
      echo -n "+"
      gmt subplot set
      #---
      gmt plot ./Wave/${station}"BHZ_obs"  -W0.8p
      gmt plot ./Wave/${station}"BHZ_syn" -W0.8p,red -N
      #---
      gmt text -Dj-${half}/0.0 -F+cLT+jTL+t"${st_list[$i]} BHZ" -N
      gmt text -Dj-${half}/0.30c -F+cLT+jTL+t"${amp_list[$i]} @~m@~m/s" -N
      gmt text -Dj-${half}/0.60c -F+cLT+jTL+t"Az. = ${azi_list[$i]}" -N
      gmt text -Dj-${half}/0.90c -F+cLT+jTL+t"Del.= ${del_list[$i]}" -N
      #---
      let i++
      #---
      if [ ${i} -eq ${#st_list[@]} ]; then gmt basemap -BS -Bx+l"Time (s)" ;fi
    done
  gmt subplot end
  # echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/20/0/0.5  -JX10/1.  -Y-0.2
gmt end show
echo -e " done !"
if [ -e ./const/wave.obs_f ]; then  rm -rf wave.obs_f ;fi
if [ -e ./const/wave.syn   ]; then  rm -rf wave.syn   ;fi
rm -rf gmt.history gmt.conf
