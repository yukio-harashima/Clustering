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
  i1=5   #snap interval
  i2=51   #total duration
  i3=0    #Average: 0
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
# 表示する図の列数
colum=1
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
interval=$(echo "scale=2;$slipmax/5" | bc)
ratio=$(echo " scale=4;(${ymax}-(${ymin}))/(${xmax}-(${xmin}))" | bc)
mapwidth=8; mapheight=$(echo "${mapwidth}*${ratio}" | bc)
wk=$(echo " scale=0;10*(${ymax}-(${ymin}))/(${xmax}-(${xmin}))" | bc)
column=5; if [ ${wk} -le  3 ]; then column=1 ;fi
row=$(echo "scale=2;${num_tw}/${colum}+0.99" | bc | awk '{printf("%d\n", $1)}')
mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.5*'${mecascale}'}' )
echo max PRD: ${slipmax} m/s,  interval: ${interval}  m/s
slipmax=$(echo "${slipmax}*1.1" |bc )

#---

info2=($(gmt gmtinfo clusteringinfo.dat -C))
labmax=$(echo "${info2[1]}+0.5" | bc)
echo ${labmax}

# 検索文字列を引数から取得
# search_string1 = $1
# search_string2 = $2
# search_string3 = $3 

search_strings=("1016" "4" "1" "5")
joined_string=$(IFS=','; echo "${search_strings[*]}")

# 検索値をカンマ区切りで定義"969"  "975"
# search_values=$1
# search_values="4,15,18"

# 変数の定義

mapheight1=$(echo "$mapheight * $row" | bc)
line_count=$(wc -l < clusteringSnap2.dat)
# 30列の値が0である行の数を取得し、null_row に格納
null_row=$(awk -F'\t' '$30 == 0 {count++} END {print count}' clusteringSnap2.dat)
# line_count から null_row を引いて data_count を算出
data_count=$((line_count - null_row))
#---
gmt begin ffi_PlaneMecaSnapCluster3_${joined_string} pdf
  gmt makecpt -Cvik -T-0.5/${labmax}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white

  # Begin main subplot with multiple columns to hold separate plots
  # gmt subplot begin 1x${#search_strings[@]}  -Fs${mapwidth}/${mapheight1} -M0.10/0.2 -A -SCb -SRl -R${r} -JX${mapwidth}/${mapheight1}
  
  
    # gmt subplot set
    
  gmt subplot begin ${row}x${#search_strings[@]} -Fs${mapwidth}/${mapheight} -M0.30/0.2 -SCb -SRl -Bx+l"Strike (km)" -By+l"Dip (km)" -Bwsne -R${r} -JX${mapwidth}/${mapheight}
  for sn in $(seq 1 ${#search_strings[@]}); do 
  # for search_string in "${#search_strings[@]}"; do 
    ss=$((sn-1))
    search_string=${search_strings[$ss]}
    search_number=$(echo "$search_string" | bc)
    for tw in $(seq 1 ${num_tw}); do
      echo -n "+"
      # echo $tw
      gmt subplot set $(($tw - 1)),$(($sn - 1))
      #---
      # awk '{if ($3=='${tw}') print $4,$5,$6}' snap2.dat | gmt xyz2grd -Gtmp.grd -I0.5 -R${r}
      # gmt grdimage tmp.grd -C  -t35
      # gmt grdcontour tmp.grd -C${interval}
      # search_string = search_strings[@]
      # awk -v search="$search_string" -F'\t' '{if ($29 ~ search && $3=='${tw}' && $6 > 0.00001) print $4,$5,$30,$7,$8,$9,$10,$11,$12,23}' clusteringSnap2.dat | sort -k 3 -n | gmt meca -Sm${mecasize}c+l -T0 -C -N
      if [ "$search_number" -le "$data_count" ]; then
        # line_count 以下の場合、$30 列を参照
        awk -v search="$search_string" -F'\t' '{if ($30 ~ search && $3=='${tw}' && $6 > 0.00001) print $4,$5,$30,$7,$8,$9,$10,$11,$12,23}' clusteringSnap2.dat | sort -k 3 -n | gmt meca -Sm${mecasize}c+l -T0 -C
      else
        # line_count より大きい場合、$29 列を参照
        awk -v search="$search_string" -F'\t' '{if ($29 ~ search && $3=='${tw}' && $6 > 0.00001) print $4,$5,$30,$7,$8,$9,$10,$11,$12,23}' clusteringSnap2.dat | sort -k 3 -n | gmt meca -Sm${mecasize}c+l -T0 -C
      fi
      #---
      echo 0 0 | gmt plot -Sa0.2c
      #---
      t1=$(echo "scale=1;(${tw}-1)*${i1}" | bc)
      t2=$(echo "scale=1;${tw}*${i1}" | bc)
      if [ ${i3} -eq 1 ]; then sec=${t2}" s"; else sec=${t1}"-"${t2}" s"; fi
      gmt text -Dj0.0/-0.35 -F+cLT+jTL+t"${sec}" -N
      #---
      if [ ${tw} -eq ${num_tw} ]; then gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba1+l"Cluster Number" -I; fi
      #---
      if [ ${tw} -eq  1 ]; then gmt text -Dj3.7/-0.35 -F+cLT+jTL+t"${search_strings[$ss]}" -N; fi




    done
   
    # gmt subplot end
  done
  
  # echo 0 0  `pwd` | gmt text -F+f6,0,black+jLB -R0/10/0/0.5 -JX20/1 -Y-1.2
  gmt subplot end
gmt end show
echo -e " done !"
rm -rf gmt.history gmt.conf
