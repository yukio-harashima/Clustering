#!/bin/bash
#---
#--- 海底地形描画版
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
echo $dur $tmax $typ | fgensnap_mf
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 15p
gmt set FONT_ANNOT_PRIMARY 15p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.8p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
gmt set FORMAT_GEO_MAP D
#---
info=($(gmt gmtinfo snap2.dat -C))

#+++++  2025/3/30 Yuji Yagi  ++++++++
DX=$(echo "${info[15]}-(${info[14]})+0.01"|bc)
DY=$(echo "${info[13]}-(${info[12]})+0.01"|bc)
echo ${DY} ${DX} 
result=`echo "${DY}/${DX} > 1.0 " | bc`
if [ $result -eq 1 ]; then
  wk_x=$(echo "(${DY}-${DX})*0.5+${DY}*0.025"|bc)
  wk_y=$(echo "${DY}*0.05"|bc)
else
  wk_x=$(echo "${DX}*0.05"|bc)
  wk_y=$(echo "(${DX}-${DY})*0.5+${DX}*0.025"|bc)
fi
xmin=$(echo "scale=2;((${info[14]}-${wk_x})/1.)" | bc); xmax=$(echo "scale=2;((${info[15]}+${wk_x})/1.)" | bc)
ymin=$(echo "scale=2;((${info[12]}-${wk_y})/1.)" | bc); ymax=$(echo "scale=2;((${info[13]}+${wk_y})/1.)" | bc)
#++++++++++++++++++++++++++++++++++++
# r=${xmin}/${xmax}/${ymin}/${ymax}
r=138.34/140.19/41.54/43.48

num_tw=${info[5]}
slipmax=${info[11]}
slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
# mid6=$(awk '{a[NR]=$6} END{n=asort(a); if(n%2==1) print a[(n+1)/2]; else print (a[n/2] + a[n/2+1])/2}' snap.dat)
interval=$(echo "scale=2;$slipmax/5" | bc)
half=$(echo "scale=2;$slipmax/2" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
column=10; row=$(echo "scale=2;${num_tw}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')


#---

#--- 1. カラーパレットの準備 (gmt beginの前) ---
# 地形用のCPTファイルを作成
gmt grdcut @earth_relief_15s -R120/150/20/50 -Gdem.nc
gmt makecpt -Cgray -T-10000/0/100 -Z > topo.cpt
# 断層すべり用のCPTファイルを作成
gmt makecpt -Chawaii -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > slip.cpt

for tw in $(seq 1 ${num_tw}); do echo -n "+"; done; echo -e " " ${num_tw}
gmt begin ffi_MapXSnap_2 png
    
    gmt subplot begin ${row}x${column} -Fs${mapwidth} -M0.15/0.3 -SCb -SRl -Bxa0.5 -Bya0.5 -BWSne -R${r} -JM${mapwidth}
    for tw in $(seq 1 ${num_tw}); do
      echo -n "+"
      gmt subplot set
      # 地形用プロット
      gmt grdgradient dem.nc -Ggrad.grd -A45 -Ne0.8
      gmt grdimage dem.nc -Igrad.grd -Ctopo.cpt -t20
      gmt coast -Df -W0.25 -Glightgray
      #---
      if "${p_s}"; then # plot strike
        # awk '{if ($3=='${tw}' &&  $6 > 0.08) print $20,$19,$6,$13,$6*5}' snap.dat |sort -k 3 -n  | gmt plot -SV0c+jc -W1.5+cl -Cslip.cpt
        awk '{if ($3=='${tw}') print $20,$19,$6,$13,$6*1.5}' snap.dat |sort -k 3 -n  | gmt plot -SV0c+jc -W1.2+cl -Cslip.cpt -t30
        
        # awk '{if ($3=='${tw}' &&  $6 > 0.1) print $20,$19,$6,$13,$6*5}' snap.dat |sort -k 3 -n  | gmt plot -SV0c+jc -W1.2+cl -Cslip.cpt -t30

        # awk -v tw_val=${tw} -v mid_val=${mid6} '{if ($3 == tw_val && $6 > mid_val) print $20,$19,$6,$13,0.5}' snap.dat | gmt plot -SV0c+jc -W1.5+cl -Cslip.cpt
        # awk '{if ($3=='${tw}') print $20,$19,$6,$14,$6/'${slipmax}'*1.4}' snap.dat | gmt plot -SV0c+jc -W1.5+cl -C
      fi
      # if "${p_p}"; then # plot P axis
      #   awk '{if ($3=='${tw}') print $20,$19,$6,$22,$6/'${slipmax}'*0.4}' snap.dat | gmt plot -SV0c+jc -W1.5+cl -C
      # fi
      # if "${p_t}"; then # plot T axis
      #   awk '{if ($3=='${tw}') print $20,$19,$6,$23,$6/'${slipmax}'*0.4}' snap.dat | gmt plot -SV0c+jc -W1.5+cl -C
      # fi
      # if "${p_b}"; then # plot B axis
      #   awk '{if ($3=='${tw}') print $20,$19,$6,$24,$6/'${slipmax}'*0.4}' snap.dat | gmt plot -SV0c+jc -W1.5+cl -C
      # fi
      #---
      
      # awk '{print $2,$1}' faultline.dat | gmt plot
      awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.2c
      #---
      t1=$(echo "scale=1;(${tw}-1)*$dur" | bc)
      t2=$(echo "scale=1;${tw}*$dur" | bc)
      if [ $typ -eq 1 ]; then sec=${t2}" s"; else sec=${t1}"-"${t2}" s"; fi
      gmt text -Dj0.0/-0.35 -F+cLT+jTL+t"${sec}" -N
      #---
      if [ ${tw} -eq ${num_tw} ]; then gmt colorbar -Cslip.cpt -DjBR+w${mapheight}/0.1+jBL+o0.15/0.0 -Ba${half}+l"Potency-rate density (m/s)"; fi
    done
  gmt subplot end
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
echo -e " done !"
rm -rf topo.cpt slip.cpt gmt.history gmt.conf
