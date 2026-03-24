#!/bin/bash
# クラスタリング　一番最初の処理
# snap_y.dat(n,vベクトルの格納されたsnap.dat)を生成する
# 実行時は必ず[TR, st_max, 0]で実行する（fort.40 のデータそのままをクラスタリングに使用するため）
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
mecascale=0.8      # To change the size of the mechanism solution, change the value of mecascale.
#---
source FileChecker.bash
Fort40Checker
RigidChecker

#-----PMcl4用のメカニズム解の回転したsnap_yr.datを生成
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
#----  Rotation 
strike=$(awk 'NR==4 {print $1}' fort.40 )
dip=$(awk 'NR==4 {print $2}' fort.40 )
str=$(echo "scale=1; 90-$strike" |bc)
echo  ${str}   -${dip}  > meca_rotation.info
#  see "RotFiveTensors" in Git-hub (SELT)
#-----
if [ $# -ne 3 ]; then
  echo -e \\n " MapMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1]" 
  fgenmrf
  info=($(gmt gmtinfo mrf.dat -C))
  tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
  interval=$(echo "scale=0;${tmax}/10." | bc)
  if [ $interval -le 0 ]; then
      interval=1
  fi
  echo -e  " Default fun: ffi_MapMecaSnap.bash ${interval}  ${tmax}  0 " \\n
  echo ${interval}  ${tmax}  0  | fgensnap_y  > .tmpout  ; rm -rf .tmpout
  i1=${interval}; i2=${tmax}; i3=0
else
  echo ${1} ${2} ${3} | fgensnap_y
  i1=${1}; i2=${2}; i3=${3}
fi

# snap_y.dat を snap_yr.dat にリネーム
if [ -f snap_y.dat ]; then
  mv snap_y.dat snap_yr.dat
  echo "snap_y.dat was renamed to snap_yr.dat"
else
  echo "snap_y.dat not found. Skipping renaming."
fi

#-----MMcl4用のメカニズム解の回転していないsnap_y.datを生成
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi
if [ $# -ne 3 ]; then
  echo -e \\n " MapMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1]" 
  fgenmrf
  info=($(gmt gmtinfo mrf.dat -C))
  tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
  interval=$(echo "scale=0;${tmax}/10." | bc)
  if [ $interval -le 0 ]; then
      interval=1
  fi
  echo -e  " Default fun: ffi_MapMecaSnap.bash ${interval}  ${tmax}  0 " \\n
  echo ${interval}  ${tmax}  0  | fgensnap_y  > .tmpout  ; rm -rf .tmpout
  i1=${interval}; i2=${tmax}; i3=0
else
  echo ${1} ${2} ${3} | fgensnap_y
  i1=${1}; i2=${2}; i3=${3}
fi

# #-----


# if [ $# -ne 4 ]; then
#   echo -e \\n " ffi_PlaneMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1] [Rotate-Focal-Mechanism : 0, No-Rotation : 1]" 
#   fgenmrf
#   info=($(gmt gmtinfo mrf.dat -C))
#   tmax=$(echo "scale=0;${info[1]}*1.05/1." | bc)
#   interval=$(echo "scale=0;${tmax}/10." | bc)
#   if [ $interval -le 0 ]; then
#       interval=1
#   fi
#   i1=${interval}
#   i2=${tmax}
#   i3=0
#   i4=0
#   echo -e  " Default fun: ffi_PlaneMecaSnap.bash ${i1} ${i2} ${i3} ${i4} " \\n
# else
#   i1=${1}; i2=${2}; i3=${3}; i4=${4}
# fi
# #----  Rotation 
# if [ ${i4}  -eq  0 ]; then
#   strike=$(awk 'NR==4 {print $1}' fort.40 )
#   dip=$(awk 'NR==4 {print $2}' fort.40 )
#   str=$(echo "scale=1; 90-$strike" |bc)
#   echo  ${str}   -${dip}  > meca_rotation.info
#   #  see "RotFiveTensors" in Git-hub (SELT)
# else
#   rm -rf meca_rotation.info
# fi
# #-----
# echo ${i1} ${i2} ${i3} | fgensnap_y