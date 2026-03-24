#!/bin/bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
# 地図上に一つの選択した指定しクラスターの時間進展を表示
# awkコマンドに直接指定
mecascale=2.0    # To change the size of the mechanism solution, change the value of mecascale.

clno=1 # Cluster number, seems fixed
#---
# Source external checker scripts if they exist
if [ -f FileChecker.bash ]; then source FileChecker.bash; fi
if [ -f Fort40Checker.bash ]; then source Fort40Checker.bash; fi # Assuming these are separate scripts or functions
if [ -f RigidChecker.bash ]; then source RigidChecker.bash; fi   # Or part of FileChecker.bash

if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi

# Default values
DEFAULT_SNAP_INTERVAL=""
DEFAULT_TOTAL_DURATION=""
DEFAULT_MODE=0
DEFAULT_ENDTIME="-1" # Using -1 to signify 'not set by user'

# Argument parsing
if [ $# -eq 0 ]; then # No arguments, try to use defaults entirely
    echo "情報: 引数が指定されませんでした。デフォルト値を使用して実行を試みます。"
    if command -v fgenmrf &> /dev/null && [ -f mrf.dat ]; then
        fgenmrf # Assuming this updates/creates mrf.dat
        # Get the second field (max time) from mrf.dat info. Original used info[1] from multi-field output.
        info_mrf_time=$(gmt gmtinfo mrf.dat -C -o1)
        if [ -n "${info_mrf_time}" ]; then
            DEFAULT_TOTAL_DURATION=$(echo "scale=0; ${info_mrf_time} * 1.05 / 1" | bc)
            DEFAULT_SNAP_INTERVAL=$(echo "scale=0; ${DEFAULT_TOTAL_DURATION} / 10" | bc)
            if [ $(echo "${DEFAULT_SNAP_INTERVAL} <= 0" | bc -l) -eq 1 ]; then
                DEFAULT_SNAP_INTERVAL=1
            fi
            echo "情報: デフォルト値を使用します: スナップ間隔=${DEFAULT_SNAP_INTERVAL}, 総期間=${DEFAULT_TOTAL_DURATION}, モード=${DEFAULT_MODE}"
            i1=${DEFAULT_SNAP_INTERVAL}; i2=${DEFAULT_TOTAL_DURATION}; i3=${DEFAULT_MODE}; i4=${DEFAULT_ENDTIME}
            echo ${i1} ${i2} ${i3} | fgensnap > .tmpout_fgensnap ; rm -rf .tmpout_fgensnap # Call fgensnap with derived defaults
            if command -v swap_mmcl5.py &> /dev/null; then swap_mmcl5.py; else echo "警告: swap_mmcl5.py が見つかりません。"; fi
        else
            echo "エラー: mrf.dat からデフォルトの総期間を導出できませんでした。終了します。"
            exit 1
        fi
    else
        echo "エラー: fgenmrf コマンドが見つからないか、mrf.dat が存在しません。デフォルト値を導出できません。終了します。"
        echo "使用法: $0 [スナップ間隔] [総期間] [Average: 0, Moment: 1] [オプション: 終了時刻]"
        exit 1
    fi
elif [ $# -eq 3 ]; then # Standard 3 arguments
    i1=${1}; i2=${2}; i3=${3}; i4=${DEFAULT_ENDTIME}
    echo "情報: 指定された引数を使用します: スナップ間隔=${i1}, 総期間=${i2}, モード=${i3}"
    echo ${i1} ${i2} ${i3} | fgensnap
    if command -v swap_mmcl5.py &> /dev/null; then swap_mmcl5.py; else echo "警告: swap_mmcl5.py が見つかりません。"; fi
elif [ $# -eq 4 ]; then # 4 arguments including endtime
    i1=${1}; i2=${2}; i3=${3}; i4=${4}
    echo "情報: 指定された引数を使用します: スナップ間隔=${i1}, 総期間=${i2}, モード=${i3}, 終了時刻=${i4}"
    # fgensnap is still called with total_duration (i2)
    echo ${i1} ${i2} ${i3} | fgensnap
    if command -v swap_mmcl5.py &> /dev/null; then swap_mmcl5.py; else echo "警告: swap_mmcl5.py が見つかりません。"; fi
else # Incorrect number of arguments (1, 2, or >4)
    echo -e \\n "使用法: $0 [スナップ間隔] [総期間] [Average: 0, Moment: 1] [オプション: 終了時刻]"
    echo "エラー: $# 個の引数が指定されました。0個（デフォルト実行）、3個、または4個の引数を指定してください。"
    # Try to provide default suggestion like original script
    if command -v fgenmrf &> /dev/null && [ -f mrf.dat ]; then
        _info_mrf_sug_time=$(gmt gmtinfo mrf.dat -C -o1)
        if [ -n "${_info_mrf_sug_time}" ]; then
            _sug_tmax=$(echo "scale=0; ${_info_mrf_sug_time} * 1.05 / 1" | bc)
            _sug_interval=$(echo "scale=0; ${_sug_tmax} / 10" | bc)
            if [ $(echo "${_sug_interval} <= 0" | bc -l) -eq 1 ]; then _sug_interval=1; fi
            echo -e  " デフォルト実行の提案: $0 ${_sug_interval} ${_sug_tmax} 0"
        fi
    fi
    exit 1
fi

#--- GMT settings
gmt set MAP_TICK_LENGTH 1p
gmt set FONT_LABEL 20p
gmt set FONT_ANNOT_PRIMARY 20p
gmt set MAP_FRAME_PEN 0.4p
gmt set MAP_TICK_PEN 0.4p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 1p
gmt set PLOT_DEGREE_FORMAT D

#--- Info from snap2.dat
if [ ! -f snap2.dat ]; then
    echo "エラー: snap2.dat が見つかりません。fgensnapが正しく実行されたか確認してください。"
    exit 1
fi
info=($(gmt gmtinfo snap2.dat -C))
if [ ${#info[@]} -lt 16 ]; then # Check if info array has enough elements
    echo "エラー: snap2.dat から情報を読み取れませんでした。ファイルが空か、形式が不正です。"
    echo "snap2.dat の内容:"
    cat snap2.dat
    echo "gmt gmtinfo snap2.dat -C の出力:"
    gmt gmtinfo snap2.dat -C
    exit 1
fi

#+++++  2025/3/30 Yuji Yagi  ++++++++
DX=$(echo "(${info[15]}-(${info[14]}))+0.01"|bc)
DY=$(echo "(${info[13]}-(${info[12]}))+0.01"|bc)
echo "DX: ${DX}, DY: ${DY}"
result=$(echo "${DY}/${DX} > 1.0" | bc) # Corrected command substitution
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
r=${xmin}/${xmax}/${ymin}/${ymax}
original_num_tw=${info[5]} # Total number of time windows available in snap2.dat
slipmax=${info[11]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc)
len=1;  if [ ${wk} -le 1  ]; then len=3; fi;  if [ ${wk} -ge 100 ]; then len=0; fi
half=$(echo "scale=${len};${slipmax}/2" | bc)
interval_slip=$(echo "scale=2;$slipmax/5" | bc) # Renamed from 'interval' to avoid clash with snap_interval (i1)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
column=5 # Fixed number of columns for subplot

# mecasize calculation (ensure slipmax is not zero)
if [ $(echo "${slipmax} == 0" | bc -l) -eq 1 ]; then
    echo "警告: slipmax が0です。mecasize をデフォルト値に設定します。"
    mecasize=0.5 # Default mecasize or handle error
else
    mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.4*'${mecascale}'}' )
fi

dip=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 )
echo "最大 PRD: ${slipmax} m/s,  スリップ間隔: ${interval_slip} m/s"
# slipmax_plot=$(echo "${slipmax}*1.1" |bc ) # This variable was defined but not used later, original slipmax is used by mecasize

#--- Paths
path_bin="$(dirname "$0")"
path_data=${path_bin/bin/data} # Assumes a specific directory structure

#--- Clustering info
if [ ! -f clusteringinfo.dat ]; then
    echo "警告: clusteringinfo.dat が見つかりません。labmaxのデフォルト値を使用します。"
    labmax=1.5 # Default labmax if file is missing (e.g. for CPT 0.5/1.5/1)
else
    info2=($(gmt gmtinfo clusteringinfo.dat -C))
    if [ ${#info2[@]} -lt 2 ]; then
        echo "警告: clusteringinfo.dat から情報を読み取れませんでした。labmaxのデフォルト値を使用します。"
        labmax=1.5
    else
        labmax=$(echo "${info2[1]}+0.5" | bc)
    fi
fi
echo "Labmax (for CPT): ${labmax}"


#--- Determine the number of time windows to plot based on endtime (i4)
num_tw_to_plot=${original_num_tw}
echo "情報: snap2.dat に基づく元の時間窓の総数: ${original_num_tw}"

if [ $(echo "${i4} > -1" | bc -l) -eq 1 ]; then # If endtime (i4) is specified (not the default -1)
    echo "情報: 指定された終了時刻 (i4): ${i4}"
    if [ $(echo "${i1} > 0" | bc -l) -eq 1 ]; then # And snap_interval (i1) is positive
        max_tw_for_endtime=$(echo "scale=0; (${i4}) / (${i1})" | bc) # Integer division for floor
        if [ -z "$max_tw_for_endtime" ]; then max_tw_for_endtime=0; fi # Safety for empty bc output

        # Ensure max_tw_for_endtime is not negative if i4 was negative (though endtime should be positive)
        if [ $(echo "${max_tw_for_endtime} < 0" | bc -l) -eq 1 ]; then
             max_tw_for_endtime=0
        fi
        echo "情報: 終了時刻に基づいて計算された最大時間窓数: ${max_tw_for_endtime}"

        if [ ${max_tw_for_endtime} -lt ${num_tw_to_plot} ]; then
            num_tw_to_plot=${max_tw_for_endtime}
        fi
        # If endtime is less than one snap_interval, max_tw_for_endtime will be 0.
        # num_tw_to_plot will become 0.
        if [ $(echo "${i4} > 0 && ${i4} < ${i1}" | bc -l) -eq 1 ]; then # If endtime is positive but less than one interval
             num_tw_to_plot=0
        elif [ $(echo "${i4} <= 0" | bc -l) -eq 1 ]; then # If endtime itself is zero or negative
             num_tw_to_plot=0
        fi
    else
        echo "警告: スナップ間隔 (i1=${i1}) が正ではありません。終了時刻は使用されません。全てのスナップショットをプロットします。"
    fi
else
    echo "情報: 終了時刻 (i4) が指定されていません。利用可能な全ての ${original_num_tw} 個のスナップショットをプロットします。"
fi

echo "情報: プロットする時間窓の数: ${num_tw_to_plot}"

# Adjust 'row' calculation for subplot based on num_tw_to_plot
if [ ${num_tw_to_plot} -gt 0 ]; then
    row=$(echo "scale=2;${num_tw_to_plot}/${column}+0.99" | bc | awk '{printf("%d\n", $1)}')
else
    row=1 # Default row to 1 if no plots, though the main plot loop won't run.
fi
echo "情報: サブプロットの行数: ${row}, 列数: ${column}"

#--- Progress bar display
if [ ${num_tw_to_plot} -gt 0 ]; then
    for (( tw_idx=1; tw_idx<=${num_tw_to_plot}; tw_idx++ )); do echo -n "+"; done; echo -e " " ${num_tw_to_plot} "スナップショット"
fi

#--- Begin GMT plotting
output_filename="ffi_MapMecaClusteringSnap6_${i1}_cl${clno}"
if [ $(echo "${i4} > -1" | bc -l) -eq 1 ]; then
    output_filename+="_end${i4}" # Add endtime to filename if specified
fi
output_filename+=".png"

# Handle case where no time windows are to be plotted
if [ ${num_tw_to_plot} -le 0 ]; then
    echo "情報: プロットする時間窓がありません（終了時刻 ${i4} が短すぎるか、0のため）。"
    gmt begin ${output_filename}
      gmt basemap -R0/10/0/1 -JX10c/1c -B+t"データなし: 終了時刻 ${i4} (s)"
    gmt end show
    echo -e "完了！ (プロットなし)"
    rm -rf gmt.history gmt.conf .gmtcommands
    exit 0
fi

gmt begin ${output_filename}
  # CPT for cluster number. Assuming cluster numbers are integers.
  gmt makecpt -Cvik -T0.5/${labmax}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
  gmt subplot begin ${row}x${column} -Fs${mapwidth}c/${mapheight}c -M0.1/1.0 -SCb -SRl -Bwsne -R${r} -JM${mapwidth}c -Ba1 # Adjusted -Fs and -JM to use 'c' for cm
    for tw in $(seq 1 ${num_tw_to_plot}); do
      echo -n "+"
      gmt subplot set # Set current subplot panel
      #--- Plotting meca symbols
      # Ensure clusteringSnap2_v2.dat exists
      if [ ! -f clusteringSnap2_v2.dat ]; then
          echo "error: subplot ${tw}: clusteringSnap2_v2.dat not found"
          gmt basemap -B+t"error: no data" # Indicate error on subplot
          continue # Skip to next subplot
      fi
      awk -F'\t' -v tw_val="${tw}" -v clnum="${clno}" '$3==tw_val && $30==clnum {print $20,$19,$30,$7,$8,$9,$10,$11,$12,23}' clusteringSnap2_v2.dat | sort -k 3 -n | gmt meca -Sm${mecasize}c+l -T0 -C
      #--- Coastlines and fault lines
      gmt coast -Df -W0.25p -Ggray
      if [ -f faultline.dat ]; then
          awk '{print $2,$1}' faultline.dat | gmt plot -W0.5p,black
          if [ ${dip} -gt 0 ]; then
              sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1
          fi
      else
          echo "警告: faultline.dat が見つかりません。"
      fi
      # Plate boundaries (optional, ensure files exist)
      # gmt set IO_SEGMENT_MARKER "***"
      # if [ -f "${path_data}/plateboundary_Bird_2003.dat" ]; then
      #   gmt psxy -W1p "${path_data}/plateboundary_Bird_2003.dat" -t20
      # fi
      # gmt set IO_SEGMENT_MARKER ">"
      # if [ -f "${path_data}/gem_active_faults.gmt.txt" ]; then
      #   gmt psxy -W1p "${path_data}/gem_active_faults.gmt.txt" -t60
      # fi
      if [ -f fort.40 ]; then
        awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.2c
      else
        echo "警告: fort.40 が見つかりません。"
      fi
      #--- Time label for the subplot
      t1=$(echo "scale=1;(${tw}-1)*${i1}" | bc)
      t2=$(echo "scale=1;${tw}*${i1}" | bc)
      if [ ${i3} -eq 1 ]; then sec=${t2}" s"; else sec=${t1}"-"${t2}" s"; fi
      gmt text -Dj0.0/0.1 -F+f20p+cLT+jBL+t"${sec}" -N # Reduced font size for subplot
      #--- Colorbar on the last plotted snapshot
      if [ ${tw} -eq ${num_tw_to_plot} ]; then
          # Using mapheight for width of vertical colorbar, or define a specific width
          # Colorbar for cluster number, annotation step should be 1
          gmt colorbar -DjBR+w${mapheight}/0.2+jBL+o0.05/0.0 -Ba1+l"Cluster Number"  # Use current CPT
      fi
    done
  gmt subplot end
  # Timestamp or other info at the bottom
  current_dir=$(pwd)
  echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
echo -e "\n完了！ 生成されたファイル: ${output_filename}"
rm -rf gmt.history gmt.conf .gmtcommands # Clean up GMT temp files
