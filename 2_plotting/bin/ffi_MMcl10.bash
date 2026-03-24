#!/bin/bash
#---
#---   The size of meca sphere was changed to be proportional to the potency [Yagi 2022.04.01]
#---
# 地図上にすべてのクラスターの時間進展を表示 メカニズム解でなく点で描画
# awkコマンドのフィルタリング条件を修正
mecascale=0.05   # To change the size of the mechanism solution, change the value of mecascale.

clno=4 # Cluster number,
#---
# Source external checker scripts if they exist
if [ -f FileChecker.bash ]; then source FileChecker.bash; fi
if [ -f Fort40Checker.bash ]; then source Fort40Checker.bash; fi
if [ -f RigidChecker.bash ]; then source RigidChecker.bash; fi
if [ -e meca_rotation.info ]; then rm -rf  meca_rotation.info  ; fi

# Default values
DEFAULT_SNAP_INTERVAL=""
DEFAULT_TOTAL_DURATION=""
DEFAULT_MODE=0
DEFAULT_STARTTIME=0
DEFAULT_ENDTIME="-1"

# Argument parsing
if [ $# -eq 0 ]; then
    echo "情報: 引数が指定されませんでした。デフォルト値を使用して実行を試みます。"
    if command -v fgenmrf &> /dev/null && [ -f mrf.dat ]; then
        fgenmrf
        info_mrf_time=$(gmt gmtinfo mrf.dat -C -o1)
        if [ -n "${info_mrf_time}" ]; then
            DEFAULT_TOTAL_DURATION=$(echo "scale=0; ${info_mrf_time} * 1.05 / 1" | bc)
            DEFAULT_SNAP_INTERVAL=$(echo "scale=0; ${DEFAULT_TOTAL_DURATION} / 10" | bc)
            if [ $(echo "${DEFAULT_SNAP_INTERVAL} <= 0" | bc -l) -eq 1 ]; then DEFAULT_SNAP_INTERVAL=1; fi
            DEFAULT_ENDTIME=${DEFAULT_TOTAL_DURATION}
            echo "情報: デフォルト値を使用します: スナップ間隔=${DEFAULT_SNAP_INTERVAL}, 総期間=${DEFAULT_TOTAL_DURATION}, モード=${DEFAULT_MODE}, 開始時刻=${DEFAULT_STARTTIME}, 終了時刻=${DEFAULT_ENDTIME}"
            i1=${DEFAULT_SNAP_INTERVAL}; i2=${DEFAULT_TOTAL_DURATION}; i3=${DEFAULT_MODE}; i4=${DEFAULT_STARTTIME}; i5=${DEFAULT_ENDTIME}
            # fgensnap is still called with original total_duration (i2) for snap2.dat generation
            echo ${i1} ${i2} ${i3} | fgensnap > .tmpout_fgensnap ; rm -rf .tmpout_fgensnap
        else
            echo "エラー: mrf.dat からデフォルトの総期間を導出できませんでした。終了します。"
            exit 1
        fi
    else
        echo "エラー: fgenmrf コマンドが見つからないか、mrf.dat が存在しません。デフォルト値を導出できません。終了します。"
        echo "使用法: $0 [スナップ間隔] [総期間] [Average: 0, Moment: 1] [オプション: 開始時刻] [オプション: 終了時刻]"
        exit 1
    fi
elif [ $# -eq 3 ]; then
    i1=${1}; i2=${2}; i3=${3}
    i4=${DEFAULT_STARTTIME}
    i5=${i2} # Default end time is total duration
    echo "情報: 指定された引数を使用します: スナップ間隔=${i1}, 総期間=${i2}, モード=${i3}, 開始時刻=${i4}, 終了時刻=${i5}"
    echo ${i1} ${i2} ${i3} | fgensnap
elif [ $# -eq 5 ]; then
    i1=${1}; i2=${2}; i3=${3}; i4=${4}; i5=${5}
    echo "情報: 指定された引数を使用します: スナップ間隔=${i1}, 総期間=${i2}, モード=${i3}, 開始時刻=${i4}, 終了時刻=${i5}"
    # fgensnap is still called with total_duration (i2) to ensure snap2.dat has data up to that point
    echo ${i1} ${i2} ${i3} | fgensnap
else
    echo -e \\n "使用法: $0 [スナップ間隔] [総期間] [Average: 0, Moment: 1] [オプション: 開始時刻] [オプション: 終了時刻]"
    echo "エラー: $# 個の引数が指定されました。0個（デフォルト実行）、3個、または5個の引数を指定してください。"
    if command -v fgenmrf &> /dev/null && [ -f mrf.dat ]; then
        _info_mrf_sug_time=$(gmt gmtinfo mrf.dat -C -o1)
        if [ -n "${_info_mrf_sug_time}" ]; then
            _sug_tmax=$(echo "scale=0; ${_info_mrf_sug_time} * 1.05 / 1" | bc)
            _sug_interval=$(echo "scale=0; ${_sug_tmax} / 10" | bc)
            if [ $(echo "${_sug_interval} <= 0" | bc -l) -eq 1 ]; then _sug_interval=1; fi
            echo -e  " デフォルト実行の提案: $0 ${_sug_interval} ${_sug_tmax} 0"
            echo -e  " 例: $0 ${_sug_interval} ${_sug_tmax} 0 10 50 (開始時刻 10s, 終了時刻 50s でプロット)"
        fi
    fi
    exit 1
fi

# Validate arguments
if [ $(echo "${i4} < 0" | bc -l) -eq 1 ] || [ $(echo "${i5} < 0" | bc -l) -eq 1 ]; then
    echo "エラー: 開始時刻 (${i4}) および終了時刻 (${i5}) は非負である必要があります。"
    exit 1
fi
if [ $(echo "${i4} >= ${i5}" | bc -l) -eq 1 ] && [ $(echo "${i5} == 0 && ${i4} == 0" | bc -l) -ne 1 ]; then
    echo "エラー: 開始時刻 (${i4}) は終了時刻 (${i5}) より小さくする必要があります（両方が0の場合を除く）。"
    exit 1
fi
if [ $(echo "${i1} <= 0" | bc -l) -eq 1 ]; then
     echo "エラー: スナップ間隔 (${i1}) は正である必要があります。"
     exit 1
fi

# --- Generate new time window definition file for swap_mmcl5_revised.py ---
new_tw_file="new_time_windows.dat"
echo "情報: 新しい時間窓定義ファイル '${new_tw_file}' を生成します。"
rm -f ${new_tw_file}
current_new_tw_idx=1
current_plot_time=${i4} # Start time for the first window

# Handle i4=0, i5=0 case: no windows
if [ $(echo "${i4} == 0 && ${i5} == 0" | bc -l) -eq 1 ]; then
    echo "情報: 開始時刻と終了時刻が両方0のため、時間窓は生成されません。"
    num_new_tw=0
else
    while (( $(echo "${current_plot_time} < ${i5}" | bc -l) )); do
        window_start_time=${current_plot_time}
        window_end_time=$(echo "${current_plot_time} + ${i1}" | bc)
        
        if (( $(echo "${window_end_time} > ${i5}" | bc -l) )); then
            window_end_time=${i5}
        fi

        # Ensure start_time is less than end_time before writing
        if (( $(echo "${window_start_time} < ${window_end_time}" | bc -l) )); then
            echo "${current_new_tw_idx} ${window_start_time} ${window_end_time}" >> ${new_tw_file}
            current_new_tw_idx=$((current_new_tw_idx + 1))
        else
             # This case (start >= end) should ideally not happen if i5 > i4 and i1 > 0
             # but as a safeguard, break if the window is zero or negative duration.
             break
        fi
        current_plot_time=$(echo "${current_plot_time} + ${i1}" | bc)
    done
    num_new_tw=$((current_new_tw_idx - 1))
fi

echo "情報: ${new_tw_file} に ${num_new_tw} 個の時間窓を定義しました。"

if [ ${num_new_tw} -le 0 ] && ! [ -f ${new_tw_file} ]; then
    echo "警告: 生成された新しい時間窓がありません。処理を継続しますが、プロットは空になる可能性があります。"
    touch ${new_tw_file} # Create an empty file to avoid swap_mmcl5_revised.py error
fi

# --- Execute swap_mmcl5_revised.py (assuming it's in the PATH or same directory) ---
export CLNO=${clno} 
echo "情報: swap_mmcl5_revised.py を実行します..."
if command -v swap_mmcl5_revised.py &> /dev/null; then
    swap_mmcl5_revised.py
    if [ $? -ne 0 ]; then
        echo "エラー: swap_mmcl5_revised.py の実行に失敗しました。ログを確認してください。"
        exit 1
    fi
else
    echo "エラー: swap_mmcl5_revised.py が見つかりません。"
    exit 1
fi
echo "情報: swap_mmcl5_revised.py の実行が完了しました。"


#--- GMT settings
gmt set MAP_TICK_LENGTH 1p
gmt set FONT_LABEL 20p
gmt set FONT_ANNOT_PRIMARY 20p
gmt set MAP_FRAME_PEN 0.4p
gmt set MAP_TICK_PEN 0.4p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 1p
gmt set PLOT_DEGREE_FORMAT D


# --- 変数設定 ---
# 断層データが格納されているフォルダ
FAULT_DIR="/Users/harashima-yukio/Desktop/JpGU2025/prz/fig/hk/output_fault_data"
# 震源情報 (スクリプトの引数から取得)
EPI_LAT=$1
EPI_LON=$2
# 出力ファイル名
OUT_FILE="HKfault_map"

# --- 入力データチェック ---
# FAULT_DIRが存在し、F*.datファイルがあるか確認します。
if [ ! -d "$FAULT_DIR" ] || [ -z "$(ls -A $FAULT_DIR/F*.dat 2>/dev/null)" ]; then
    echo "エラー: フォルダ '$FAULT_DIR' が見つからないか、'F*.dat' ファイルが含まれていません。"
    exit 1
fi

# --- データ範囲の決定 ---
# フォルダ内の全断層データの情報（緯度経度、深度）を一度に取得します。
# これにより、全ての断層を包含する地図範囲やカラーパレットのレンジを決定できます。
# -Cオプションにより、w e s n zmin zmax の順で値が返されます。
all_info=($(gmt gmtinfo $FAULT_DIR/F*.dat -C))

#--- Info from snap2.dat (still needed for map region, etc.)
if [ ! -f snap2.dat ]; then
    echo "エラー: snap2.dat が見つかりません。fgensnapが正しく実行されたか確認してください。"
    exit 1
fi
info=($(gmt gmtinfo snap2.dat -C))
if [ ${#info[@]} -lt 16 ]; then
    echo "エラー: snap2.dat から情報を読み取れませんでした。ファイルが空か、形式が不正です。"
    exit 1
fi

#+++++  2025/3/30 Yuji Yagi  ++++++++
DX=$(echo "${info[15]}-(${info[14]})+0.01"|bc | awk '{printf "%.5f", $0}')
DY=$(echo "${info[13]}-(${info[12]})+0.01"|bc | awk '{printf "%.5f", $0}')
echo ${DY} ${DX} 
result=$(echo "${DY} > ${DX}" | bc -l)
if [ $result -eq 1 ]; then
  wk_x=$(echo "(${DY}-${DX})*0.5+${DY}*0.025"|bc | awk '{printf "%.5f", $0}')
  wk_y=$(echo "${DY}*0.05"|bc | awk '{printf "%.5f", $0}')
else
  wk_x=$(echo "${DX}*0.05"|bc | awk '{printf "%.5f", $0}')
  wk_y=$(echo "(${DX}-${DY})*0.5+${DX}*0.025"|bc | awk '{printf "%.5f", $0}')
fi
echo ${wk_y}
echo ${info[12]}
xmin=$(echo "scale=2;((${info[14]} - ${wk_x})/1.)" | bc | awk '{printf "%.5f", $0}') 
xmax=$(echo "scale=2;((${info[15]} + ${wk_x})/1.)" | bc | awk '{printf "%.5f", $0}')
ymin=$(echo "scale=2;((${info[12]} - ${wk_y})/1.)" | bc | awk '{printf "%.5f", $0}')
ymax=$(echo "scale=2;((${info[13]} + ${wk_y})/1.)" | bc | awk '{printf "%.5f", $0}')
#++++++++++++++++++++++++++++++++++++
echo ${ymin}
r=${xmin}/${xmax}/${ymin}/${ymax}
# r=137.62000/140.24000/40.94000/43.69000
echo ${r}
slipmax=${info[11]}
wk=$(echo " scale=0;${slipmax}*10/1" | bc | awk '{printf "%.5f", $0}')
len=1;  if [ $(echo "${wk} <= 1" | bc -l) -eq 1 ]; then len=3; fi;  if [ $(echo "${wk} >= 100" | bc -l) -eq 1 ]; then len=0; fi
# half=$(echo "scale=${len};${slipmax}/2" | bc) # Unused
interval_slip=$(echo "scale=2;$slipmax/5" | bc)
mapwidth=6; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth}c -R${r} | awk '{print $2}')
column=5

if [ $(echo "${slipmax} == 0" | bc -l) -eq 1 ]; then
    echo "警告: slipmax が0です。mecasize をデフォルト値に設定します。"
    mecasize=0.5
else
    mecasize=$(echo ${slipmax} | awk '{print 5/$1*0.4*'${mecascale}'}' )
fi

dip_val_fort40=$(awk 'NR==4 {printf("%d\n",$2)}' fort.40 ) # Renamed to avoid conflict
echo "最大 PRD: ${slipmax} m/s,  スリップ間隔: ${interval_slip} m/s"

#--- Clustering info
if [ ! -f clusteringinfo.dat ]; then
    echo "警告: clusteringinfo.dat が見つかりません。labmaxのデフォルト値を使用します。"
    labmax=1.5
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

#--- Determine the number of new time windows to plot from the generated file ---
if [ -f ${new_tw_file} ]; then
    num_tw_to_plot=$(wc -l < ${new_tw_file} | awk '{print $1}') # Get actual line count
else
    num_tw_to_plot=0
fi
echo "情報: プロットする新しい時間窓の数 (${new_tw_file} より): ${num_tw_to_plot}"

if [ ${num_tw_to_plot} -gt 0 ]; then
    row=$(echo "scale=0;(${num_tw_to_plot} + ${column} - 1) / ${column}" | bc)
else
    row=1
fi
echo "情報: サブプロットの行数: ${row}, 列数: ${column}"

#--- Progress bar display
if [ ${num_tw_to_plot} -gt 0 ]; then
    echo -n "プロット進行状況: ["
    for (( p_idx=1; p_idx<=${num_tw_to_plot}; p_idx++ )); do echo -n "-"; done; echo -e "]"
    echo -n "                  " # Initial spacing for '+'
fi

#--- Begin GMT plotting
# ===== 変更点: 出力ファイル名を変更 =====
output_filename="ffi_MapMecaClusteringSnap10_${i1}_cl_all"
if [ $# -eq 5 ]; then
    output_filename+="_start${i4}_end${i5}"
fi
# output_filename+=".png"
processed_plot_data_file="clusteringSnap2_v2.dat" # swap_mmcl5_revised.py output

if [ ! -f ${processed_plot_data_file} ]; then
    echo "エラー: ${processed_plot_data_file} が見つかりません。swap_mmcl5_revised.pyが正しく実行されたか確認してください。"
    # Create a dummy plot indicating missing data
    gmt begin ${output_filename}
      gmt basemap -R0/10/0/1 -JX10c/1c -B+t"エラー: ${processed_plot_data_file} なし"
    gmt end show
    exit 1
fi


if [ ${num_tw_to_plot} -le 0 ]; then
    echo "情報: プロットする時間窓がありません（指定された時間範囲 [${i4}, ${i5}] s に対応する有効なスナップショットがないため）。"
    gmt begin ${output_filename}
      gmt basemap -R0/10/0/1 -JX10c/1c -B+t"データなし: 時間範囲 [${i4}, ${i5}] s"
    gmt end show
    echo -e "完了！ (プロットなし)"
    rm -f gmt.history gmt.conf .gmtcommands ${new_tw_file}
    exit 0
fi


# r=138.1/140.42/41.17/43.81
gmt begin ${output_filename} png
  gmt makecpt -Cvik -T0.5/${labmax}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
#   gmt subplot begin ${row}x${column} -Fs${mapwidth}c/${mapheight}c -M0.1c/1.0c -SCb -SRl -Bwsne -R${r} -JM${mapwidth}c -Ba1f0.5
  gmt subplot begin ${row}x${column} -Fs${mapwidth}+c2 -M0.1/1.0 -SCb -SRl -Bwsne -R${r} -JM${mapwidth} -Ba1
    plot_panel_idx=0
    # Read new_time_windows.dat line by line
    # Format: new_tw_index window_start_time window_end_time
    while IFS=" " read -r current_tw_for_plot start_time_label end_time_label; do
        plot_panel_idx=$((plot_panel_idx + 1))
        echo -n "+" # Progress indicator
        gmt subplot set # Set current subplot panel
        gmt coast -Df -W0.25p -Ggray

        ####
        #  echo "以下の断層データをプロットします:"
        # for fault_file in $FAULT_DIR/F*.dat; do
        #     # echo "  -> ${fault_file}"

        #     # 断層のトレース（矩形全体）を赤線でプロット
        #     # awkで経度($2), 緯度($1)の順に列を入れ替えてGMTに渡します。
        #     awk '{print $2, $1}' "${fault_file}" | gmt plot -W0.7p,red

        #     # 各断層の最浅深度を取得
        #     fault_info=($(gmt gmtinfo "${fault_file}" -C))
        #     top_depth=${fault_info[4]}

        #     # 断層の上端（地表に最も近い部分）を太い黒線でプロット
        #     awk -v depth="$top_depth" '{if($3 == depth) print $2, $1}' "${fault_file}" | gmt plot -W1.5p,black
        # done


        ####
      

        # Filter data for the current_tw_for_plot (which is the new_tw from the file)
        # The 30th column in clusteringSnap2_v2.dat is 'no' (cluster number)
        # The 3rd column is 'tw' (new time window index)
        awk -F'\t' '{OFS="\t"; print $0}' ${processed_plot_data_file} | sort -k 6 -n > tmp.dat
        
        # ===== 変更点: すべてのデータをプロットするようにフィルタ条件を修正 =====
        # 元のコード: awk -F'\t' -v tw_val="${current_tw_for_plot}" -v clnum="${clno}" '$3==tw_val && $30==clnum && $6> 0.00001{...}'
        # $30==clnum の条件を削除
        awk -F'\t' -v tw_val="${current_tw_for_plot}" \
            '$3==tw_val && $30 > 0 {print $20,$19,$30}' tmp.dat \
            | gmt plot -Sc0.12c -W0.025c -C -N
        
        
        if [ -f faultline.dat ]; then
            awk '{print $2,$1}' faultline.dat | gmt plot -W0.5p,black
            if [ ${dip_val_fort40} -gt 0 ]; then # Use renamed variable
                sort -k 3 -n faultline.dat | awk 'NR <= 2 {print $2,$1}' | gmt plot -W1p
            fi
        else
            echo "警告: subplot ${plot_panel_idx}: faultline.dat が見つかりません。"
        fi
        if [ -f fort.40 ]; then
          awk '{if(NR==2) print $5,$4}' fort.40 | gmt plot -Sa0.2c
        else
          echo "警告: subplot ${plot_panel_idx}: fort.40 が見つかりません。"
        fi

        # --- Time label for the subplot ---
        # Use start_time_label and end_time_label from new_time_windows.dat
        # Ensure labels are formatted correctly, bc output might not be ideal for direct use
        start_label_fmt=$(printf "%.0f" "${start_time_label}")
        end_label_fmt=$(printf "%.0f" "${end_time_label}")

        if [ ${i3} -eq 1 ]; then # Moment snapshot, use end time
             sec_label="${end_label_fmt} s"
        else # Average snapshot, use range
             sec_label="${start_label_fmt}-${end_label_fmt} s"
        fi
        gmt text -Dj0.0c/0.1c -F+f20p+cLT+jBL+t"${sec_label}" -N

        # --- Colorbar on the last plotted snapshot ---
        if [ ${plot_panel_idx} -eq ${num_tw_to_plot} ]; then
            gmt colorbar -DjBR+w${mapheight}c/0.2c+jBL+o0.05c/0.0c -Ba1+l"Cluster Number"
        fi
    done < "${new_tw_file}" # Read from the generated file
  gmt subplot end
  current_dir=$(pwd)
#   gmt basemap -R0/10/0/1 -JX20c/1.5c -Y-1.5c -B0 # Dummy basemap for text positioning
  echo 0 0.1 "${current_dir}" | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX20/1.  -Y-1.2
gmt end show
echo -e "\n完了！ 生成されたファイル: ${output_filename}"
rm -f gmt.history gmt.conf tmp.dat .gmtcommands "${new_tw_file}" # Clean up new_time_windows.dat