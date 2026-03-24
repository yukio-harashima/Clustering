#!/bin/bash
#---
# ファイル生成や共通チェックは元のスクリプトを踏襲
source FileChecker.bash
Fort40Checker
RigidChecker
#---
# データの生成処理
if [ $# -ne 1 ] || [ $1 -le 0 ]; then
    echo -e \\n " RupStrDip.bash  [duration (duration > 0.)]" 
    echo -e     "   default run" \\n
    vr=$(awk 'NR==2 {print $7}' fort.40 )
    duration=$(awk 'NR==6 {print ($1*$3)/'${vr}'+$7*$8}' fort.40 )
    echo ${duration} | fgenstrdiprup_mf
else
    echo $1 | fgenstrdiprup_mf
fi

# GMTの共通設定
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 11p
gmt set FONT_ANNOT_PRIMARY 11p
gmt set MAP_FRAME_PEN 0.5p
gmt set MAP_TICK_PEN 0.5p
gmt set PAGE_COLOR 254/254/254

#---
# 1. ファイル検索とソート処理
# vr_str_*.dat を探し、1列目(Strike)の最小値に基づいてソートする
# 同時に、全ファイルを通した3列目(値)の最大値(global_max)も探索する

# 一時ファイル作成
tmp_list="tmp_file_list.txt"
rm -f "$tmp_list"

# ファイルリスト取得
files=(vr_str_*.dat)
if [ ! -e "${files[0]}" ]; then
    if [ -e "vr_str.dat" ]; then
        files=("vr_str.dat")
    else
        echo "Error: No vr_str_*.dat files found."
        exit 1
    fi
fi

echo "Sorting files by Strike minimum value and finding global max value..."

global_slipmax=0

# 各ファイルの情報を取得してリスト化
for f in "${files[@]}"; do
    info=($(gmt gmtinfo "$f" -C))
    min_strike=${info[0]}
    local_max=${info[5]}
    
    if (( $(echo "$local_max > $global_slipmax" | bc -l) )); then
        global_slipmax=$local_max
    fi

    echo "$min_strike $f" >> "$tmp_list"
done

# 全体の最大値に少し余裕を持たせる
global_slipmax=$(echo " $global_slipmax * 1.01" | bc)
# 等高線間隔
global_interval=$(echo "scale=2;$global_slipmax/10" | bc)

echo "Global Max Value: $global_slipmax (Interval: $global_interval)"

# ソート実行
sorted_files=($(sort -n -k1 "$tmp_list" | awk '{print $2}'))
total_files=${#sorted_files[@]}

rm -f "$tmp_list"

echo "Processing order: ${sorted_files[*]}"

#---
# 2. 描画開始
out_name="ffi_Strike_Stack"

# 設定
margin=0.2
is_first=1
count=0

# 高さ追跡用変数
total_shift=0
first_height=0

gmt begin "$out_name" png
    
    # ソートされた順序でループ
    for file_str in "${sorted_files[@]}"; do
        count=$((count + 1))
        
        # ID抽出
        basename=$(basename "$file_str" .dat)
        id=${basename#vr_str}
        id=${id#_}
        if [ -z "$id" ]; then disp_id="Single"; else disp_id="$id"; fi

        echo "Drawing $file_str (ID: $disp_id)..."

        # --- サイズ計算 ---
        aspect=0.2
        strinfo=($(gmt gmtinfo "$file_str" -C))
        strmin=${strinfo[0]}; strmax=${strinfo[1]}
        tmax=${strinfo[3]}
        strr=0/${tmax}/${strmin}/${strmax}

        str_diff=$(echo "${strmax} - ${strmin}" | bc)
        if [ "$(echo "$str_diff == 0" | bc)" -eq 1 ]; then str_diff=1; fi
        ratio=$(echo " scale=4;${str_diff}/(${tmax})" | bc)
        
        strwidth=10
        mapheight=$(echo "${strwidth}*${ratio}*${aspect}" | bc)

        # --- 描画位置の調整  ---
        if [ $is_first -eq 1 ]; then
            # 最初の図
            move_cmd=""
            is_first=0
            total_shift=0
            # 最初の図の高さを保存（カラーバー計算用）
            first_height=${mapheight}
        else
            # 2つ目以降: 「今回の図の高さ + マージン」分だけ下に移動 (積み上げロジック維持)
            # ※ mapheight は現在の図の高さ
            shift_val=$(echo "${margin} + ${mapheight}" | bc)
            move_cmd="-Y-${shift_val}"
            
            # 総移動量を更新
            total_shift=$(echo "${total_shift} + ${shift_val}" | bc)
        fi
        
        # --- 軸ラベル制御 ---
        extra_opts=()
        if [ $count -eq $total_files ]; then
            axes="-BWSne"
            extra_opts=(-Bx+l"Time (s)")
        else
            axes="-BWsne"
        fi

        # --- 描画実行 ---
        # カラーパレット作成 (共通スケール)
        gmt makecpt -Chot -T0/${global_slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
        
        # Basemap
        gmt basemap -R${strr} -JX${strwidth}/-${mapheight} ${axes} "${extra_opts[@]}" -By+l"Strike (km)" ${move_cmd}
        
        # Plot
        awk '{print $1,$2,$3}' "$file_str" | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${strr}
        gmt grdimage tmp.grd -C

        # Strike = -25km の点線を描画
        # X軸: 0 to tmax, Y軸: -25 固定
        echo -e "0 -25\n${tmax} -25" | gmt plot -W0.5p,black,dashed -t40

        gmt grdcontour tmp.grd -C${global_interval}
        
        
        
        gmt basemap -Bnswe
        
        # ID text
        echo "0 1 Fault: ${disp_id}" | gmt text -F+f10p,Helvetica-Bold,black+jTL -R0/1/0/1 -JX${strwidth}/${mapheight} -Dj0.1/0.1

        # prev_height変数は不使用

    done

    # --- 共通カラーバーの描画 ---
    
    # 1. 全体の高さ = 累積移動量 + 最初の図の高さ
    # (積み上げロジックに基づき、最初の図のTopから最後の図のBottomまでの距離を算出)
    total_height=$(echo "${total_shift} + ${first_height}" | bc)
    
    # 2. カラーバーの中心位置オフセット
    # 基準点(最後の図の右中央)から、全体の中央までの距離
    # オフセット = (最初の図の高さ + 累積移動量 - 最後の図の高さ) / 2
    # ※ mapheight はループ終了時なので「最後の図の高さ」が入っています
    cbar_offset_y=$(echo "scale=4; (${first_height} + ${total_shift} - ${mapheight}) / 2" | bc)
    
    # カラーバーを描画
    # -DjMR: 最後の図の右中央基準
    # +w: 長さ=total_height
    # +o: Xオフセットを2.0cに設定 (重なり防止のため正の値へ変更)
    # +m: ラベルと数値を右側に移動
    # +e: 両端に三角形
    
    gmt colorbar -DjMR+w${total_height}/0.3c+o-0.5c/${cbar_offset_y}c+m \
        -Ba${global_interval}+l"Potency-rate Density (m/s)"

gmt end show

# 後始末
rm -rf gmt.history gmt.conf tmp tmp.grd tmp_file_list.txt
echo "Done. Created $out_name.png"