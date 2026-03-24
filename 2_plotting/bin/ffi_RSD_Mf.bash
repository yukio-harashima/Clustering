#!/bin/bash
#---
source FileChecker.bash
Fort40Checker
RigidChecker
#---
# データの生成処理
# 引数がなければ自動計算、あれば指定時間で fgenstrdiprup_mf を実行
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
# ファイル検索とループ処理
# vr_str_*.dat にマッチするファイルをすべて処理する
# マッチするファイルがない場合のエラー回避のため nullglob を設定する手もあるが、
# ここでは [ -e file ] でチェックする。

# 処理対象のファイルリストを取得 (vr_str_*.dat)
# もし vr_str.dat (番号なし) も対象にしたい場合は files=(vr_str*.dat) とする
files=(vr_str_*.dat)

# ファイルが見つからなかった場合のチェック
if [ ! -e "${files[0]}" ]; then
    echo "Warning: No vr_str_*.dat files found."
    # 単一ファイル vr_str.dat がある場合はそれをリストに入れるフォールバック
    if [ -e "vr_str.dat" ]; then
        echo "Found vr_str.dat. Processing single file."
        files=("vr_str.dat")
    else
        exit 1
    fi
fi

# ループ開始
for file_str in "${files[@]}"; do
    # ファイルが存在するか確認（念のため）
    [ -e "$file_str" ] || continue

    # ファイル名からID（番号部分）を抽出する
    # 例: vr_str_01.dat -> 01, vr_str.dat -> (空文字)
    basename=$(basename "$file_str" .dat) # vr_str_01
    id=${basename#vr_str}               # _01
    id=${id#_}                          # 01 (先頭の_があれば削除)

    # IDが空の場合は出力ファイル名に番号を付けない、ある場合は付ける
    if [ -z "$id" ]; then
        file_dip="vr_dip.dat"
        out_name="ffi_RSD_Mf"
        echo "Processing single file..."
    else
        file_dip="vr_dip_${id}.dat"
        out_name="ffi_RSD_Mf_${id}"
        echo "Processing ID: ${id} ..."
    fi

    # 対応する dip ファイルが存在するか確認
    if [ ! -f "$file_dip" ]; then
        echo "Warning: Pair file $file_dip not found. Skipping."
        continue
    fi

    # --- ここから図の描画処理 (変数を file_str, file_dip に置き換え) ---

    aspect=2.0 # larger value means higher
    
    # gmtinfoで各ファイルの範囲を取得
    strinfo=($(gmt gmtinfo "$file_str" -C))
    dipinfo=($(gmt gmtinfo "$file_dip" -C))
    
    strmin=${strinfo[0]}; strmax=${strinfo[1]}
    dipmin=${dipinfo[0]}; dipmax=${dipinfo[1]}
    tmax=${strinfo[3]}
    
    strr=0/${tmax}/${strmin}/${strmax}
    dipr=0/${tmax}/${dipmin}/${dipmax}

    # サイズ計算（ゼロ除算回避のため分母チェックなどを入れるとなお良いが、元のロジックを踏襲）
    # 分母が0になる（最大=最小）場合の簡易ガード
    str_diff=$(echo "${strmax} - ${strmin}" | bc)
    dip_diff=$(echo "${dipmax} - ${dipmin}" | bc)
    
    if [ "$(echo "$str_diff == 0" | bc)" -eq 1 ]; then str_diff=1; fi
    if [ "$(echo "$dip_diff == 0" | bc)" -eq 1 ]; then dip_diff=1; fi

    ratio=$(echo " scale=4;${tmax}/(${str_diff})" | bc)
    strwidth=10
    mapheight=$(echo "${strwidth}*${ratio}*${aspect}" | bc)
    
    ratio=$(echo " scale=4;${tmax}/(${dip_diff})" | bc)
    dipheight=$(echo "(${mapheight}/${ratio})/${aspect}" | bc)
    
    echo ${strinfo[5]} > tmp; echo ${dipinfo[5]} >> tmp; info=($(gmt gmtinfo tmp -C))
    slipmax=${info[1]}
    slipmax=$(echo " $slipmax * 1.01" | bc)  # To avoid Non-val
    interval=$(echo "scale=2;$slipmax/5" | bc)
    half=$(echo "scale=2;$slipmax/2" | bc)
    xmv=$(echo "${strwidth}+0.05" | bc)

    #---
    cororheight=$(echo "${mapheight}+${dipheight}*1.2" | bc)
    # echo cororheight${cororheight}, mapheight${mapheight}, dipheight ${dipheight}

    # gold ratio
    # widthratio=1.618
    # shilver ratio
    widthratio=1.414
    rat=$(echo "scale=5;${strwidth}/(${widthratio}*${cororheight})"|bc)
    cororheight=$(echo "${rat}*${cororheight}" | bc)
    mapheight=$(echo "${rat}*${mapheight}" | bc)
    dipheight=$(echo "${rat}*${dipheight}" | bc)
    echo "ID:${id} -> cororheight:${cororheight}, mapheight:${mapheight}, dipheight:${dipheight}"

    # 【修正】yshiftの計算を変更
    # 元のコード(重複): yshift=$(echo "${dipheight}*1.2" | bc)
    # 元のコード(過剰): yshift=$(echo "${mapheight} + 2.0" | bc)
    # 修正後: Dip図の高さ(dipheight) + 余白(1.5cm) 分だけ下に移動する
    yshift=$(echo "${dipheight} + 1.5" | bc)

    #---
    # gmt begin の出力ファイル名を out_name に変更
    gmt begin "$out_name" png
      # strike
      gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
      gmt basemap -R${strr} -JX${strwidth}/-${mapheight} -BWsne -Bx+l"Time (s)" -By+l"Strike (km)"
      
      # file_str を使用
      awk '{print $1,$2,$3}' "$file_str" | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${strr}
      gmt grdimage tmp.grd -C
      gmt grdcontour tmp.grd -C${interval}
      gmt basemap -Bnswe
      
      # dip
      # Strike図の下に移動 (-Y-${yshift})
      gmt basemap -R${dipr} -JX${strwidth}/${dipheight} -BWSne -Bx+l"Time (s)" -By+l"Dip (km)" -Y-${yshift}
      
      # file_dip を使用
      awk '{print $1,$2,$3}' "$file_dip" | awk '{print $2,$1,$3}' | gmt xyz2grd -Gtmp.grd -I0.1 -R${dipr}
      gmt grdimage tmp.grd -C
      gmt grdcontour tmp.grd -C${interval}
      gmt colorbar -DjBR+w${cororheight}/0.2+jBL+o0.05/0.0 -Ba${half}+l"Potency-rate Density (m/s)"
      gmt basemap -Bnswe
      echo 0 0  `pwd`      | gmt text -F+f6,0,black+jLB   -R0/10/0/0.5  -JX10/1.  -Y-1.2
    gmt end show

done

# ループ終了後の後始末
rm -rf gmt.history gmt.conf tmp tmp.grd
echo "Done."