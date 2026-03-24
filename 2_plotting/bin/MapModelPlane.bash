#!/bin/bash
if [ $# -ne 11 ]; then
  echo -e \\n " MapModelPlane.bash [strike] [dip] [lat] [lon] [depth] [xx] [yy] [mn] [nn] [m0] [n0]" \\n
  # 北海道非矩形
  # MapModelPlane.bash 179.00 0.00 42.851 139.197 16.5 5 5 38 11 12 5
  exit 1
fi
echo ${1} ${2} 45. ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} 1. | fgenmpinfo # 1.: spatial knot / 2.: subfault
#---
gmt set MAP_TICK_LENGTH 2p
gmt set FONT_LABEL 11p
gmt set FONT_ANNOT_PRIMARY 11p
gmt set MAP_FRAME_PEN 0.8p
gmt set MAP_TICK_PEN 0.9p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_ANNOT_OFFSET 2p
#---
info=($(gmt gmtinfo knotedge.dat -C))
if [ `echo "${info[4]}<0." | bc` == 1 ]; then echo " !!! WARNING !!! Model plane is beyond the ground surface"; fi
xmin=$(echo "${info[2]}-0.1" | bc); xmax=$(echo "${info[3]}+0.1" | bc)
ymin=$(echo "${info[0]}-0.1" | bc); ymax=$(echo "${info[1]}+0.1" | bc)
r=${xmin}/${xmax}/${ymin}/${ymax}
mapwidth=9; mapheight=$(echo ${xmax} ${ymax} | gmt mapproject -JM${mapwidth} -R${r} | awk '{print $2}')
gmt begin MapModelPlane pdf
  gmt coast -R${r} -JM${mapwidth} -BWSne -Ba1 -Df -Ggray

  if [ ! -e "zeroone.txt" ]; then

    #--- faultline 1
    if [ `echo "${2}!=0." | bc` == 1 ]; then
      awk '{print $2,$1}' knotedge.dat | gmt plot -Wdotted
    else
      awk '{print $2,$1}' knotedge.dat | gmt plot
    fi
    #--- circle
    # awk '{print $4,$3}' knotplot.dat | gmt plot -Sp0.08 -W1,white
    #--- gradation circle
    gmt makecpt -Cviridis -T${info[4]}/${info[5]} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
    awk '{print $4,$3,$5}' knotplot.dat | gmt plot -Sp0.08 -W1,white -C
    gmt colorbar -DjBR+w${mapheight}/0.1+jBL+o0.05/0.0 -Ba+l"Depth (km)"
    #--- spatial knot square
    #     for l in $(seq 1 $(cat knotcorner.dat | wc -l )); do
    #       gmt plot << END
    #     $(awk '{if(NR=='$l') print $4,$3}' knotcorner.dat)
    #     $(awk '{if(NR=='$l') print $7,$6}' knotcorner.dat)
    #     $(awk '{if(NR=='$l') print $10,$9}' knotcorner.dat)
    #     $(awk '{if(NR=='$l') print $13,$12}' knotcorner.dat)
    #     $(awk '{if(NR=='$l') print $4,$3}' knotcorner.dat)
    # END
    #     done
  if [ -f usgs_events.csv ]; then
      awk -F, 'NR==1{for(i=1;i<=NF;i++){if($i=="longitude")lon=i;if($i=="latitude")lat=i;if($i=="depth")dep=i;if($i=="mag")mag=i}} NR>1{print $lon, $lat, $dep, 0.04*$mag}' usgs_events.csv  | gmt plot -Sc -W0.2,black -C  -t30
  fi
   
    #--- faultline 2
    if [ `echo "${2}!=0." | bc` == 1 ]; then
      awk '{{if($3=='${info[4]}') print $2,$1}}' knotedge.dat | gmt plot -W1
    fi

  else
    # --- 準備：zeroone.txtからzeroone_meta.txtを生成 ---
    mn=$8
    nn=$9
    meta_mn=$((mn + 2))
    meta_nn=$((nn + 2))
    top_border=$(for i in $(seq 1 $meta_mn); do echo -n "0 "; done)
    echo "$top_border" > zeroone_meta.txt
    while read -r line; do
        echo "0 $line 0" >> zeroone_meta.txt
    done < zeroone.txt
    bottom_border=$top_border
    echo "$bottom_border" >> zeroone_meta.txt

    # --- ステップ1: 条件に合うセルを '3' に置き換えた zeroone_modified.txt を作成 ---
    awk '
        { for (i=1; i<=NF; i++) { a[NR, i] = $i; } num_cols = NF; }
        END {
            num_rows = NR;
            for (r=1; r<=num_rows; r++) { for (c=1; c<=num_cols; c++) {
                if (r > 1 && r < num_rows && c > 1 && c < num_cols) {
                    if (a[r, c] == 1) {
                        neighbor_sum = a[r-1, c-1] + a[r-1, c] + a[r-1, c+1] + a[r, c-1] + a[r, c+1] + a[r+1, c-1] + a[r+1, c] + a[r+1, c+1];
                        if (neighbor_sum == 7) { b[r, c] = 3; } else { b[r, c] = a[r, c]; }
                    } else { b[r, c] = a[r, c]; }
                } else { b[r, c] = a[r, c]; }
            } }
            for (r=1; r<=num_rows; r++) {
                for (c=1; c<=num_cols; c++) { printf "%s%s", b[r, c], (c==num_cols ? "" : " "); }
                printf "\n";
            }
        }
    ' zeroone_meta.txt > zeroone_modified.txt

    # --- 【修正点】ステップ1.5: 除外マスクファイルを効率的に作成 ---
    EXCLUSION_MASK="exclusion_mask.txt"
    awk -v nn=${nn} -v mn=${mn} '
        # 最初に zeroone_modified.txt を読み込み、値が3のセルの中心座標(n,m)を配列に記録する
        FNR==NR {
            n_orig = nn - FNR + 2; # ファイル行番号(FNR)を元のグリッドのn座標に変換
            for (m_orig=1; m_orig<=mn; m_orig++) {
                if ($(m_orig+1) == 3) {
                    center_count++;
                    n_centers[center_count] = n_orig;
                    m_centers[center_count] = m_orig;
                }
            }
            next;
        }
        # 2回目の読み込みの前に、マスクを生成して出力する (BEGINブロックの代わり)
        FNR==1 && FNR!=NR {
            # まずマスク配列を全て0で初期化
            for (r=1; r<=nn; r++) { for (c=1; c<=mn; c++) mask[r,c] = 0; }
            # 記録しておいた中心座標を元に、3x3領域のマスクを1に設定
            for (i=1; i<=center_count; i++) {
                nc = n_centers[i];
                mc = m_centers[i];
                for (dn=-1; dn<=1; dn++) {
                    for (dm=-1; dm<=1; dm++) {
                        r = nc + dn; c = mc + dm;
                        if (r>=1 && r<=nn && c>=1 && c<=mn) {
                            mask[r, c] = 1;
                        }
                    }
                }
            }
            # 生成したマスクファイルを出力
            for (r=1; r<=nn; r++) {
                for (c=1; c<=mn; c++) {
                    printf "%s%s", mask[r,c], (c==mn ? "" : " ");
                }
                printf "\n";
            }
        }
    ' zeroone_modified.txt zeroone_modified.txt > ${EXCLUSION_MASK}

    # ===============================================================
    # --- 既存の線画処理 (ここのロジックは変更なし) ---
    # ===============================================================
    for n in $(seq 1 $9 | sort -nr); do
      for m in $(seq 1 $8); do
        if [ $(awk '{if (NR=='$(($9+1-$n))') print $'$m'}' zeroone.txt) -eq 1 ]; then
          awk '{if($1=='$n' && $2=='$m') print $4,$3}' knotplot.dat | gmt plot -Sp0.08 -W1,white
          if [ $m -eq 1 ] || [ $(awk '{if (NR=='$(($9+1-$n))') print $'$(($m-1))'}' zeroone.txt) -eq 0 ]; then # left
            mask_val_current=$(awk -v r=${n} -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            mask_val_neighbor=$(awk -v r=${n} -v c=$((m-1)) 'NR==r{print $c}' ${EXCLUSION_MASK})
            if ! ( [ "${mask_val_current:-0}" -eq 1 ] && [ "${mask_val_neighbor:-0}" -eq 1 ] ); then
              gmt plot << END
    $(awk '{if($1=='$n' && $2=='$m') print $4,$3}' knotcorner.dat)
    $(awk '{if($1=='$n' && $2=='$m') print $7,$6}' knotcorner.dat)
END
            fi
          fi
          if [ $m -eq $8 ] || [ $(awk '{if (NR=='$(($9+1-$n))') print $'$(($m+1))'}' zeroone.txt) -eq 0 ]; then # right
            mask_val_current=$(awk -v r=${n} -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            mask_val_neighbor=$(awk -v r=${n} -v c=$((m+1)) 'NR==r{print $c}' ${EXCLUSION_MASK})
            if ! ( [ "${mask_val_current:-0}" -eq 1 ] && [ "${mask_val_neighbor:-0}" -eq 1 ] ); then
              gmt plot << END
    $(awk '{if($1=='$n' && $2=='$m') print $10,$9}' knotcorner.dat)
    $(awk '{if($1=='$n' && $2=='$m') print $13,$12}' knotcorner.dat)
END
            fi
          fi
          if [ $n -eq $9 ] || [ $(awk '{if (NR=='$(($9-$n))') print $'$m'}' zeroone.txt) -eq 0 ]; then # top
            mask_val_current=$(awk -v r=${n} -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            mask_val_neighbor=$(awk -v r=$((n+1)) -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            if ! ( [ "${mask_val_current:-0}" -eq 1 ] && [ "${mask_val_neighbor:-0}" -eq 1 ] ); then
              gmt plot << END
    $(awk '{if($1=='$n' && $2=='$m') print $7,$6}' knotcorner.dat)
    $(awk '{if($1=='$n' && $2=='$m') print $10,$9}' knotcorner.dat)
END
            fi
          fi
          if [ $n -eq 1 ] || [ $(awk '{if (NR=='$(($9+2-$n))') print $'$m'}' zeroone.txt) -eq 0 ]; then # bottom
            mask_val_current=$(awk -v r=${n} -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            mask_val_neighbor=$(awk -v r=$((n-1)) -v c=${m} 'NR==r{print $c}' ${EXCLUSION_MASK})
            if ! ( [ "${mask_val_current:-0}" -eq 1 ] && [ "${mask_val_neighbor:-0}" -eq 1 ] ); then
              gmt plot << END
    $(awk '{if($1=='$n' && $2=='$m') print $4,$3}' knotcorner.dat)
    $(awk '{if($1=='$n' && $2=='$m') print $13,$12}' knotcorner.dat)
END
            fi
          fi
        fi
      done
    done

    # ===============================================================
    # --- 赤い線の描画処理 (変更なし) ---
    # ===============================================================
    TMP_RED_LINES="tmp_red_lines.txt"
    rm -f ${TMP_RED_LINES}
    for n in $(seq 1 ${nn}); do
      for m in $(seq 1 ${mn}); do
        meta_n=$(( nn - n + 2 )); meta_m=$(( m + 1 ));
        cell_value=$(awk -v row="${meta_n}" -v col="${meta_m}" 'NR==row {print $col}' zeroone_modified.txt)
        if [ "${cell_value:-0}" -eq 3 ]; then
          n_tl=$((n-1)); m_tl=$((m-1)); n_t=$((n-1)); m_t=$m; n_tr=$((n-1)); m_tr=$((m+1));
          n_r=$n; m_r=$((m+1)); n_br=$((n+1)); m_br=$((m+1)); n_b=$((n+1)); m_b=$m;
          n_bl=$((n+1)); m_bl=$((m-1)); n_l=$n; m_l=$((m-1));
          awk \
            -v nTL="$n_tl" -v mTL="$m_tl" -v nT="$n_t" -v mT="$m_t" -v nTR="$n_tr" -v mTR="$m_tr" \
            -v nR="$n_r" -v mR="$m_r" -v nBR="$n_br" -v mBR="$m_br" -v nB="$n_b" -v mB="$m_b" \
            -v nBL="$n_bl" -v mBL="$m_bl" -v nL="$n_l" -v mL="$m_l" '
            { coords[$1, $2] = $4 " " $3; }
            END {
                if ( (nTL, mTL) in coords && (nT, mT) in coords && (nTR, mTR) in coords && (nR, mR) in coords && (nBR, mBR) in coords && (nB, mB) in coords && (nBL, mBL) in coords && (nL, mL) in coords ) {
                    print ">"; print coords[nTL, mTL]; print coords[nT, mT]; print coords[nTR, mTR]; print coords[nR, mR]; print coords[nBR, mBR]; print coords[nB, mB]; print coords[nBL, mBL]; print coords[nL, mL]; print coords[nTL, mTL];
                }
            }
          ' knotplot.dat >> ${TMP_RED_LINES}
        fi
      done
    done
    if [ -s ${TMP_RED_LINES} ]; then
      gmt plot ${TMP_RED_LINES} -W1p,red
    fi

    # 一時ファイルを削除
    rm -f ${TMP_RED_LINES} zeroone_modified.txt ${EXCLUSION_MASK}
  fi

  #--- epicenter
  echo ${4} ${3} | gmt plot -Sa0.2c
gmt end show
rm -rf gmt.history gmt.conf