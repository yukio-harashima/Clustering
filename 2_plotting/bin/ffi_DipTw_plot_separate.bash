#!/bin/bash
# dipの時間変化をクラスターごとにプロット

# --- 設定 ---
# 時間間隔 (秒/スナップショット) -> Pythonスクリプトに渡す
tw_interval=0.5

# 【プロットするクラスター番号の指定】
# "ALL" : すべてのクラスターをプロット
# 数値  : 指定した番号のクラスターのみプロット (例: 1)
# TARGET_CLUSTER="ALL"
TARGET_CLUSTER=4
echo "Info: Target Cluster No: $TARGET_CLUSTER"
LABEL="cluster_${TARGET_CLUSTER}"

# 入力データファイル名
INFILE="clusteringSnap2_mod.dat"
# 出力画像ファイル名
OUTFILE="ffi_dip_time_plot_sep_cl${TARGET_CLUSTER}"
# データ処理用Pythonスクリプト (カレントディレクトリにある前提)
PY_SCRIPT="sub_DipTw.py"

# 中間データファイル名
# PROCESSED_DATA: 合算データ (範囲計算用)
# DATA0: dip0用, DATA1: dip1用
PROCESSED_DATA="processed_dip_data.tmp"
DATA0="processed_dip0.tmp"
DATA1="processed_dip1.tmp"


# --- 1. ファイルの存在確認 ---
if [ ! -f "$INFILE" ]; then
  echo "エラー: データファイル $INFILE が見つかりません。"
  exit 1
fi

# if [ ! -f "$PY_SCRIPT" ]; then
#   echo "エラー: Pythonスクリプト $PY_SCRIPT が見つかりません。"
#   echo "同じディレクトリに sub_DipTw.py を配置し、実行権限(chmod +x)を与えてください。"
#   exit 1
# fi

# # 実行権限のチェック (なければ付与を試みる)
# if [ ! -x "$PY_SCRIPT" ]; then
#   echo "警告: $PY_SCRIPT に実行権限がありません。chmod +x を実行します。"
#   chmod +x "$PY_SCRIPT"
# fi

# --- 2. データ処理 (Pythonスクリプト実行) ---
echo "Running python script to process data..."

# Pythonスクリプトを直接実行し、全データを一時ファイルに保存
# シェバンがあるため python コマンドは不要
"$PY_SCRIPT" "$INFILE" --interval "$tw_interval" > "$PROCESSED_DATA"

# 処理結果が空でないか確認
if [ ! -s "$PROCESSED_DATA" ]; then
  echo "エラー: データ処理に失敗しました（出力が空です）。"
  rm -f "$PROCESSED_DATA"
  exit 1
fi

# データを dip0 (ID=0) と dip1 (ID=1) に分割
# 出力列: 1:Time, 2:Dip, 3:ClusterNo
awk '$4 == 0 {print $1, $2, $3}' "$PROCESSED_DATA" > "$DATA0"
awk '$4 == 1 {print $1, $2, $3}' "$PROCESSED_DATA" > "$DATA1"

echo "Data separated into:"
echo "  - $DATA0 (dip0)"
echo "  - $DATA1 (dip1)"

# --- 3. 範囲と最大値の自動取得 ---
echo "Calculating data ranges..."

# 時間軸 (1列目) の範囲取得
info_t=($(gmt info "$PROCESSED_DATA" -i0 -C))
tmin=0                  # 開始時間を0秒固定
tmax=${info_t[1]}

# クラスター番号 (3列目) の最大値取得
info_no=($(gmt info "$PROCESSED_DATA" -i2 -C))
nomax=${info_no[1]}

# 安全策
if [ -z "$nomax" ]; then
    nomax=1
    echo "Info: nomaxをデフォルト値(1)に設定しました。"
fi

echo "Time range (sec): $tmin - $tmax"
echo "Cluster No max: $nomax"

# --- 4. カラーパレット (CPT) の作成 ---
cpt_max=$(echo "$nomax + 0.5" | bc)
gmt makecpt -Cvik -T0.5/${cpt_max}/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > vector_color.cpt

# --- 5. プロット開始 ---
gmt begin $OUTFILE png

    # --- 節面１ ---
    # 横軸: 時間
    # 縦軸: Dip (0～90)
    gmt basemap -R${tmin}/${tmax}/0/90 -JX15c/-10c \
        -BWSne \
        -Bxaf \
        -Byaf+l"Dip1" \
        --MAP_GRID_PEN_PRIMARY=0.2p,gray,dash

    # --- データプロット (dip0) ---
    if [ "$TARGET_CLUSTER" == "ALL" ]; then
        gmt plot "$DATA0" -Sc0.2c -Cvector_color.cpt -W0.1p,black
    else
        awk -v target="$TARGET_CLUSTER" '$3 == target {print $0}' "$DATA0" | \
        gmt plot -Sc0.1c -Cvector_color.cpt -t90
        # gmt plot -Sc0.1c -Ggray -t90
    fi

    # --- 節面２ ---
    # 横軸: 時間
    # 縦軸: Dip (0～90)
    gmt basemap -R${tmin}/${tmax}/0/90 -JX15c/-10c \
        -BWSne \
        -Bxaf+l"Time (s)" \
        -Byaf+l"Dip2" \
        --MAP_GRID_PEN_PRIMARY=0.2p,gray,dash -Y-11

    # --- データプロット (dip1) ---
    if [ "$TARGET_CLUSTER" == "ALL" ]; then
        gmt plot "$DATA1" -Sc0.2c -Cvector_color.cpt -W0.1p,black
    else
        awk -v target="$TARGET_CLUSTER" '$3 == target {print $0}' "$DATA1" | \
        gmt plot -Sc0.1c -Cvector_color.cpt -t90
        # gmt plot -Sc0.1c -Ggray -t90
    fi


    # --- カラーバー ---
    gmt colorbar -Cvector_color.cpt -DjBR+w5c/0.3c+o0.5c/0.5c+h -Ba1f1 -B+l"Cluster Number" -Y-1.5
    
    # --- ラベル ---
    echo 0 1 ${LABEL} | gmt text -F+f15,0,black+jLB -JX10/1.  

gmt end show

# 一時ファイルの削除
rm -f vector_color.cpt "$PROCESSED_DATA" "$DATA0" "$DATA1" gmt.conf gmt.history