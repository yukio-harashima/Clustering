#!/bin/bash
# 東西方向の走向方向のクラスターごとのプロット（東ー北ー西）　

# --- 設定 ---
# 時間間隔 (秒/スナップショット) -> Pythonスクリプトに渡す
tw_interval=0.5

# 【追加】プロットするクラスター番号の指定
# "ALL" : すべてのクラスターをプロット
# 数値  : 指定した番号のクラスターのみプロット (例: 1)
# TARGET_CLUSTER="ALL"
TARGET_CLUSTER=1

# 入力データファイル名  calc_fault_v3.pyからの出力
INFILE="clusteringSnap2_mod.dat"
# 出力画像ファイル名
OUTFILE="ffi_strike_time_plot_v2_cl${TARGET_CLUSTER}"
# データ処理用Pythonスクリプト
PY_SCRIPT="sub_StrTw.py"

# 中間データファイル名
# PROCESSED_DATA: 合算データ (範囲計算用)
# DATA0: stk0用, DATA1: stk1用
PROCESSED_DATA="processed_data.tmp"
DATA0="processed_data0.tmp"
DATA1="processed_data1.tmp"


# --- 1. ファイルの存在確認 ---
if [ ! -f "$INFILE" ]; then
  echo "エラー: データファイル $INFILE が見つかりません。"
  exit 1
fi

# if [ ! -f "$PY_SCRIPT" ]; then
#   echo "エラー: Pythonスクリプト $PY_SCRIPT が見つかりません。"
#   echo "同じディレクトリに sub_StTw.py を配置してください。"
#   exit 1
# # fi

# # 実行権限の確認と付与
# if [ ! -x "$PY_SCRIPT" ]; then
#   chmod +x "$PY_SCRIPT"
# fi

# --- 2. データ処理 (Pythonスクリプト実行と分割) ---
echo "Running python script to process data..."

# Pythonスクリプトを実行し、全データを一時ファイルに保存
"$PY_SCRIPT" "$INFILE" --interval "$tw_interval" > "$PROCESSED_DATA"

# 処理結果が空でないか確認
if [ ! -s "$PROCESSED_DATA" ]; then
  echo "エラー: データ処理に失敗しました（出力が空です）。"
  rm -f "$PROCESSED_DATA"
  exit 1
fi

# データを stk0 (ID=0) と stk1 (ID=1) に分割
# 出力列: 1:Time, 2:Strike, 3:ClusterNo
awk '$4 == 0 {print $1, $2, $3}' "$PROCESSED_DATA" > "$DATA0"
awk '$4 == 1 {print $1, $2, $3}' "$PROCESSED_DATA" > "$DATA1"

echo "Data separated into:"
echo "  - $DATA0 (stk0)"
echo "  - $DATA1 (stk1)"

# --- 3. 範囲と最大値の自動取得 ---
echo "Calculating data ranges..."

# 時間軸 (1列目) の範囲取得 (合算データから取得)
info_t=($(gmt info "$PROCESSED_DATA" -i0 -C))
tmin=0                  # 開始時間を0秒固定にする場合
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

# --- 5. Y軸カスタムアノテーションファイルの作成 ---
cat << EOF > y_custom_annot.txt
-90 ag 90W
-45 af 45W
0   af N
45  af 45E
90  ag 90E
EOF

# --- 6. プロット開始 ---
gmt begin $OUTFILE png

    # --- ベースマップ ---
    # 横軸: 時間
    # 縦軸: 走向 (-90～90)
    gmt basemap -R${tmin}/${tmax}/-90/90 -JX15c/-10c \
        -BWSne \
        -Bxaf+l"Time (s)" \
        -By+l"Strike" \
        -Byc"y_custom_annot.txt" \
        --MAP_GRID_PEN_PRIMARY=0.2p,gray,dash

    # --- データプロット (stk0) ---
    if [ "$TARGET_CLUSTER" == "ALL" ]; then
        gmt plot "$DATA0" -Sc0.1c -Cvector_color.cpt -W0.1p,black
    else
        awk -v target="$TARGET_CLUSTER" '$3 == target {print $0}' "$DATA0" | \
        gmt plot -Sc0.1c -Cvector_color.cpt -t90
    fi

    # --- データプロット (stk1) ---
    # 必要であればシンボルを変える (例: -Ss0.2c)
    if [ "$TARGET_CLUSTER" == "ALL" ]; then
        gmt plot "$DATA1" -Sc0.1c -Cvector_color.cpt -W0.1p,black
    else
        awk -v target="$TARGET_CLUSTER" '$3 == target {print $0}' "$DATA1" | \
        gmt plot -Sc0.1c -Gblack -t90
    fi

    # --- カラーバー ---
    gmt colorbar -Cvector_color.cpt -DjBR+w5c/0.3c+o0.5c/0.5c+h -Ba1f1 -B+l"Cluster Number" -Y-1.5

gmt end show

# 一時ファイルの削除
rm -f vector_color.cpt y_custom_annot.txt "$PROCESSED_DATA" "$DATA0" "$DATA1" gmt.conf gmt.history