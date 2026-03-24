#!/bin/bash
# n_vector.datのシュミット投影（非クラスタリング結果）

# --- 設定 ---
# 入力データファイル名
INFILE="schmidt_vector.dat"
# 出力画像ファイル名
OUTFILE="ffi_schmidt_plot_nv"
# 図の縮尺 (1単位あたりのcm)
SCALE=5

# --- 1. データファイルの存在確認 ---
if [ ! -f "$INFILE" ]; then
  echo "エラー: データファイル $INFILE が見つかりません。先にschimdt_cal.pyを実行します."
  schimdt_cal.py
fi

# --- 2. sliprate (8列目) の最大値自動取得 ---
# -i7 は 0から数えて7番目 = 8列目を意味します
info=($(gmt info "$INFILE" -i7 -C))
slipmax=${info[1]}

# 安全策: slipmaxが取得できない、あるいは0以下の場合のデフォルト値
if [ -z "$slipmax" ] || [ $(echo "$slipmax <= 0" | bc -l) -eq 1 ]; then
    slipmax=1.0
    echo "Info: slipmaxをデフォルト値(1.0)に設定しました。"
else
    echo "Sliprate Max: ${slipmax}"
fi

# --- 3. カラーパレット (CPT) の作成 ---
# hotパレットを使用。-I で色を反転 (白->赤->黒) させ、背景白で見やすくします。
gmt makecpt -Chot -T0/${slipmax} -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white > slip_color.cpt

# --- 4. プロット開始 ---
gmt begin $OUTFILE
  gmt figure $OUTFILE png A+m0.5c

  # サブプロット開始 (1行2列)
  # -Fs: パネルサイズ指定 (12cm x 12cm)
  gmt subplot begin 1x2 -Fs12c/12c -M0.1c -R-1.2/1.2/-1.2/1.2 -B+n

    # === 左側 (0): n vector ===
    gmt subplot set 0
    
    # タイトルを表示
    echo 0 1.25 "n vector" | gmt text -F+f18p,Helvetica-Bold,black -N

    # 外枠の円
    echo 0 0 | gmt plot -Sc$((${SCALE}*2))c -W1p,black

    # 方位記号
    echo 0 1.1 N | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 1.1 0 E | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 0 -1.1 S | gmt text -F+f12p,Helvetica-Bold,black -N
    echo -1.1 0 W | gmt text -F+f12p,Helvetica-Bold,black -N

    # nベクトル (丸) のプロット
    # 1. sort -k8 -n : 8列目(sliprate)で昇順(小さい順)にソート
    # 2. awk : n1_X($4), n1_Y($5), sliprate($8) を使用
    sort -k8 -n "$INFILE" | \
    awk -F'\t' '{print $4, $5, $8}' | \
    gmt plot -Sc0.2c -W0.2p,black -Cslip_color.cpt

    # === 右側 (1): v vector ===
    gmt subplot set 1
    
    # タイトルを表示
    echo 0 1.25 "v vector" | gmt text -F+f18p,Helvetica-Bold,black -N

    # 外枠の円
    echo 0 0 | gmt plot -Sc$((${SCALE}*2))c -W1p,black

    # 方位記号
    echo 0 1.1 N | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 1.1 0 E | gmt text -F+f12p,Helvetica-Bold,black -N
    echo 0 -1.1 S | gmt text -F+f12p,Helvetica-Bold,black -N
    echo -1.1 0 W | gmt text -F+f12p,Helvetica-Bold,black -N

    # vベクトル (四角) のプロット
    # 1. sort -k8 -n : 8列目(sliprate)で昇順(小さい順)にソート
    # 2. awk : n2_X($6), n2_Y($7), sliprate($8) を使用
    sort -k8 -n "$INFILE" | \
    awk -F'\t' '{print $6, $7, $8}' | \
    gmt plot -Ss0.2c -W0.2p,black -Cslip_color.cpt

  gmt subplot end

  # カラーバー (共通)
  # 図の下中央に配置
  gmt colorbar -Cslip_color.cpt -DJBC+w15c/0.3c+o0c/0c+h -Baf -B+l"Slip rate"

gmt end show

# 一時ファイルの削除
rm -f slip_color.cpt gmt.conf gmt.history