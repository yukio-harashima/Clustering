#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# ffi_DipTw_plot_xxxx.bashのデータ成形用サブルーチン
import sys
import argparse

def process_data(infile, interval):
    data_list = []
    try:
        with open(infile, 'r') as f:
            for line in f:
                # 空行スキップ
                if not line.strip(): continue
                parts = line.strip().split()
                
                # 最低限必要な列数があるか確認
                # no列(最後)まで含めると最低でも16列程度必要
                # 列構成目安: 1:n, 2:m, 3:tw, ..., 15:dip0, 16:dip1, ..., last:no
                if len(parts) < 16: continue
                
                try:
                    # tw: 3列目(idx 2)
                    # dip0: 15列目(idx 14)
                    # dip1: 16列目(idx 15)
                    # no: 最後(idx -1)
                    tw_val = float(parts[2])
                    dip0_val = float(parts[14])
                    dip1_val = float(parts[15])
                    no_val = float(parts[-1])
                    
                    time_scaled = tw_val * interval
                    
                    # plane_id: 0 for dip0, 1 for dip1
                    data_list.append({
                        'time': time_scaled,
                        'dip': dip0_val,
                        'cluster_no': no_val,
                        'plane_id': 0
                    })
                    data_list.append({
                        'time': time_scaled,
                        'dip': dip1_val,
                        'cluster_no': no_val,
                        'plane_id': 1
                    })
                    
                except ValueError:
                    # ヘッダー行や数値変換できない行はスキップ
                    continue

    except FileNotFoundError:
        print(f"Error: File '{infile}' not found.", file=sys.stderr)
        sys.exit(1)

    # 時間順にソート
    data_list.sort(key=lambda x: x['time'])
    
    for d in data_list:
        # 出力フォーマット: Time Dip ClusterNo PlaneID
        print(f"{d['time']:.6f} {d['dip']:.6f} {d['cluster_no']:.0f} {d['plane_id']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process dip data and sort by time.")
    parser.add_argument("infile", help="Input data file path")
    parser.add_argument("--interval", type=float, default=1.0, help="Time interval multiplier")
    args = parser.parse_args()
    
    process_data(args.infile, args.interval)