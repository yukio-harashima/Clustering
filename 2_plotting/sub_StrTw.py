#!/Users/harashima-yukio/.pyenv/versions/3.12.3/bin/python3
# ffi_StrTW_plot_xxx.bash のサブルーチン
import sys
import argparse

def convert_strike(strike):
    """
    走向(Strike)を -90 (N90W) ～ 90 (N90E) の範囲に変換する。
    S と S+180 は等価とみなし、北(0)を中心とした表現にする。
    
    変換ロジック:
      - 270 <= S < 360  -> S - 360
      - 90 < S < 270    -> S - 180
      - 0 <= S <= 90    -> そのまま
    """
    s = strike % 360
    if s >= 270:
        return s - 360
    elif s > 90:
        return s - 180
    return s

def process_data(infile, interval):
    data_list = []
    try:
        with open(infile, 'r') as f:
            for line in f:
                if not line.strip(): continue
                parts = line.strip().split()
                # no列(最後)まで含めると最低でも15列程度必要
                if len(parts) < 15: continue
                
                try:
                    # tw: 3列目(idx 2), stk0: 13列目(idx 12), stk1: 14列目(idx 13), no: 最後(idx -1)
                    tw_val = float(parts[2])
                    stk0_val = float(parts[12])
                    stk1_val = float(parts[13])
                    no_val = float(parts[-1])
                    
                    time_scaled = tw_val * interval
                    
                    # plane_id: 0 for stk0, 1 for stk1
                    data_list.append({
                        'time': time_scaled,
                        'strike': convert_strike(stk0_val),
                        'cluster_no': no_val,
                        'plane_id': 0
                    })
                    data_list.append({
                        'time': time_scaled,
                        'strike': convert_strike(stk1_val),
                        'cluster_no': no_val,
                        'plane_id': 1
                    })
                    
                except ValueError:
                    continue

    except FileNotFoundError:
        sys.exit(1)

    # 時間順にソート
    data_list.sort(key=lambda x: x['time'])
    
    for d in data_list:
        # 4列目に plane_id (0 or 1) を追加して出力
        print(f"{d['time']:.6f} {d['strike']:.6f} {d['cluster_no']:.0f} {d['plane_id']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("--interval", type=float, default=1.0)
    args = parser.parse_args()
    process_data(args.infile, args.interval)