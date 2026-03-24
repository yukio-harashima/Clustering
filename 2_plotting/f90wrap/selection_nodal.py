# selection_nodal.py
#
# This Python module provides a function equivalent to the Fortran
# subroutine `selection_nodal`. It selects one of two nodal planes
# based on its proximity to a reference fault plane.

import numpy as np
import math

def calculate_normal_vector(strike_deg, dip_deg, rad_conversion_factor):
    """
    断層面の走向(strike)と傾斜(dip)から法線ベクトルを計算します。
    法線ベクトルは North-East-Down (NED) 座標系で返されます。
    走向は北から時計回り、傾斜は水平面からの角度です。
    法線は断層の下向きを指します。

    Args:
        strike_deg (float): 走向 (度).
        dip_deg (float): 傾斜 (度).
        rad_conversion_factor (float): 度からラジアンへの変換係数.

    Returns:
        numpy.ndarray: 3要素の法線ベクトル [nx, ny, nz] (N, E, D).
    """
    strike_rad = strike_deg * rad_conversion_factor
    dip_rad = dip_deg * rad_conversion_factor

    # 法線ベクトルの成分 (Aki & Richards convention for fault normal pointing downwards)
    # n_z = -cos(dip)  (Down component)
    # n_x = -sin(dip)sin(strike) (North component)
    # n_y =  sin(dip)cos(strike) (East component)
    
    nx = -math.sin(dip_rad) * math.sin(strike_rad)  # North component
    ny =  math.sin(dip_rad) * math.cos(strike_rad)  # East component
    nz = -math.cos(dip_rad)                         # Down component
    
    return np.array([nx, ny, nz], dtype=np.float64)

def selection_nodal_py(
    str1_deg, dip1_deg, 
    str2_deg, dip2_deg, 
    str0_deg, dip0_deg, 
    rad_conversion_factor
    ):
    """
    Fortranのselection_nodalサブルーチンのPython版。
    与えられた2つの節面のうち、参照断層面(str0, dip0)により法線ベクトルが近い方を
    主断層面として選択します。

    Args:
        str1_deg (float): 候補1の走向 (度).
        dip1_deg (float): 候補1の傾斜 (度).
        str2_deg (float): 候補2の走向 (度).
        dip2_deg (float): 候補2の傾斜 (度).
        str0_deg (float): 参照断層面の走向 (度).
        dip0_deg (float): 参照断層面の傾斜 (度).
        rad_conversion_factor (float): 度からラジアンへの変換係数.

    Returns:
        int: 選択された節面を示す番号 (1 または 2).
    """
    # 各断層面の法線ベクトルを計算
    n1 = calculate_normal_vector(str1_deg, dip1_deg, rad_conversion_factor)
    n2 = calculate_normal_vector(str2_deg, dip2_deg, rad_conversion_factor)
    n0 = calculate_normal_vector(str0_deg, dip0_deg, rad_conversion_factor)

    # 法線ベクトル同士の内積を計算
    # 内積の絶対値が大きいほど、法線ベクトルの向きが近い（または真逆で平行）ことを示す
    # (法線は向きを持つため、完全に同じ向きか真逆の向きが最も「近い」と解釈される)
    dot_n1_n0 = np.dot(n1, n0)
    dot_n2_n0 = np.dot(n2, n0)

    # 内積の絶対値が大きい方を選択
    if abs(dot_n1_n0) >= abs(dot_n2_n0):
        nodal_p = 1
    else:
        nodal_p = 2
        
    return nodal_p

if __name__ == '__main__':
    print("--- Testing selection_nodal.py ---")
    
    # 度からラジアンへの変換係数
    RAD = math.pi / 180.0

    # テストケース1
    s1, d1 = 30.0, 60.0
    s2, d2 = 210.0, 30.0 # s1,d1の共役面に近い例 (走向が180度違い、傾斜が異なる)
    s0, d0 = 45.0, 55.0  # 参照面
    
    print(f"\nTest Case 1:")
    print(f"  Plane 1: Strike={s1}, Dip={d1}")
    print(f"  Plane 2: Strike={s2}, Dip={d2}")
    print(f"  Reference Plane: Strike={s0}, Dip={d0}")
    
    selected_plane = selection_nodal_py(s1, d1, s2, d2, s0, d0, RAD)
    print(f"  Selected plane: {selected_plane}")

    n1_vec = calculate_normal_vector(s1, d1, RAD)
    n2_vec = calculate_normal_vector(s2, d2, RAD)
    n0_vec = calculate_normal_vector(s0, d0, RAD)
    print(f"  Normal vector N1: {n1_vec}")
    print(f"  Normal vector N2: {n2_vec}")
    print(f"  Normal vector N0: {n0_vec}")
    print(f"  Dot(N1, N0): {np.dot(n1_vec, n0_vec):.4f} (abs: {abs(np.dot(n1_vec, n0_vec)):.4f})")
    print(f"  Dot(N2, N0): {np.dot(n2_vec, n0_vec):.4f} (abs: {abs(np.dot(n2_vec, n0_vec)):.4f})")


    # テストケース2: 参照面が平面2に近い場合
    s1_t2, d1_t2 = 10.0, 80.0
    s2_t2, d2_t2 = 100.0, 20.0 
    s0_t2, d0_t2 = 95.0, 25.0  # 参照面 (平面2に近い)

    print(f"\nTest Case 2:")
    print(f"  Plane 1: Strike={s1_t2}, Dip={d1_t2}")
    print(f"  Plane 2: Strike={s2_t2}, Dip={d2_t2}")
    print(f"  Reference Plane: Strike={s0_t2}, Dip={d0_t2}")

    selected_plane_t2 = selection_nodal_py(s1_t2, d1_t2, s2_t2, d2_t2, s0_t2, d0_t2, RAD)
    print(f"  Selected plane: {selected_plane_t2}")

    n1_vec_t2 = calculate_normal_vector(s1_t2, d1_t2, RAD)
    n2_vec_t2 = calculate_normal_vector(s2_t2, d2_t2, RAD)
    n0_vec_t2 = calculate_normal_vector(s0_t2, d0_t2, RAD)
    print(f"  Normal vector N1: {n1_vec_t2}")
    print(f"  Normal vector N2: {n2_vec_t2}")
    print(f"  Normal vector N0: {n0_vec_t2}")
    print(f"  Dot(N1, N0): {np.dot(n1_vec_t2, n0_vec_t2):.4f} (abs: {abs(np.dot(n1_vec_t2, n0_vec_t2)):.4f})")
    print(f"  Dot(N2, N0): {np.dot(n2_vec_t2, n0_vec_t2):.4f} (abs: {abs(np.dot(n2_vec_t2, n0_vec_t2)):.4f})")

