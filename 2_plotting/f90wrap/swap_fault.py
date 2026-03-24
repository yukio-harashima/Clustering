# swap_fault.py
#
# This Python module provides a function equivalent to the Fortran
# subroutine `swap_fault`. It swaps two sets of fault parameters
# (strike, dip, rake) based on a selection flag.

def swap_fault_py(str1, dip1, rake1, str2, dip2, rake2, nodal_p):
    """
    Fortranのswap_faultサブルーチンのPython版。
    nodal_pの値に基づき、2つの断層面のパラメータセットを入れ替えます。

    Args:
        str1 (float): 候補1の走向 (度).
        dip1 (float): 候補1の傾斜 (度).
        rake1 (float): 候補1のすべり角 (度).
        str2 (float): 候補2の走向 (度).
        dip2 (float): 候補2の傾斜 (度).
        rake2 (float): 候補2のすべり角 (度).
        nodal_p (int): 選択フラグ。2の場合にパラメータを入れ替えます。

    Returns:
        tuple: (str1_out, dip1_out, rake1_out, str2_out, dip2_out, rake2_out)
               入れ替え後の断層パラメータセット。
    """
    # 受け取った値をfloat型に変換しておく（念のため）
    s1_out, d1_out, r1_out = float(str1), float(dip1), float(rake1)
    s2_out, d2_out, r2_out = float(str2), float(dip2), float(rake2)

    if nodal_p == 2:
        # パラメータセット1とパラメータセット2を入れ替える
        s_temp = s1_out
        d_temp = d1_out
        r_temp = r1_out
        
        s1_out = s2_out
        d1_out = d2_out
        r1_out = r2_out
        
        s2_out = s_temp
        d2_out = d_temp
        r2_out = r_temp
        
    return s1_out, d1_out, r1_out, s2_out, d2_out, r2_out

if __name__ == '__main__':
    print("--- Testing swap_fault.py ---")

    # テストケース1: nodal_p = 1 (入れ替えなし)
    s1_in, d1_in, r1_in = 30.0, 60.0, 90.0
    s2_in, d2_in, r2_in = 210.0, 30.0, -90.0
    nodal_p_test = 1
    
    print(f"\nTest Case 1: nodal_p = {nodal_p_test}")
    print(f"  Input Plane 1: Strike={s1_in}, Dip={d1_in}, Rake={r1_in}")
    print(f"  Input Plane 2: Strike={s2_in}, Dip={d2_in}, Rake={r2_in}")
    
    res_s1, res_d1, res_r1, res_s2, res_d2, res_r2 = swap_fault_py(
        s1_in, d1_in, r1_in, s2_in, d2_in, r2_in, nodal_p_test
    )
    
    print(f"  Output Plane 1: Strike={res_s1}, Dip={res_d1}, Rake={res_r1}")
    print(f"  Output Plane 2: Strike={res_s2}, Dip={res_d2}, Rake={res_r2}")
    
    assert res_s1 == s1_in and res_d1 == d1_in and res_r1 == r1_in, "Test Case 1 Failed (Plane 1)"
    assert res_s2 == s2_in and res_d2 == d2_in and res_r2 == r2_in, "Test Case 1 Failed (Plane 2)"
    print("  Test Case 1 Passed.")

    # テストケース2: nodal_p = 2 (入れ替えあり)
    nodal_p_test = 2
    
    print(f"\nTest Case 2: nodal_p = {nodal_p_test}")
    print(f"  Input Plane 1: Strike={s1_in}, Dip={d1_in}, Rake={r1_in}")
    print(f"  Input Plane 2: Strike={s2_in}, Dip={d2_in}, Rake={r2_in}")
    
    res_s1, res_d1, res_r1, res_s2, res_d2, res_r2 = swap_fault_py(
        s1_in, d1_in, r1_in, s2_in, d2_in, r2_in, nodal_p_test
    )
    
    print(f"  Output Plane 1: Strike={res_s1}, Dip={res_d1}, Rake={res_r1}")
    print(f"  Output Plane 2: Strike={res_s2}, Dip={res_d2}, Rake={res_r2}")

    assert res_s1 == s2_in and res_d1 == d2_in and res_r1 == r2_in, "Test Case 2 Failed (Plane 1)"
    assert res_s2 == s1_in and res_d2 == d1_in and res_r1 == r1_in, "Test Case 2 Failed (Plane 2)"
    print("  Test Case 2 Passed.")

    # テストケース3: 異なる値で nodal_p = 2
    s1_in_3, d1_in_3, r1_in_3 = 10.0, 20.0, 30.0
    s2_in_3, d2_in_3, r2_in_3 = 40.0, 50.0, 60.0
    nodal_p_test_3 = 2

    print(f"\nTest Case 3: nodal_p = {nodal_p_test_3}")
    print(f"  Input Plane 1: Strike={s1_in_3}, Dip={d1_in_3}, Rake={r1_in_3}")
    print(f"  Input Plane 2: Strike={s2_in_3}, Dip={d2_in_3}, Rake={r2_in_3}")

    res_s1_3, res_d1_3, res_r1_3, res_s2_3, res_d2_3, res_r2_3 = swap_fault_py(
        s1_in_3, d1_in_3, r1_in_3, s2_in_3, d2_in_3, r2_in_3, nodal_p_test_3
    )

    print(f"  Output Plane 1: Strike={res_s1_3}, Dip={res_d1_3}, Rake={res_r1_3}")
    print(f"  Output Plane 2: Strike={res_s2_3}, Dip={res_d2_3}, Rake={res_r2_3}")
    
    assert res_s1_3 == s2_in_3 and res_d1_3 == d2_in_3 and res_r1_3 == r2_in_3, "Test Case 3 Failed (Plane 1)"
    assert res_s2_3 == s1_in_3 and res_d2_3 == d1_in_3 and res_r1_3 == r1_in_3, "Test Case 3 Failed (Plane 2)"
    print("  Test Case 3 Passed.")

