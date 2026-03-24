# fpsset.py
#
# This Python module replicates the functionality of the Fortran
# subroutine `fpsset` from the FPSPACK library. It defines
# common constants used by other seismological routines.
#
# In Python, these constants are typically defined at the module level
# and are available once the module is imported. The `fpsset_py` function
# is provided for nominal compatibility with code expecting a setup call,
# but the constants are directly accessible via the FPSCOM_CONSTANTS dictionary.

import math

# 定数 (FortranのCOMMON /fpscom/ に対応)
# これらの値はFPSPACK.FOR内のfpssetサブルーチンから取得
FPSCOM_CONSTANTS = {
    'amistr': -360.0,       # strike lower limit
    'amastr': 360.0,        # strike upper limit
    'amidip': 0.0,          # dip lower limit
    'amadip': 90.0,         # dip upper limit
    'amirak': -360.0,       # rake lower limit
    'amarak': 360.0,        # rake upper limit
    'amitre': -360.0,       # trend lower limit
    'amatre': 360.0,        # trend upper limit
    'amiplu': 0.0,          # plunge lower limit
    'amaplu': 90.0,         # plunge upper limit
    'orttol': 2.0,          # orthogonality tolerance (degrees)
    'ovrtol': 0.001,        # dip overtaking tolerance
    'tentol': 0.0001,       # moment tensor symmetry tolerance
    'dtor': math.pi / 180.0, # degree to radians conversion factor (0.017453292519943296)
    'c360': 360.0,
    'c90': 90.0,
    'c0': 0.0,
    'c1': 1.0,
    'c2': 2.0,
    'c3': 3.0,
    'io': 6                  # error messages file unit (Fortran unit number)
}

# FortranのSAVE ifl; DATA ifl/0/; if(ifl.eq.0) then ... ifl=1; endif
# のロジックは、Pythonモジュールでは不要です。
# モジュールレベルの変数はインポート時に一度だけ初期化されます。

def fpsset_py():
    """
    Python版のfpsset関数。
    Fortran版のfpssetサブルーチン呼び出しのインターフェースを模倣します。
    実際の定数設定はモジュールレベルで行われているため、
    この関数は主に定数が利用可能であることを確認する役割や、
    ログ出力などに使用できます。
    
    この実装では、定数は既に FPSCOM_CONSTANTS 辞書として定義されているため、
    この関数自体は実質的な処理を行いません。
    """
    # print("fpsset.py: Constants are defined in FPSCOM_CONSTANTS dictionary.")
    # 必要であれば、ここで何らかの初期化確認やログ出力を追加できます。
    pass

if __name__ == '__main__':
    print("--- Testing fpsset.py ---")
    
    # fpsset_py() を呼び出す (実質的には何もしないが、呼び出しの確認)
    fpsset_py()
    
    print("FPSCOM_CONSTANTS dictionary content:")
    for key, value in FPSCOM_CONSTANTS.items():
        print(f"  {key}: {value}")
        
    # 特定の定数の使用例
    print(f"\nExample usage:")
    print(f"  Degrees to radians (dtor): {FPSCOM_CONSTANTS['dtor']:.10f}")
    print(f"  Orthogonality tolerance (orttol): {FPSCOM_CONSTANTS['orttol']} degrees")

    # Fortranのdtorと比較
    fortran_dtor = 0.017453292519943296
    python_dtor = math.pi / 180.0
    print(f"\nComparison of dtor values:")
    print(f"  Fortran dtor value used in FPSPACK: {fortran_dtor}")
    print(f"  Python math.pi / 180.0:         {python_dtor}")
    assert math.isclose(FPSCOM_CONSTANTS['dtor'], fortran_dtor), \
        "dtor value mismatch with typical Fortran value"
    print("  dtor value in FPSCOM_CONSTANTS matches expected value.")

