# ar2plp.py
#
# This Python module aims to replicate the functionality of the ar2plp subroutine
# and its core dependencies from the FPSPACK Fortran library.
#
# IMPORTANT NOTES:
# 1.  This is a Python-based re-implementation. While efforts are made to follow
#     the logic of the original Fortran code (e.g., Aki & Richards, 1980, FPSPACK conventions),
#     exact numerical equivalence with compiled Fortran code is not guaranteed due
#     to potential differences in floating-point arithmetic, algorithm details,
#     and library functions (e.g., eigenvalue solvers).
# 2.  The calculation of scalar moments (M0, M1, M0b) and CLVD component (eta)
#     attempts to follow the logic in FPSPACK's ar2pt. However, subtle details
#     might differ. Full fidelity would require a line-by-line port of ar2pt and avec.
# 3.  Error handling (ierr) is basic and may not cover all cases handled by FPSPACK.
# 4.  Thorough validation against trusted results or the original Fortran
#     code is highly recommended if precision and exact compatibility are critical.
# 5.  This script uses NumPy for numerical operations.

import numpy as np
import math

# --- Constants ---
RAD_CONST = np.pi / 180.0
DEG_CONST = 180.0 / np.pi
PI_CONST = np.pi
SQRT2 = np.sqrt(2.0)

# FPSPACK's COMMON /fpscom/ like constants
FPSCOM_CONSTANTS = {
    'amistr': -360.0, 'amastr': 360.0,
    'amidip': 0.0,   'amadip': 90.0,
    'amirak': -360.0, 'amarak': 360.0,
    'amitre': -360.0, 'amatre': 360.0,
    'amiplu': 0.0,   'amaplu': 90.0,
    'orttol': 2.0,
    'ovrtol': 0.001,
    'tentol': 0.0001, # Moment tensor symmetry tolerance
    'dtor': RAD_CONST,
    'c360': 360.0, 'c90': 90.0, 'c0': 0.0, 'c1': 1.0, 'c2': 2.0, 'c3': 3.0,
    'io': 6 
}

# --- Helper Functions ---

def _normalize_vector(vec):
    """正規化されたベクトルを返す。ゼロベクトルの場合はゼロベクトルを返す。"""
    norm = np.linalg.norm(vec)
    if norm < 1e-15:
        return np.zeros_like(vec)
    return vec / norm

def _vector_to_trend_plunge(vector_ned):
    """
    3Dベクトル(North, East, Down)からトレンドとプランジを計算。
    トレンド: 北から時計回りの角度 (0-360度)
    プランジ: 水平面からの下向きの角度 (0-90度)
    """
    vn, ve, vd = vector_ned[0], vector_ned[1], vector_ned[2]
    vh = np.sqrt(vn**2 + ve**2)

    if vh < 1e-9: 
        plunge_deg = 90.0 if vd >= 0 else -90.0 
        trend_deg = 0.0 
    else:
        plunge_rad = np.arctan2(vd, vh)
        plunge_deg = plunge_rad * DEG_CONST
        trend_rad = np.arctan2(ve, vn) 
        trend_deg = trend_rad * DEG_CONST
        if trend_deg < 0:
            trend_deg += 360.0
            
    if plunge_deg < 0:
        plunge_deg = -plunge_deg
        trend_deg = (trend_deg + 180.0) % 360.0
        
    return trend_deg, plunge_deg

def _sdr_from_normal_slip(normal_vec_ned, slip_vec_ned):
    """
    断層面の法線ベクトル(n)とすべりベクトル(u) (N, E, Down 座標系)
    から走向(strike), 傾斜(dip), すべり角(rake)を計算する。
    (Aki & Richards, 1980, pp. 114-116 および FPSPACK の規約を参照)
    """
    n_orig = np.array(normal_vec_ned, dtype=np.float64)
    u_orig = np.array(slip_vec_ned, dtype=np.float64)

    n = _normalize_vector(n_orig)
    u = _normalize_vector(u_orig)

    # 法線ベクトルが上向き(D<0)なら反転 (FPSPACKのnd2plは入力法線を下向きにする処理を含む)
    if n[2] < 0:
        n = -n
        u = -u 

    # Dip (0-90 deg)
    # FPSPACK nd2pl: delta = acos(-anz) where anz is UP-component of normal.
    # If our n is NED (N,E,Down), then anz = -n[2]. So delta = acos(n[2]).
    # This delta is the dip.
    val_for_acos_dip = n[2]
    if val_for_acos_dip > 1.0: val_for_acos_dip = 1.0
    if val_for_acos_dip < -1.0: val_for_acos_dip = -1.0 # Should be positive if n[2]>=0
    dip_rad = np.arccos(val_for_acos_dip)
    dip_deg = dip_rad * DEG_CONST

    # Strike (0-360 deg from North, clockwise)
    sin_dip = np.sin(dip_rad)

    if abs(sin_dip) < 1e-7: # Dip is 0 or 180 (effectively 0 or 90 if normal points down)
                            # If dip_rad is from arccos(n[2]), then n[2]~1 means dip~0.
        strike_deg = 0.0 # Strike undefined, conventionally 0
        # Rake for horizontal fault (dip=0): angle of slip vector in horizontal plane
        # from strike direction (North if strike=0).
        # Rake = atan2(slip_E, slip_N)
        rake_rad = np.arctan2(u[1], u[0])
    else:
        # FPSPACK nd2pl: phi = atan2(-anx, any) where anx,any are N,E components of normal.
        strike_rad = np.arctan2(-n[0], n[1])
        strike_deg = strike_rad * DEG_CONST
        if strike_deg < 0:
            strike_deg += 360.0
        
        # Rake (-180 to 180 deg)
        # FPSPACK nd2pl: alam = atan2(-dz/sin(delta), dx*cos(phi)+dy*sin(phi))
        # where dz is UP-component of slip, dx,dy are N,E components of slip.
        # delta is dip, phi is strike.
        # Our u is (N,E,Down). So dz_fps = -u[2]. dx_fps = u[0]. dy_fps = u[1].
        
        sin_lambda_num = -(-u[2]) / sin_dip # u[2]/sin_dip
        cos_lambda_num = u[0] * np.cos(strike_rad) + u[1] * np.sin(strike_rad)
        rake_rad = np.arctan2(sin_lambda_num, cos_lambda_num)

    rake_deg = rake_rad * DEG_CONST
    
    while rake_deg > 180.0: rake_deg -= 360.0
    while rake_deg < -180.0: rake_deg += 360.0
        
    return strike_deg, dip_deg, rake_deg

def ar2plp_py(am_tensor_aki_richards):
    """
    Python implementation for ar2plp.
    Calculates P/T/B axes, their trend/plunge, and two nodal planes.
    
    Input: am_tensor_aki_richards (3x3 NumPy array, Aki & Richards: N,E,D)
    Returns a 20-element tuple:
    (M0, M1, e_iso, M0b, 
     s1,d1,r1,dd1, s2,d2,r2,dd2,
     tp,pp, tt,pt, tb,pb, 
     eta, ierr)
    """
    ierr = 0
    mt = np.array(am_tensor_aki_richards, dtype=np.float64)

    # 0. Check symmetry (FPSPACK ar2pt does this)
    if not np.allclose(mt, mt.T, atol=FPSCOM_CONSTANTS['tentol']):
        ierr = ierr | 1 # Bitwise OR for errors
        mt = (mt + mt.T) / 2.0 # Force symmetry

    # 1. Eigenvalue decomposition (mimics FPSPACK's avec -> jacobi, eigsrt)
    # NumPy's eig sorts eigenvalues differently (usually ascending or no specific order for value)
    # and eigenvectors are columns.
    try:
        eig_vals_raw, eig_vecs_raw = np.linalg.eig(mt)
    except np.linalg.LinAlgError:
        return (0.,)*19 + (-100,) # Indicate fatal error in eigenvalue decomp

    # Sort eigenvalues as M_T >= M_B >= M_P (algebraically, for P,T,B axis definition)
    # This is standard geophysical convention for principal stresses/strains.
    # eig_vals_raw are the principal moments.
    idx_sort_descending = np.argsort(eig_vals_raw)[::-1]
    principal_moments = eig_vals_raw[idx_sort_descending] # M_T, M_B, M_P
    
    # Corresponding eigenvectors form the T, B, P axes
    t_axis_vec = _normalize_vector(eig_vecs_raw[:, idx_sort_descending[0]])
    b_axis_vec = _normalize_vector(eig_vecs_raw[:, idx_sort_descending[1]])
    p_axis_vec = _normalize_vector(eig_vecs_raw[:, idx_sort_descending[2]])

    # Conventionally orient P and T axes to point into the lower hemisphere (D component > 0)
    # B-axis orientation is then set for a right-handed system (e.g., B = P x T, then ensure D > 0)
    if p_axis_vec[2] < 0.0: p_axis_vec = -p_axis_vec
    if t_axis_vec[2] < 0.0: t_axis_vec = -t_axis_vec
    
    # Ensure B = P x T (or T x P depending on convention) and then adjust B to point down
    # Let's try B = T x P for consistency with some conventions, then check direction
    b_cross = np.cross(t_axis_vec, p_axis_vec)
    if np.dot(b_axis_vec, b_cross) < 0: # If original b_axis is opposite to T x P
        b_axis_vec = -b_axis_vec # Flip it to align with T x P sense
    if b_axis_vec[2] < 0.0: b_axis_vec = -b_axis_vec # Ensure B also points down


    # 2. Calculate Trend and Plunge for P, T, B axes
    trendt, plungt = _vector_to_trend_plunge(t_axis_vec)
    trendb, plungb = _vector_to_trend_plunge(b_axis_vec)
    trendp, plungp = _vector_to_trend_plunge(p_axis_vec)

    # 3. Calculate Nodal Planes from P and T axes
    # (Aki & Richards, 1980, p. 110, eq. 4.85, or FPSPACK's pt2nd)
    # Normal vector n_fault = (t_axis +/- p_axis)/sqrt(2)
    # Slip vector u_fault = (t_axis -/+ p_axis)/sqrt(2)
    
    n_plane1 = _normalize_vector(t_axis_vec + p_axis_vec) # Normal for plane 1
    u_plane1 = _normalize_vector(t_axis_vec - p_axis_vec) # Slip for plane 1
    s1, d1, r1 = _sdr_from_normal_slip(n_plane1, u_plane1)
    
    n_plane2 = _normalize_vector(t_axis_vec - p_axis_vec) # Normal for plane 2
    u_plane2 = _normalize_vector(t_axis_vec + p_axis_vec) # Slip for plane 2
    s2, d2, r2 = _sdr_from_normal_slip(n_plane2, u_plane2)

    dd1 = (s1 + 90.0) % 360.0
    dd2 = (s2 + 90.0) % 360.0

    # 4. Scalar moments and CLVD component (eta) - mimicking FPSPACK's ar2pt logic
    e_iso = np.trace(mt) / 3.0
    
    # Deviatoric eigenvalues (already sorted as T_dev, B_dev, P_dev if principal_moments are T,B,P)
    dev_eig_vals = principal_moments - e_iso
    
    # FPSPACK ar2pt sorts deviatoric eigenvalues by *inverse order of their modulus magnitude*
    # Let these be val_fps[0], val_fps[1], val_fps[2]
    # where |val_fps[0]| >= |val_fps[1]| >= |val_fps[2]|
    abs_dev_eig_vals = np.abs(dev_eig_vals)
    fps_sort_indices = np.argsort(abs_dev_eig_vals)[::-1] 
    
    val_fps_0_dev = dev_eig_vals[fps_sort_indices[0]] # Corresponds to FPSPACK's val(1) after deviatoric & sort
    val_fps_1_dev = dev_eig_vals[fps_sort_indices[1]] # Corresponds to FPSPACK's val(2)
    val_fps_2_dev = dev_eig_vals[fps_sort_indices[2]] # Corresponds to FPSPACK's val(3)

    # In FPSPACK ar2pt, am0 is val_fps_0_dev (can be negative)
    # Then P,T axes are assigned based on sign of this am0.
    # Our P,T axes are already from sorted principal moments.
    # Let's use the principal moments directly for M0, etc.
    # M0 = (M_T_dev - M_P_dev) / 2 for pure DC (Jost & Herrmann, 1989)
    # Or, M0 = max_deviatoric_eigenvalue (FPSPACK's am0 = val(1))
    
    m0_fps_signed_from_dev = val_fps_0_dev
    
    if abs(m0_fps_signed_from_dev) < 1e-30: 
        eta_clvd = 0.0
        if abs(val_fps_2_dev) > 1e-15 : ierr = ierr | 8 
    else:
        eta_clvd = -val_fps_2_dev / (2.0 * m0_fps_signed_from_dev)
    
    m0_out = abs(m0_fps_signed_from_dev) # Scalar Seismic Moment M0 (always positive by this def)
    m1_fps = abs(val_fps_2_dev)          # Secondary DC moment
    m0b_fps = (abs(val_fps_0_dev) + abs(val_fps_1_dev)) / 2.0 # Best DC moment

    # FPSPACK's error codes are specific. ierr=1 for non-symmetric.
    # Other errors (5-10) in ar2plp are internal from ca2ax, pt2nd, nd2pl.
    # We don't have direct equivalents for those internal errors here.

    return (m0_out, m1_fps, e_iso, m0b_fps,
            s1, d1, r1, dd1,
            s2, d2, r2, dd2,
            trendp, plungp, trendt, plungt, trendb, plungb,
            eta_clvd, int(ierr))

# --- Functions for selection and swapping (from previous python_subroutines.py) ---
def fpsset_py():
    pass # Constants are defined at module level

def selection_nodal_py(str1_deg, dip1_deg, str2_deg, dip2_deg, str0_deg, dip0_deg):
    n1 = np.zeros(3); n2 = np.zeros(3); n0 = np.zeros(3)
    s1_rad, d1_rad = str1_deg * RAD_CONST, dip1_deg * RAD_CONST
    s2_rad, d2_rad = str2_deg * RAD_CONST, dip2_deg * RAD_CONST
    s0_rad, d0_rad = str0_deg * RAD_CONST, dip0_deg * RAD_CONST

    # Normal vector components (N, E, D), assuming Z is UP for this comparison
    # as in original selection_nodal.
    n1[0] = -np.sin(d1_rad) * np.sin(s1_rad); n1[1] = np.sin(d1_rad) * np.cos(s1_rad); n1[2] = -np.cos(d1_rad)
    n2[0] = -np.sin(d2_rad) * np.sin(s2_rad); n2[1] = np.sin(d2_rad) * np.cos(s2_rad); n2[2] = -np.cos(d2_rad)
    n0[0] = -np.sin(d0_rad) * np.sin(s0_rad); n0[1] = np.sin(d0_rad) * np.cos(s0_rad); n0[2] = -np.cos(d0_rad)
    
    dot_n1_n0 = np.dot(n1, n0)
    dot_n2_n0 = np.dot(n2, n0)

    if abs(dot_n1_n0) >= abs(dot_n2_n0):
        return 1
    else:
        return 2

def swap_fault_py(str1, dip1, rake1, str2, dip2, rake2, nodal_p):
    if nodal_p == 2:
        return float(str2), float(dip2), float(rake2), float(str1), float(dip1), float(rake1)
    return float(str1), float(dip1), float(rake1), float(str2), float(dip2), float(rake2)


if __name__ == '__main__':
    print("--- Testing ar2plp.py ---")
    fpsset_py() # Call to ensure constants are notionally "set"

    # Example Moment Tensor from GlobalCMT: 2004 Parkfield M6.0
    # Harvard (Up, South, East) components, scaled by 1e18 Nm
    mrr_h = 0.233; mtt_h = 0.101; mpp_h = -0.334; # Mtt is South-South, Mpp is East-East
    mrt_h = 0.830; mrp_h = -0.395; mtp_h = 0.061; # Mrt is Up-South, Mrp is Up-East, Mtp is South-East
    
    # Convert Harvard (r, theta(S), phi(E)) to Aki-Richards (N,E,D) for ar2plp_py input
    # M_NN = M_SS (Mtt)
    # M_EE = M_EE (Mpp)
    # M_DD = M_RR (Mrr)
    # M_NE = -M_SE (-Mtp)
    # M_ND = -M_SR (-Mrt)  <-- Note: Mrt is Up-South, so M_SR = Mrt. M_ND = -Mrt
    # M_ED =  M_ER (Mrp)   <-- Note: Mrp is Up-East, so M_ER = Mrp. M_ED = Mrp
    
    mt_test_ned = np.array([
        [mtt_h,  -mtp_h, -mrt_h],  # M_NN, M_NE, M_ND
        [-mtp_h, mpp_h,   mrp_h],  # M_EN, M_EE, M_ED
        [-mrt_h, mrp_h,   mrr_h]   # M_DN, M_DE, M_DD
    ]) * 1e18 # Scale factor

    print("Input Moment Tensor (Aki & Richards N,E,D convention):\n", mt_test_ned)
    
    results = ar2plp_py(mt_test_ned)
    
    print("\nPython ar2plp_py Results (20 outputs):")
    labels = [
        "M0 (DC)", "M1 (secondary DC)", "e_iso (isotropic)", "M0b (best DC)",
        "Strike1", "Dip1", "Rake1", "DipDir1",
        "Strike2", "Dip2", "Rake2", "DipDir2",
        "TrendP", "PlungeP", "TrendT", "PlungeT", "TrendB", "PlungeB",
        "Eta (CLVD)", "Error Code"
    ]
    for label, value in zip(labels, results):
        if "Error Code" in label:
            print(f"{label:>20s}: {int(value)}")
        elif "Eta" in label: # Eta is a ratio
             print(f"{label:>20s}: {value:.4f}")
        elif any(x in label for x in ["Strike", "Dip", "Rake", "DipDir", "Trend", "Plunge"]):
            print(f"{label:>20s}: {value:.2f} deg")
        else: # Moments
            print(f"{label:>20s}: {value:.3e}")

    # GlobalCMT solution for Parkfield 2004/09/28:
    # NP1: strike=148, dip=80, rake=-174
    # NP2: strike=241, dip=84, rake=-10
    # P-axis: plunge=9, trend=102
    # T-axis: plunge=5, trend=194
    # B-axis: plunge=80, trend=8
    # Scalar Moment = 3.5e+17 Nm (M0 from ar2plp should be ~ this / 2 if it's deviatoric M0)
    # FPSPACK M0 is often the largest deviatoric eigenvalue.
    # The M0 from ar2plp_py (abs(val_fps_0_dev)) should be compared carefully.
    # GlobalCMT M0 is total scalar moment, often (M_T_dev - M_P_dev)/2 or sqrt(sum(Mij_dev^2)/2).
    
    print("\nCompare with GlobalCMT for Parkfield 2004/09/28 (approx):")
    print("NP1: s=148, d=80, r=-174 | NP2: s=241, d=84, r=-10")
    print("P-axis: T=102, P=9 | T-axis: T=194, P=5 | B-axis: T=8, P=80")
    
    # Test selection and swap
    s1,d1,r1,dd1 = results[4:8]
    s2,d2,r2,dd2 = results[8:12]
    ref_s, ref_d = 148, 80 # Example reference plane (from GCMT NP1)

    nodal_p = selection_nodal_py(s1,d1,s2,d2, ref_s, ref_d)
    print(f"\nSelected nodal plane (1 or 2) based on ref ({ref_s},{ref_d}): {nodal_p}")

    s1_final,d1_final,r1_final, s2_final,d2_final,r2_final = swap_fault_py(s1,d1,r1, s2,d2,r2, nodal_p)
    print("Final Plane 1 (after swap if needed):")
    print(f"  Strike={s1_final:.2f}, Dip={d1_final:.2f}, Rake={r1_final:.2f}")
    print("Final Plane 2 (after swap if needed):")
    print(f"  Strike={s2_final:.2f}, Dip={d2_final:.2f}, Rake={r2_final:.2f}")

