#!/usr/bin/env python3
"""
Combination-based calibration script (paper-faithful).

Expect 6 CSV files named:
    idle1_M0.csv, idle2_M0.csv, ..., idle6_M0.csv

Each CSV must contain columns: Time, Ax, Ay, Az, Gx, Gy, Gz
We only use Ax, Ay, Az (mean per file).

Orientation order (row i of V and G):
 1 -> -Y   : (0, -g, 0)
 2 -> +X   : ( g,  0, 0)
 3 -> +Y   : (0,  g, 0)
 4 -> -X   : (-g, 0, 0)
 5 -> +Z   : (0,  0, g)
 6 -> -Z   : (0,  0, -g)

This script follows Eq. (5)-(7) from the paper and uses all combinations:
 Z: choose 2 (6C2 = 15)
 Y: choose 3 (6C3 = 20)
 X: choose 4 (6C4 = 15)

Outputs final K_a, M_a, b_a and diagnostics (means, stds, residuals).
"""

import os
import itertools
import numpy as np
import pandas as pd

# -----------------------------
# User parameters
# -----------------------------
BASE_PATH = "/home/kjw2kor/shared_folder/Calibration Metrices/"                # path where idle*.csv are located
FILE_PATTERN = "idle{}_M0.csv" # file name pattern, i = 1..6
N_ORIENT = 6
G_MAG = 9.8                    # gravity magnitude (m/s^2). Change to 9.80665 if desired.

# Tolerance for considering a combo singular or unstable
SINGULAR_TOL = 1e-9
K_ZERO_TOL = 1e-9

# -----------------------------
# Helper: load and compute V
# -----------------------------
def compute_V_matrix(base_path=BASE_PATH, pattern=FILE_PATTERN, n=6):
    rows = []
    for i in range(1, n+1):
        fname = os.path.join(base_path, pattern.format(i))
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"Expected file not found: {fname}")
        df = pd.read_csv(fname)
        # check columns exist
        for c in ("Ax","Ay","Az"):
            if c not in df.columns:
                raise KeyError(f"Column '{c}' not found in {fname}")
        mean_ax = float(df["Ax"].mean())
        mean_ay = float(df["Ay"].mean())
        mean_az = float(df["Az"].mean())
        rows.append([mean_ax, mean_ay, mean_az])
        print(f"Loaded {fname}: mean Ax,Ay,Az = [{mean_ax:.6g}, {mean_ay:.6g}, {mean_az:.6g}]")
    V = np.array(rows, dtype=float)  # shape (6,3)
    return V

# -----------------------------
# Helper: G matrix in requested order
# -----------------------------
def get_G_matrix(g=G_MAG):
    # order: -Y, +X, +Y, -X, +Z, -Z
    return np.array([
        [0.0, -g,  0.0],   # 1: -Y
        [ g,  0.0, 0.0],   # 2: +X
        [0.0,  g,  0.0],   # 3: +Y
        [-g, 0.0, 0.0],    # 4: -X
        [0.0, 0.0,  g],    # 5: +Z
        [0.0, 0.0, -g],    # 6: -Z
    ], dtype=float)

# -----------------------------
# Solve small linear system safely
# -----------------------------
def solve_linear(A, y):
    """
    Solve A x = y for x.
    Returns (x, ok_flag). ok_flag False if singular/ill-conditioned.
    Uses numpy.linalg.solve when square and well-conditioned; otherwise least squares.
    """
    A = np.asarray(A, dtype=float)
    y = np.asarray(y, dtype=float)
    m, n = A.shape
    try:
        if m == n:
            # try direct solve
            # check condition (avoid singular)
            cond = np.linalg.cond(A)
            if cond < 1.0 / SINGULAR_TOL:
                x = np.linalg.solve(A, y)
                return x, True
            else:
                # ill-conditioned, fallback to lstsq
                x, *_ = np.linalg.lstsq(A, y, rcond=None)
                return x, True
        else:
            # over/under determined: use least-squares
            x, *_ = np.linalg.lstsq(A, y, rcond=None)
            return x, True
    except Exception as e:
        return None, False

# -----------------------------
# Combination-based solvers
# -----------------------------
def combos_for_axis(V, G):
    """
    For a given V (6x3) and G (6x3) compute parameter estimates for each axis
    using all valid combinations. Returns dictionaries of per-combination results and stats.
    """
    n = V.shape[0]
    if n != 6:
        raise ValueError("This routine expects exactly 6 orientations (rows).")

    results = {
        "Z": [],  # list of dicts { 'kz':..., 'bz':..., 'combo':(i,j) }
        "Y": [],  # list of dicts { 'ky':..., 'by':..., 'alpha_zx':..., 'combo':(i,j,k) }
        "X": [],  # list of dicts { 'kx':..., 'bx':..., 'alpha_yz':..., 'alpha_zy':..., 'combo':(i,j,k,l) }
        "skipped": {"Z":0, "Y":0, "X":0}
    }

    # --- Z axis: choose 2 (6C2) ---
    for (i,j) in itertools.combinations(range(n), 2):
        A = np.column_stack([np.ones(2), np.array([G[i,2], G[j,2]])]).reshape(2,2)
        y = np.array([V[i,2], V[j,2]])
        x, ok = solve_linear(A, y)
        if not ok:
            results["skipped"]["Z"] += 1
            continue
        bz, kz = x[0], x[1]
        results["Z"].append({"kz": float(kz), "bz": float(bz), "combo": (i+1,j+1)})
    # --- Y axis: choose 3 (6C3) ---
    for (i,j,k) in itertools.combinations(range(n), 3):
        A = np.column_stack([
            np.ones(3),
            np.array([G[i,1], G[j,1], G[k,1]]),  # gy
            np.array([G[i,0], G[j,0], G[k,0]])   # gx
        ])
        y = np.array([V[i,1], V[j,1], V[k,1]])
        x, ok = solve_linear(A, y)
        if not ok:
            results["skipped"]["Y"] += 1
            continue
        by = float(x[0])
        ky = float(x[1])
        ky_alphazx = float(x[2])
        if abs(ky) < K_ZERO_TOL:
            # unstable; skip
            results["skipped"]["Y"] += 1
            continue
        alphazx = ky_alphazx / ky
        results["Y"].append({"ky": ky, "by": by, "alpha_zx": alphazx, "combo": (i+1,j+1,k+1)})
    # --- X axis: choose 4 (6C4) ---
    for (i,j,k,l) in itertools.combinations(range(n), 4):
        A = np.column_stack([
            np.ones(4),
            np.array([G[i,0], G[j,0], G[k,0], G[l,0]]),  # gx
            np.array([G[i,2], G[j,2], G[k,2], G[l,2]]),  # gz
            np.array([G[i,1], G[j,1], G[k,1], G[l,1]])   # gy
        ])
        y = np.array([V[i,0], V[j,0], V[k,0], V[l,0]])
        x, ok = solve_linear(A, y)
        if not ok:
            results["skipped"]["X"] += 1
            continue
        bx = float(x[0])
        kx = float(x[1])
        kx_alphayz = float(x[2])
        kx_alphazy = float(x[3])
        if abs(kx) < K_ZERO_TOL:
            results["skipped"]["X"] += 1
            continue
        alphayz = kx_alphayz / kx
        alphazy = kx_alphazy / kx
        results["X"].append({
            "kx": kx, "bx": bx,
            "alpha_yz": alphayz, "alpha_zy": alphazy,
            "combo": (i+1,j+1,k+1,l+1)
        })

    return results

# -----------------------------
# Aggregation: mean and std
# -----------------------------
def aggregate_results(results):
    """
    Compute mean and std for each parameter from combination results.
    Returns dict with aggregated parameter stats and counts.
    """
    agg = {}
    # Z
    kz_list = [r["kz"] for r in results["Z"]]
    bz_list = [r["bz"] for r in results["Z"]]
    agg["Z_count"] = len(kz_list)
    agg["kz_mean"] = np.mean(kz_list) if kz_list else None
    agg["kz_std"] = np.std(kz_list, ddof=1) if len(kz_list) > 1 else 0.0
    agg["bz_mean"] = np.mean(bz_list) if bz_list else None
    agg["bz_std"] = np.std(bz_list, ddof=1) if len(bz_list) > 1 else 0.0

    # Y
    ky_list = [r["ky"] for r in results["Y"]]
    by_list = [r["by"] for r in results["Y"]]
    azx_list = [r["alpha_zx"] for r in results["Y"]]
    agg["Y_count"] = len(ky_list)
    agg["ky_mean"] = np.mean(ky_list) if ky_list else None
    agg["ky_std"] = np.std(ky_list, ddof=1) if len(ky_list) > 1 else 0.0
    agg["by_mean"] = np.mean(by_list) if by_list else None
    agg["by_std"] = np.std(by_list, ddof=1) if len(by_list) > 1 else 0.0
    agg["alpha_zx_mean"] = np.mean(azx_list) if azx_list else None
    agg["alpha_zx_std"] = np.std(azx_list, ddof=1) if len(azx_list) > 1 else 0.0

    # X
    kx_list = [r["kx"] for r in results["X"]]
    bx_list = [r["bx"] for r in results["X"]]
    ayz_list = [r["alpha_yz"] for r in results["X"]]
    azy_list = [r["alpha_zy"] for r in results["X"]]
    agg["X_count"] = len(kx_list)
    agg["kx_mean"] = np.mean(kx_list) if kx_list else None
    agg["kx_std"] = np.std(kx_list, ddof=1) if len(kx_list) > 1 else 0.0
    agg["bx_mean"] = np.mean(bx_list) if bx_list else None
    agg["bx_std"] = np.std(bx_list, ddof=1) if len(bx_list) > 1 else 0.0
    agg["alpha_yz_mean"] = np.mean(ayz_list) if ayz_list else None
    agg["alpha_yz_std"] = np.std(ayz_list, ddof=1) if len(ayz_list) > 1 else 0.0
    agg["alpha_zy_mean"] = np.mean(azy_list) if azy_list else None
    agg["alpha_zy_std"] = np.std(azy_list, ddof=1) if len(azy_list) > 1 else 0.0

    # skipped counts
    agg["skipped"] = results["skipped"]

    return agg

# -----------------------------
# Build final matrices from aggregated means
# -----------------------------
def build_final_matrices(agg):
    # check values exist
    for key in ("kx_mean","ky_mean","kz_mean","bx_mean","by_mean","bz_mean",
                "alpha_yz_mean","alpha_zy_mean","alpha_zx_mean"):
        if key not in agg or agg[key] is None:
            raise RuntimeError(f"Not enough valid combinations to compute {key}. Aborting.")

    K = np.diag([agg["kx_mean"], agg["ky_mean"], agg["kz_mean"]])
    M = np.array([
        [1.0,           0.0,            0.0],
        [0.0,           1.0,    agg["alpha_yz_mean"]],
        [agg["alpha_zx_mean"], agg["alpha_zy_mean"], 1.0]
    ])
    b = np.array([agg["bx_mean"], agg["by_mean"], agg["bz_mean"]])
    S = K.dot(M)
    return K, M, b, S

# -----------------------------
# Diagnostics: residuals
# -----------------------------
def compute_residuals(V, G, K, M, b):
    """
    Compute residuals per orientation: r_i = v_i - (K M g_i + b)
    Returns residual matrix (6x3) and RMS.
    """
    S = K.dot(M)                       # 3x3
    V_est = G.dot(S.T) + np.ones((G.shape[0],1)).dot(b.reshape(1,3))
    res = V - V_est
    rms = np.sqrt(np.mean(res**2))
    return res, rms

# -----------------------------
# Main orchestrator
# -----------------------------
def main():
    print("=== Combination-based calibration (paper-faithful) ===")
    V = compute_V_matrix()
    G = get_G_matrix()

    print("\nV (means) matrix:\n", V)
    print("\nG matrix (order -Y, +X, +Y, -X, +Z, -Z):\n", G)

    print("\n--- Solving all combinations per axis ---")
    results = combos_for_axis(V, G)
    agg = aggregate_results(results)

    # print stats
    print("\n--- Combination results summary ---")
    print(f"Z-axis: {agg['Z_count']} valid combos, skipped {results['skipped']['Z']}")
    print(f"  k_z mean = {agg['kz_mean']:.6g}  std = {agg['kz_std']:.6g}")
    print(f"  b_z mean = {agg['bz_mean']:.6g}  std = {agg['bz_std']:.6g}")

    print(f"\nY-axis: {agg['Y_count']} valid combos, skipped {results['skipped']['Y']}")
    print(f"  k_y mean = {agg['ky_mean']:.6g}  std = {agg['ky_std']:.6g}")
    print(f"  b_y mean = {agg['by_mean']:.6g}  std = {agg['by_std']:.6g}")
    print(f"  alpha_zx mean = {agg['alpha_zx_mean']:.6g}  std = {agg['alpha_zx_std']:.6g}")

    print(f"\nX-axis: {agg['X_count']} valid combos, skipped {results['skipped']['X']}")
    print(f"  k_x mean = {agg['kx_mean']:.6g}  std = {agg['kx_std']:.6g}")
    print(f"  b_x mean = {agg['bx_mean']:.6g}  std = {agg['bx_std']:.6g}")
    print(f"  alpha_yz mean = {agg['alpha_yz_mean']:.6g}  std = {agg['alpha_yz_std']:.6g}")
    print(f"  alpha_zy mean = {agg['alpha_zy_mean']:.6g}  std = {agg['alpha_zy_std']:.6g}")

    # build final matrices
    K, M, b, S = build_final_matrices(agg)

    print("\n=== Final calibration matrices ===")
    np.set_printoptions(precision=8, suppress=True)
    print("K_a (scale factors):\n", K)
    print("M_a (misalignment):\n", M)
    print("b_a (bias vector):\n", b)
    print("S = K_a * M_a:\n", S)

    # residuals
    res, rms = compute_residuals(V, G, K, M, b)
    print("\nResiduals (V - (K M G + b)) per orientation (rows = orientations):\n", res)
    print(f"RMS residual (all axes, all orientations) = {rms:.6g}")

    # Optionally return results
    return {"V":V, "G":G, "results":results, "agg":agg, "K":K, "M":M, "b":b, "S":S, "res":res, "rms":rms}

if __name__ == "__main__":
    out = main()
