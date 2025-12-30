#!/usr/bin/env python3
"""
FIRE-2 ANISOTROPY CHECKER
=========================
Determines the orbital structure of the Dark Matter Halo to identify
the mechanism behind the 1.878 velocity slope.

The Critical Test:
- Is the halo Radial (Plunging)? (Beta > 0)
- Is the halo Isotropic (Random)? (Beta = 0)
- Is the halo Tangential (Shells)? (Beta < 0) <-- PREDICTION: Beta approx -0.3

Method:
1. Load Dark Matter (PartType1) in the 35-50 kpc shell.
2. Decompose velocity vectors into Radial (v_r) and Tangential (v_t) components.
3. Calculate Velocity Dispersion Tensor.
4. Compute Anisotropy Parameter Beta = 1 - (sigma_t^2 / 2*sigma_r^2).
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os

# ============================================================================
# CONFIGURATION
# ============================================================================

# Center of the Galaxy (Updated from your previous logs)
COM = np.array([29335.6, 30990.8, 32488.4])

# Target Shell
R_MIN, R_MAX = 35, 50

# Velocity Window for slope check (optional, but good for context)
V_MIN, V_MAX = 40, 130

# ============================================================================
# LOGIC
# ============================================================================

def load_and_process_shell(filenames):
    print(f"1. Loading FIRE-2 Dark Matter (Shell {R_MIN}-{R_MAX} kpc)...")
    
    v_rad_list = []
    v_tan_list = []
    v_tot_list = []
    
    total_loaded = 0
    
    for fname in filenames:
        if os.path.exists(fname):
            try:
                with h5py.File(fname, 'r') as f:
                    # Load Positions & Velocities
                    pos = f['PartType1/Coordinates'][:]
                    vel = f['PartType1/Velocities'][:]
                    
                    # 1. Center the Halo
                    r_vec = pos - COM
                    dist = np.sqrt(np.sum(r_vec**2, axis=1))
                    
                    # 2. Filter for Shell
                    mask = (dist >= R_MIN) & (dist <= R_MAX)
                    n_shell = np.sum(mask)
                    
                    if n_shell > 0:
                        r_shell = r_vec[mask]
                        v_shell = vel[mask]
                        d_shell = dist[mask]
                        
                        # 3. Vector Decomposition
                        # Radial Unit Vector: r_hat = r / |r|
                        # Avoid div by zero (dist is constrained > 35, so safe)
                        r_hat = r_shell / d_shell[:, np.newaxis]
                        
                        # Radial Velocity: v_r = v dot r_hat
                        # (N,3) * (N,3) -> sum -> (N,)
                        vr = np.sum(v_shell * r_hat, axis=1)
                        
                        # Tangential Velocity Vector: v_t_vec = v - vr * r_hat
                        vt_vec = v_shell - vr[:, np.newaxis] * r_hat
                        
                        # Tangential Speed: |v_t|
                        vt = np.sqrt(np.sum(vt_vec**2, axis=1))
                        
                        # Total Speed
                        vtot = np.sqrt(np.sum(v_shell**2, axis=1))
                        
                        v_rad_list.append(vr)
                        v_tan_list.append(vt)
                        v_tot_list.append(vtot)
                        
                        total_loaded += n_shell
                        # print(f"   -> Found {n_shell} particles in {fname}")
                        
            except Exception as e:
                print(f"   x Error {fname}: {e}")

    if total_loaded == 0:
        return None, None, None

    return np.concatenate(v_rad_list), np.concatenate(v_tan_list), np.concatenate(v_tot_list)

def main():
    # Detect Files
    files = [f for f in os.listdir('.') if f.startswith('snapshot') and f.endswith('.hdf5')]
    if not files:
        print("No snapshot .hdf5 files found!")
        return

    # 1. Get Data
    vr, vt, vtot = load_and_process_shell(files)
    if vr is None: return
    
    print(f"\n2. Statistics for {len(vr):,} Particles")
    
    # 2. Calculate Dispersions (Sigma)
    # Sigma^2 = Mean(v^2) - Mean(v)^2
    # For Radial, mean might be non-zero (infall/outflow), but for anisotropy usually we take variance about mean.
    
    sigma_r = np.std(vr)
    sigma_t = np.sqrt(np.mean(vt**2) / 2) # Convention: sigma_t is 1D dispersion. vt^2 has 2 components (theta, phi).
                                          # Usually beta = 1 - (sigma_theta^2 + sigma_phi^2) / (2 sigma_r^2)
                                          # vt^2 = v_theta^2 + v_phi^2.
                                          # So sigma_tan_total^2 = Mean(vt^2).
                                          # Beta = 1 - Mean(vt^2) / (2 * sigma_r^2)
    
    sigma_tan_sq = np.mean(vt**2)
    sigma_rad_sq = np.var(vr) # Variance matches sigma_r^2
    
    beta = 1.0 - (sigma_tan_sq / (2.0 * sigma_rad_sq))
    
    print("-" * 40)
    print(f"   Radial Dispersion (sigma_r):     {np.sqrt(sigma_rad_sq):.2f} km/s")
    print(f"   Tangential Dispersion (sigma_t): {np.sqrt(sigma_tan_sq/2):.2f} km/s (1D equiv)")
    print(f"   ANISOTROPY PARAMETER (Beta):     {beta:.4f}")
    print("-" * 40)
    
    # 3. Interpret Beta
    print("\nINTERPRETATION:")
    if beta > 0.1:
        print(">> RADIAL BIAS (Plunging Orbits)")
        print("   Status: Contradicts the Shell/1.74 Hypothesis.")
    elif beta < -0.1:
        print(">> TANGENTIAL BIAS (Shell Structure)")
        print("   Status: CONFIRMS the Shell/1.74 Hypothesis!")
    else:
        print(">> ISOTROPIC (Random/Thermal)")
        print("   Status: Requires stiff gravity (n > 2) to explain 1.878.")

    # 4. Plot Distributions
    plt.figure(figsize=(10, 6))
    
    # Normalize histograms
    plt.hist(np.abs(vr), bins=100, density=True, alpha=0.5, color='red', label=f'Radial |Vr| (sig={np.std(vr):.0f})')
    # For Vt, theoretical max is different, visualize the magnitude
    plt.hist(vt, bins=100, density=True, alpha=0.5, color='blue', label=f'Tangential Vt (sig={np.sqrt(np.mean(vt**2)):.0f})')
    
    plt.xlabel('Velocity Component (km/s)')
    plt.ylabel('Probability Density')
    plt.title(f'Orbital Architecture of the Halo (Beta = {beta:.3f})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig('fire2_anisotropy.png')
    print("Saved plot to 'fire2_anisotropy.png'")

if __name__ == "__main__":
    main()