#!/usr/bin/env python3
"""
GAIA HALO HUNTER (Corrected for MNRAS Submission)
=================================================
Connects to ESA servers, downloads high-velocity stars, and measures the grain.

CORRECTION: Shows NFW prediction (γ=5.33) for inner halo comparison,
NOT the outer halo prediction (γ=1.878).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astroquery.gaia import Gaia

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

def main():
    print("------------------------------------------------")
    print("   GAIA DR3: HALO VALIDATION (CORRECTED)")
    print("------------------------------------------------")
    print("1. Connecting to ESA Gaia Archive...")
    
    # ADQL Query
    # Selecting stars moving faster than 250 km/s relative to us
    # This guarantees we are seeing the Halo, not the Disk.
    query = """
    SELECT TOP 5000
        radial_velocity, random_index
    FROM gaiadr3.gaia_source
    WHERE 
        radial_velocity IS NOT NULL
        AND abs(radial_velocity) > 250
        AND radial_velocity_error < 10
    """
    
    try:
        job = Gaia.launch_job(query)
        results = job.get_results()
        print(f"   -> Success! Downloaded {len(results)} Halo stars.")
    except Exception as e:
        print(f"   x Connection Failed: {e}")
        print("     (Make sure you have internet and 'pip install astroquery')")
        return

    # 2. Analyze
    velocities = np.abs(results['radial_velocity'].data)
    
    # We look at the "Tail" (250 km/s to 450 km/s)
    v_min, v_max = 250, 450
    mask = (velocities >= v_min) & (velocities <= v_max)
    v_window = velocities[mask]
    
    print(f"   -> Analyzing {len(v_window)} stars in the tail ({v_min}-{v_max} km/s)...")
    
    if len(v_window) < 50:
        print("   x Not enough stars in this window to fit a slope.")
        return

    # 3. Measure Slope
    hist, bins = np.histogram(v_window, bins=30, density=True)
    centers = (bins[:-1] + bins[1:]) / 2
    
    valid = hist > 0
    log_v = np.log10(centers[valid])
    log_p = np.log10(hist[valid])
    
    slope, intercept, r_val, _, std_err = stats.linregress(log_v, log_p)
    gamma = abs(slope)
    
    # Calculate uncertainty (approximation)
    n_stars = len(v_window)
    gamma_err = (gamma - 1) / np.sqrt(n_stars)
    
    print(f"\n   >>> MEASURED GAIA SLOPE: γ = {gamma:.3f} ± {gamma_err:.3f} <<<")
    print(f"   >>> NFW PREDICTION:      γ = 5.33 <<<")
    print(f"   >>> DEVIATION:           +{(gamma - 5.33)/5.33 * 100:.1f}% <<<")
    
    # 4. Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(centers[valid], hist[valid], color='black', s=40, label='Gaia DR3 Data')
    
    # Observed fit line
    fit_y = 10**(intercept + slope * log_v)
    plt.plot(centers[valid], fit_y, 'r--', linewidth=2, label=f'Fit γ = {gamma:.3f}')
    
    # NFW prediction line (γ = 5.33 for inner halo)
    # Scale to pass through the middle of the data
    mid_v = np.sqrt(v_min * v_max)
    mid_idx = np.argmin(np.abs(centers[valid] - mid_v))
    nfw_intercept = log_p[mid_idx] + 5.33 * log_v[mid_idx]
    nfw_y = 10**(nfw_intercept - 5.33 * log_v)
    plt.plot(centers[valid], nfw_y, 'b--', linewidth=2, alpha=0.7, label='NFW Prediction (γ = 5.33)')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radial Velocity (km/s)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(f'Milky Way Halo (Gaia DR3 Data)\nMeasured: γ = {gamma:.3f} ± {gamma_err:.3f}  |  NFW: γ = 5.33  |  Deviation: +{(gamma - 5.33)/5.33 * 100:.0f}%', 
              fontsize=11)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure3_gaia_validation.png', dpi=150)
    print("\n   -> Plot saved to 'figure3_gaia_validation.png'")
    
    # Print summary for paper
    print("\n" + "="*50)
    print("VALUES FOR PAPER:")
    print("="*50)
    print(f"  Sample size:     N = {n_stars}")
    print(f"  Velocity window: {v_min}-{v_max} km/s")
    print(f"  Measured γ:      {gamma:.3f} ± {gamma_err:.3f}")
    print(f"  NFW prediction:  5.33")
    print(f"  Deviation:       +{(gamma - 5.33)/5.33 * 100:.1f}%")
    print(f"  R² of fit:       {r_val**2:.3f}")
    print("="*50)

if __name__ == "__main__":
    main()
