#!/usr/bin/env python3
"""
GAIA GRAIN SCANNER (CORRECTED)
==============================
Proper halo selection:
- v_rad > 250 km/s (excludes thick disk contamination)
- parallax 0.33-1.0 mas (1-3 kpc shell)
- Sorted by sample purity, not R²
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astroquery.gaia import Gaia
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

def measure_grain(velocities, v_min, v_max):
    """Measure power-law slope in velocity window."""
    mask = (velocities >= v_min) & (velocities <= v_max)
    v_window = velocities[mask]
    
    if len(v_window) < 30:
        return None, None, None, len(v_window)
    
    hist, bins = np.histogram(v_window, bins=25, density=True)
    centers = (bins[:-1] + bins[1:]) / 2
    
    valid = hist > 0
    if np.sum(valid) < 5:
        return None, None, None, len(v_window)
    
    log_v = np.log10(centers[valid])
    log_p = np.log10(hist[valid])
    
    slope, intercept, r_val, _, std_err = stats.linregress(log_v, log_p)
    gamma = abs(slope)
    r2 = r_val**2
    
    return gamma, r2, std_err, len(v_window)

def main():
    print("=" * 60)
    print("GAIA DR3 GRAIN SCANNER (HALO-PURE)")
    print("=" * 60)
    
    print("\n1. Querying Gaia Archive (HALO ONLY)...")
    print("   Cuts: |v_rad| > 250 km/s, distance 1-3 kpc")
    
    # CORRECTED QUERY:
    # - v_rad > 250 km/s excludes thick disk
    # - parallax 0.33-1.0 mas = distance 1-3 kpc (local halo shell)
    query = """
    SELECT TOP 10000
        radial_velocity, radial_velocity_error, parallax
    FROM gaiadr3.gaia_source
    WHERE 
        radial_velocity IS NOT NULL
        AND abs(radial_velocity) > 250
        AND parallax BETWEEN 0.33 AND 1.0
        AND radial_velocity_error < 15
    """
    
    try:
        job = Gaia.launch_job(query)
        results = job.get_results()
        print(f"   -> Downloaded {len(results)} PURE HALO stars")
    except Exception as e:
        print(f"   x Connection Failed: {e}")
        return
    
    velocities = np.abs(results['radial_velocity'].data)
    
    # Scan windows starting at 250 (no disk contamination)
    v_mins = [250, 260, 270, 280, 290, 300, 320]
    v_maxs = [400, 450, 500, 550]
    
    print(f"\n2. Scanning {len(v_mins) * len(v_maxs)} velocity windows...")
    
    results_list = []
    
    for v_min in v_mins:
        for v_max in v_maxs:
            if v_max <= v_min + 50:
                continue
                
            gamma, r2, std_err, n_stars = measure_grain(velocities, v_min, v_max)
            
            if gamma is not None and r2 is not None and r2 > 0.85:
                results_list.append({
                    'v_min': v_min,
                    'v_max': v_max,
                    'gamma': gamma,
                    'r2': r2,
                    'std_err': std_err,
                    'n_stars': n_stars
                })
    
    df = pd.DataFrame(results_list)
    
    if len(df) == 0:
        print("   x No valid measurements. Sample may be too small.")
        return
    
    # Sort by v_min (purity) then by n_stars (statistics)
    df = df.sort_values(['v_min', 'n_stars'], ascending=[True, False])
    
    df.to_csv('gaia_grain_scan_results.csv', index=False)
    print(f"   -> Saved {len(df)} measurements to gaia_grain_scan_results.csv")
    
    # Statistics
    print("\n" + "=" * 60)
    print("SCAN RESULTS (HALO-PURE SAMPLE)")
    print("=" * 60)
    
    print(f"\nTotal valid measurements: {len(df)}")
    print(f"Gamma range: {df['gamma'].min():.3f} - {df['gamma'].max():.3f}")
    print(f"Mean gamma: {df['gamma'].mean():.3f} ± {df['gamma'].std():.3f}")
    print(f"NFW prediction: 5.33")
    
    # High-purity subset (v_min >= 270)
    pure = df[df['v_min'] >= 270]
    if len(pure) > 0:
        print(f"\nHIGH-PURITY SUBSET (v_min >= 270 km/s):")
        print(f"  Mean gamma: {pure['gamma'].mean():.3f} ± {pure['gamma'].std():.3f}")
        print(f"  N measurements: {len(pure)}")
    
    # Top 10 by v_min (purest samples first)
    print("\n--- TOP 10 BY SAMPLE PURITY ---")
    print(f"{'V-Window':<12} {'Gamma':<8} {'R²':<8} {'N Stars':<10}")
    print("-" * 40)
    for _, row in df.head(10).iterrows():
        print(f"{row['v_min']:.0f}-{row['v_max']:.0f} km/s   {row['gamma']:.3f}    {row['r2']:.3f}    {row['n_stars']:.0f}")
    
    # Deviation from NFW
    df['nfw_deviation'] = (df['gamma'] - 5.33) / 5.33 * 100
    
    print(f"\n--- DEVIATION FROM NFW (γ = 5.33) ---")
    print(f"Mean deviation: +{df['nfw_deviation'].mean():.1f}%")
    print(f"Measurements above NFW: {(df['gamma'] > 5.33).sum()}/{len(df)}")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Gamma distribution
    ax1 = axes[0, 0]
    ax1.hist(df['gamma'], bins=15, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(5.33, color='red', linestyle='--', linewidth=2, label='NFW = 5.33')
    ax1.axvline(df['gamma'].mean(), color='green', linestyle='-', linewidth=2, label=f'Mean = {df["gamma"].mean():.2f}')
    ax1.set_xlabel('Gamma (γ)')
    ax1.set_ylabel('Count')
    ax1.set_title('Distribution of Measured γ (Halo-Pure)')
    ax1.legend()
    
    # 2. Gamma vs v_min
    ax2 = axes[0, 1]
    scatter = ax2.scatter(df['v_min'], df['gamma'], c=df['r2'], cmap='viridis', s=60)
    ax2.axhline(5.33, color='red', linestyle='--', label='NFW = 5.33')
    ax2.set_xlabel('V_min (km/s)')
    ax2.set_ylabel('Gamma (γ)')
    ax2.set_title('γ vs Minimum Velocity (Higher = Purer)')
    plt.colorbar(scatter, ax=ax2, label='R²')
    
    # 3. R² distribution
    ax3 = axes[1, 0]
    ax3.hist(df['r2'], bins=15, color='coral', edgecolor='black', alpha=0.7)
    ax3.axvline(0.9, color='green', linestyle='--', linewidth=2, label='R² = 0.9')
    ax3.set_xlabel('R² (Fit Quality)')
    ax3.set_ylabel('Count')
    ax3.set_title('Distribution of Fit Quality')
    ax3.legend()
    
    # 4. Gamma vs window width
    df['width'] = df['v_max'] - df['v_min']
    ax4 = axes[1, 1]
    scatter2 = ax4.scatter(df['width'], df['gamma'], c=df['n_stars'], cmap='plasma', s=60)
    ax4.axhline(5.33, color='red', linestyle='--', label='NFW = 5.33')
    ax4.set_xlabel('Window Width (km/s)')
    ax4.set_ylabel('Gamma (γ)')
    ax4.set_title('γ vs Window Width')
    plt.colorbar(scatter2, ax=ax4, label='N Stars')
    
    plt.suptitle(f'Gaia DR3 HALO-PURE Scan: {len(df)} Measurements\nMean γ = {df["gamma"].mean():.2f} ± {df["gamma"].std():.2f} | NFW = 5.33 | Distance: 1-3 kpc', 
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig('gaia_grain_scan_results.png', dpi=150)
    print("\n   -> Saved plot to gaia_grain_scan_results.png")
    
    print("\n" + "=" * 60)
    print("SCAN COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    main()