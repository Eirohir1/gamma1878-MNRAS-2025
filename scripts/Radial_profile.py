#!/usr/bin/env python3
"""
FIRE-2 GRAIN PARAMETER SCANNER
==============================

Systematically scans across:
- Velocity windows (v_min, v_max combinations)
- Radial shells (various widths and centers)

Goal: Find optimal parameters and characterize the oscillatory structure

Author: Claude + Gemini methodology
Date: December 25, 2025
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

SNAPSHOT_FILES = [
    'snapshot_600.0.hdf5',
    'snapshot_600.1.hdf5',
    'snapshot_600.2.hdf5',
    'snapshot_600.3.hdf5'
]

# Target grain
GRAIN_TARGET = 1.878

# Velocity windows to test
VELOCITY_WINDOWS = [
    (20, 100),
    (20, 120),
    (20, 150),
    (30, 100),
    (30, 120),
    (30, 150),
    (40, 100),
    (40, 120),
    (40, 130),   # Original
    (40, 150),
    (50, 120),
    (50, 150),
    (60, 140),
]

# Radial shells - focus on regions where oscillations were seen
RADIAL_SHELLS = [
    # Inner disk/halo
    (5, 10), (10, 15), (15, 20), (20, 25), (25, 30),
    (5, 15), (10, 20), (15, 25), (20, 30),
    # Mid halo - strong signal region
    (30, 40), (35, 45), (40, 50), (45, 55), (50, 60),
    (30, 45), (35, 50), (40, 55), (45, 60),
    (30, 50), (40, 60), (35, 55),
    # Outer halo
    (60, 80), (70, 90), (80, 100), (90, 110), (100, 120),
    (60, 90), (70, 100), (80, 110), (90, 120),
    # Transition zone (where you found crossings)
    (100, 130), (110, 140), (120, 150), (130, 160),
    (140, 170), (150, 180), (160, 190), (170, 200),
    # Far field
    (200, 250), (250, 300), (300, 350), (350, 400),
    (200, 300), (250, 350), (300, 400),
    # Fine grid around best matches (40-50, 170-180)
    (38, 48), (42, 52), (45, 55),
    (165, 175), (168, 178), (170, 180), (172, 182), (175, 185),
]

# Minimum particles
MIN_PARTICLES = 50000

# Null test iterations
N_NULL_ITERATIONS = 100

# ============================================================================
# DATA LOADING
# ============================================================================

def load_fire2_data(file_list):
    """Load FIRE-2 data."""
    all_pos = []
    all_vel = []
    
    for filename in file_list:
        try:
            with h5py.File(filename, 'r') as f:
                print(f"  Loading {filename}...")
                
                if 'PartType1' in f:
                    pos = f['PartType1/Coordinates'][:]
                    vel = f['PartType1/Velocities'][:]
                else:
                    continue
                
                all_pos.append(pos)
                all_vel.append(vel)
                print(f"    → {len(pos):,} particles")
                
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    if not all_pos:
        return None, None, None
    
    pos = np.concatenate(all_pos)
    vel = np.concatenate(all_vel)
    center = np.median(pos, axis=0)
    
    return pos, vel, center


# ============================================================================
# GRAIN MEASUREMENT
# ============================================================================

def measure_grain(velocities, v_min, v_max, n_bins=100):
    """Measure power-law exponent."""
    if len(velocities) < 1000:
        return np.nan, np.nan, 0
    
    counts, bin_edges = np.histogram(velocities, bins=n_bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    mask = (bin_centers >= v_min) & (bin_centers <= v_max) & (counts > 0)
    
    if np.sum(mask) < 5:
        return np.nan, np.nan, 0
    
    log_v = np.log10(bin_centers[mask])
    log_p = np.log10(counts[mask])
    
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_v, log_p)
        return abs(slope), r_value**2, np.sum(mask)
    except:
        return np.nan, np.nan, 0


def maxwell_boltzmann_grain(velocities, v_min, v_max, n_iterations=N_NULL_ITERATIONS):
    """Get expected grain from Maxwell-Boltzmann."""
    v_mean = np.mean(velocities)
    a = v_mean / np.sqrt(8 / np.pi)
    
    mb_grains = []
    for _ in range(n_iterations):
        samples = stats.chi.rvs(df=3, scale=a, size=len(velocities))
        grain, _, _ = measure_grain(samples, v_min, v_max)
        if not np.isnan(grain):
            mb_grains.append(grain)
    
    if mb_grains:
        return np.mean(mb_grains), np.std(mb_grains)
    return np.nan, np.nan


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("FIRE-2 GRAIN PARAMETER SCANNER")
    print("=" * 70)
    print()
    print(f"Target grain: {GRAIN_TARGET}")
    print(f"Velocity windows: {len(VELOCITY_WINDOWS)}")
    print(f"Radial shells: {len(RADIAL_SHELLS)}")
    print()
    
    # Load data
    print("Loading FIRE-2 data...")
    pos, vel, center = load_fire2_data(SNAPSHOT_FILES)
    
    if pos is None:
        print("ERROR: Could not load data")
        return
    
    print(f"\nTotal particles: {len(pos):,}")
    print(f"Center: [{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}]")
    
    # Compute distances and velocity magnitudes
    dist = np.sqrt(np.sum((pos - center)**2, axis=1))
    v_mag = np.sqrt(np.sum(vel**2, axis=1))
    
    print()
    
    # ========================================================================
    # SCAN
    # ========================================================================
    
    print("=" * 70)
    print("SCANNING PARAMETER SPACE")
    print("=" * 70)
    print()
    
    results = []
    total = len(VELOCITY_WINDOWS) * len(RADIAL_SHELLS)
    count = 0
    
    for v_min, v_max in VELOCITY_WINDOWS:
        for r_min, r_max in RADIAL_SHELLS:
            count += 1
            
            # Select particles
            mask = (dist >= r_min) & (dist < r_max)
            n_parts = np.sum(mask)
            
            if n_parts < MIN_PARTICLES:
                continue
            
            velocities = v_mag[mask]
            
            # Measure grain
            grain, r_sq, n_bins = measure_grain(velocities, v_min, v_max)
            
            if np.isnan(grain):
                continue
            
            deviation = abs(grain - GRAIN_TARGET) / GRAIN_TARGET * 100
            
            results.append({
                'v_min': v_min,
                'v_max': v_max,
                'r_min': r_min,
                'r_max': r_max,
                'r_mid': (r_min + r_max) / 2,
                'grain': grain,
                'r_squared': r_sq,
                'n_particles': n_parts,
                'deviation': deviation
            })
            
            if count % 50 == 0:
                print(f"  Progress: {count}/{total}")
    
    print(f"\nValid measurements: {len(results)}")
    print()
    
    # ========================================================================
    # TOP MATCHES
    # ========================================================================
    
    print("=" * 70)
    print("TOP 30 CLOSEST MATCHES TO γ = 1.878")
    print("=" * 70)
    print()
    
    results_sorted = sorted(results, key=lambda x: x['deviation'])
    
    print(f"{'Rank':<5} {'V-Window':<12} {'R-Shell':<14} {'Grain':>8} {'Δ%':>7} {'R²':>6} {'N Parts':>12}")
    print("-" * 75)
    
    for i, r in enumerate(results_sorted[:30]):
        print(f"{i+1:<5} {r['v_min']}-{r['v_max']:<7} "
              f"{r['r_min']:>3}-{r['r_max']:<3} kpc    "
              f"{r['grain']:>7.4f} {r['deviation']:>6.2f}% "
              f"{r['r_squared']:>5.3f} {r['n_particles']:>12,}")
    
    print()
    
    # ========================================================================
    # BEST WITH QUALITY CONSTRAINTS
    # ========================================================================
    
    print("=" * 70)
    print("BEST MATCHES WITH QUALITY CONSTRAINTS")
    print("=" * 70)
    print()
    
    for min_rsq in [0.99, 0.98, 0.95, 0.90]:
        filtered = [r for r in results if r['r_squared'] >= min_rsq]
        if filtered:
            best = min(filtered, key=lambda x: x['deviation'])
            print(f"R² ≥ {min_rsq}:")
            print(f"  V-window: {best['v_min']}-{best['v_max']} km/s")
            print(f"  R-shell:  {best['r_min']}-{best['r_max']} kpc")
            print(f"  Grain:    {best['grain']:.4f} (Δ = {best['deviation']:.2f}%)")
            print(f"  R²:       {best['r_squared']:.4f}")
            print(f"  N:        {best['n_particles']:,}")
            print()
    
    # ========================================================================
    # OSCILLATION ANALYSIS
    # ========================================================================
    
    print("=" * 70)
    print("OSCILLATION STRUCTURE ANALYSIS")
    print("=" * 70)
    print()
    
    # Use best velocity window
    best_v = (40, 130)  # Original that showed oscillations
    
    # Get all shells with this v-window
    v_filtered = [r for r in results if r['v_min'] == best_v[0] and r['v_max'] == best_v[1]]
    v_filtered = sorted(v_filtered, key=lambda x: x['r_mid'])
    
    print(f"Using V-window: {best_v[0]}-{best_v[1]} km/s")
    print()
    
    # Find crossings of 1.878
    crossings = []
    for i in range(len(v_filtered) - 1):
        g1 = v_filtered[i]['grain']
        g2 = v_filtered[i+1]['grain']
        r1 = v_filtered[i]['r_mid']
        r2 = v_filtered[i+1]['r_mid']
        
        if (g1 - GRAIN_TARGET) * (g2 - GRAIN_TARGET) < 0:
            frac = (GRAIN_TARGET - g1) / (g2 - g1)
            r_cross = r1 + frac * (r2 - r1)
            crossings.append(r_cross)
    
    print(f"Crossings of γ = 1.878:")
    for i, r in enumerate(crossings):
        print(f"  {i+1}. ~{r:.1f} kpc")
    
    if len(crossings) >= 2:
        spacings = np.diff(crossings)
        print(f"\nSpacing between crossings:")
        for i, s in enumerate(spacings):
            print(f"  {crossings[i]:.1f} → {crossings[i+1]:.1f}: {s:.1f} kpc")
        print(f"\nMean spacing: {np.mean(spacings):.1f} kpc")
        print(f"Std spacing:  {np.std(spacings):.1f} kpc")
    
    print()
    
    # Shells within 5% of target
    close_matches = [r for r in v_filtered if r['deviation'] < 5]
    print(f"Shells within 5% of target: {len(close_matches)}")
    for r in close_matches:
        print(f"  {r['r_min']}-{r['r_max']} kpc: γ = {r['grain']:.4f} (Δ = {r['deviation']:.2f}%)")
    
    print()
    
    # ========================================================================
    # MB COMPARISON FOR TOP MATCH
    # ========================================================================
    
    print("=" * 70)
    print("MAXWELL-BOLTZMANN COMPARISON")
    print("=" * 70)
    print()
    
    best = results_sorted[0]
    mask = (dist >= best['r_min']) & (dist < best['r_max'])
    velocities = v_mag[mask]
    
    mb_mean, mb_std = maxwell_boltzmann_grain(velocities, best['v_min'], best['v_max'])
    
    print(f"Best match: V={best['v_min']}-{best['v_max']}, R={best['r_min']}-{best['r_max']} kpc")
    print(f"  Observed grain:     {best['grain']:.4f}")
    print(f"  Maxwell-Boltzmann:  {mb_mean:.4f} ± {mb_std:.4f}")
    print(f"  Target:             {GRAIN_TARGET}")
    
    if mb_std > 0:
        z_obs_vs_mb = abs(best['grain'] - mb_mean) / mb_std
        print(f"  Deviation from MB:  {z_obs_vs_mb:.1f}σ")
    
    obs_to_target = abs(best['grain'] - GRAIN_TARGET)
    mb_to_target = abs(mb_mean - GRAIN_TARGET)
    
    if obs_to_target < mb_to_target:
        print(f"  ✓ Observed is CLOSER to 1.878 than thermal expectation")
    
    print()
    
    # ========================================================================
    # VISUALIZATIONS
    # ========================================================================
    
    print("Generating visualizations...")
    create_visualizations(results, v_filtered, crossings)
    print("✓ Saved: fire2_grain_scan_results.png")
    
    # Save CSV
    import pandas as pd
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('deviation')
    results_df.to_csv('fire2_grain_scan_results.csv', index=False)
    print("✓ Saved: fire2_grain_scan_results.csv")
    
    print()
    print("=" * 70)
    print("SCAN COMPLETE")
    print("=" * 70)


def create_visualizations(results, v_filtered, crossings):
    """Create visualizations."""
    import pandas as pd
    results_df = pd.DataFrame(results)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Grain vs Radius for best V-window
    ax1 = axes[0, 0]
    if v_filtered:
        radii = [r['r_mid'] for r in v_filtered]
        grains = [r['grain'] for r in v_filtered]
        deviations = [r['deviation'] for r in v_filtered]
        
        sc = ax1.scatter(radii, grains, c=deviations, cmap='RdYlGn_r', 
                        s=50, edgecolors='black', linewidth=0.5)
        ax1.plot(radii, grains, 'b-', alpha=0.3)
        
        for r_cross in crossings:
            ax1.axvline(r_cross, color='red', linestyle=':', alpha=0.5)
        
        ax1.axhline(GRAIN_TARGET, color='green', linestyle='--', linewidth=2)
        ax1.axhspan(GRAIN_TARGET * 0.95, GRAIN_TARGET * 1.05, alpha=0.2, color='green')
        
        plt.colorbar(sc, ax=ax1, label='Deviation (%)')
    
    ax1.set_xlabel('Radius (kpc)', fontsize=11)
    ax1.set_ylabel('Grain', fontsize=11)
    ax1.set_title('Grain vs Radius (V=40-130 km/s)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Top 15 matches
    ax2 = axes[0, 1]
    top15 = results_df.nsmallest(15, 'deviation')
    labels = [f"R:{int(r['r_min'])}-{int(r['r_max'])}" for _, r in top15.iterrows()]
    colors = ['green' if d < 1 else 'limegreen' if d < 2 else 'yellow' if d < 5 else 'orange'
              for d in top15['deviation']]
    ax2.barh(range(len(top15)), top15['deviation'], color=colors, edgecolor='black')
    ax2.set_yticks(range(len(top15)))
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.axvline(1, color='green', linestyle='--', label='1%')
    ax2.axvline(5, color='orange', linestyle='--', label='5%')
    ax2.set_xlabel('Deviation from target (%)', fontsize=11)
    ax2.set_title('Top 15 Matches', fontsize=12, fontweight='bold')
    ax2.legend()
    
    # Panel 3: R² distribution
    ax3 = axes[1, 0]
    ax3.hist(results_df['r_squared'], bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax3.axvline(0.99, color='green', linestyle='--', label='R²=0.99')
    ax3.axvline(0.95, color='orange', linestyle='--', label='R²=0.95')
    ax3.set_xlabel('R² (fit quality)', fontsize=11)
    ax3.set_ylabel('Count', fontsize=11)
    ax3.set_title('Distribution of Fit Quality', fontsize=12, fontweight='bold')
    ax3.legend()
    
    # Panel 4: Grain distribution
    ax4 = axes[1, 1]
    ax4.hist(results_df['grain'], bins=40, color='purple', edgecolor='black', alpha=0.7)
    ax4.axvline(GRAIN_TARGET, color='green', linestyle='--', linewidth=2, label=f'Target: {GRAIN_TARGET}')
    ax4.axvspan(GRAIN_TARGET * 0.95, GRAIN_TARGET * 1.05, alpha=0.2, color='green')
    ax4.set_xlabel('Grain', fontsize=11)
    ax4.set_ylabel('Count', fontsize=11)
    ax4.set_title('Distribution of Measured Grains', fontsize=12, fontweight='bold')
    ax4.legend()
    
    plt.suptitle(f'FIRE-2 Grain Parameter Scan\nTarget: γ = {GRAIN_TARGET}', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('fire2_grain_scan_results.png', dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()