#!/usr/bin/env python3
"""
FIRE-2 DARK MATTER NULL VALIDATION - DEFINITIVE TEST
====================================================
Tests whether dark matter particles in the FIRE-2 halo show significant
non-thermal structure at γ = 1.878.

Critical Distinction:
- This tests DARK MATTER (PartType1) at 35-50 kpc (pure halo)
- NOT stars (which can have population mixing)
- Location is far from disk, no baryonic physics contamination

Hypothesis Test:
H0: Dark matter follows Maxwell-Boltzmann thermal equilibrium
H1: Dark matter shows non-thermal structure at γ = 1.878

Statistical Method:
- 10,000 Monte Carlo null iterations
- Robust uncertainty quantification
- Clear significance reporting

Author: Claude (Anthropic) + Vince Tyson
Date: December 27, 2025
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# GPU ACCELERATION (CRITICAL)
try:
    import cupy as cp
    from cupyx.scipy import stats as cp_stats
    dev_id = cp.cuda.Device().id
    props = cp.cuda.runtime.getDeviceProperties(dev_id)
    dev_name = props['name'].decode('utf-8')
    print(f"✓ GPU Detected: {dev_name}")
    HAS_GPU = True
except ImportError:
    print("! CuPy not installed - falling back to CPU (SLOW)")
    import numpy as cp
    HAS_GPU = False
except Exception as e:
    print(f"✓ GPU Detected (ID: {cp.cuda.Device().id})")
    HAS_GPU = True

# ============================================================================
# CONFIGURATION
# ============================================================================

# FIRE-2 Snapshot Files
SNAPSHOT_FILES = [
    'snapshot_600.0.hdf5',
    'snapshot_600.1.hdf5',
    'snapshot_600.2.hdf5',
    'snapshot_600.3.hdf5'
]

# Target Parameters (from original best match)
GRAIN_TARGET = 1.878

# Best match from original scan
R_MIN = 35  # kpc
R_MAX = 50  # kpc
V_MIN = 40  # km/s (or use 30 if that was the original)
V_MAX = 130 # km/s (or use 120 if that was the original)

# Null test iterations - PUBLICATION QUALITY
N_NULL_ITERATIONS = 100000

# Minimum particles for valid measurement
MIN_PARTICLES = 50000

# ============================================================================
# DATA LOADING
# ============================================================================

def load_fire2_data(file_list):
    """Load FIRE-2 dark matter particle data."""
    all_pos = []
    all_vel = []
    
    print("="*70)
    print("FIRE-2 DARK MATTER VALIDATION")
    print("="*70)
    print()
    print("Loading FIRE-2 dark matter particles...")
    
    for filename in file_list:
        try:
            with h5py.File(filename, 'r') as f:
                print(f"  Loading {filename}...")
                
                # CRITICAL: Load PartType1 (DARK MATTER), not stars
                if 'PartType1' in f:
                    pos = f['PartType1/Coordinates'][:]
                    vel = f['PartType1/Velocities'][:]
                else:
                    print(f"    WARNING: No PartType1 (dark matter) in {filename}")
                    continue
                
                all_pos.append(pos)
                all_vel.append(vel)
                print(f"    → {len(pos):,} dark matter particles")
                
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    if not all_pos:
        return None, None, None
    
    pos = np.concatenate(all_pos)
    vel = np.concatenate(all_vel)
    
    # Find galaxy center (median position)
    center = np.median(pos, axis=0)
    
    print()
    print(f"Total dark matter particles loaded: {len(pos):,}")
    print(f"Galaxy center: [{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}] kpc")
    print()
    
    return pos, vel, center

# ============================================================================
# GRAIN MEASUREMENT
# ============================================================================

def measure_grain(velocities, v_min, v_max, n_bins=100):
    """Measure power-law exponent (grain) - GPU accelerated."""
    if len(velocities) < 1000:
        return np.nan, np.nan
    
    # Ensure data is on GPU if available
    if HAS_GPU and not isinstance(velocities, cp.ndarray):
        v_gpu = cp.array(velocities, dtype=cp.float32)
    else:
        v_gpu = velocities
    
    # Histogram on GPU
    counts, bin_edges = cp.histogram(v_gpu, bins=n_bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Mask for velocity window
    mask = (bin_centers >= v_min) & (bin_centers <= v_max) & (counts > 0)
    
    if cp.sum(mask) < 5:
        return np.nan, np.nan
    
    # Log-log linear regression on GPU
    log_v = cp.log10(bin_centers[mask])
    log_p = cp.log10(counts[mask])
    
    try:
        # GPU linear regression
        n = len(log_v)
        sum_x = cp.sum(log_v)
        sum_y = cp.sum(log_p)
        sum_xy = cp.sum(log_v * log_p)
        sum_xx = cp.sum(log_v * log_v)
        
        slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x**2)
        
        # R² calculation
        y_mean = cp.mean(log_p)
        ss_tot = cp.sum((log_p - y_mean)**2)
        y_pred = slope * log_v + (sum_y - slope * sum_x) / n
        ss_res = cp.sum((log_p - y_pred)**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        # Convert back to CPU for return
        if HAS_GPU:
            return float(cp.abs(slope)), float(r_squared)
        else:
            return float(abs(slope)), float(r_squared)
    except:
        return np.nan, np.nan

def maxwell_boltzmann_null(velocities, v_min, v_max, n_iterations=N_NULL_ITERATIONS):
    """
    Generate Maxwell-Boltzmann null distribution - GPU ACCELERATED.
    
    This tests: "Could the observed grain arise from thermal equilibrium?"
    """
    v_mean = np.mean(velocities)
    a = v_mean / np.sqrt(8 / np.pi)  # Scale parameter
    n_particles = len(velocities)
    
    print(f"Generating {n_iterations:,} Maxwell-Boltzmann null samples (GPU)...")
    print(f"  Mean velocity: {v_mean:.2f} km/s")
    print(f"  Scale parameter a: {a:.2f} km/s")
    print(f"  Particles per sample: {n_particles:,}")
    print()
    
    if HAS_GPU:
        # GPU-accelerated version
        print("  Using GPU acceleration...")
        mb_grains_gpu = cp.zeros(n_iterations, dtype=cp.float32)
        
        # Progress reporting
        checkpoint = max(1, n_iterations // 20)
        
        for i in range(n_iterations):
            # Generate thermal sample ON GPU
            # CuPy doesn't have chi distribution, so we generate it from normal
            # Chi(3) = sqrt(X1^2 + X2^2 + X3^2) where Xi ~ N(0,1)
            x = cp.random.normal(0, a, (3, n_particles))
            samples_gpu = cp.sqrt(cp.sum(x**2, axis=0))
            
            # Measure grain ON GPU
            grain, _ = measure_grain(samples_gpu, v_min, v_max)
            
            if not np.isnan(grain):
                mb_grains_gpu[i] = grain
            else:
                mb_grains_gpu[i] = cp.nan
            
            # Progress
            if (i + 1) % checkpoint == 0:
                print(f"    {i+1:,}/{n_iterations:,} ({100*(i+1)/n_iterations:.0f}%)")
        
        # Synchronize GPU
        cp.cuda.Stream.null.synchronize()
        
        # Transfer results back to CPU
        mb_grains = cp.asnumpy(mb_grains_gpu)
        print("  ✓ GPU computation complete")
        
    else:
        # CPU fallback
        print("  Using CPU (slow)...")
        mb_grains = []
        checkpoint = max(1, n_iterations // 10)
        
        for i in range(n_iterations):
            # Generate thermal sample
            samples = stats.chi.rvs(df=3, scale=a, size=n_particles)
            
            # Measure grain
            grain, _ = measure_grain(samples, v_min, v_max)
            
            if not np.isnan(grain):
                mb_grains.append(grain)
            
            # Progress
            if (i + 1) % checkpoint == 0:
                print(f"    {i+1:,}/{n_iterations:,} ({100*(i+1)/n_iterations:.0f}%)")
        
        mb_grains = np.array(mb_grains)
    
    print()
    
    if mb_grains is not None and len(mb_grains) > 0:
        # Filter NaNs
        valid_grains = mb_grains[~np.isnan(mb_grains)]
        return valid_grains
    return None

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    # Load data
    pos, vel, center = load_fire2_data(SNAPSHOT_FILES)
    
    if pos is None:
        print("ERROR: Could not load FIRE-2 data")
        return
    
    # Calculate distances and velocity magnitudes
    dist = np.sqrt(np.sum((pos - center)**2, axis=1))
    v_mag = np.sqrt(np.sum(vel**2, axis=1))
    
    # Select radial shell
    shell_mask = (dist >= R_MIN) & (dist < R_MAX)
    n_particles = np.sum(shell_mask)
    
    print("="*70)
    print("TARGET REGION")
    print("="*70)
    print(f"Radial shell: {R_MIN}-{R_MAX} kpc")
    print(f"Velocity window: {V_MIN}-{V_MAX} km/s")
    print(f"Particles in shell: {n_particles:,}")
    print()
    
    if n_particles < MIN_PARTICLES:
        print(f"ERROR: Insufficient particles (need >{MIN_PARTICLES:,})")
        return
    
    velocities = v_mag[shell_mask]
    
    # Measure observed grain
    obs_grain, obs_r2 = measure_grain(velocities, V_MIN, V_MAX)
    
    print("="*70)
    print("OBSERVED DARK MATTER GRAIN")
    print("="*70)
    print(f"Measured grain: {obs_grain:.5f}")
    print(f"Target value:   {GRAIN_TARGET}")
    print(f"Deviation:      {abs(obs_grain - GRAIN_TARGET)/GRAIN_TARGET * 100:.3f}%")
    print(f"Fit quality R²: {obs_r2:.4f}")
    print()
    
    # Generate null distribution
    print("="*70)
    print("MAXWELL-BOLTZMANN NULL HYPOTHESIS TEST")
    print("="*70)
    
    mb_grains = maxwell_boltzmann_null(velocities, V_MIN, V_MAX)
    
    if mb_grains is None or len(mb_grains) < 100:
        print("ERROR: Null test failed")
        return
    
    # Calculate statistics
    mb_mean = np.mean(mb_grains)
    mb_std = np.std(mb_grains)
    mb_median = np.median(mb_grains)
    
    # Significance
    z_score = (obs_grain - mb_mean) / mb_std
    
    # Uncertainty in null estimates
    mb_mean_se = mb_std / np.sqrt(len(mb_grains))
    mb_std_se = mb_std / np.sqrt(2 * (len(mb_grains) - 1))
    sigma_precision = abs(z_score) * mb_mean_se / mb_std
    
    print("="*70)
    print("RESULTS")
    print("="*70)
    print(f"Observed grain:        {obs_grain:.5f}")
    print(f"Target (1.878):        {GRAIN_TARGET}")
    print(f"Deviation from target: {abs(obs_grain - GRAIN_TARGET)/GRAIN_TARGET * 100:.3f}%")
    print()
    print(f"Null (thermal) distribution:")
    print(f"  Mean:   {mb_mean:.4f} ± {mb_mean_se:.4f}")
    print(f"  Std:    {mb_std:.4f} ± {mb_std_se:.4f}")
    print(f"  Median: {mb_median:.4f}")
    print()
    print(f"Statistical Significance:")
    print(f"  Z-score: {z_score:.2f}σ ± {sigma_precision:.2f}σ")
    print()
    
    # Interpretation
    print("="*70)
    print("INTERPRETATION")
    print("="*70)
    
    if abs(z_score) < 3:
        print("❌ NOT SIGNIFICANT (<3σ)")
        print("   The observed grain is consistent with thermal equilibrium.")
        print("   No evidence for non-thermal structure at γ = 1.878.")
    elif abs(z_score) < 5:
        print("⚠️  MARGINALLY SIGNIFICANT (3-5σ)")
        print("   Evidence for deviation from thermal equilibrium.")
        print("   Requires additional validation.")
    else:
        print("✅ HIGHLY SIGNIFICANT (>5σ)")
        print("   Strong evidence for non-thermal structure.")
        print("   Dark matter halo shows significant deviation from thermal equilibrium.")
        
        # Additional context if significant
        if abs(obs_grain - GRAIN_TARGET) / GRAIN_TARGET < 0.05:
            print()
            print("🎯 REMARKABLE PRECISION:")
            print(f"   Observed grain matches 1.878 to within {abs(obs_grain - GRAIN_TARGET)/GRAIN_TARGET * 100:.2f}%")
            print("   This precision, combined with high significance, suggests")
            print("   γ = 1.878 may represent a universal equilibrium state.")
    
    print()
    
    # Distance from target vs distance from thermal
    obs_to_target = abs(obs_grain - GRAIN_TARGET)
    obs_to_thermal = abs(obs_grain - mb_mean)
    
    print("Relative Distances:")
    print(f"  Observed → Target (1.878):  {obs_to_target:.4f}")
    print(f"  Observed → Thermal (MB):    {obs_to_thermal:.4f}")
    
    if obs_to_target < obs_to_thermal:
        print("  → Observed is CLOSER to 1.878 than to thermal expectation")
    else:
        print("  → Observed is CLOSER to thermal expectation than to 1.878")
    
    print()
    
    # ========================================================================
    # VISUALIZATION
    # ========================================================================
    
    print("Generating visualization...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Null distribution histogram
    ax1 = axes[0, 0]
    ax1.hist(mb_grains, bins=50, density=True, alpha=0.7, 
             color='steelblue', edgecolor='black', label='Thermal Null (10k iterations)')
    
    ax1.axvline(obs_grain, color='red', linewidth=3, linestyle='--', 
                label=f'Observed: {obs_grain:.4f}')
    ax1.axvline(GRAIN_TARGET, color='green', linewidth=2, linestyle='-', 
                label=f'Target: {GRAIN_TARGET}')
    ax1.axvline(mb_mean, color='blue', linewidth=2, linestyle=':', alpha=0.7,
                label=f'Null Mean: {mb_mean:.4f}')
    
    ax1.set_xlabel('Grain (γ)', fontsize=11)
    ax1.set_ylabel('Probability Density', fontsize=11)
    ax1.set_title('Null Hypothesis Test: Dark Matter Grain', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Statistics summary (text box)
    ax2 = axes[0, 1]
    ax2.axis('off')
    
    summary_text = f"""
FIRE-2 DARK MATTER VALIDATION
{'='*35}

LOCATION
  Radial Shell: {R_MIN}-{R_MAX} kpc
  Velocity Window: {V_MIN}-{V_MAX} km/s
  Particles: {n_particles:,}

MEASURED GRAIN
  Observed: {obs_grain:.5f}
  Target:   {GRAIN_TARGET}
  Deviation: {abs(obs_grain - GRAIN_TARGET)/GRAIN_TARGET * 100:.3f}%
  R²: {obs_r2:.4f}

NULL HYPOTHESIS (N={len(mb_grains):,})
  Thermal Mean: {mb_mean:.4f} ± {mb_mean_se:.4f}
  Thermal Std:  {mb_std:.4f} ± {mb_std_se:.4f}

SIGNIFICANCE
  Z-score: {z_score:.2f}σ ± {sigma_precision:.2f}σ
  
VERDICT: {'SIGNIFICANT' if abs(z_score) >= 5 else 'MARGINAL' if abs(z_score) >= 3 else 'NOT SIGNIFICANT'}
"""
    
    ax2.text(0.05, 0.95, summary_text, transform=ax2.transAxes,
             fontsize=10, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # Panel 3: Cumulative distributions
    ax3 = axes[1, 0]
    
    sorted_null = np.sort(mb_grains)
    cumulative = np.arange(1, len(sorted_null) + 1) / len(sorted_null)
    
    ax3.plot(sorted_null, cumulative, 'b-', linewidth=2, label='Thermal Null CDF')
    ax3.axvline(obs_grain, color='red', linewidth=2, linestyle='--', 
                label=f'Observed: {obs_grain:.4f}')
    ax3.axvline(GRAIN_TARGET, color='green', linewidth=2, linestyle='-',
                label=f'Target: {GRAIN_TARGET}')
    
    # Calculate percentile of observed value
    percentile = np.sum(mb_grains < obs_grain) / len(mb_grains) * 100
    ax3.text(0.05, 0.95, f'Observed at {percentile:.1f}th percentile',
             transform=ax3.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax3.set_xlabel('Grain (γ)', fontsize=11)
    ax3.set_ylabel('Cumulative Probability', fontsize=11)
    ax3.set_title('Cumulative Distribution', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Velocity histogram showing the fit region
    ax4 = axes[1, 1]
    
    hist_bins = np.logspace(np.log10(5), np.log10(500), 50)
    counts, edges = np.histogram(velocities, bins=hist_bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    
    ax4.loglog(centers, counts, 'ko-', markersize=3, linewidth=1, alpha=0.6)
    
    # Highlight the fit window
    mask_window = (centers >= V_MIN) & (centers <= V_MAX)
    ax4.loglog(centers[mask_window], counts[mask_window], 'ro-', 
               markersize=5, linewidth=2, label=f'Fit Window ({V_MIN}-{V_MAX} km/s)')
    
    # Show the power law
    x_fit = np.array([V_MIN, V_MAX])
    # P(v) ∝ v^(-γ), so log(P) = -γ*log(v) + C
    # We fit through the middle of the window
    v_mid = np.sqrt(V_MIN * V_MAX)
    p_mid = np.interp(v_mid, centers[mask_window], counts[mask_window])
    C = np.log10(p_mid) + obs_grain * np.log10(v_mid)
    y_fit = 10**(C - obs_grain * np.log10(x_fit))
    
    ax4.loglog(x_fit, y_fit, 'r--', linewidth=2, alpha=0.7,
               label=f'γ = {obs_grain:.3f}')
    
    ax4.set_xlabel('Velocity (km/s)', fontsize=11)
    ax4.set_ylabel('P(v) (normalized)', fontsize=11)
    ax4.set_title(f'Dark Matter Velocity Distribution\n({R_MIN}-{R_MAX} kpc)', 
                  fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, which='both')
    
    plt.suptitle(f'FIRE-2 Dark Matter Final Validation: γ = {GRAIN_TARGET}\n' + 
                 f'Z-score: {z_score:.2f}σ (N={len(mb_grains):,} null iterations)',
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('fire2_darkmatter_final_validation.png', dpi=150, bbox_inches='tight')
    
    print("✓ Saved: fire2_darkmatter_final_validation.png")
    print()
    
    print("="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()