#!/usr/bin/env python3
"""
GAIA HALO HUNTER - HIGH POWER VERSION
======================================
Maximum sample size for maximum statistical power.

Strategy:
- Remove TOP limit (get ALL qualifying stars)
- Slightly loosened cuts to maximize sample
- Expected: 100k-500k stars (vs original 5k)
- Expected significance: 10-25σ (vs original 2.44σ)

SMOKING GUN MODE: If Gaia validates, it should VALIDATE.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astroquery.gaia import Gaia
import warnings
warnings.filterwarnings("ignore")

def main():
    print("=" * 60)
    print("   GAIA DR3: HIGH-POWER VALIDATION")
    print("   'BB gun → Howitzer' Mode")
    print("=" * 60)
    
    # STRATEGY 1: Try full sample first
    print("\n[1] Attempting FULL SAMPLE query (no TOP limit)...")
    print("    Expected: 100k-500k stars (vs 5k original)")
    
    # Query without TOP limit
    # Loosened cuts slightly:
    # - radial_velocity_error < 15 km/s (was 10) - more stars, still good quality
    # - abs(radial_velocity) > 200 km/s (was 250) - captures more halo
    query_full = """
    SELECT 
        radial_velocity, 
        radial_velocity_error,
        parallax,
        parallax_error
    FROM gaiadr3.gaia_source
    WHERE 
        radial_velocity IS NOT NULL
        AND abs(radial_velocity) > 200
        AND radial_velocity_error < 15
        AND parallax IS NOT NULL
        AND parallax_error IS NOT NULL
    """
    
    try:
        print("    Launching query (this may take 2-5 minutes)...")
        job = Gaia.launch_job_async(query_full)  # async for large queries
        results = job.get_results()
        n_total = len(results)
        print(f"    ✓ SUCCESS! Downloaded {n_total:,} stars")
        print(f"      (That's {n_total/5000:.1f}× more than original!)")
        
    except Exception as e:
        print(f"    ✗ Full query failed: {e}")
        print("    Falling back to TOP 500000...")
        
        # STRATEGY 2: Fallback to large TOP limit
        query_fallback = query_full.replace("SELECT ", "SELECT TOP 500000 ")
        try:
            job = Gaia.launch_job(query_fallback)
            results = job.get_results()
            n_total = len(results)
            print(f"    ✓ Got {n_total:,} stars with TOP limit")
        except Exception as e:
            print(f"    ✗ Even fallback failed: {e}")
            return
    
    # Extract data
    velocities = np.abs(results['radial_velocity'].data)
    v_errors = results['radial_velocity_error'].data
    
    # Quality filter: Only keep high-quality measurements
    quality_mask = (v_errors < 10)  # Tighten after download
    velocities_hq = velocities[quality_mask]
    n_hq = len(velocities_hq)
    
    print(f"\n[2] Quality filtering:")
    print(f"    Total downloaded: {n_total:,}")
    print(f"    High quality (error < 10 km/s): {n_hq:,}")
    print(f"    Quality fraction: {n_hq/n_total*100:.1f}%")
    
    # Analyze high-velocity tail
    v_min, v_max = 250, 450  # Standard window
    mask = (velocities_hq >= v_min) & (velocities_hq <= v_max)
    v_window = velocities_hq[mask]
    n_window = len(v_window)
    
    print(f"\n[3] Analysis window ({v_min}-{v_max} km/s):")
    print(f"    Stars in window: {n_window:,}")
    print(f"    Original script: 5,000")
    print(f"    IMPROVEMENT: {n_window/5000:.1f}× more stars")
    
    if n_window < 100:
        print("    ✗ Insufficient stars even with expanded query")
        return
    
    # Measure slope (same method as original)
    hist, bins = np.histogram(v_window, bins=50, density=True)  # More bins for more data
    centers = (bins[:-1] + bins[1:]) / 2
    
    valid = hist > 0
    log_v = np.log10(centers[valid])
    log_p = np.log10(hist[valid])
    
    slope, intercept, r_val, _, std_err = stats.linregress(log_v, log_p)
    gamma_measured = abs(slope)
    
    # Calculate uncertainty (statistical only)
    gamma_err_stat = std_err  # Regression uncertainty
    
    # Calculate Z-score vs NFW prediction
    gamma_nfw = 5.33
    z_score = (gamma_measured - gamma_nfw) / gamma_err_stat
    
    print(f"\n" + "=" * 60)
    print("   RESULTS (HIGH-POWER)")
    print("=" * 60)
    print(f"  Measured γ:      {gamma_measured:.3f} ± {gamma_err_stat:.3f}")
    print(f"  NFW prediction:  {gamma_nfw:.2f}")
    print(f"  Deviation:       +{(gamma_measured - gamma_nfw):.3f} (+{(gamma_measured - gamma_nfw)/gamma_nfw * 100:.1f}%)")
    print(f"  Z-score:         {z_score:.2f}σ")
    print(f"  R²:              {r_val**2:.4f}")
    print(f"  Sample size:     N = {n_window:,}")
    print("=" * 60)
    
    # Compare to original
    print(f"\n[4] POWER IMPROVEMENT:")
    print(f"    Original: 5,000 stars → 2.44σ")
    print(f"    Now:      {n_window:,} stars → {z_score:.2f}σ")
    print(f"    Factor:   {z_score/2.44:.1f}× more significant")
    
    significance_level = "DISCOVERY" if z_score >= 5 else "EVIDENCE" if z_score >= 3 else "SUGGESTIVE"
    print(f"    Status:   {significance_level} LEVEL")
    
    # Calculate what we'd need for 5σ
    if z_score < 5:
        n_needed = (5 / z_score)**2 * n_window
        print(f"\n    Note: Would need ~{n_needed:,.0f} stars for 5σ")
        print(f"          (Gaia DR4 expected 2026 will have more)")
    
    # Plot (same as original but with updated numbers)
    plt.figure(figsize=(12, 7))
    
    # Data points
    plt.scatter(centers[valid], hist[valid], 
                color='black', s=50, alpha=0.6, 
                label=f'Gaia DR3 (N={n_window:,})', zorder=3)
    
    # Observed fit
    fit_y = 10**(intercept + slope * log_v)
    plt.plot(centers[valid], fit_y, 'r-', linewidth=3, 
             label=f'Measured: γ = {gamma_measured:.3f}±{gamma_err_stat:.3f}', zorder=4)
    
    # NFW prediction
    mid_v = np.sqrt(v_min * v_max)
    mid_idx = np.argmin(np.abs(centers[valid] - mid_v))
    nfw_intercept = log_p[mid_idx] + gamma_nfw * log_v[mid_idx]
    nfw_y = 10**(nfw_intercept - gamma_nfw * log_v)
    plt.plot(centers[valid], nfw_y, 'b--', linewidth=3, alpha=0.7,
             label=f'NFW: γ = {gamma_nfw:.2f}', zorder=2)
    
    # Shaded uncertainty region
    gamma_upper = gamma_measured + gamma_err_stat
    gamma_lower = gamma_measured - gamma_err_stat
    fit_y_upper = 10**(intercept + (-gamma_upper) * log_v)
    fit_y_lower = 10**(intercept + (-gamma_lower) * log_v)
    plt.fill_between(centers[valid], fit_y_lower, fit_y_upper, 
                     color='red', alpha=0.2, zorder=1)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radial Velocity (km/s)', fontsize=14, fontweight='bold')
    plt.ylabel('Probability Density', fontsize=14, fontweight='bold')
    
    title = f'Milky Way Halo Validation (Gaia DR3)\n'
    title += f'γ = {gamma_measured:.3f}±{gamma_err_stat:.3f}  |  '
    title += f'NFW = {gamma_nfw:.2f}  |  '
    title += f'Deviation: +{(gamma_measured - gamma_nfw)/gamma_nfw * 100:.1f}%  |  '
    title += f'Significance: {z_score:.2f}σ ({significance_level})'
    plt.title(title, fontsize=12, fontweight='bold')
    
    plt.legend(fontsize=11, loc='best', framealpha=0.95)
    plt.grid(True, which="both", alpha=0.3, linestyle='--')
    
    # Add text box with statistics
    textstr = f'N = {n_window:,} stars\n'
    textstr += f'R² = {r_val**2:.4f}\n'
    textstr += f'{significance_level} LEVEL'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.02, textstr, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='bottom', bbox=props)
    
    plt.tight_layout()
    plt.savefig('figure3_gaia_validation_HIGH_POWER.png', dpi=200)
    print(f"\n[5] Plot saved: figure3_gaia_validation_HIGH_POWER.png")
    
    # Save results to file (UTF-8 encoding for Greek letters)
    with open('gaia_high_power_results.txt', 'w', encoding='utf-8') as f:
        f.write("GAIA DR3 HIGH-POWER VALIDATION RESULTS\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sample size:        {n_window:,} stars\n")
        f.write(f"Velocity window:    {v_min}-{v_max} km/s\n")
        f.write(f"Measured gamma:     {gamma_measured:.3f} ± {gamma_err_stat:.3f}\n")
        f.write(f"NFW prediction:     {gamma_nfw:.2f}\n")
        f.write(f"Deviation:          +{(gamma_measured - gamma_nfw):.3f} (+{(gamma_measured - gamma_nfw)/gamma_nfw * 100:.1f}%)\n")
        f.write(f"Statistical sig:    {z_score:.2f} sigma ({significance_level} LEVEL)\n")
        f.write(f"R² of fit:          {r_val**2:.4f}\n")
        f.write(f"\nImprovement over original:\n")
        f.write(f"  Stars: 5,000 -> {n_window:,} ({n_window/5000:.1f}x more)\n")
        f.write(f"  Significance: 2.44 sigma -> {z_score:.2f} sigma ({z_score/2.44:.1f}x stronger)\n")
        f.write(f"  Status: SUGGESTIVE -> {significance_level}\n")
    
    print(f"[6] Results saved: gaia_high_power_results.txt")
    
    print("\n" + "=" * 60)
    print("   HIGH-POWER VALIDATION COMPLETE")
    print("=" * 60)
    if z_score >= 5:
        print("   🎯 DISCOVERY-LEVEL SIGNIFICANCE ACHIEVED!")
        print("   Gaia is now a SMOKING GUN, not a BB gun.")
    elif z_score >= 3:
        print("   ✓ EVIDENCE-LEVEL SIGNIFICANCE ACHIEVED!")
        print("   Gaia validation is now STRONG.")
    else:
        print("   ~ Suggestive evidence, stronger than before")
        print("   Gaia DR4 (2026) will provide even more stars")
    print("=" * 60)

if __name__ == "__main__":
    main()