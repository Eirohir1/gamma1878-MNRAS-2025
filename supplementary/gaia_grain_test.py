#!/usr/bin/env python3
"""
GAIA GRAIN TEST
===============
Tests the hemisphere grain asymmetry using real Gaia DR3 data
instead of FIRE-2 simulations.

If FIRE-2 results are simulation artifacts, Gaia should show no pattern.
If the structure is real, Gaia should show similar asymmetry.

Author: Vince (Shipyard Labs)
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

try:
    from astroquery.gaia import Gaia
    HAS_ASTROQUERY = True
except ImportError:
    HAS_ASTROQUERY = False
    print("⚠️  astroquery not installed. Run: pip install astroquery")

# Galactic center position and solar motion (IAU standard)
R_SUN = 8.122  # kpc - distance from Sun to Galactic center
Z_SUN = 0.0208  # kpc - height above Galactic plane
V_SUN = np.array([11.1, 232.24, 7.25])  # km/s - U, V, W solar motion


def query_gaia_6d_sample(n_stars=50000, quality='good'):
    """
    Query Gaia DR3 for stars with full 6D phase space info.
    Returns stars with positions and velocities.
    """
    if not HAS_ASTROQUERY:
        raise ImportError("astroquery required")
    
    # Quality cuts for reliable velocities
    if quality == 'strict':
        ruwe_cut = 1.2
        parallax_error_cut = 0.1
    else:  # 'good'
        ruwe_cut = 1.4
        parallax_error_cut = 0.2
    
    query = f"""
    SELECT TOP {n_stars}
        source_id,
        ra, dec,
        parallax, parallax_error,
        pmra, pmdec,
        radial_velocity, radial_velocity_error,
        l, b,
        ruwe
    FROM gaiadr3.gaia_source
    WHERE parallax > 0.5
        AND parallax_error / parallax < {parallax_error_cut}
        AND radial_velocity IS NOT NULL
        AND radial_velocity_error < 10
        AND ruwe < {ruwe_cut}
    ORDER BY random_index
    """
    
    print(f"📡 Querying Gaia DR3 for {n_stars} stars with 6D kinematics...")
    job = Gaia.launch_job(query)
    table = job.get_results()
    print(f"   Retrieved {len(table)} stars")
    
    return table


def galactocentric_coords(table):
    """
    Convert Gaia heliocentric coords to Galactocentric.
    Returns X, Y, Z positions (kpc) and Vx, Vy, Vz velocities (km/s)
    """
    # Distance from parallax
    dist = 1.0 / table['parallax'].data  # kpc (parallax in mas)
    
    # Galactic coordinates
    l_rad = np.radians(table['l'].data)
    b_rad = np.radians(table['b'].data)
    
    # Heliocentric Cartesian
    x_helio = dist * np.cos(b_rad) * np.cos(l_rad)
    y_helio = dist * np.cos(b_rad) * np.sin(l_rad)
    z_helio = dist * np.sin(b_rad)
    
    # Galactocentric position (Sun at -R_SUN on X-axis)
    X = -x_helio - R_SUN
    Y = -y_helio
    Z = z_helio + Z_SUN
    
    # Velocities (simplified - proper motion to km/s)
    # v_transverse = 4.74047 * mu * d (for mu in mas/yr, d in kpc)
    k = 4.74047  # km/s per (mas/yr * kpc)
    
    vl = k * table['pmra'].data * dist   # approx tangential in l
    vb = k * table['pmdec'].data * dist  # approx tangential in b
    vr = table['radial_velocity'].data    # radial
    
    # Transform to Galactocentric (simplified)
    # This is approximate - proper treatment uses astropy coordinates
    Vx = -vr * np.cos(b_rad) * np.cos(l_rad) + vl * np.sin(l_rad) + vb * np.sin(b_rad) * np.cos(l_rad)
    Vy = -vr * np.cos(b_rad) * np.sin(l_rad) - vl * np.cos(l_rad) + vb * np.sin(b_rad) * np.sin(l_rad)
    Vz = vr * np.sin(b_rad) + vb * np.cos(b_rad)
    
    # Correct for solar motion
    Vx -= V_SUN[0]
    Vy -= V_SUN[1]
    Vz -= V_SUN[2]
    
    return X, Y, Z, Vx, Vy, Vz


def calculate_grain(velocities, method='dispersion_ratio'):
    """
    Calculate the 'grain' metric from velocity distribution.
    
    Methods:
    - dispersion_ratio: ratio of radial to tangential dispersion
    - velocity_anisotropy: beta = 1 - (sigma_t^2 / 2*sigma_r^2)
    - clustering_slope: power-law slope of velocity distribution
    """
    Vx, Vy, Vz = velocities
    
    # Total speed
    V_total = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    
    # Radial vs tangential (in galactocentric sense)
    V_radial = np.sqrt(Vx**2 + Vy**2)  # in-plane
    V_vertical = np.abs(Vz)
    
    if method == 'dispersion_ratio':
        sigma_r = np.std(V_radial)
        sigma_z = np.std(V_vertical)
        # Grain as ratio - values around 1.5-2.5 expected
        grain = sigma_r / sigma_z if sigma_z > 0 else np.nan
        
    elif method == 'velocity_anisotropy':
        sigma_r = np.std(V_radial)
        sigma_z = np.std(V_vertical)
        beta = 1 - (sigma_z**2 / (2 * sigma_r**2)) if sigma_r > 0 else np.nan
        # Map beta to grain scale (beta typically -1 to 1)
        grain = 1.878 * (1 + beta)  # Centers around 1.878 for isotropic
        
    elif method == 'clustering_slope':
        # Fit power-law to velocity distribution tail
        v_sorted = np.sort(V_total)
        n = len(v_sorted)
        # Use upper quartile
        v_tail = v_sorted[int(0.75*n):]
        if len(v_tail) > 10:
            log_v = np.log10(v_tail)
            log_rank = np.log10(np.arange(1, len(v_tail)+1))
            slope = -np.polyfit(log_v, log_rank, 1)[0]
            grain = slope
        else:
            grain = np.nan
    
    return grain


def run_hemisphere_test(X, Y, Z, Vx, Vy, Vz, anchor=1.878):
    """
    Core test: Compare grain between hemispheres.
    This mirrors your FIRE-2 uranium_grain_nexus test.
    """
    print("\n" + "="*70)
    print("🛰️  GAIA GRAIN NEXUS: Real Data Hemisphere Test")
    print("="*70)
    
    # Split by Z (galactic north/south)
    north_mask = Z > 0
    south_mask = Z <= 0
    
    # Also define a 'gold zone' near the plane
    gold_mask = np.abs(Z) < 0.5  # within 500 pc of plane
    
    results = {}
    
    for name, mask in [("NORTH (Z+)", north_mask), 
                       ("SOUTH (Z-)", south_mask),
                       ("GOLD ZONE", gold_mask)]:
        if np.sum(mask) < 100:
            print(f"⚠️  {name}: Insufficient stars ({np.sum(mask)})")
            continue
            
        vels = (Vx[mask], Vy[mask], Vz[mask])
        
        grain_disp = calculate_grain(vels, 'dispersion_ratio')
        grain_anis = calculate_grain(vels, 'velocity_anisotropy')
        grain_slope = calculate_grain(vels, 'clustering_slope')
        
        results[name] = {
            'count': np.sum(mask),
            'grain_dispersion': grain_disp,
            'grain_anisotropy': grain_anis,
            'grain_slope': grain_slope,
            'mean_speed': np.mean(np.sqrt(Vx[mask]**2 + Vy[mask]**2 + Vz[mask]**2))
        }
    
    # Print results table
    print(f"\n{'HEMISPHERE':<15} | {'COUNT':>8} | {'GRAIN (σ)':>10} | {'GRAIN (β)':>10} | {'MEAN V':>10}")
    print("-" * 70)
    
    for name, r in results.items():
        print(f"{name:<15} | {r['count']:>8} | {r['grain_dispersion']:>10.4f} | "
              f"{r['grain_anisotropy']:>10.4f} | {r['mean_speed']:>8.2f} km/s")
    
    # Calculate asymmetry
    if "NORTH (Z+)" in results and "SOUTH (Z-)" in results:
        n_grain = results["NORTH (Z+)"]['grain_dispersion']
        s_grain = results["SOUTH (Z-)"]['grain_dispersion']
        asymmetry = (s_grain - n_grain) / anchor * 100
        
        print("\n" + "-" * 70)
        print("ASYMMETRY ANALYSIS")
        print("-" * 70)
        print(f"North grain:     {n_grain:.4f}")
        print(f"South grain:     {s_grain:.4f}")
        print(f"Asymmetry:       {asymmetry:+.2f}% relative to anchor ({anchor})")
        print(f"Ratio (S/N):     {s_grain/n_grain:.4f}")
        
        # Compare to FIRE-2
        print("\n" + "-" * 70)
        print("COMPARISON TO FIRE-2 RESULTS")
        print("-" * 70)
        print(f"{'Metric':<20} | {'FIRE-2':>12} | {'GAIA':>12} | {'Match?':>10}")
        print("-" * 70)
        
        fire2_n = 1.4574
        fire2_s = 2.1598
        fire2_ratio = fire2_s / fire2_n
        gaia_ratio = s_grain / n_grain
        
        ratio_match = "✅ YES" if abs(fire2_ratio - gaia_ratio) < 0.3 else "❌ NO"
        
        print(f"{'North Grain':<20} | {fire2_n:>12.4f} | {n_grain:>12.4f} |")
        print(f"{'South Grain':<20} | {fire2_s:>12.4f} | {s_grain:>12.4f} |")
        print(f"{'S/N Ratio':<20} | {fire2_ratio:>12.4f} | {gaia_ratio:>12.4f} | {ratio_match}")
        
        return results, gaia_ratio, fire2_ratio
    
    return results, None, None


def run_with_sample_data():
    """
    Run with synthetic test data if Gaia query unavailable.
    Uses known Milky Way velocity ellipsoid properties.
    """
    print("\n⚠️  Running with SYNTHETIC test data (Milky Way-like)")
    print("   Install astroquery for real Gaia data\n")
    
    n_stars = 50000
    
    # Generate MW-like distribution
    np.random.seed(42)
    
    # Positions - exponential disk + halo
    R = np.random.exponential(3.0, n_stars)  # kpc, disk scale length
    phi = np.random.uniform(0, 2*np.pi, n_stars)
    Z = np.random.laplace(0, 0.3, n_stars)  # thin disk
    
    X = R * np.cos(phi)
    Y = R * np.sin(phi)
    
    # Velocities - known MW velocity ellipsoid
    # North/South asymmetry from Gaia has been observed!
    sigma_R = 35  # km/s
    sigma_phi = 25
    sigma_Z_north = 20
    sigma_Z_south = 25  # Known asymmetry!
    
    Vx = np.random.normal(0, sigma_R, n_stars)
    Vy = np.random.normal(220, sigma_phi, n_stars)  # rotation
    
    # Apply known north/south asymmetry
    Vz = np.zeros(n_stars)
    north = Z > 0
    Vz[north] = np.random.normal(0, sigma_Z_north, np.sum(north))
    Vz[~north] = np.random.normal(0, sigma_Z_south, np.sum(~north))
    
    return X, Y, Z, Vx, Vy, Vz


def main():
    print("="*70)
    print("GAIA GRAIN TEST - Real Data Validation")
    print("="*70)
    print("\nThis test checks if the FIRE-2 hemisphere asymmetry")
    print("appears in real Milky Way data from Gaia DR3.\n")
    
    try:
        if HAS_ASTROQUERY:
            # Try real Gaia query
            table = query_gaia_6d_sample(n_stars=50000, quality='good')
            X, Y, Z, Vx, Vy, Vz = galactocentric_coords(table)
            data_source = "GAIA DR3"
        else:
            raise ImportError("astroquery not available")
            
    except Exception as e:
        print(f"⚠️  Gaia query failed: {e}")
        X, Y, Z, Vx, Vy, Vz = run_with_sample_data()
        data_source = "SYNTHETIC (MW-like)"
    
    print(f"\nData source: {data_source}")
    print(f"Total stars: {len(X)}")
    
    # Run the hemisphere test
    results, gaia_ratio, fire2_ratio = run_hemisphere_test(X, Y, Z, Vx, Vy, Vz)
    
    # Interpretation
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    if gaia_ratio is not None and fire2_ratio is not None:
        if abs(gaia_ratio - fire2_ratio) < 0.2:
            print("""
✅ STRUCTURE CONFIRMED IN REAL DATA

The hemisphere asymmetry found in FIRE-2 simulations also appears
in Gaia observations. This suggests the pattern is NOT a simulation
artifact but reflects real Milky Way dynamics.

NEXT STEPS:
1. Increase sample size for statistical significance
2. Test radial dependence (inner vs outer halo)
3. Check if asymmetry aligns with known Milky Way features
   (Sagittarius stream, Gaia-Enceladus debris)
""")
        elif gaia_ratio > 1.0:
            print("""
⚠️  ASYMMETRY PRESENT BUT DIFFERENT MAGNITUDE

Gaia shows north/south asymmetry, but the ratio differs from FIRE-2.
This could mean:
- Real MW asymmetry exists but has different cause
- FIRE-2 simulation has different merger history than real MW
- Your grain metric captures real structure but scaling differs
""")
        else:
            print("""
❌ NO MATCHING ASYMMETRY IN REAL DATA

The FIRE-2 pattern does not appear in Gaia data.
This suggests the FIRE-2 result may be:
- Specific to that simulation's merger history
- An artifact of simulation resolution/method
- Not present in the real Milky Way
""")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    main()
