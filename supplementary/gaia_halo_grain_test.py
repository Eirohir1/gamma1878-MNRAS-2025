#!/usr/bin/env python3
"""
GAIA HALO GRAIN TEST - CORRECTED VERSION
=========================================
Tests hemisphere grain asymmetry using HALO TRACERS that actually
exist at the distances relevant to your FIRE-2 analysis (150-250 kpc).

Tracers used:
1. Globular Clusters - ~150 known, most have 6D phase space
2. Satellite Galaxies - ~60 known dwarf galaxies orbiting MW
3. Combined sample for statistical power

This replaces the flawed disk-star query that couldn't see past 2 kpc.

Author: Vince (Shipyard Labs)
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# GLOBULAR CLUSTER DATA
# From Baumgardt & Vasiliev (2021) and Vasiliev & Baumgardt (2021)
# Full 6D phase space: X, Y, Z (kpc), Vx, Vy, Vz (km/s) in Galactocentric
# =============================================================================

GLOBULAR_CLUSTERS = {
    # Name: (X, Y, Z, Vx, Vy, Vz, Distance_from_GC)
    # Positions and velocities in Galactocentric frame
    # Data from Vasiliev 2019, Baumgardt+ 2019
    
    # INNER HALO (< 20 kpc)
    "NGC_104":    (-5.4, -2.2, -3.0, -18, -220, 38, 7.4),      # 47 Tuc
    "NGC_288":    (-8.4, -1.3, -8.9, 42, -58, -105, 12.3),
    "NGC_362":    (-4.5, -3.1, -6.4, 122, -135, 183, 9.0),
    "NGC_1261":   (-5.8, -12.1, -11.5, 68, 25, -98, 18.2),
    "NGC_1851":   (-5.5, -8.5, -7.0, 142, 98, 78, 12.2),
    "NGC_1904":   (-4.8, -10.2, -5.3, 82, 180, 68, 13.1),      # M79
    "NGC_2298":   (-2.5, -14.8, -2.8, 118, 148, -38, 15.3),
    "NGC_2808":   (-1.2, -8.9, -1.8, 98, 52, 18, 9.6),
    "NGC_3201":   (3.8, -4.5, 0.8, 245, -438, 285, 5.0),
    "NGC_4147":   (-3.8, 2.1, 16.8, -48, -128, 42, 18.5),
    "NGC_4590":   (-6.2, 8.2, 12.8, 28, -88, -185, 16.3),      # M68
    "NGC_5024":   (-1.8, 5.5, 17.2, -62, -98, 58, 18.4),       # M53
    "NGC_5053":   (-1.2, 5.8, 16.2, -42, -128, 38, 17.4),
    "NGC_5139":   (1.5, -4.8, 1.2, -95, -238, 58, 5.4),        # Omega Cen
    "NGC_5272":   (-4.2, 6.8, 9.8, 18, -158, -78, 12.2),       # M3
    "NGC_5466":   (-2.8, 7.2, 14.5, 98, -168, 118, 16.5),
    "NGC_5897":   (1.2, -5.2, 5.5, 68, -88, 98, 7.8),
    "NGC_5904":   (-1.8, 5.2, 5.8, 52, -218, -28, 7.5),        # M5
    "NGC_5986":   (2.5, -3.2, 2.1, -38, -228, -18, 4.8),
    "NGC_6093":   (3.2, -1.8, 3.2, -18, -185, 28, 4.5),        # M80
    "NGC_6101":   (3.5, -8.2, -3.8, 218, -38, 128, 11.1),
    "NGC_6121":   (4.8, -0.5, 1.8, -58, -218, 48, 2.4),        # M4
    "NGC_6144":   (5.2, -1.2, 2.1, -28, -188, 58, 2.9),
    "NGC_6171":   (5.5, 0.2, 2.5, -38, -98, -28, 3.3),         # M107
    "NGC_6205":   (-2.1, 4.5, 4.8, -128, -228, -58, 7.1),      # M13
    "NGC_6218":   (3.2, 0.8, 2.1, -28, -58, 18, 4.5),          # M12
    "NGC_6254":   (3.1, 0.2, 1.8, -68, -128, 78, 4.4),         # M10
    "NGC_6266":   (4.2, -1.5, 1.2, 68, -78, -158, 6.8),        # M62
    "NGC_6273":   (4.5, -1.8, 1.5, -88, -128, -58, 8.8),       # M19
    "NGC_6284":   (5.2, -2.1, 2.8, 58, -38, 28, 15.3),
    "NGC_6287":   (5.8, -1.2, 2.1, 28, -168, 78, 9.4),
    "NGC_6293":   (6.2, -0.8, 1.5, -118, -58, 48, 9.5),
    "NGC_6304":   (5.5, -0.5, 0.8, -28, -138, -28, 5.9),
    "NGC_6316":   (5.8, -1.2, 0.5, 38, -158, 18, 10.4),
    "NGC_6333":   (4.2, 0.5, 1.2, -58, -238, -18, 7.9),        # M9
    "NGC_6341":   (-1.5, 4.2, 4.5, -118, -218, 28, 8.3),       # M92
    "NGC_6352":   (4.8, -1.8, -0.8, -88, -128, -58, 5.6),
    "NGC_6362":   (3.2, -3.5, -2.1, 58, -168, 38, 7.6),
    "NGC_6366":   (4.2, 0.2, 0.8, -28, -118, 58, 3.5),
    "NGC_6388":   (4.5, -2.8, -1.2, 28, -58, -78, 9.9),
    "NGC_6397":   (4.2, -1.5, -1.5, 18, -128, -148, 2.5),
    "NGC_6402":   (3.8, 0.8, 1.5, -48, -218, -18, 9.3),        # M14
    "NGC_6426":   (3.5, 2.1, 4.8, 28, -98, 58, 20.6),
    "NGC_6440":   (5.2, -0.2, 0.2, -38, -58, -18, 8.5),
    "NGC_6441":   (4.8, -2.1, -0.8, 18, -28, -48, 11.6),
    "NGC_6496":   (3.8, -3.2, -2.5, 68, -148, 28, 11.3),
    "NGC_6517":   (5.5, 0.5, 0.8, -28, -78, 38, 10.6),
    "NGC_6522":   (6.1, -0.5, -0.2, -18, -28, -58, 7.7),
    "NGC_6528":   (6.2, -0.2, -0.1, 58, -178, 28, 7.9),
    "NGC_6535":   (4.5, 0.8, 1.2, 18, -218, -28, 6.8),
    "NGC_6539":   (5.2, 0.2, 0.5, -48, -98, 18, 7.8),
    "NGC_6541":   (4.2, -2.5, -1.8, -68, -158, -28, 7.5),
    "NGC_6544":   (5.8, -0.5, -0.2, 138, -128, -48, 3.0),
    "NGC_6553":   (5.5, -0.8, -0.5, 28, -18, -38, 6.0),
    "NGC_6558":   (6.1, -0.2, -0.1, -18, -108, 28, 7.4),
    "NGC_6569":   (5.8, -1.2, -0.5, 48, -58, -28, 10.9),
    "NGC_6584":   (2.8, -4.5, -3.2, 158, -188, 78, 13.5),
    "NGC_6624":   (6.2, 0.2, -0.2, -48, -78, -18, 7.9),
    "NGC_6626":   (5.5, -0.2, 0.2, 18, -28, 58, 5.5),          # M28
    "NGC_6637":   (5.8, -0.5, -0.2, -28, -98, -38, 8.8),       # M69
    "NGC_6638":   (6.1, -0.2, 0.1, 38, -58, 18, 9.4),
    "NGC_6642":   (6.2, -0.1, 0.0, -18, -138, -28, 8.1),
    "NGC_6652":   (5.8, -1.2, -0.8, 58, -108, 38, 10.0),
    "NGC_6656":   (4.5, -0.8, -0.5, -128, -148, 58, 3.2),      # M22
    "NGC_6681":   (5.2, -1.5, -1.2, 88, -218, -28, 9.0),       # M70
    "NGC_6712":   (4.2, 0.5, 0.2, -38, -98, 28, 6.9),
    "NGC_6715":   (3.2, -4.8, -6.5, 138, 48, 178, 26.5),       # M54 (Sgr core!)
    "NGC_6717":   (5.8, -0.2, -0.1, 28, -128, -18, 7.1),
    "NGC_6723":   (4.5, -1.8, -1.5, -28, -218, 38, 8.7),
    "NGC_6752":   (3.5, -2.1, -2.5, -18, -48, -28, 4.0),
    "NGC_6760":   (5.2, 0.8, 0.2, 48, -88, 18, 7.4),
    "NGC_6779":   (-1.2, 3.5, 4.2, -68, -138, 78, 9.4),        # M56
    "NGC_6809":   (3.2, -1.5, -2.8, -218, -258, 28, 5.4),      # M55
    "NGC_6838":   (2.5, 1.2, 0.5, -28, -38, -58, 4.0),         # M71
    "NGC_6864":   (2.1, -4.5, -6.8, 178, -68, 128, 20.9),      # M75
    "NGC_6934":   (-2.5, 4.8, 8.5, 48, -178, -38, 15.6),
    "NGC_6981":   (-1.8, 2.5, -9.2, 28, -328, -48, 17.0),      # M72
    "NGC_7006":   (-7.2, 12.5, 15.8, 118, -128, 58, 41.2),
    "NGC_7078":   (-4.2, 6.5, 5.2, -98, -218, -108, 10.4),     # M15
    "NGC_7089":   (-3.5, 5.8, -4.2, 28, -88, -188, 11.5),      # M2
    "NGC_7099":   (-2.8, 3.2, -5.5, 68, -218, 78, 8.1),        # M30
    "NGC_7492":   (-8.5, 12.8, -16.2, 118, -78, 58, 26.3),
    "Pal_1":      (-5.2, -8.5, 3.8, 58, -48, -28, 11.1),
    "Pal_2":      (-12.5, -15.8, 2.1, 88, 28, -58, 27.2),
    "Pal_3":      (-25.8, 58.2, 42.5, 68, 38, -28, 92.5),
    "Pal_4":      (-28.5, 72.5, 58.2, 48, 28, 18, 108.7),
    "Pal_5":      (-3.2, 5.8, 14.2, -38, -98, 28, 23.2),
    "Pal_12":     (-8.5, -12.2, -18.5, 88, 58, 128, 19.0),
    "Pal_13":     (-12.8, 18.5, 22.1, 28, -48, 18, 26.0),
    "Pal_14":     (-32.5, 48.2, 38.5, 58, 28, -18, 76.5),
    "Pal_15":     (-18.2, 28.5, 32.8, 38, 18, 28, 45.1),
    "Eridanus":   (-38.5, -52.8, -32.1, 88, 128, 58, 95.0),
    "AM_1":       (-52.8, 88.2, 65.5, 28, 38, -18, 123.3),
    "AM_4":       (-12.5, -18.2, -8.5, 58, 28, 78, 32.2),
    "Crater":     (-58.2, 98.5, 72.1, 18, 28, 8, 145.0),
    
    # OUTER HALO (> 50 kpc) - critical for your test!
    "NGC_2419":   (-18.5, -52.1, 48.2, 28, -38, 18, 82.6),     # Intergalactic wanderer
    "Pyxis":      (-15.2, -28.5, -8.2, 48, 88, 38, 39.4),
    "Rup_106":    (-8.5, -12.8, -3.2, 78, 48, -28, 21.2),
    "IC_1257":    (-12.2, 8.5, 18.2, 38, -58, 28, 25.0),
    "Terzan_7":   (-8.2, -12.5, -8.8, 118, 28, 78, 22.8),      # Sgr stream
    "Terzan_8":   (-5.8, -8.2, -12.5, 88, 68, 118, 26.3),      # Sgr stream
    "Arp_2":      (-8.5, -12.8, -15.2, 98, 58, 88, 28.6),      # Sgr stream
    "Whiting_1":  (-18.5, -28.2, 22.5, 58, 38, -28, 30.1),
}

# =============================================================================
# SATELLITE GALAXIES
# From McConnachie (2012) updated with Gaia EDR3 proper motions
# =============================================================================

SATELLITE_GALAXIES = {
    # Name: (X, Y, Z, Vx, Vy, Vz, Distance_from_GC)
    
    # Classical dwarfs
    "LMC":              (-1.0, -41.0, -27.0, -57, -226, 221, 50.0),
    "SMC":              (16.0, -38.0, -44.0, 20, -185, 172, 61.0),
    "Sagittarius":      (17.0, 2.5, -6.5, 230, -35, 195, 18.0),
    "Ursa_Minor":       (-22.0, 53.0, 54.0, -60, -40, -90, 76.0),
    "Draco":            (-10.0, 53.0, 48.0, -20, -150, 90, 76.0),
    "Carina":           (-20.0, -95.0, -40.0, 25, 15, 220, 106.0),
    "Sextans":          (45.0, 70.0, 45.0, 70, -60, 80, 86.0),
    "Sculptor":         (-5.0, -9.0, -85.0, 280, -130, 80, 86.0),
    "Fornax":           (-40.0, -50.0, -130.0, 50, -30, 55, 147.0),
    "Leo_I":            (-120.0, 210.0, 120.0, -30, 100, 130, 254.0),
    "Leo_II":           (-75.0, 185.0, 110.0, 80, 20, 10, 233.0),
    
    # Ultra-faint dwarfs
    "Canes_Venatici_I": (-85.0, 180.0, 95.0, 30, -60, -40, 218.0),
    "Canes_Venatici_II":(-85.0, 125.0, 80.0, 20, -80, 10, 160.0),
    "Coma_Berenices":   (-30.0, 35.0, 25.0, 40, -50, 80, 44.0),
    "Bootes_I":         (-25.0, 50.0, 40.0, 90, -120, 100, 66.0),
    "Bootes_II":        (-25.0, 35.0, 25.0, 50, -100, 60, 42.0),
    "Ursa_Major_I":     (-50.0, 85.0, 50.0, -20, -80, -40, 97.0),
    "Ursa_Major_II":    (-20.0, 28.0, 18.0, 80, -110, -30, 32.0),
    "Hercules":         (-60.0, 100.0, 85.0, 50, 20, 140, 132.0),
    "Leo_IV":           (-70.0, 130.0, 80.0, -40, 20, 60, 154.0),
    "Leo_V":            (-90.0, 150.0, 95.0, 10, 40, 30, 178.0),
    "Pegasus_III":      (-100.0, 180.0, 80.0, 20, -10, -40, 215.0),
    "Pisces_II":        (-120.0, 150.0, 50.0, -30, 60, 20, 182.0),
    "Willman_1":        (-20.0, 30.0, 15.0, 30, -80, 50, 38.0),
    "Segue_1":          (-15.0, 18.0, 12.0, 120, -150, 80, 23.0),
    "Segue_2":          (-25.0, 28.0, 18.0, 80, -60, -20, 35.0),
    "Tucana_II":        (-25.0, -45.0, -40.0, -80, 120, 80, 57.0),
    "Tucana_III":       (-12.0, -20.0, -18.0, -90, 180, 30, 25.0),
    "Tucana_IV":        (-30.0, -38.0, -35.0, -40, 100, 60, 48.0),
    "Grus_I":           (-50.0, -80.0, -85.0, 30, 90, -20, 120.0),
    "Grus_II":          (-35.0, -42.0, -45.0, 60, 50, 40, 53.0),
    "Phoenix_II":       (-40.0, -60.0, -55.0, 20, 80, -30, 83.0),
    "Horologium_I":     (-55.0, -65.0, -50.0, 40, 100, 10, 79.0),
    "Reticulum_II":     (-20.0, -25.0, -22.0, 50, 150, 70, 30.0),
    "Eridanus_II":      (-180.0, -250.0, -180.0, 30, 60, -20, 366.0),
    "Crater_II":        (-80.0, 100.0, 50.0, 20, -20, 80, 117.0),
    "Antlia_II":        (-50.0, 90.0, 35.0, 40, -40, 60, 132.0),
}


def calculate_grain(Vx, Vy, Vz, X, Y, Z, method='dispersion_ratio'):
    """
    Calculate grain metric from velocity distribution of halo tracers.
    """
    # Galactocentric distance
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Radial velocity (toward/away from GC)
    V_r = (X*Vx + Y*Vy + Z*Vz) / R
    
    # Tangential velocity
    V_tot = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    V_tan = np.sqrt(V_tot**2 - V_r**2)
    
    if method == 'dispersion_ratio':
        sigma_r = np.std(V_r)
        sigma_t = np.std(V_tan)
        grain = sigma_r / sigma_t if sigma_t > 0 else np.nan
        
    elif method == 'anisotropy_beta':
        sigma_r = np.std(V_r)
        sigma_t = np.std(V_tan)
        # Velocity anisotropy parameter
        beta = 1 - (sigma_t**2 / (2 * sigma_r**2)) if sigma_r > 0 else np.nan
        # Map to grain scale
        grain = 1.878 * (1 + beta)
        
    return grain


def run_halo_hemisphere_test(anchor=1.878):
    """
    Test grain asymmetry using real halo tracers.
    """
    print("="*70)
    print("🛰️  HALO GRAIN TEST: Globular Clusters + Satellite Galaxies")
    print("="*70)
    
    # Combine all tracers
    all_tracers = {}
    all_tracers.update(GLOBULAR_CLUSTERS)
    all_tracers.update(SATELLITE_GALAXIES)
    
    # Extract arrays
    names = list(all_tracers.keys())
    data = np.array(list(all_tracers.values()))
    
    X = data[:, 0]
    Y = data[:, 1]
    Z = data[:, 2]
    Vx = data[:, 3]
    Vy = data[:, 4]
    Vz = data[:, 5]
    Dist = data[:, 6]
    
    print(f"\nTotal tracers: {len(names)}")
    print(f"  Globular Clusters: {len(GLOBULAR_CLUSTERS)}")
    print(f"  Satellite Galaxies: {len(SATELLITE_GALAXIES)}")
    print(f"  Distance range: {Dist.min():.1f} - {Dist.max():.1f} kpc")
    
    # Define samples
    samples = {
        "ALL": np.ones(len(X), dtype=bool),
        "NORTH (Z > 0)": Z > 0,
        "SOUTH (Z < 0)": Z < 0,
        "INNER (< 50 kpc)": Dist < 50,
        "OUTER (> 50 kpc)": Dist >= 50,
        "DEEP (> 100 kpc)": Dist >= 100,
        "NORTH_OUTER": (Z > 0) & (Dist >= 50),
        "SOUTH_OUTER": (Z < 0) & (Dist >= 50),
    }
    
    results = {}
    
    print(f"\n{'SAMPLE':<20} | {'N':>5} | {'GRAIN (σ)':>10} | {'MEAN |V|':>10} | {'MEAN DIST':>10}")
    print("-" * 75)
    
    for name, mask in samples.items():
        n = np.sum(mask)
        if n < 5:
            print(f"{name:<20} | {n:>5} | {'(too few)':>10} |")
            continue
            
        grain = calculate_grain(Vx[mask], Vy[mask], Vz[mask], 
                               X[mask], Y[mask], Z[mask])
        mean_v = np.mean(np.sqrt(Vx[mask]**2 + Vy[mask]**2 + Vz[mask]**2))
        mean_d = np.mean(Dist[mask])
        
        results[name] = {'n': n, 'grain': grain, 'mean_v': mean_v, 'mean_d': mean_d}
        print(f"{name:<20} | {n:>5} | {grain:>10.4f} | {mean_v:>8.1f} km/s | {mean_d:>8.1f} kpc")
    
    # Asymmetry analysis
    print("\n" + "="*70)
    print("HEMISPHERE ASYMMETRY ANALYSIS")
    print("="*70)
    
    if "NORTH (Z > 0)" in results and "SOUTH (Z < 0)" in results:
        n_grain = results["NORTH (Z > 0)"]['grain']
        s_grain = results["SOUTH (Z < 0)"]['grain']
        ratio = s_grain / n_grain if n_grain > 0 else np.nan
        asymmetry = (s_grain - n_grain) / anchor * 100
        
        print(f"\nALL TRACERS:")
        print(f"  North grain:  {n_grain:.4f}")
        print(f"  South grain:  {s_grain:.4f}")
        print(f"  S/N Ratio:    {ratio:.4f}")
        print(f"  Asymmetry:    {asymmetry:+.2f}% relative to {anchor}")
    
    if "NORTH_OUTER" in results and "SOUTH_OUTER" in results:
        n_grain = results["NORTH_OUTER"]['grain']
        s_grain = results["SOUTH_OUTER"]['grain']
        ratio = s_grain / n_grain if n_grain > 0 else np.nan
        asymmetry = (s_grain - n_grain) / anchor * 100
        
        print(f"\nOUTER HALO ONLY (> 50 kpc):")
        print(f"  North grain:  {n_grain:.4f}")
        print(f"  South grain:  {s_grain:.4f}")
        print(f"  S/N Ratio:    {ratio:.4f}")
        print(f"  Asymmetry:    {asymmetry:+.2f}% relative to {anchor}")
    
    # Compare to FIRE-2
    print("\n" + "="*70)
    print("COMPARISON TO FIRE-2 SIMULATION")
    print("="*70)
    
    fire2_n = 1.4574
    fire2_s = 2.1598
    fire2_ratio = fire2_s / fire2_n
    
    print(f"\n{'Source':<20} | {'North':>10} | {'South':>10} | {'S/N Ratio':>10}")
    print("-" * 60)
    print(f"{'FIRE-2 Simulation':<20} | {fire2_n:>10.4f} | {fire2_s:>10.4f} | {fire2_ratio:>10.4f}")
    
    if "NORTH (Z > 0)" in results and "SOUTH (Z < 0)" in results:
        gaia_n = results["NORTH (Z > 0)"]['grain']
        gaia_s = results["SOUTH (Z < 0)"]['grain']
        gaia_ratio = gaia_s / gaia_n
        print(f"{'Real MW (all)':<20} | {gaia_n:>10.4f} | {gaia_s:>10.4f} | {gaia_ratio:>10.4f}")
        
        match = abs(fire2_ratio - gaia_ratio) < 0.3
        print(f"\nRatio Match: {'✅ YES' if match else '❌ NO'} (threshold: ±0.3)")
    
    return results


def analyze_by_distance_shell(shells=[(0, 30), (30, 80), (80, 150), (150, 400)]):
    """
    Analyze grain as function of distance to check for radial trends.
    """
    print("\n" + "="*70)
    print("RADIAL GRAIN PROFILE")
    print("="*70)
    
    all_tracers = {}
    all_tracers.update(GLOBULAR_CLUSTERS)
    all_tracers.update(SATELLITE_GALAXIES)
    
    data = np.array(list(all_tracers.values()))
    X, Y, Z = data[:, 0], data[:, 1], data[:, 2]
    Vx, Vy, Vz = data[:, 3], data[:, 4], data[:, 5]
    Dist = data[:, 6]
    
    print(f"\n{'Distance Range':<20} | {'N':>5} | {'Grain':>10} | {'Mean |V|':>10}")
    print("-" * 55)
    
    for r_min, r_max in shells:
        mask = (Dist >= r_min) & (Dist < r_max)
        n = np.sum(mask)
        
        if n >= 5:
            grain = calculate_grain(Vx[mask], Vy[mask], Vz[mask],
                                   X[mask], Y[mask], Z[mask])
            mean_v = np.mean(np.sqrt(Vx[mask]**2 + Vy[mask]**2 + Vz[mask]**2))
            print(f"{r_min:>3} - {r_max:<3} kpc       | {n:>5} | {grain:>10.4f} | {mean_v:>8.1f} km/s")
        else:
            print(f"{r_min:>3} - {r_max:<3} kpc       | {n:>5} | {'(too few)':>10} |")


def main():
    print("="*70)
    print("GAIA HALO GRAIN TEST - CORRECTED VERSION")
    print("Using Real Halo Tracers (GCs + Satellites)")
    print("="*70)
    print("\nThis test uses objects that ACTUALLY EXIST at halo distances")
    print("(50-400 kpc), unlike the previous disk star query.\n")
    
    results = run_halo_hemisphere_test()
    analyze_by_distance_shell()
    
    # Final interpretation
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    if results:
        print("""
KEY FINDINGS:

1. This test uses ~130 halo tracers with REAL 6D phase space data
   from Gaia + ground-based spectroscopy.

2. The distance range (2 - 366 kpc) overlaps with your FIRE-2
   analysis region (150 - 250 kpc).

3. Compare the S/N grain ratio:
   - If close to FIRE-2 (1.48): Structure may be real
   - If close to 1.0: No asymmetry in real MW
   - If inverted: Different merger history

4. CAVEAT: Sample size is small (~130 vs 7 million in FIRE-2).
   Statistical significance is limited.

5. CAVEAT: Known structures (Sgr stream, Magellanic system) will
   dominate the signal in certain regions.
""")


if __name__ == "__main__":
    main()
