import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def solve_lane_emden(n):
    """
    Solves the Lane-Emden equation for polytropic index n.
    Standard stellar structure equation for self-gravitating fluids.
    """
    # Boundary conditions near center (xi -> 0)
    xi_0 = 1e-4
    theta_0 = 1.0 - (xi_0**2)/6.0
    dtheta_0 = -xi_0/3.0
    
    def system(xi, y):
        theta, z = y
        if theta < 0: return [0, 0] # Vacuum
        
        # Singularity handling
        if xi < 1e-6: dz = -1
        else: dz = -np.abs(theta)**n - (2/xi)*z
        
        return [z, dz]

    sol = solve_ivp(system, [xi_0, 100], [theta_0, dtheta_0], 
                    method='RK45', dense_output=True, rtol=1e-6)
    return sol.t, sol.y[0], sol.y[1]

def analyze_tsallis_gravity():
    print("🌌 TSALLIS-EFT GRAVITY TEST: q = 1.878")
    print("=======================================")
    
    # 1. DERIVE POLYTROPIC INDEX (n)
    # For a Tsallis gas, the Polytropic Index n is related to q by:
    # P ~ rho^q  =>  gamma = q  =>  n = 1 / (q - 1)
    q = 1.878
    n = 1.0 / (q - 1.0)
    
    print(f"Input q-parameter:    {q}")
    print(f"Derived Index (n):    {n:.4f}")
    
    if n < 5:
        print("   Note: n < 5 implies a finite system (has an edge).")
    else:
        print("   Note: n >= 5 implies an infinite system (fits halo).")

    # 2. SOLVE GRAVITY (Lane-Emden)
    xi, theta, dtheta = solve_lane_emden(n)
    
    # Filter vacuum
    mask = theta > 0
    xi = xi[mask]
    theta = theta[mask]
    dtheta = dtheta[mask]
    
    # 3. COMPUTE ROTATION CURVE
    # v^2 propto M(r)/r
    # M(xi) = -xi^2 * dtheta/dxi
    # v^2 ~ (-xi^2 * dtheta) / xi = -xi * dtheta
    v_sq = -xi * dtheta
    v_circ = np.sqrt(np.maximum(0, v_sq))
    
    # Normalize
    v_norm = v_circ / np.max(v_circ)
    r_norm = xi
    
    # 4. COMPUTE EFT (RADIAL ACCELERATION RELATION)
    # The RAR relates g_obs (Total Gravity) to g_bar (Newtonian Gravity).
    # In this simulation:
    # g_obs = v^2 / r  (The actual acceleration required to hold the fluid)
    # g_bar = Newtonian gravity of the mass distribution calculated so far.
    # For a self-consistent polytrope, g_obs IS g_bar (Newtonian).
    # BUT, we want to see if the SHAPE mimics the Dark Matter "Boost".
    
    # Let's compare the "Tsallis Potential" to a "Keplerian Potential"
    # g_obs (from Tsallis) vs g_newton (Point Mass expectation)
    g_obs = v_sq / xi 
    g_newton = 1.0 / xi**2 # Point source falloff
    
    # Normalize them to match at some inner radius (r=1)
    scale_idx = np.argmin(np.abs(xi - 1.0))
    scale_factor = g_obs[scale_idx] / g_newton[scale_idx]
    g_newton *= scale_factor
    
    # The "Dark Matter Boost" Factor
    boost = g_obs / g_newton

    # 5. PLOT 1: ROTATION CURVE
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    ax1.plot(r_norm, v_norm, 'c-', lw=3, label=f'Tsallis q={q} (EFT Model)')
    
    # Plot Reference Lines
    # Kepler (Point Mass)
    v_kepler = 1.0/np.sqrt(r_norm)
    v_kepler = v_kepler * (v_norm[scale_idx] / v_kepler[scale_idx])
    ax1.plot(r_norm, v_kepler, 'r--', alpha=0.5, label='Newtonian (No DM)')
    
    # Flat (Ideal DM)
    ax1.axhline(v_norm[scale_idx], color='lime', linestyle=':', label='Flat Rotation (Ideal DM)')
    
    ax1.set_title("Rotation Curve prediction", fontsize=14)
    ax1.set_xlabel("Radius (xi)")
    ax1.set_ylabel("Velocity")
    ax1.set_ylim(0, 1.2)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 6. PLOT 2: THE EFT / RAR RELATION
    # This plots the "Boost" vs Acceleration
    # Real galaxies show a boost at low acceleration.
    
    ax2.plot(np.log10(g_newton), np.log10(boost), 'm-', lw=3, label='Gravity Boost Factor')
    ax2.axhline(0, color='black', linestyle='--', label='Newtonian Baseline')
    
    ax2.set_title("EFT Check: The Dark Matter 'Boost'", fontsize=14)
    ax2.set_xlabel("Log(Gravitational Acceleration)")
    ax2.set_ylabel("Log(Boost Factor: g_obs / g_bar)")
    ax2.invert_xaxis() # Astronomy convention: Low acceleration on right
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('tsallis_gravity_result.png')
    print("✅ PLOTS GENERATED: 'tsallis_gravity_result.png'")
    
    # 7. VERDICT
    # Does the curve flatten?
    v_peak = np.max(v_norm)
    v_outer = v_norm[-1]
    flatness = v_outer / v_peak
    
    print(f"\n--- PHYSICS VERDICT ---")
    print(f"Velocity Retention at Edge: {flatness:.1%}")
    
    if flatness > 0.9:
        print("✅ SUCCESS: The q=1.878 model generates a FLAT Rotation Curve.")
        print("   This parameter effectively replaces Dark Matter.")
    elif flatness > 0.5:
        print("⚠️ PARTIAL: The curve is extended, but not perfectly flat.")
        print("   It creates a 'Heavy Halo', but strictly finite.")
    else:
        print("❌ FAILURE: The curve drops like a stone.")
        print("   q=1.878 creates a compact star, not a galaxy halo.")

if __name__ == "__main__":
    analyze_tsallis_gravity()