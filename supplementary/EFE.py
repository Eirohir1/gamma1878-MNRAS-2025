import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def solve_tov_tsallis(q_index):
    print(f"🌌 RELATIVISTIC GRAVITY CHECK (TOV Equation)")
    print(f"   Input Parameter q: {q_index}")
    
    # Constants (dimensionless units for stability)
    # G = c = 1
    # Equation of State: P = K * rho^Gamma
    gamma = q_index 
    K = 0.1 # Stiffness
    
    # Initial Conditions (Center)
    r0 = 1e-3
    rho_c = 1.0 
    P_c = K * rho_c**gamma
    m_0 = (4/3) * np.pi * r0**3 * rho_c
    
    def system(r, y):
        P, m, phi = y
        
        if P <= 0: return [0, 0, 0] # Vacuum boundary
        
        rho = (P / K)**(1/gamma)
        
        # TOV Equations
        factor = -1.0 / (r * (r - 2*m))
        term1 = (rho + P)
        term2 = (m + 4 * np.pi * r**3 * P)
        
        dp_dr = factor * term1 * term2
        dm_dr = 4 * np.pi * r**2 * rho
        dphi_dr = -1.0 * dp_dr / (rho + P) 
        
        return [dp_dr, dm_dr, dphi_dr]

    # Integrate
    sol = solve_ivp(system, [r0, 50], [P_c, m_0, 0.0], 
                    method='LSODA', rtol=1e-6)
    
    r_raw = sol.t
    P_raw = sol.y[0]
    m_raw = sol.y[1]
    phi_raw = sol.y[2]
    
    # --- FIX: Synchronized Filtering ---
    # We create a mask for valid pressure (non-vacuum)
    mask = P_raw > 1e-6
    
    # Apply mask to ALL arrays so they are the same length
    r = r_raw[mask]
    P = P_raw[mask]
    m = m_raw[mask]
    phi = phi_raw[mask]
    
    # --- Recalculate Velocity with clean arrays ---
    # rho must be derived from the FILTERED P
    rho = (P / K)**(1/gamma)
    
    # dphi/dr uses the FILTERED r, m, P
    # Metric factor: 1 / (1 - 2M/r)
    # This prevents the "Singularity" error if 2M ~ r
    metric_factor = 1.0 / (r * (r - 2*m))
    dphi_dr = (m + 4*np.pi*r**3*P) * metric_factor
    
    # Orbital Velocity: v/c = sqrt( r * dphi/dr )
    v_sq = r * dphi_dr
    
    # Safety check for negative squares ( shouldn't happen in stable region)
    v_sq = np.maximum(0, v_sq)
    v_circ = np.sqrt(v_sq) 
    
    return r, v_circ, m

def visualize_results(r, v, q):
    plt.figure(figsize=(10, 6))
    
    # Normalize
    v_norm = v / np.max(v)
    
    plt.plot(r, v_norm, 'c-', lw=3, label=f'Relativistic Halo (q={q})')
    
    # Keplerian Reference (Newton)
    idx_peak = np.argmax(v_norm)
    v_newton = 1.0 / np.sqrt(r / r[idx_peak])
    v_newton[0:idx_peak] = np.nan
    
    plt.plot(r, v_newton, 'r--', alpha=0.5, label='Einstein Vacuum (Keplerian)')
    plt.axhline(v_norm[idx_peak], color='lime', linestyle=':', label='Dark Matter Target (Flat)')
    
    plt.xlabel('Radius (spacetime units)')
    plt.ylabel('Orbital Velocity (v/c)')
    plt.title(f'EFE Solution for Polylensic Matter (q={q})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1.2)
    
    plt.savefig('efe_tsallis_result.png')
    print("✅ PLOT GENERATED: 'efe_tsallis_result.png'")
    
    # VERDICT
    v_outer = v_norm[-1]
    drop = v_outer / 1.0
    
    print(f"\n--- RELATIVISTIC VERDICT ---")
    print(f"Velocity at edge: {drop:.1%} of peak")
    
    if drop > 0.9:
        print("✅ FLAT ROTATION: The q-parameter effectively creates Dark Matter geometry.")
    elif drop > 0.6:
        print("⚠️ EXTENDED: It helps, but it still falls off eventually.")
    else:
        print("❌ KEPLERIAN: No effect. Spacetime behaves normally.")

if __name__ == "__main__":
    GRAIN_PARAMETER = 1.878 
    r_res, v_res, m_res = solve_tov_tsallis(GRAIN_PARAMETER)
    visualize_results(r_res, v_res, GRAIN_PARAMETER)