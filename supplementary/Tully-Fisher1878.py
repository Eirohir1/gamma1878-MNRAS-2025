import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

def solve_galaxy_velocity(mass_scale, radius_scale, q_val=1.878):
    """
    Solves the Polylensic Gravity (TOV) for a galaxy of a given size/density.
    Returns: (Total Mass, Flat Velocity)
    """
    gamma = q_val
    K = 0.1 # Stiffness constant
    
    # We vary the central density based on mass_scale
    # Physics note: Galaxies have roughly constant central surface brightness
    # rho_c * r_scale ~ constant? Let's try varying r0 and keeping rho_c fixed
    # or varying both. Let's vary the 'size' (r0) significantly.
    
    r0 = radius_scale
    rho_c = 1.0 # Normalized central density
    
    P_c = K * rho_c**gamma
    m_0 = (4/3) * np.pi * r0**3 * rho_c
    
    def system(r, y):
        P, m, phi = y
        if P <= 0: return [0, 0, 0]
        rho = (P / K)**(1/gamma)
        
        factor = -1.0 / (r * (r - 2*m)) if r > 2*m else 0
        term1 = (rho + P)
        term2 = (m + 4 * np.pi * r**3 * P)
        
        dp_dr = factor * term1 * term2
        dm_dr = 4 * np.pi * r**2 * rho
        dphi_dr = -1.0 * dp_dr / (rho + P)
        return [dp_dr, dm_dr, dphi_dr]

    # Integrate to edge of galaxy (e.g., 50 * scale radius)
    r_max = 50 * r0
    try:
        sol = solve_ivp(system, [r0, r_max], [P_c, m_0, 0.0], method='LSODA', rtol=1e-5)
    except:
        return None
        
    # Extract properties at the "edge" (Simulating observation)
    r = sol.t
    P = sol.y[0]
    m = sol.y[1]
    
    # Filter valid
    mask = P > 1e-6
    if np.sum(mask) < 10: return None
    
    m_final = m[mask][-1]
    
    # Calculate V_flat
    # v^2 approx M/r in Newtonian limit, but we use GR metric deriv
    # Simplified: V_flat is the max/asymptotic velocity
    r_val = r[mask]
    P_val = P[mask]
    m_val = m[mask]
    metric_factor = 1.0 / (r_val * (r_val - 2*m_val))
    dphi_dr = (m_val + 4*np.pi*r_val**3*P_val) * metric_factor
    v_curve = np.sqrt(r_val * dphi_dr)
    v_flat = np.max(v_curve) # Peak/Flat velocity
    
    return m_final, v_flat

def run_tully_fisher_test():
    print("📏 TULLY-FISHER TEST: Checking Universal Scaling...")
    
    masses = []
    velocities = []
    
    # Simulate 30 Galaxies of different sizes (Dwarfs to Giants)
    # We vary the scale radius r0 from 0.001 to 0.1 (dimensionless units)
    scales = np.logspace(-3, -1, 30)
    
    print(f"   Simulating {len(scales)} galaxies with q=1.878...")
    
    for r_scale in scales:
        res = solve_galaxy_velocity(1.0, r_scale)
        if res:
            m, v = res
            masses.append(m)
            velocities.append(v)
            
    # Log-Log Plot
    log_M = np.log10(masses)
    log_V = np.log10(velocities)
    
    # Linear Fit
    slope, intercept = np.polyfit(log_V, log_M, 1)
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(log_V, log_M, c='cyan', edgecolors='k', s=100, label='Synthetic Galaxies')
    
    x_fit = np.linspace(min(log_V), max(log_V), 100)
    y_fit = slope * x_fit + intercept
    
    plt.plot(x_fit, y_fit, 'r--', lw=2, label=f'Fit Slope = {slope:.2f}')
    
    # Reference Slope (Tully-Fisher is ~4)
    # We plot a reference line passing through the center of data
    y_ref = 4.0 * (x_fit - np.mean(log_V)) + np.mean(log_M)
    plt.plot(x_fit, y_ref, 'g:', lw=2, label='Tully-Fisher Law (Slope=4)')
    
    plt.xlabel('Log(Rotation Speed)')
    plt.ylabel('Log(Galaxy Mass)')
    plt.title(f"Tully-Fisher Emergence (q=1.878) -> Slope: {slope:.2f}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('tully_fisher_test.png')
    
    print(f"✅ PLOT: 'tully_fisher_test.png'")
    print(f"\n--- SCALING LAW VERDICT ---")
    print(f"Derived Slope: {slope:.2f}")
    print(f"Expected Slope: ~4.0")
    
    if 3.5 < slope < 4.5:
        print("✅ SUCCESS: The Polylensic Vacuum naturally reproduces the Tully-Fisher relation!")
    else:
        print("⚠️ DIVERGENCE: The scaling doesn't match standard observations.")

if __name__ == "__main__":
    run_tully_fisher_test()