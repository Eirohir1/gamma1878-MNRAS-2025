import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def solve_for_flatness(q_val):
    # Constants
    gamma = q_val
    if gamma <= 1: return 0.0 # Physics breaks
    K = 0.1
    
    # Initial Conditions
    r0 = 1e-3
    rho_c = 1.0 
    P_c = K * rho_c**gamma
    m_0 = (4/3) * np.pi * r0**3 * rho_c
    
    def system(r, y):
        P, m, phi = y
        if P <= 0: return [0, 0, 0]
        rho = (P / K)**(1/gamma)
        
        factor = -1.0 / (r * (r - 2*m))
        term1 = (rho + P)
        term2 = (m + 4 * np.pi * r**3 * P)
        dp_dr = factor * term1 * term2
        dm_dr = 4 * np.pi * r**2 * rho
        dphi_dr = -1.0 * dp_dr / (rho + P)
        return [dp_dr, dm_dr, dphi_dr]

    try:
        sol = solve_ivp(system, [r0, 50], [P_c, m_0, 0.0], method='LSODA', rtol=1e-5)
    except:
        return 0.0
        
    # Analyze Curve
    P = sol.y[0]
    mask = P > 1e-5
    r = sol.t[mask]
    m = sol.y[1][mask]
    P = P[mask]
    
    if len(r) < 10: return 0.0
    
    dphi_dr = (m + 4*np.pi*r**3*P) / (r * (r - 2*m))
    v = np.sqrt(r * dphi_dr)
    
    # Metric: Retention at Edge (How flat is it?)
    # We want V_edge / V_peak
    if np.max(v) == 0: return 0.0
    return v[-1] / np.max(v)

def run_sensitivity_test():
    print("🔬 SENSITIVITY CHECK: Is 1.878 special?")
    
    q_values = np.linspace(1.2, 3.0, 50)
    flatness = []
    
    print(f"   Testing range q=[1.2 ... 3.0]")
    
    for q in q_values:
        f = solve_for_flatness(q)
        flatness.append(f)
        
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(q_values, flatness, 'b-', lw=2)
    plt.axvline(1.878, color='r', linestyle='--', label='Your Value (1.878)')
    
    # Mark the Peak
    best_idx = np.argmax(flatness)
    best_q = q_values[best_idx]
    best_val = flatness[best_idx]
    
    plt.plot([best_q], [best_val], 'ro', label=f'Peak Flatness (q={best_q:.2f})')
    
    plt.xlabel('Tsallis Parameter (q)')
    plt.ylabel('Rotation Curve Flatness (V_edge / V_peak)')
    plt.title('Why 1.878? Searching for the Dark Matter Sweet Spot')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('sensitivity_check.png')
    
    print(f"✅ PLOT: 'sensitivity_check.png'")
    print(f"   Peak Flatness found at q = {best_q:.3f}")
    
    # Rigor Check
    if abs(best_q - 1.878) < 0.1:
        print("✅ RIGOR CONFIRMED: 1.878 is physically special (Optimizes Flatness).")
    else:
        print(f"❌ RIGOR FAILED: The best fit is actually q={best_q:.3f}.")
        print("   Your value works, but it isn't the unique solution.")

if __name__ == "__main__":
    run_sensitivity_test()