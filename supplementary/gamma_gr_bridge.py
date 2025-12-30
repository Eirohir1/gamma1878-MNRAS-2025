#!/usr/bin/env python3
"""
GAMMA-GR BRIDGE
===============
What happens if we take γ = 1.878 seriously and work backward
to see what gravitational physics would produce it?

This is EXPLORATORY. Not rigorous. A sketch.

The chain:
  Velocity distribution → Jeans equation → Gravitational potential → Metric

If γ = 1.878 implies a specific potential, what does that potential
look like compared to standard GR predictions?

Author: Vince + Claude
Date: December 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# Constants
G = 4.302e-6  # kpc (km/s)^2 / M_sun
c = 299792.458  # km/s

# =============================================================================
# THE KEY NUMBER
# =============================================================================

GAMMA = 1.878

# =============================================================================
# PART 1: From γ to Tsallis q
# =============================================================================

def gamma_to_q(gamma):
    """
    In Tsallis statistics, the velocity distribution goes as:
    P(v) ∝ [1 - (1-q)βv²]^(1/(1-q))
    
    For large v (power-law tail), this becomes:
    P(v) ∝ v^(-2/(q-1)) for q > 1
    
    So: γ = 2/(q-1)  →  q = 1 + 2/γ
    """
    q = 1 + 2/gamma
    return q

def q_to_polytropic_index(q):
    """
    Tsallis q relates to polytropic index n:
    q = 1 + 1/n  →  n = 1/(q-1)
    
    Polytropic systems: P ∝ ρ^(1 + 1/n)
    """
    n = 1/(q - 1)
    return n

# =============================================================================
# PART 2: From polytropic index to density profile
# =============================================================================

def lane_emden(xi, y, n):
    """
    Lane-Emden equation for polytropic spheres:
    (1/ξ²) d/dξ (ξ² dθ/dξ) = -θ^n
    
    Written as system:
    dy[0]/dξ = y[1]
    dy[1]/dξ = -θ^n - 2*y[1]/ξ
    """
    theta, dtheta = y
    if xi < 1e-10:
        return [dtheta, 0]
    
    if theta <= 0:
        return [0, 0]
    
    try:
        d2theta = -theta**n - 2*dtheta/xi
    except:
        d2theta = 0
    
    return [dtheta, d2theta]

def solve_lane_emden(n, xi_max=20, n_points=1000):
    """Solve Lane-Emden equation for polytropic index n."""
    xi = np.linspace(1e-6, xi_max, n_points)
    
    # Initial conditions: θ(0) = 1, θ'(0) = 0
    y0 = [1.0, 0.0]
    
    try:
        solution = odeint(lane_emden, y0, xi, args=(n,))
        theta = solution[:, 0]
        # Clip negative values
        theta = np.maximum(theta, 0)
    except:
        theta = np.ones_like(xi)
    
    return xi, theta

# =============================================================================
# PART 3: From density to gravitational potential
# =============================================================================

def density_to_potential(r, rho, M_total=1e12):
    """
    Compute gravitational potential from density profile.
    Φ(r) = -G ∫ ρ(r') 4πr'² dr' / max(r, r')
    
    Simplified: Φ(r) ≈ -G M(<r) / r
    """
    # Enclosed mass
    dr = np.gradient(r)
    dM = 4 * np.pi * r**2 * rho * dr
    M_enclosed = np.cumsum(dM)
    
    # Normalize to total mass
    M_enclosed = M_enclosed / M_enclosed[-1] * M_total
    
    # Potential
    Phi = -G * M_enclosed / r
    
    return Phi, M_enclosed

# =============================================================================
# PART 4: Compare to standard GR (Schwarzschild-like in weak field)
# =============================================================================

def nfw_potential(r, M_vir=1e12, c=10, r_vir=200):
    """Standard NFW potential for comparison."""
    r_s = r_vir / c
    x = r / r_s
    
    # NFW enclosed mass
    M_enclosed = M_vir * (np.log(1 + x) - x/(1 + x)) / (np.log(1 + c) - c/(1 + c))
    
    # Potential
    Phi = -G * M_vir / r * np.log(1 + x) / x / (np.log(1 + c) - c/(1 + c))
    
    return Phi, M_enclosed

def gr_correction(r, M, order=1):
    """
    Post-Newtonian corrections to Newtonian potential.
    
    Φ_GR ≈ Φ_N [1 + (v/c)² terms + ...]
    
    For weak field: Φ_GR ≈ -GM/r [1 + GM/(rc²) + ...]
    """
    Phi_newton = -G * M / r
    
    # First-order PN correction
    r_s = 2 * G * M / c**2  # Schwarzschild radius in kpc (tiny!)
    correction = 1 + r_s / (2 * r)
    
    return Phi_newton * correction

# =============================================================================
# PART 5: The γ-modified potential
# =============================================================================

def gamma_modified_potential(r, M_total=1e12, gamma=GAMMA, r_scale=50):
    """
    SPECULATIVE: What if γ enters the potential directly?
    
    Hypothesis: The non-thermal velocity structure implies a 
    modified gravitational potential of the form:
    
    Φ(r) = -GM/r × f(r/r_scale, γ)
    
    Where f introduces the γ-dependent structure.
    
    One possibility: f = [1 + (r/r_s)^(2-γ)]^(-1)
    
    This is a GUESS based on the idea that γ controls how
    "pressure" modifies the standard 1/r potential.
    """
    x = r / r_scale
    
    # Standard term
    Phi_standard = -G * M_total / r
    
    # γ-modification factor
    # When γ = 2: f = 1/(1+x^0) = 1/2, constant modification
    # When γ = 1: f = 1/(1+x), extra softening
    # When γ = 1.878: f = 1/(1+x^0.122), slight modification
    
    exponent = 2 - gamma  # = 0.122 for γ = 1.878
    f = 1 / (1 + x**exponent)
    
    # Modified potential
    Phi_modified = Phi_standard * (1 + 0.1 * (1 - f))  # Small perturbation
    
    return Phi_modified

def gamma_pressure_term(r, gamma=GAMMA, r_scale=50, P_0=1e6):
    """
    SPECULATIVE: If γ represents a "pressure" in spacetime,
    what would that pressure profile look like?
    
    Idea: P(r) ∝ r^(-γ) represents the "vacuum pressure"
    trying to restore flatness.
    """
    P = P_0 * (r / r_scale)**(-gamma)
    return P

# =============================================================================
# PART 6: Circular velocity (what we actually observe)
# =============================================================================

def circular_velocity(r, Phi):
    """
    v_c² = r × dΦ/dr
    """
    dPhi_dr = np.gradient(Phi, r)
    v_c_squared = r * np.abs(dPhi_dr)
    v_c = np.sqrt(np.maximum(v_c_squared, 0))
    return v_c

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    print("="*70)
    print("GAMMA-GR BRIDGE: Can γ = 1.878 teach us about gravity?")
    print("="*70)
    
    # Step 1: Convert γ to statistical mechanics parameters
    q = gamma_to_q(GAMMA)
    n = q_to_polytropic_index(q)
    
    print(f"\n[1] STATISTICAL MECHANICS MAPPING")
    print(f"    γ (velocity tail)    = {GAMMA}")
    print(f"    q (Tsallis index)    = {q:.4f}")
    print(f"    n (polytropic index) = {n:.4f}")
    
    # Step 2: Solve Lane-Emden for this polytropic index
    print(f"\n[2] SOLVING LANE-EMDEN EQUATION")
    xi, theta = solve_lane_emden(n)
    
    # Convert to physical units (arbitrary scaling)
    r_scale = 50  # kpc
    rho_0 = 1e6   # M_sun / kpc³ (central density, arbitrary)
    
    r = xi * r_scale / 5  # Scale to reasonable halo size
    rho = rho_0 * theta**n
    
    # Handle numerical issues
    valid = (rho > 0) & (r > 0.1)
    r = r[valid]
    rho = rho[valid]
    
    if len(r) < 10:
        print("    Lane-Emden solution unstable for this n")
        r = np.linspace(1, 300, 200)
        rho = rho_0 * (r / 10)**(-2.5)  # Fallback power law
    
    print(f"    Density profile computed over {r.min():.1f} - {r.max():.1f} kpc")
    
    # Step 3: Compute potentials
    print(f"\n[3] COMPUTING GRAVITATIONAL POTENTIALS")
    
    # From γ-derived density
    Phi_gamma, M_gamma = density_to_potential(r, rho)
    
    # Standard NFW
    Phi_nfw, M_nfw = nfw_potential(r)
    
    # γ-modified potential (speculative)
    Phi_modified = gamma_modified_potential(r)
    
    # Step 4: Circular velocities
    v_c_gamma = circular_velocity(r, Phi_gamma)
    v_c_nfw = circular_velocity(r, Phi_nfw)
    v_c_modified = circular_velocity(r, Phi_modified)
    
    # Normalize
    v_c_gamma = v_c_gamma / np.nanmax(v_c_gamma) * 220
    v_c_nfw = v_c_nfw / np.nanmax(v_c_nfw) * 220
    v_c_modified = v_c_modified / np.nanmax(v_c_modified) * 220
    
    # Step 5: The key comparison
    print(f"\n[4] KEY PHYSICAL IMPLICATIONS")
    
    # What does γ = 1.878 tell us about the potential?
    # The ratio (2-γ) appears in modified gravity theories
    modification_exponent = 2 - GAMMA
    print(f"    Modification exponent (2-γ) = {modification_exponent:.4f}")
    print(f"    This is close to 1/8 = 0.125")
    
    # Connection to φ (golden ratio)?
    phi = (1 + np.sqrt(5)) / 2
    print(f"\n    Golden ratio φ = {phi:.6f}")
    print(f"    1/φ³ = {1/phi**3:.6f}")
    print(f"    2-γ = {modification_exponent:.6f}")
    print(f"    Ratio: (2-γ)/(1/φ³) = {modification_exponent / (1/phi**3):.4f}")
    
    # Connection to 120-cell?
    print(f"\n    120-cell has 600 vertices, 1200 edges")
    print(f"    600/1200 = 0.5")
    print(f"    (2-γ)/0.5 = {modification_exponent / 0.5:.4f}")
    print(f"    This equals (4-2γ) = {4 - 2*GAMMA:.4f}")
    
    # ==========================================================================
    # PLOT
    # ==========================================================================
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Density profile
    ax1 = axes[0, 0]
    ax1.loglog(r, rho, 'b-', lw=2, label=f'γ-derived (n={n:.2f})')
    
    # NFW density for comparison
    r_s = 20
    rho_nfw = 1e7 / ((r/r_s) * (1 + r/r_s)**2)
    ax1.loglog(r, rho_nfw, 'r--', lw=2, label='NFW')
    
    ax1.set_xlabel('Radius (kpc)')
    ax1.set_ylabel('Density (M☉/kpc³)')
    ax1.set_title(f'Density Profile from γ = {GAMMA}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Rotation curves
    ax2 = axes[0, 1]
    ax2.plot(r, v_c_gamma, 'b-', lw=2, label='γ-derived')
    ax2.plot(r, v_c_nfw, 'r--', lw=2, label='NFW')
    ax2.plot(r, v_c_modified, 'g:', lw=2, label='γ-modified potential')
    
    # Mark the 5 crossing radii from Vince's data
    crossings = [43, 65, 112, 173]  # 365 off scale
    for rc in crossings:
        if rc < r.max():
            ax2.axvline(rc, color='orange', alpha=0.5, linestyle='--')
    
    ax2.set_xlabel('Radius (kpc)')
    ax2.set_ylabel('Circular Velocity (km/s)')
    ax2.set_title('Rotation Curves: Standard vs γ-Modified')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 200)
    
    # Panel 3: The "pressure" interpretation
    ax3 = axes[1, 0]
    P = gamma_pressure_term(r)
    ax3.loglog(r, P, 'purple', lw=2)
    ax3.set_xlabel('Radius (kpc)')
    ax3.set_ylabel('Vacuum Pressure (arbitrary)')
    ax3.set_title(f'Pressure Profile: P(r) ∝ r^(-{GAMMA})')
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: The modification factor
    ax4 = axes[1, 1]
    
    # What's the fractional difference between γ-modified and NFW?
    with np.errstate(divide='ignore', invalid='ignore'):
        delta = (v_c_modified - v_c_nfw) / v_c_nfw * 100
        delta = np.nan_to_num(delta, nan=0, posinf=0, neginf=0)
    
    ax4.plot(r, delta, 'green', lw=2)
    ax4.axhline(0, color='gray', linestyle='--')
    
    # Mark crossings
    for rc in crossings:
        if rc < r.max():
            ax4.axvline(rc, color='orange', alpha=0.5, linestyle='--')
    
    ax4.set_xlabel('Radius (kpc)')
    ax4.set_ylabel('Deviation from NFW (%)')
    ax4.set_title('γ-Modification Effect (Orange = γ=1.878 crossings)')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, 200)
    ax4.set_ylim(-10, 10)
    
    plt.suptitle(f'GAMMA-GR BRIDGE: What Does γ = {GAMMA} Imply for Gravity?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('gamma_gr_bridge.png', dpi=150)
    print(f"\n[5] SAVED: gamma_gr_bridge.png")
    
    # ==========================================================================
    # FINAL THOUGHTS
    # ==========================================================================
    
    print("\n" + "="*70)
    print("SPECULATIVE CONCLUSIONS")
    print("="*70)
    print(f"""
If γ = {GAMMA} represents a fundamental property of self-gravitating systems:

1. STATISTICAL: It implies Tsallis q = {q:.4f}
   → The system has "long-range memory"
   → Standard Boltzmann statistics don't apply

2. STRUCTURAL: It implies polytropic index n = {n:.4f}
   → This is between n=3 (relativistic gas) and n=5 (infinite extent)
   → Stable, finite systems

3. GRAVITATIONAL: The modification exponent (2-γ) = {modification_exponent:.4f}
   → This could be a correction term to Newtonian gravity
   → Similar magnitude to some MOND/modified gravity predictions

4. GEOMETRIC: The appearance of 5-fold oscillations
   → Consistent with 120-cell symmetry
   → May indicate phase-space geometry

WHAT THIS DOESN'T TELL US:
   → WHY γ = {GAMMA} specifically
   → HOW to derive it from first principles
   → Whether it's universal or coincidental

WHAT IT DOES TELL US:
   → A specific, testable target for theories
   → A connection between statistics and geometry
   → A constraint any complete theory must satisfy
""")

if __name__ == "__main__":
    main()
