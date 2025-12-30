import numpy as np
import matplotlib.pyplot as plt

def generate_skymap():
    print("Generating QROCODILE Prediction Map...")

    # --- 1. DEFINE COORDINATES (ICRS / J2000) ---
    
    # Target A: Standard Dark Matter "Wind"
    # Direction: Solar Apex (Motion of Sun relative to LSR)
    # Approx: RA 18h 03m, Dec +30d (Near Vega/Hercules/Cygnus border)
    # RA = 270.75 deg, Dec = +30.0 deg
    name_A = "Standard Model\n(Solar Wind)"
    ra_A_deg = 270.75
    dec_A_deg = 30.0

    # Target B: Polylensic "Grain"
    # Direction: Supergalactic North Pole (SGP)
    # Defined by de Vaucouleurs: RA 286.32, Dec +15.7
    name_B = "Polylensic Grain\n(Supergalactic Plane)"
    ra_B_deg = 286.32
    dec_B_deg = 15.7

    # --- 2. CALCULATE SEPARATION ANGLE (Spherical Trig) ---
    # Convert to Radians
    ra_A = np.radians(ra_A_deg)
    dec_A = np.radians(dec_A_deg)
    ra_B = np.radians(ra_B_deg)
    dec_B = np.radians(dec_B_deg)

    # Great Circle Distance Formula (Haversine or Dot Product)
    # cos(theta) = sin(d1)sin(d2) + cos(d1)cos(d2)cos(ra1-ra2)
    cos_theta = np.sin(dec_A)*np.sin(dec_B) + \
                np.cos(dec_A)*np.cos(dec_B)*np.cos(ra_A - ra_B)
    
    theta_rad = np.arccos(cos_theta)
    theta_deg = np.degrees(theta_rad)

    print(f"--- GEOMETRIC PREDICTION ---")
    print(f"Standard Vector: RA {ra_A_deg:.1f}, Dec {dec_A_deg:.1f}")
    print(f"Fractal Vector:  RA {ra_B_deg:.1f}, Dec {dec_B_deg:.1f}")
    print(f"---------------------------")
    print(f"SEPARATION ANGLE: {theta_deg:.2f} degrees")
    
    # Check the "Golden Angle of Quadrant" Hypothesis
    # 90 degrees * (1 - 1/phi) approx 34.4 degrees? 
    # Or just raw degrees. 
    # Let's just output the raw angle for the paper.

    # --- 3. PLOTTING (Mollweide Projection) ---
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection='mollweide')

    # Matplotlib Mollweide expects:
    # x in radians [-pi, pi] (Longitude)
    # y in radians [-pi/2, pi/2] (Latitude)
    # Note: Astronomy sky maps usually have RA increasing to the LEFT.
    # We will map RA [0, 360] to x [-pi, pi] carefully.
    # Standard: x = - (RA - 180) * pi/180  (flips direction and centers 180)
    
    def rad_to_plot(ra_radians):
        # Shift so 180 (pi) is at center (0)
        # And multiply by -1 to make East left (Astronomy convention)
        return -(ra_radians - np.pi)

    x_A = rad_to_plot(ra_A)
    y_A = dec_A
    x_B = rad_to_plot(ra_B)
    y_B = dec_B

    # Plot Background Grid
    ax.grid(True, color='gray', alpha=0.3)
    ax.set_facecolor('#f0f0f5')

    # Plot Points
    ax.scatter([x_A], [y_A], c='blue', s=250, marker='*', label='Standard (Cygnus)')
    ax.scatter([x_B], [y_B], c='red', s=250, marker='o', label='Polylensic (SGP)')

    # Draw Connecting Line (Geodesic)
    # We interpolate points along the great circle for a smooth curved line
    t = np.linspace(0, 1, 100)
    # Vector Slerp or just linear interp in angular space?
    # Simple approximation for visual line is fine for small separation, 
    # but let's do it properly-ish by just plotting the segment.
    # Since they are close, a straight line in lat/lon is "okay" visually, 
    # but let's just draw a dashed line.
    ax.plot([x_A, x_B], [y_A, y_B], 'k--', linewidth=2, alpha=0.6)

    # Annotate
    mid_x = (x_A + x_B) / 2
    mid_y = (y_A + y_B) / 2
    ax.text(mid_x, mid_y + 0.1, f"Gap: {theta_deg:.1f}°", 
            ha='center', fontsize=12, fontweight='bold', color='purple')

    ax.text(x_A, y_A - 0.15, name_A, ha='center', va='top', color='blue', fontsize=10)
    ax.text(x_B, y_B - 0.15, name_B, ha='center', va='top', color='red', fontsize=10)

    # Titles
    plt.title(f"The QROCODILE Trap: Directional Anomaly Prediction\nSeparation Angle = {theta_deg:.1f}°", fontsize=14, fontweight='bold')
    
    # Legend location
    plt.legend(loc='lower right')
    
    # Save
    plt.tight_layout()
    plt.savefig('directionality_skymap.png', dpi=300)
    print(">> Plot saved to 'directionality_skymap.png'")

    # --- 4. TEXT REPORT ---
    with open("qrocodile_prediction.txt", "w") as f:
        f.write("--- QROCODILE DIRECTIONALITY PREDICTION ---\n")
        f.write(f"Standard Target (Solar Wind): RA {ra_A_deg:.2f}, Dec {dec_A_deg:.2f}\n")
        f.write(f"Polylensic Target (SGP):      RA {ra_B_deg:.2f}, Dec {dec_B_deg:.2f}\n")
        f.write(f"-------------------------------------------\n")
        f.write(f"PREDICTED ANGULAR OFFSET:     {theta_deg:.2f} degrees\n")
        f.write(f"-------------------------------------------\n")
        f.write("Implication: If QROCODILE detects a signal offset by ~20-25 deg\n")
        f.write("from Cygnus towards the Equator/SGP, it confirms the Fractal Grain.\n")

if __name__ == "__main__":
    generate_skymap()