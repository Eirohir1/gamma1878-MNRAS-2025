import h5py
import numpy as np
import sys

# THE USS: GEMINCLAUDE - VACUUM STATE PROBE
# Mission: Resolve the 3.7% tension by escaping the Galactic CGM

def vacuum_grain_check(filenames):
    all_mags = []
    com = np.array([29355.95, 31028.46, 32500.65])
    
    print(f"🛰️ VACUUM PROBE: Scanning the Far Field (350-500 kpc)...")
    
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            pos, vel = f['PartType1/Coordinates'][:], f['PartType1/Velocities'][:]
            dist = np.sqrt(np.sum((pos - com)**2, axis=1))
            
            # THE VACUUM FILTER: Beyond the 'Dirty Water' of the CGM
            mask = (dist > 350) & (dist < 500)
            mag = np.sqrt(np.sum(vel[mask]**2, axis=1))
            all_mags.extend(mag)
    
    mags = np.array(all_mags)
    counts, bins = np.histogram(mags, bins=100, density=True)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Analyze the Grain at the lower velocity boundary
    v_mask = (bin_centers > 30) & (bin_centers < 120) & (counts > 0)
    log_x, log_y = np.log10(bin_centers[v_mask]), np.log10(counts[v_mask])
    slope, _ = np.polyfit(log_x, log_y, 1)
    
    detected = abs(slope)
    target = 1.878
    accuracy = 100 - abs((detected - target)/target)*100
    
    print(f"\n--- 🗺️ FAR-FIELD VACUUM REPORT ---")
    print(f"PARTICLES: {len(mags):,}")
    print(f"DEEP SPACE GRAIN: {detected:.4f}")
    print(f"CGM GRAIN (PREVIOUS): 1.8089")
    print(f"TARGET GRAIN: 1.8780")
    print(f"ACCURACY: {accuracy:.2f}%")
    
    if detected > 1.81:
        print("RESULT: GRAIN IS RECOVERING TOWARD 1.878. TENSION IS SYSTEMATIC.")
    else:
        print("RESULT: GRAIN IS CONSTANT. 1.81 IS THE PHYSICAL LIMIT.")

if __name__ == "__main__":
    vacuum_grain_check(sys.argv[1:])