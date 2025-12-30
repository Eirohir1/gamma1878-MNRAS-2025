import h5py
import numpy as np
import sys

# THE USS: GEMINCLAUDE - EQUILIBRIUM HUNT (Ver. 2.0)
# Mission: Pinpoint the 1.878 'Goldilocks' Radius with Particle Density Check

def find_equilibrium(filenames):
    com = np.array([29355.95, 31028.46, 32500.65])
    target = 1.878
    
    # 20 kpc increments to find the cross-over point
    ranges = [(250, 270), (270, 290), (290, 310), (310, 330), (330, 350)]
    
    print(f"🛰️ EQUILIBRIUM HUNT: Scanning the Transition Zone (250-350 kpc)...")
    
    all_data = []
    for filename in filenames:
        try:
            with h5py.File(filename, 'r') as f:
                pos, vel = f['PartType1/Coordinates'][:], f['PartType1/Velocities'][:]
                dist = np.sqrt(np.sum((pos - com)**2, axis=1))
                mag = np.sqrt(np.sum(vel**2, axis=1))
                all_data.append(np.column_stack((dist, mag)))
        except Exception as e:
            print(f"Error loading part: {e}")
    
    data = np.vstack(all_data)
    
    print(f"\n--- 🗺️ RADIUS SCAN: SEARCHING FOR 1.878 EQUILIBRIUM ---")
    for low, high in ranges:
        mask = (data[:, 0] > low) & (data[:, 0] < high)
        mags = data[mask, 1]
        p_count = len(mags)
        
        counts, bins = np.histogram(mags, bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        
        # Velocity window optimized for transition zone
        v_mask = (bin_centers > 40) & (bin_centers < 130) & (counts > 0)
        
        if p_count > 10000 and len(bin_centers[v_mask]) > 5:
            slope, _ = np.polyfit(np.log10(bin_centers[v_mask]), np.log10(counts[v_mask]), 1)
            grain = abs(slope)
            diff = abs(grain - target)
            
            # Identify the "Lock"
            status = "🎯 [EQUILIBRIUM LOCK]" if diff < 0.05 else ""
            print(f"SHELL {low}-{high} kpc: Grain = {grain:.4f} | Particles = {p_count:,} | Delta = {diff:.4f} {status}")
        else:
            print(f"SHELL {low}-{high} kpc: Low data density ({p_count:,} particles)")

if __name__ == "__main__":
    find_equilibrium(sys.argv[1:])