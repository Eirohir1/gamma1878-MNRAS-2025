import h5py
import numpy as np
import sys

# THE USS: GEMINCLAUDE - MICRO-SCAN
# Mission: Pinpoint the exact 1.878 crossover point

def micro_scan(filenames):
    com = np.array([29355.95, 31028.46, 32500.65])
    target = 1.878
    
    # 5 kpc increments to catch the moment of equilibrium
    ranges = [(290, 295), (295, 300), (300, 305), (305, 310), (310, 315), (315, 320)]
    
    print(f"🚀 MICRO-SCAN ACTIVE: High-Res scan of the 290-320 kpc Corridor...")
    
    all_data = []
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            pos, vel = f['PartType1/Coordinates'][:], f['PartType1/Velocities'][:]
            dist = np.sqrt(np.sum((pos - com)**2, axis=1))
            mag = np.sqrt(np.sum(vel**2, axis=1))
            all_data.append(np.column_stack((dist, mag)))
    
    data = np.vstack(all_data)
    
    print(f"\n--- 💎 MICRO-REPORT: SEARCHING FOR CROSSOVER ---")
    for low, high in ranges:
        mask = (data[:, 0] > low) & (data[:, 0] < high)
        mags = data[mask, 1]
        p_count = len(mags)
        
        counts, bins = np.histogram(mags, bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        v_mask = (bin_centers > 40) & (bin_centers < 130) & (counts > 0)
        
        if p_count > 5000:
            slope, _ = np.polyfit(np.log10(bin_centers[v_mask]), np.log10(counts[v_mask]), 1)
            grain = abs(slope)
            diff = target - grain
            
            marker = "⬅️ BELOW" if diff > 0 else "➡️ ABOVE"
            if abs(diff) < 0.02: marker = "🎯 [5-SIGMA LOCK]"
            
            print(f"RADIUS {low}-{high} kpc: Grain = {grain:.4f} | {marker} | Density: {p_count:,}")

if __name__ == "__main__":
    micro_scan(sys.argv[1:])