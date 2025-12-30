import h5py
import numpy as np
import sys

def find_directional_grain(filenames):
    data = {'X': [], 'Y': [], 'Z': []}
    com = np.array([29355.95, 31028.46, 32500.65])
    
    print(f"🛰️ DIRECTIONAL PROBE ACTIVE: Scanning for Manifold Seams...")
    
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            vel = f['PartType1/Velocities'][:]
            pos = f['PartType1/Coordinates'][:]
            dist = np.sqrt(np.sum((pos - com)**2, axis=1))
            
            # Filter for our High-Res "Sweet Spot"
            mask = (dist > 100) & (dist < 200)
            valid_vel = vel[mask]
            
            # Split velocities into dominant directional components
            # (Which axis is the particle pushing against most?)
            for i, axis in enumerate(['X', 'Y', 'Z']):
                # We look at the magnitude of velocity along each specific axis
                mag_axis = np.abs(valid_vel[:, i])
                data[axis].extend(mag_axis)

    target_grain = 1.878
    print(f"\n--- 🗺️ AXIS ANALYSIS: DIRECTIONAL STIFFNESS ---")
    
    for axis in ['X', 'Y', 'Z']:
        mags = np.array(data[axis])
        counts, bins = np.histogram(mags, bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        
        # Grain Filter: 50-150 km/s (The Interaction Zone)
        mask = (bin_centers > 50) & (bin_centers < 150) & (counts > 0)
        log_x = np.log10(bin_centers[mask])
        log_y = np.log10(counts[mask])
        
        slope, _ = np.polyfit(log_x, log_y, 1)
        detected_grain = abs(slope)
        accuracy = 100 - abs((detected_grain - target_grain)/target_grain)*100
        
        print(f"AXIS {axis}: Grain = {detected_grain:.4f} | Accuracy = {accuracy:.2f}%")

if __name__ == "__main__":
    find_directional_grain(sys.argv[1:])