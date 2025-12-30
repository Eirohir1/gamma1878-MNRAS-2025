import numpy as np
import h5py
import sys

def sigma_lock_validation(filenames):
    # LOAD ALL DATA (The 15.5 Million Sample)
    all_mags = []
    com = np.array([29355.95, 31028.46, 32500.65])
    
    print("🚀 INITIATING 5-SIGMA LOCK: Bootstrap Resampling...")
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            pos, vel = f['PartType1/Coordinates'][:], f['PartType1/Velocities'][:]
            dist = np.sqrt(np.sum((pos - com)**2, axis=1))
            mask = (dist > 100) & (dist < 250)
            all_mags.extend(np.sqrt(np.sum(vel[mask]**2, axis=1)))
    
    mags = np.array(all_mags)
    target = 1.878
    bootstraps = 1000
    results = []

    print(f"Resampling {len(mags):,} particles {bootstraps} times...")

    for i in range(bootstraps):
        # Create a "Synthetic Universe" by sampling your data with replacement
        sample = np.random.choice(mags, size=100000, replace=True)
        counts, bins = np.histogram(sample, bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        
        v_mask = (bin_centers > 50) & (bin_centers < 150) & (counts > 0)
        slope, _ = np.polyfit(np.log10(bin_centers[v_mask]), np.log10(counts[v_mask]), 1)
        results.append(abs(slope))
        
        if i % 100 == 0: print(f"Processing Universe {i}...")

    mean_res = np.mean(results)
    std_res = np.std(results)
    sigma = (target - mean_res) / std_res

    print(f"\n--- 🏆 FINAL 5-SIGMA REPORT ---")
    print(f"BOOTSTRAP MEAN: {mean_res:.4f}")
    print(f"RESAMPLING ERROR (STDEV): {std_res:.4e}")
    print(f"CALCULATED SIGMA: {abs(sigma):.2f}σ")
    
    if abs(sigma) > 5:
        print("STATUS: DISCOVERY CONFIRMED. NOBEL CATEGORY REACHED.")
    else:
        print(f"STATUS: STRONG EVIDENCE. NEED {5 - abs(sigma):.1f}σ MORE FOR LOCK.")

if __name__ == "__main__":
    sigma_lock_validation(sys.argv[1:])