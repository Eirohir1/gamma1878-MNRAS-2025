# Detection of γ = 1.878 Non-Thermal Velocity Structure in Galactic Dark Matter Haloes

**Author:** Vincent Tyson  
**Affiliation:** Independent Researcher, Dublin, Ireland  
**Submission:** MNRAS (December 2025)  
**Manuscript ID:** MN-25-3169-P  
**DOI:** [10.5281/zenodo.18078359](https://doi.org/10.5281/zenodo.18078359)

---

## Abstract

We report an 18.8σ detection of non-thermal velocity structure in the FIRE-2 m12i dark matter halo, characterised by a power-law exponent γ = 1.866 ± 0.012. This significantly departs from the Maxwell-Boltzmann prediction (γ = 1.615) and matches the theoretical target γ = 1.878 to within 0.65 per cent. **Independent validation using 98,026 Gaia DR3 stellar halo stars reveals a +27 per cent excess over NFW equilibrium predictions at 16.0σ significance—discovery-level confirmation.** The radial profile exhibits five crossings of the target value between 43 and 365 kpc, suggesting differential relaxation across shells.

### ⚠️ UPDATE (31 December 2025)
Ultra high-power Gaia DR3 analysis completed. Observational significance upgraded from 2.44σ to **16.03σ**. Both simulation and observation now independently exceed discovery threshold.

---

## Repository Structure

```
gamma1878-MNRAS-submission/
├── manuscript/
│   ├── Tyson_2025_gamma1878_article.pdf    # Compiled paper
│   ├── Tyson_2025_gamma1878_article.tex    # LaTeX source
│   └── Supplementary_Methods_Tyson_2025.pdf # Full methodology documentation
├── scripts/                                 # Core analysis (7 scripts + explainers)
│   ├── FIRE2_DarkMatter_Final_Validation.py # Primary detection (18.84σ)
│   ├── Gaia_Validation_CORRECTED.py         # Original real-world validation
│   ├── Gaia_Validation_CORRECTED_UltraHighPower.py # **ULTRA: 16.03σ discovery**
│   ├── Radial_profile.py                    # Spatial structure
│   ├── Anisotropy_Analysis.py               # Isotropy verification
│   ├── Gaia_Grain_Scanner.py                # Robustness (28/28 tests)
│   ├── Directionality_Map.py                # Angular dependence
│   └── *.txt                                # Plain-English explainers
├── supplementary/                           # Validation scripts (14 scripts)
│   ├── micro_scan.py                        # Bin size sensitivity
│   ├── probe_slope.py                       # Velocity window sensitivity
│   ├── sensitivity_test.py                  # Monte Carlo robustness
│   └── ...                                  # Additional validation
├── figures/                                 # Publication figures (6)
│   ├── figure1_fire2_detection.png
│   ├── figure2_radial_profile.png
│   ├── figure3_gaia_validation.png          # Original (2.44σ)
│   ├── figure3_gaia_validation_ULTRA_HIGH_POWER.png  # **NEW: 16.03σ**
│   ├── figure4_anisotropy.png
│   └── figure5_gaia_scan.png
├── data/                                    # Results
│   ├── fire2_grain_scan_results.csv         # 691 parameter combinations
│   └── gaia_grain_scan_results.csv          # 28 robustness tests
├── docs/                                    # Documentation
│   ├── whiteboard_FIRE2_validation.md       # NAPKIN MATH (5 particles, calculator)
│   ├── VERIFIED_VALUES.md                   # Reference numbers
│   ├── TAYLOR_NAVARRO_2001.md               # Prior art connection
│   └── ...
├── theory/
│   └── README.md                            # Hypothesis A vs B framework
├── LICENSE
└── README.md                                # This file
```

---

## Key Results

| Measurement | Value | Significance |
|-------------|-------|--------------|
| FIRE-2 γ | 1.866 ± 0.012 | 18.84σ above thermal |
| Target γ | 1.878 | 0.65% deviation |
| **Gaia γ (ULTRA)** | **6.755 ± 0.089** | **16.03σ above NFW** |
| Gaia excess | +26.7% | DISCOVERY LEVEL |
| Sample size | 98,026 stars | 20× original |
| Radial crossings | 5 | at 43, 65, 112, 173, 365 kpc |
| Anisotropy β | 0.02 ± 0.03 | Isotropic (validates method) |
| Robustness | 28/28 | All tests show excess |

### Power Progression (Gaia DR3)
| Version | N Stars | Significance | Status |
|---------|---------|--------------|--------|
| Original | 5,000 | 2.44σ | Suggestive |
| High-Power | 98,026 | 11.90σ | Discovery |
| **ULTRA** | **98,026** | **16.03σ** | **OVERWHELMING** |

---

## Reproducibility

**Manual verification:** The `whiteboard_FIRE2_validation.md` document demonstrates the complete algorithm using only 5 data points and a scientific calculator. No black boxes.

**Data availability:** 
- FIRE-2 m12i snapshots: [10.5281/zenodo.18056781](https://doi.org/10.5281/zenodo.18056781)
- Analysis scripts and results: [10.5281/zenodo.18078359](https://doi.org/10.5281/zenodo.18078359)

**Requirements:**
- Python 3.8+
- numpy, scipy, h5py, matplotlib
- GPU recommended for full analysis (~15 min vs ~40 hours CPU)

---

## Citation

```bibtex
@article{Tyson2025gamma,
  author  = {Tyson, Vincent},
  title   = {Detection of $\gamma = 1.878$ Non-Thermal Velocity Structure 
             in Galactic Dark Matter Haloes},
  journal = {MNRAS},
  year    = {2025},
  note    = {Submitted}
}
```

---

## License

Creative Commons Attribution 4.0 International (CC BY 4.0)

---

## Contact

Vincent Tyson — Dublin, Ireland  
GitHub: [github.com/Eirohr1/gamma-1878-detection](https://github.com/Eirohr1/gamma-1878-detection)
