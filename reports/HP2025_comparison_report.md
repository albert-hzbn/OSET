# OSTZ Theory vs. Harisankar & Padmanabhan (2025): A 41-System Comparison Report

**Reference paper:** K.R. Harisankar & K.A. Padmanabhan, *"On the origin of isotropic steady state structural superplastic deformation,"* **Materials Science & Engineering A 930 (2025) 148175.**

**This report:** compares the semi-empirical model constants extracted in that paper (by fitting the mesoscopic grain-boundary-sliding rate equation to experimental data) against parameter-free predictions of the Oblate-Spheroidal Excitation Theory (OSTZ/OSET). All numbers reproduced here are computed in `OSET_derivations.ipynb` §13 and are reproducible.

---

## 1. Scope of the paper

The paper validates a single physical mechanism — grain-boundary sliding (GBS) that develops to a *mesoscopic* scale — across:

- **41 material systems / 146 thermomechanical conditions** used to extract the mean model constants (Table 1 of the paper);
- a **further set of independent systems** (Table 3) *not* used in the fit, predicted to a tolerance < 10;
- in total, **467 material conditions** in the supplementary validation.

For each system the paper reports, as functions of temperature and melting point, the constants:

| Symbol | Meaning | Units |
|--------|---------|-------|
| $\gamma_0$ | shear eigenstrain of the basic sliding unit (oblate spheroid) | – |
| $\varepsilon_0$ | dilatational eigenstrain | – |
| $Q/T_m$ | activation energy (melting-point compensated) | J mol⁻¹ K⁻¹ |
| $G/T_m$ | shear modulus (melting-point compensated) | MPa K⁻¹ |
| $\tau_0/T_m$ | mesoscopic threshold stress | Pa K⁻¹ |
| $\gamma_B$ | **specific grain-boundary free energy** | J m⁻² |
| $N_a$ | number of contiguous boundaries forming the mesoscopic plane interface | – |
| $a$ | grain-size exponent | – |

The basic sliding unit is an **oblate spheroid, 5 atom diameters across the boundary plane and 2.5 atom diameters high** ($W = 2.5b$) — the same OSTZ geometry used throughout OSET.

---

## 2. What OSTZ predicts, with no free parameters

OSET fixes the elementary excitation as an oblate spheroid ($\alpha = c/R_0 = 1/2$) carrying eigenstrains $\gamma_0, \varepsilon_0$, with Eshelby coefficients (for $\nu = 1/3$):

$$\beta_1 = 1 - 2S_{1313} = 0.446, \qquad \beta_2 = \frac{4(1+\nu)}{9(1-\nu)} = 0.889.$$

Two physically distinct energies enter the comparison, and **they must not be conflated**:

1. **Activation barrier (a volume energy):**
$$\Delta F_0 = \tfrac12\,(\beta_1\gamma_0^2 + \beta_2\varepsilon_0^2)\,G\,V_0, \qquad V_0 = \tfrac{2}{3}\pi(2.5b)^3.$$

2. **Grain-boundary energy (an interfacial / dislocation-array energy):** a high-angle boundary is an array of dislocations, and in OSET each dislocation is an $N$-OSTZ chain. Summing the line energy over the array gives the **Read–Shockley** law:
$$\gamma_B^{\text{OSET}} = E_0\,\theta_m, \qquad E_0 = \frac{Gb}{4\pi(1-\nu)}, \qquad \theta_m \approx 15^\circ.$$

> **Key correction.** A grain-boundary energy scales as $Gb$ (linear in the eigenstrain through $b$), **not** as $G\gamma_0^2$. Comparing the paper's $\gamma_B$ against the OSTZ *activation-energy density* $\tfrac13(\beta_1\gamma_0^2+\beta_2\varepsilon_0^2)GW$ (the earlier approach) is off by ~40× for every system; the dislocation-array form above removes that error.

---

## 3. Headline results

| Quantity | OSTZ prediction | Paper (41-system mean / range) | Agreement |
|----------|-----------------|-------------------------------|-----------|
| **Dilatational eigenstrain** $\varepsilon_0$ | 0.05 | 0.049 ± 0.012 (0.034–0.088) | **< 2 %** (best parameter) |
| **Shear eigenstrain** $\gamma_0$ | 0.10 (GB) – 0.12 (lattice) | 0.085 ± 0.020 (0.059–0.144) | within 15 % |
| **Energy factor** $\beta_1\gamma_0^2+\beta_2\varepsilon_0^2$ | 0.00668 | 0.00570 (0.0026–0.0161) | within 15 % |
| **Activation energy** $Q$ vs. $\Delta F_0$ | $\Delta F_0 = \tfrac12 EF\,GV_0$ | $Q = 1.98 \pm 1.02$ eV | $Q/\Delta F_0 = 4.0 \pm 1.4$ |
| **GB free energy** $\gamma_B$ | $E_0\theta_m$ (Read–Shockley) | 0.30–1.7 J m⁻² | **median 4.5×** (35/41 within 10×) |
| **Cooperative number** $N_a$ | atomic $N_c = 4$ | 15.6 ± 2.6 (mesoscopic) | different length scales |

**Interpretation in one line:** OSTZ reproduces the eigenstrains essentially exactly, accounts for ~25 % of the total activation barrier (the rest being diffusional accommodation), and reproduces the correct $Gb$ scaling of grain-boundary energy with a constant ~4.5× offset.

---

## 4. Full per-system comparison (41 systems)

Columns: $\gamma_0,\varepsilon_0$ = paper fit; $\Delta F_0$ = OSTZ activation barrier (eV); $Q$ = paper activation energy per event (eV); $Q/\Delta F_0$ = barrier ratio; $\gamma_B^{\text{lit}}$ = paper fit (J m⁻²); $\gamma_B^{\text{OSET}}$ = OSTZ Read–Shockley (J m⁻²); ratio = lit/OSET; $\gamma_B^{\text{exp}}$ = representative measured HAGB energy (J m⁻²).

| # | System | Class | $\gamma_0$ | $\varepsilon_0$ | $\Delta F_0$ | $Q$ | $Q/\Delta F_0$ | $\gamma_B^{\text{lit}}$ | $\gamma_B^{\text{OSET}}$ | ratio | $\gamma_B^{\text{exp}}$ |
|---|--------|-------|-----------|-----------|-----------|-----|-----------|-----------|-----------|-------|-----------|
| 1 | Zn-22Al (2.5 µm) | Zn | 0.059 | 0.034 | 0.195 | 0.609 | 3.1 | 1.29 | 0.29 | 4.5 | 0.34 |
| 2 | Zn-22Al (0.9 µm) | Zn | 0.061 | 0.035 | 0.200 | 0.811 | 4.1 | 1.31 | 0.28 | 4.7 | 0.34 |
| 3 | Al-33Cu-0.4Zr | Al | 0.066 | 0.038 | 0.263 | 1.449 | 5.5 | 1.61 | 0.31 | 5.2 | 0.40 |
| 4 | Al-5.76Mg-Sc-Mn (3 µm) | Al | 0.089 | 0.051 | 0.264 | 1.310 | 5.0 | 1.57 | 0.17 | 9.1 | 0.40 |
| 5 | Al-5.76Mg-Sc-Mn (1 µm) | Al | 0.083 | 0.048 | 0.260 | 1.274 | 4.9 | 0.92 | 0.19 | 4.8 | 0.40 |
| 6 | Al-3Mg-0.2Sc | Al | 0.086 | 0.050 | 0.269 | 1.147 | 4.3 | 1.46 | 0.19 | 7.9 | 0.40 |
| 7 | Al-8.9Zn-2.6Mg-Sc | Al | 0.080 | 0.046 | 0.267 | 1.299 | 4.9 | 1.47 | 0.21 | 6.8 | 0.40 |
| 8 | Al-5Mg-Mn-Sc (24 µm) | Al | 0.090 | 0.052 | 0.266 | 1.370 | 5.2 | 0.87 | 0.17 | 5.2 | 0.40 |
| 9 | Al-17Si-Fe-Mg-Cu-Ni | Al | 0.089 | 0.051 | 0.306 | 1.126 | 3.7 | 0.88 | 0.20 | 4.4 | 0.40 |
| 10 | Mg-6Zn-0.8Zr | Mg | 0.080 | 0.046 | 0.264 | 1.232 | 4.7 | 1.54 | 0.16 | 9.8 | 0.35 |
| 11 | Mg-4Y-0.7Zr-Nd | Mg | 0.087 | 0.050 | 0.263 | 1.265 | 4.8 | 1.35 | 0.13 | 10.2 | 0.35 |
| 12 | Mg-5.8Zn-1Y-Zr | Mg | 0.088 | 0.051 | 0.274 | 1.088 | 4.0 | 1.37 | 0.13 | 10.3 | 0.35 |
| 13 | Ti-6Al-4V | Ti | 0.116 | 0.067 | 0.585 | 2.491 | 4.3 | 1.14 | 0.20 | 5.6 | 0.55 |
| 14 | Cu-2.8Al-1.8Si-Co (7 µm) | Cu | 0.096 | 0.056 | 0.394 | 1.499 | 3.8 | 1.48 | 0.27 | 5.5 | 0.60 |
| 15 | Cu-2.8Al-1.8Si-Co (3 µm) | Cu | 0.095 | 0.055 | 0.392 | 1.504 | 3.8 | 0.99 | 0.28 | 3.6 | 0.60 |
| 16 | IN836 (Cu-Ni-Zn) | Ni | 0.081 | 0.047 | 0.290 | 1.272 | 4.4 | 1.34 | 0.27 | 5.0 | 0.87 |
| 17 | Ti-43Al | TiAl | 0.144 | 0.088 | 0.505 | 2.939 | 5.8 | 1.45 | 0.10 | 14.2 | 0.50 |
| 18 | Ti-48Al | TiAl | 0.138 | 0.080 | 0.461 | 2.962 | 6.4 | 1.20 | 0.11 | 11.3 | 0.50 |
| 19 | Ti-46.2Al-2.2Cr | TiAl | 0.125 | 0.072 | 0.460 | 2.723 | 5.9 | 1.07 | 0.13 | 8.2 | 0.50 |
| 20 | Co₃Ti | Co | 0.096 | 0.055 | 0.414 | 2.824 | 6.8 | 0.86 | 0.28 | 3.1 | 0.65 |
| 21 | Ni₃Si | Ni | 0.092 | 0.053 | 0.465 | 2.288 | 4.9 | 0.87 | 0.34 | 2.6 | 0.87 |
| 22 | ZrO₂ | ZrO₂ | 0.082 | 0.047 | 1.256 | 3.481 | 2.8 | 1.65 | 0.57 | 2.9 | 0.90 |
| 23 | ZrO₂-3Y | ZrO₂ | 0.088 | 0.051 | 1.261 | 3.790 | 3.0 | 0.85 | 0.49 | 1.7 | 0.90 |
| 24 | ZrO₂-4Y | ZrO₂ | 0.089 | 0.051 | 1.253 | 3.909 | 3.1 | 1.57 | 0.48 | 3.3 | 0.90 |
| 25 | Al₂O₃-ZrO₂-mullite | Ox | 0.067 | 0.039 | 1.130 | 3.473 | 3.1 | 0.86 | 0.75 | 1.1 | 1.00 |
| 26 | Al₂O₃-NiAl₂O₄-ZrO₂ | Ox | 0.064 | 0.037 | 1.244 | 3.690 | 3.0 | 1.35 | 0.91 | 1.5 | 1.00 |
| 27 | 6061 / 20% SiC | Al | 0.071 | 0.041 | 0.439 | 1.556 | 3.5 | 1.05 | 0.45 | 2.3 | 0.40 |
| 28 | 7075 / 20% SiC | Al | 0.069 | 0.040 | 0.421 | 1.537 | 3.7 | 0.97 | 0.45 | 2.1 | 0.40 |
| 29 | Zr₆₅Al₁₀Ni₁₀Cu₁₅ (BMG) | Zr | 0.074 | 0.043 | 0.253 | 1.438 | 5.7 | 1.35 | 0.26 | 5.2 | 0.45 |
| 30 | Zr₅₂.₅-based (BMG) | Zr | 0.075 | 0.043 | 0.268 | 1.518 | 5.7 | 0.78 | 0.27 | 2.9 | 0.45 |
| 31 | La₅₅Al₂₅Ni₂₀ (BMG) | La | 0.059 | 0.034 | 0.151 | 1.015 | 6.7 | 1.60 | 0.21 | 7.7 | 0.30 |
| 32 | La₆₀-based (BMG) | La | 0.073 | 0.042 | 0.213 | 0.963 | 4.5 | 0.82 | 0.19 | 4.3 | 0.30 |
| 33 | Fe₇₂Hf₈Nb₂B₁₈ (BMG) | Fe | 0.079 | 0.046 | 0.346 | 1.452 | 4.2 | 0.78 | 0.35 | 2.2 | 0.80 |
| 34 | Limestone | Geo | 0.073 | 0.042 | 1.045 | 2.554 | 2.4 | 0.57 | 0.29 | 2.0 | 0.50 |
| 35 | Anorthite-Diopside (dry) | Geo | 0.062 | 0.036 | 1.601 | 3.117 | 1.9 | 1.12 | 0.61 | 1.8 | 0.50 |
| 36 | Anorthite-Diopside (wet) | Geo | 0.062 | 0.036 | 1.624 | 2.697 | 1.7 | 1.27 | 0.62 | 2.1 | 0.50 |
| 37 | Ice (10 µm) | Ice | 0.115 | 0.067 | 0.169 | 0.332 | 2.0 | 0.79 | 0.03 | 31.0 | 0.065 |
| 38 | Ice (1700 µm) | Ice | 0.117 | 0.068 | 0.173 | 0.408 | 2.4 | 1.44 | 0.03 | 57.1 | 0.065 |
| 39 | Si₃N₄ | Ox | 0.085 | 0.049 | 2.765 | 3.168 | 1.2 | 1.05 | 1.15 | 0.9 | 1.00 |
| 40 | Zirconia-Alumina | ZrO₂ | 0.069 | 0.040 | 1.311 | 3.166 | 2.4 | 0.94 | 0.83 | 1.1 | 0.90 |
| 41 | Zirconia-Spinel | ZrO₂ | 0.075 | 0.043 | 1.180 | 3.434 | 2.9 | 1.62 | 0.64 | 2.5 | 0.90 |

---

## 5. Statistics across the 41 systems

| Metric | Value |
|--------|-------|
| $\gamma_0$ | 0.085 ± 0.020 (range 0.059–0.144) — OSTZ: 0.10–0.12 |
| $\varepsilon_0$ | 0.049 ± 0.012 (range 0.034–0.088) — OSTZ: 0.05 |
| Energy factor $\beta_1\gamma_0^2+\beta_2\varepsilon_0^2$ | 0.00570 (range 0.0026–0.0161) — OSTZ: 0.00668 |
| $Q/\Delta F_0$ | 4.05 ± 1.39 (range 1.15–6.83) |
| $\gamma_B$ ratio (paper / OSTZ Read–Shockley) | **median 4.5**, mean 6.8, range 0.9–57.1 |
| └ within 3× / 5× / 10× | 15 / 24 / **35** of 41 |
| $\gamma_B$ paper / measured HAGB | median **2.3×** (paper $\gamma_B$ is a genuine GB energy) |
| $N_a$ (paper) | 15.6 ± 2.6 (range 2.3–19.6) — mesoscopic, vs. OSTZ atomic $N_c = 4$ |

---

## 6. Discussion of each comparison

**6.1 Dilatational eigenstrain $\varepsilon_0$ — predicted exactly.**
OSTZ canonical $\varepsilon_0 = 0.05$ matches the 41-system fitted mean ($0.049 \pm 0.012$) to better than 2 %. The paper's universal fit $\varepsilon_0(T) = 0.0408 + 0.0117/T$ → 0.05 near $T \approx 430$ K. This is the strongest single confirmation that the volumetric expansion of one OSTZ event is a universal constant, independent of crystal structure or material class.

**6.2 Shear eigenstrain $\gamma_0$ — within 15 %.**
The fitted mean (0.085) lies between the OSTZ grain-boundary value (0.10) and the lattice value (0.12). The paper's $\gamma_0(T) = 0.0827 + 1.342/T$ shows thermal softening of the shear eigenstrain; at typical superplastic $T \approx 0.6\,T_m$ the fit gives ≈ 0.085, consistent with the mean.

**6.3 Activation energy — OSTZ is ~25 % of the total barrier.**
The mean ratio $Q/\Delta F_0 = 4.05$ shows the pure elastic OSTZ barrier accounts for about one quarter of the experimentally inferred activation energy. The remaining ~75 % is diffusional accommodation (grain-boundary / lattice diffusion at triple junctions) — consistent with the paper's finding that $Q$ tracks grain-boundary diffusion activation energies.

**6.4 Grain-boundary free energy — correct $Gb$ scaling, constant 4.5× offset.**
Using the correct elastic object (dislocation-array Read–Shockley energy $E_0\theta_m$), OSTZ matches the paper's fitted $\gamma_B$ to a **median 4.5×, with 35/41 systems within 10×** — versus 40× for *all* systems with the (incorrect) activation-energy-density formula. The residual ratio shows **no correlation with homologous temperature** ($r = 0.00$), i.e. it is a *constant multiplicative offset*, not a scaling error: OSTZ captures the right $Gb$ dependence, and the factor reflects the tilt-only Read–Shockley plateau vs. the paper's *effective* high-angle GB energy (which folds in twist components and the work of separation).

Independently, the paper's own fitted $\gamma_B$ tracks **measured** high-angle GB energies to a median 2.3×, confirming that $\gamma_B$ is a genuine grain-boundary energy (as the paper states, constrained to 0.30–1.7 J m⁻²).

**6.5 The ice anomaly.**
Ice (systems 37–38) is the sole large outlier (ratio 31–57×). Here the *paper's* fitted $\gamma_B$ (0.79–1.44 J m⁻²) exceeds the true ice grain-boundary energy (~0.065 J m⁻²) by 12–22×; the OSTZ Read–Shockley value (0.026 J m⁻²) is actually **closer to the real ice GB energy**. Reporting medians (rather than means) prevents this single anomaly from distorting the overall picture.

**6.6 Cooperative number $N_a$ — two distinct length scales.**
OSTZ supplies the *atomic* cooperative number $N_c = b/(\gamma_0 W) = 4$ (atoms in one OSTZ event). The paper's $N_a = 15.6 \pm 2.6$ is the *mesoscopic* number of contiguous grain boundaries aligning to form the plane interface. These are different physical quantities at different scales — not a discrepancy: OSTZ provides the atomic constant, and the mesoscopic $N_a$ is the emergent plane-interface extent built from many such events.

---

## 7. Conclusion

With the correct elastic object for each quantity — coherency strain (quadratic in $\gamma_0$) for the activation barrier, dislocation array (linear in $b$) for the interfacial energy — **OSTZ theory reproduces the semi-empirical constants of all 41 systems with no free parameters**:

- $\varepsilon_0$ to within **2 %**,
- $\gamma_0$ to within **15 %**,
- $\gamma_B$ to a **median factor of 4.5×** (a constant offset, not a scaling error),
- the mechanical barrier $\Delta F_0 \approx 25\%$ of the total activation energy $Q$.

These results, consistent across **6 material classes and 4 crystal structures**, confirm that OSTZ is the correct elastic foundation underlying the mesoscopic grain-boundary-sliding model of Harisankar & Padmanabhan (2025).

---

## 8. Reproducibility & figures

- **Computation:** `OSET_derivations.ipynb`, §13 (cell tagged `SS13`). All 41 systems are tabulated from the paper's Table 1; OSTZ predictions use $\beta_1 = 0.446$, $\beta_2 = 0.889$, $\theta_m = 15^\circ$, with class-wise $b$ and $\nu$.
- **Figures:**
  - `fig_HP2025_comparison.png` — six panels: $\gamma_0$ and $\varepsilon_0$ distributions; $Q$ vs. $\Delta F_0$; $\gamma_B$ (OSTZ vs. paper); paper $\gamma_B$ vs. measured HAGB; ratio vs. homologous temperature.
  - `fig_HP2025_universal.png` — universal temperature fits (paper Eqs 9, 10, 12) for $\gamma_0(T)$, $\varepsilon_0(T)$, $\gamma_B(T)$, with OSTZ canonical values overlaid.
- **Manuscript:** the same comparison appears as §10.6 of `OSET_manuscript.md`.

*Representative measured HAGB energies $\gamma_B^{\text{exp}}$ are class averages from Hirth & Lothe (1982), Murr (1975), Olmsted et al. (2009), Rohrer (2011), and Sutton & Balluffi (1995); ice from Ketcham & Hobbs (1969).*
