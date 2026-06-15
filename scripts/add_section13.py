"""Insert §13: Comparison with Harisankar & Padmanabhan (2025) into notebook."""
import json

sec13_md = r"""---
## 13. Comparison with Harisankar & Padmanabhan (2025) — 41-System Validation

**Source:** K.R. Harisankar & K.A. Padmanabhan, *Mater. Sci. Eng. A* **930** (2025) 148175.

The paper fits the OSTZ/GBS model to **41 material systems** (146 data conditions) spanning
metals, intermetallics, ceramics, BMGs, and geological materials.  This section compares
the semi-empirical constants extracted from experiment with OSTZ *ab initio* predictions
(no free parameters beyond $\beta_1$, $\beta_2$ from Eshelby theory).

Key parameters compared:

| Symbol | OSTZ canonical | Paper range |
|--------|---------------|-------------|
| $\gamma_0$ | 0.12 (crystal) / 0.10 (GB) | 0.059–0.144 |
| $\varepsilon_0$ | 0.05 | 0.034–0.088 |
| $\beta_1\gamma_0^2+\beta_2\varepsilon_0^2$ | 0.00864 | 0.0022–0.0149 |
| $\gamma_B$ (J/m²) | §12 OSET formula | 0.38–1.68 |
| $\Delta F_0$ (eV) | $(1/2)(\beta_1\gamma_0^2+\beta_2\varepsilon_0^2)GV_0$ | from $Q/T_m$ |"""

sec13_code = r"""# ═══════════════════════════════════════════════════════════════════════════════
# SS13  Comparison with Harisankar & Padmanabhan (2025)  MSEA 930:148175
# ═══════════════════════════════════════════════════════════════════════════════
# 41 material systems from Table 1 of the paper.
# Each tuple: (label, class, d_um, T_K, Thom, g0, e0, QoTm, GoTm, t0oTm, gB, a, Na)
# Units: d [µm], T [K], Q/Tm [J/(mol K)], G/Tm [MPa/K], tau0/Tm [Pa/K], gB [J/m²]

HP25 = [
# ── Metals & Alloys ──────────────────────────────────────────────────────────
("Zn22Al-2.5",    "Zn",  2.5,   423, 0.558, 0.059, 0.034,  77.49, 51.38,  2902.33, 1.29, 0.262, 14.21),
("Zn22Al-0.9",    "Zn",  0.9,   453, 0.597, 0.061, 0.035, 103.18, 49.39,  3614.76, 1.31, 0.248, 14.05),
("Al33Cu0.4Zr",   "Al",  7.6,   713, 0.869, 0.066, 0.038, 170.35, 41.65,   146.16, 1.61, 0.332, 15.76),
("Al-MgScMn-3",   "Al",  3.0,   723, 0.775, 0.089, 0.051, 135.44, 20.29,  2154.32, 1.57, 0.293, 15.57),
("Al-MgScMn-1",   "Al",  1.0,   523, 0.561, 0.083, 0.048, 131.86, 22.78, 13365.50, 0.92, 0.255, 14.53),
("Al3Mg0.2Sc",    "Al",  0.2,   573, 0.614, 0.086, 0.050, 118.56, 21.87,  5798.52, 1.46, 0.250, 14.77),
("Al-ZnMgSc",     "Al",  0.7,   493, 0.528, 0.080, 0.046, 134.19, 25.27, 26323.70, 1.47, 0.249, 14.29),
("Al5Mg-24",      "Al",  24.0,  748, 0.802, 0.090, 0.052, 141.75, 19.87,  1350.44, 0.87, 0.347, 19.57),
("Al17Si",        "Al",  1.4,   763, 0.818, 0.089, 0.051, 116.48, 23.53,  9024.68, 0.88, 0.278, 15.78),
("Mg6Zn0.8Zr",    "Mg",  0.7,   448, 0.485, 0.080, 0.046, 128.70, 18.04,  9296.55, 1.54, 0.245, 13.98),
("Mg4Y0.7Zr",     "Mg",  2.0,   597, 0.647, 0.087, 0.050, 132.28, 15.21,  5021.66, 1.35, 0.273, 15.04),
("Mg5.8Zn1Y",     "Mg",  17.5,  673, 0.728, 0.088, 0.051, 113.57, 15.38,  3528.13, 1.37, 0.358, 18.11),  # ε₀=0.051 (PDF 0.510 is a typo)
("Ti6Al4V",       "Ti",  0.9,  1023, 0.529, 0.116, 0.067, 124.29, 11.56,  4066.21, 1.14, 0.286, 16.73),
("Cu2.8Al-7",     "Cu",  7.0,   723, 0.533, 0.096, 0.056, 106.59, 24.61,  8134.25, 1.48, 0.328, 15.70),
("Cu2.8Al-3",     "Cu",  3.0,   673, 0.497, 0.095, 0.055, 107.17, 25.19,  6939.53, 0.99, 0.289, 15.41),
("IN836",         "Ni",  2.5,   735, 0.652, 0.081, 0.047, 108.91, 33.01, 19591.86, 1.34, 0.289, 15.61),
# ── Intermetallic Compounds ──────────────────────────────────────────────────
("Ti43Al",        "TiAl",5.0,  1273, 0.695, 0.144, 0.088, 154.81,  7.23,  7394.89, 1.45, 0.341, 15.07),
("Ti48Al",        "TiAl",0.9,  1163, 0.654, 0.138, 0.080, 160.71,  7.73, 23813.29, 1.20, 0.291, 17.00),
("Ti46Al2Cr",     "TiAl",0.8,  1073, 0.599, 0.125, 0.072, 146.66,  9.39, 22138.44, 1.07, 0.287, 16.88),
("Co3Ti",         "Co",  24.0, 1173, 0.776, 0.096, 0.055, 180.25, 23.78, 26514.51, 0.86, 0.389, 12.15),
("Ni3Si",         "Ni",  15.0, 1323, 0.882, 0.092, 0.053, 147.16, 31.01,  6179.96, 0.87, 0.406, 10.97),
# ── Ceramics / Composites ────────────────────────────────────────────────────
("ZrO2",          "ZrO2",0.07, 1273, 0.426, 0.082, 0.047, 112.41, 17.79, 13330.00, 1.65, 0.282, 17.56),
("ZrO2-3Y",       "ZrO2",0.51, 1523, 0.510, 0.088, 0.051, 122.47, 15.38,  1572.93, 0.85, 0.288, 17.30),
("ZrO2-4Y",       "ZrO2",0.75, 1573, 0.526, 0.089, 0.051, 126.12, 15.05,   592.37, 1.57, 0.291, 17.05),
("Al2O3-ZrO2-Si", "Ox",  0.4,  1673, 0.723, 0.067, 0.039, 144.83, 33.26,  1919.62, 0.86, 0.284, 17.16),
("Al2O3-NiAl-Zr", "Ox",  1.3,  1623, 0.702, 0.064, 0.037, 153.98, 40.41,  4690.84, 1.35, 0.298, 16.49),
("6061-20SiC",    "Al",  0.8,   773, 0.590, 0.071, 0.041, 114.57, 37.50,  3954.22, 1.05, 0.272, 15.87),
("7075-20SiC",    "Al",  5.0,   753, 0.596, 0.069, 0.040, 117.39, 39.33,  7220.88, 0.97, 0.314, 15.62),
# ── Bulk Metallic Glasses ────────────────────────────────────────────────────
("Zr65BMG",       "Zr",  10.,   653, 0.587, 0.074, 0.043, 124.71, 26.27, 13821.94, 1.35, 0.340, 16.22),
("Zr52BMG",       "Zr",  10.,   683, 0.621, 0.075, 0.043, 133.19, 27.60, 19572.74, 0.78, 0.343, 16.14),
("La55BMG",       "La",  10.,   483, 0.687, 0.059, 0.034, 139.26, 30.18, 14096.72, 1.60, 0.322, 16.57),
("La60BMG",       "La",  10.,   460, 0.713, 0.073, 0.042, 144.00, 30.35,  2170.54, 0.82, 0.320, 16.60),
("Fe72BMG",       "Fe",  10.,   863, 0.694, 0.079, 0.046, 112.65, 38.33, 29436.88, 0.78, 0.358, 15.48),
# ── Geological / Ice ─────────────────────────────────────────────────────────
("Limestone",     "Geo", 4.2,   973, 0.604, 0.073, 0.042, 152.98, 12.89, 13761.67, 0.57, 0.321, 15.80),
("AnDi-dry",      "Geo", 3.1,  1323, 0.726, 0.062, 0.036, 165.03, 24.03, 10762.47, 1.12, 0.321, 15.89),
("AnDi-wet",      "Geo", 3.1,  1223, 0.671, 0.062, 0.036, 142.75, 24.36, 22018.63, 1.27, 0.319, 16.03),
("Ice-10",        "Ice", 10.,   220, 0.806, 0.115, 0.067, 117.39,  6.65,  6227.10, 0.79, 0.289, 16.67),
("Ice-1700",      "Ice",1700.,  241, 0.883, 0.117, 0.068, 144.41,  6.57,  1391.93, 1.44, 0.414,  2.30),
# ── Nanostructured ───────────────────────────────────────────────────────────
("Si3N4",         "Ox",  0.07, 1723, 0.793, 0.085, 0.049, 140.67, 54.29,  5945.67, 1.05, 0.278, 17.36),
("ZrO2-Al2O3",    "ZrO2",0.06, 1623, 0.584, 0.069, 0.040, 109.91, 27.94,  2321.02, 0.94, 0.280, 17.55),
("ZrO2-Spinel",   "ZrO2",0.05, 1573, 0.574, 0.075, 0.043, 120.89, 21.78,  5912.42, 1.62, 0.281, 17.62),
]

# Burgers vector (m) by material class
b_nm = {"Zn":0.267,"Al":0.286,"Mg":0.320,"Ti":0.295,"TiAl":0.285,
        "Cu":0.256,"Ni":0.250,"Co":0.255,"Fe":0.248,"Zr":0.275,
        "La":0.300,"ZrO2":0.360,"Ox":0.350,"Geo":0.500,"Ice":0.452}

NAv = 6.0221e23;  eV_J = 1.6022e-19
beta1 = 1 - 2*S1313v(1/3)     # 0.4456  (α=0.5, ν=1/3 Eshelby)
beta2 = 4*(1+1/3)/(9*(1-1/3)) # 0.8889

# --- Compute per system ---
res13 = []
for lbl, cls, d_um, T, Thom, g0, e0, QoTm, GoTm, t0oTm, gB_lit, a_fit, Na_lit in HP25:
    b  = b_nm.get(cls, 0.28) * 1e-9          # m
    Tm = T / Thom                              # K
    G  = GoTm * Tm * 1e6                      # Pa
    V0 = (2/3)*np.pi * (2.5*b)**3             # m³  (W=2.5b, GB context)
    A0 = np.pi * (2.5*b)**2                   # m²

    EF = beta1*g0**2 + beta2*e0**2            # dimensionless energy factor
    dF0_J   = 0.5 * EF * G * V0               # J per OSTZ
    dF0_eV  = dF0_J / eV_J                    # eV per OSTZ
    Q_eV    = QoTm * Tm / NAv / eV_J          # eV per event (Q/Tm × Tm / NA)
    gB_oset = (1/3) * EF * G * 2.5*b         # J/m²  (§12 formula)
    tau0_oset_MPa = np.sqrt(2*G*gB_oset)*1e-6 # MPa   (order-of-magnitude)
    tau0_lit_MPa  = t0oTm * Tm * 1e-6         # Pa/K × K = Pa → MPa
    Na_oset = 1 / (2.5 * g0)                   # theoretical Na from OSTZ

    res13.append(dict(
        lbl=lbl, cls=cls, d=d_um, T=T, Thom=Thom,
        g0=g0, e0=e0, EF=EF, G_GPa=G/1e9,
        dF0=dF0_eV, Q_eV=Q_eV,
        gB_lit=gB_lit, gB_oset=gB_oset,
        tau0_oset=tau0_oset_MPa, tau0_lit=tau0_lit_MPa,
        Na_lit=Na_lit, Na_oset=Na_oset, a=a_fit,
    ))

# --- Print comparison table ---
print("SS13  OSTZ Theory vs. Harisankar & Padmanabhan (2025)  — 41 Systems")
print("="*120)
hdr = (f"{'System':<18s} {'Cl':3s} {'g0':6s} {'e0':6s} {'EF_oset':8s} "
       f"{'ΔF0(eV)':8s} {'Q(eV)':8s} {'Q/ΔF0':7s} "
       f"{'gB_lit':8s} {'gB_oset':8s} {'gB_rat':7s} "
       f"{'Na_lit':7s} {'Na_th':6s}")
print(hdr)
print('-'*120)
for r in res13:
    ratio_F = r['Q_eV']/r['dF0'] if r['dF0']>0 else float('nan')
    ratio_g = r['gB_lit']/r['gB_oset'] if r['gB_oset']>0 else float('nan')
    print(f"{r['lbl']:<18s} {r['cls']:3s} "
          f"{r['g0']:6.3f} {r['e0']:6.3f} {r['EF']:8.5f} "
          f"{r['dF0']:8.3f} {r['Q_eV']:8.3f} {ratio_F:7.2f} "
          f"{r['gB_lit']:8.3f} {r['gB_oset']:8.3f} {ratio_g:7.2f} "
          f"{r['Na_lit']:7.2f} {r['Na_oset']:6.2f}")
print()

# --- Statistics ---
EFs   = np.array([r['EF']    for r in res13])
g0s   = np.array([r['g0']    for r in res13])
e0s   = np.array([r['e0']    for r in res13])
ratF  = np.array([r['Q_eV']/r['dF0'] for r in res13])
ratG  = np.array([r['gB_lit']/r['gB_oset'] for r in res13])
NaL   = np.array([r['Na_lit'] for r in res13])
NaT   = np.array([r['Na_oset']for r in res13])

print(f"γ₀:           range {g0s.min():.3f}–{g0s.max():.3f},  mean {g0s.mean():.3f}±{g0s.std():.3f}  |  OSTZ canonical: 0.12")
print(f"ε₀:           range {e0s.min():.3f}–{e0s.max():.3f},  mean {e0s.mean():.3f}±{e0s.std():.3f}  |  OSTZ canonical: 0.05")
print(f"Energy factor: range {EFs.min():.5f}–{EFs.max():.5f}, mean {EFs.mean():.5f}  |  OSTZ canon: {beta1*0.10**2+beta2*0.05**2:.5f}")
print(f"Q/ΔF₀ ratio:  range {ratF.min():.2f}–{ratF.max():.2f},  mean {ratF.mean():.2f}±{ratF.std():.2f}")
print(f"γ_B ratio:    range {ratG.min():.2f}–{ratG.max():.2f},  mean {ratG.mean():.2f}±{ratG.std():.2f}")
print(f"Na (paper):   range {NaL.min():.1f}–{NaL.max():.1f},    mean {NaL.mean():.1f}±{NaL.std():.1f}")
print(f"Na (theory):  range {NaT.min():.1f}–{NaT.max():.1f},    mean {NaT.mean():.1f}±{NaT.std():.1f}  (1/(2.5γ₀))")

# --- Spot checks ---
# 1. Mean γ₀ should be close to 0.10 (GB canonical, within 30%)
check('mean γ₀ (41 systems) within 30% of 0.10', g0s.mean(), 0.10, tol=0.30)
# 2. Mean ε₀ within 5% of 0.05  (after PDF typo fix: 0.510→0.051 for Mg5.8Zn1Y)
check('mean ε₀ (41 systems) within 5% of 0.05', e0s.mean(), 0.05, tol=0.05)
# 3. Q/ΔF₀ ratio is positive and finite for all
assert all(ratF > 0) and all(np.isfinite(ratF)), "Q/ΔF₀ has non-positive or infinite values"
print("✓ PASS  Q/ΔF₀ > 0 for all 41 systems")
# 4. OSTZ systematically underestimates γ_B  (mean ratio > 1)
assert ratG.mean() > 1, "Expected OSTZ to underestimate γ_B"
print(f"✓ PASS  OSTZ underestimates γ_B: mean ratio = {ratG.mean():.1f}×")
# 5. Na_theory < Na_paper  (OSTZ gives lower bound)
assert NaT.mean() < NaL.mean()
print(f"✓ PASS  Na_theory ({NaT.mean():.1f}) < Na_paper ({NaL.mean():.1f}): OSTZ gives lower bound")

# --- Figures ---
fig13, axes13 = plt.subplots(2, 3, figsize=(16, 10))
fig13.suptitle('OSTZ Theory vs. Harisankar & Padmanabhan (2025) — 41 Systems', fontsize=13, fontweight='bold')

# Panel A: γ₀ distribution
ax = axes13[0, 0]
ax.hist(g0s, bins=14, color='steelblue', alpha=0.75, edgecolor='k')
ax.axvline(0.12, color='red', lw=2, ls='--', label='OSTZ crystal (0.12)')
ax.axvline(0.10, color='darkorange', lw=2, ls='--', label='OSTZ GB (0.10)')
ax.axvline(g0s.mean(), color='navy', lw=2, label=f'paper mean ({g0s.mean():.3f})')
ax.set_xlabel(r'$\gamma_0$ (shear eigenstrain)');  ax.set_ylabel('Count')
ax.set_title(r'$\gamma_0$ distribution');  ax.legend(fontsize=8)

# Panel B: ε₀ distribution
ax = axes13[0, 1]
ax.hist(e0s, bins=12, color='teal', alpha=0.75, edgecolor='k')
ax.axvline(0.05, color='red', lw=2, ls='--', label='OSTZ canonical (0.05)')
ax.axvline(e0s.mean(), color='navy', lw=2, label=f'paper mean ({e0s.mean():.3f})')
ax.set_xlabel(r'$\varepsilon_0$ (dilatational eigenstrain)');  ax.set_ylabel('Count')
ax.set_title(r'$\varepsilon_0$ distribution');  ax.legend(fontsize=8)

# Panel C: Energy factor
ax = axes13[0, 2]
EF_canon_gb = beta1*0.10**2 + beta2*0.05**2
EF_canon_cr = beta1*0.12**2 + beta2*0.05**2
ax.hist(EFs, bins=14, color='purple', alpha=0.75, edgecolor='k')
ax.axvline(EF_canon_cr, color='red', lw=2, ls='--', label=f'OSTZ crystal ({EF_canon_cr:.5f})')
ax.axvline(EF_canon_gb, color='darkorange', lw=2, ls='--', label=f'OSTZ GB ({EF_canon_gb:.5f})')
ax.set_xlabel(r'$\beta_1\gamma_0^2+\beta_2\varepsilon_0^2$');  ax.set_ylabel('Count')
ax.set_title('Energy factor'); ax.legend(fontsize=8)

# Panel D: Q_paper vs ΔF₀_OSTZ (scatter)
ax = axes13[1, 0]
dF0s = np.array([r['dF0'] for r in res13])
Q_evs= np.array([r['Q_eV'] for r in res13])
for cls_v, mk, co in [("Al",'o','royalblue'),("Mg",'s','teal'),("Ti",'^','orange'),
                       ("Cu",'D','brown'),("Ni",'v','olive'),("TiAl",'p','magenta'),
                       ("ZrO2",'*','green'),("Ox",'h','darkgreen'),
                       ("Zr",'<','gray'),("La",'>','pink'),("Fe",'8','red'),
                       ("Geo",'P','saddlebrown'),("Ice",'X','cornflowerblue'),("Co",'x','navy')]:
    idx = [i for i,r in enumerate(res13) if r['cls']==cls_v]
    if idx:
        ax.scatter(dF0s[idx], Q_evs[idx], marker=mk, color=co, label=cls_v, s=50)
lim = max(dF0s.max(), Q_evs.max())*1.1
ax.plot([0,lim],[0,lim],'k--',lw=1,label='1:1')
ax.set_xlabel(r'$\Delta F_0$ OSTZ (eV)'); ax.set_ylabel(r'$Q$ paper (eV/event)')
ax.set_title(r'Activation energy: OSTZ $\Delta F_0$ vs. paper $Q$')
ax.legend(fontsize=7, ncol=2)

# Panel E: γ_B comparison
ax = axes13[1, 1]
gB_lit_arr  = np.array([r['gB_lit']  for r in res13])
gB_oset_arr = np.array([r['gB_oset'] for r in res13])
ax.scatter(gB_oset_arr, gB_lit_arr, c='purple', alpha=0.7, s=50, edgecolors='k', lw=0.5)
lim = max(gB_lit_arr.max(), gB_oset_arr.max())*1.1
ax.plot([0,lim],[0,lim],'k--',lw=1,label='1:1')
for fac, ls_v in [(5,'--'),(10,':')]:
    ax.plot([0,lim],[0,fac*lim],'gray',lw=0.8,ls=ls_v,label=f'{fac}× line')
ax.set_xlabel(r'$\gamma_B$ OSTZ (J/m²)');  ax.set_ylabel(r'$\gamma_B$ paper (J/m²)')
ax.set_title(r'Grain boundary free energy $\gamma_B$');  ax.legend(fontsize=8)

# Panel F: Na comparison
ax = axes13[1, 2]
ax.scatter(NaT, NaL, c='darkorange', alpha=0.7, s=50, edgecolors='k', lw=0.5)
lim = max(NaT.max(), NaL.max())*1.1
ax.plot([0,lim],[0,lim],'k--',lw=1,label='1:1')
ax.set_xlabel(r'$N_a$ theory: $1/(2.5\gamma_0)$');  ax.set_ylabel(r'$N_a$ paper (fitted)')
ax.set_title(r'Cooperative grain count $N_a$');  ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('fig_HP2025_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: fig_HP2025_comparison.png")

# --- Universal temperature fits ---
T_range = np.linspace(300, 2000, 300)
g0_fit   = 0.0827 + 1.3422 / T_range    # Eq.(9)
e0_fit   = 0.0408 + 0.0117 / T_range    # Eq.(10)
GoTm_fit = 26.45  - 928.98  / T_range   # Eq.(8)
t0oTm_fit= 4.593e-3 + 0.4511 / T_range  # Eq.(11)

T_data = np.array([r['T'] for r in res13])
g0_data = np.array([r['g0'] for r in res13])
e0_data = np.array([r['e0'] for r in res13])

fig13b, axes13b = plt.subplots(1, 2, figsize=(12, 5))
fig13b.suptitle('Universal Temperature Dependences  (Eqs 8–11)', fontsize=12, fontweight='bold')

ax = axes13b[0]
ax.plot(T_range, g0_fit, 'r-', lw=2, label=r'$\gamma_0 = 0.0827+1.3422/T$ (Eq. 9)')
ax.axhline(0.12, color='blue', ls='--', lw=1.5, label='OSTZ canonical 0.12')
ax.axhline(0.10, color='navy', ls=':', lw=1.5, label='OSTZ GB 0.10')
ax.scatter(T_data, g0_data, c='gray', s=30, zorder=5, label='Paper Table 1')
ax.set_xlabel('T (K)');  ax.set_ylabel(r'$\gamma_0$')
ax.set_title(r'Shear eigenstrain $\gamma_0$');  ax.legend(fontsize=9)

ax = axes13b[1]
ax.plot(T_range, e0_fit, 'g-', lw=2, label=r'$\varepsilon_0 = 0.0408+0.0117/T$ (Eq. 10)')
ax.axhline(0.05, color='blue', ls='--', lw=1.5, label='OSTZ canonical 0.05')
ax.scatter(T_data, e0_data, c='gray', s=30, zorder=5, label='Paper Table 1')
ax.set_xlabel('T (K)');  ax.set_ylabel(r'$\varepsilon_0$')
ax.set_title(r'Dilatational eigenstrain $\varepsilon_0$');  ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig('fig_HP2025_universal.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: fig_HP2025_universal.png")
print()
print("SS13 Summary:")
print(f"  OSTZ predicts ΔF₀ = {dF0s.mean():.3f} ± {dF0s.std():.3f} eV; paper Q = {Q_evs.mean():.3f} ± {Q_evs.std():.3f} eV")
print(f"  Mean Q/ΔF₀ = {ratF.mean():.2f} → OSTZ elastic energy is ~{1/ratF.mean()*100:.0f}% of total activation barrier")
print(f"  Mean γ_B ratio = {ratG.mean():.1f}× → OSTZ gives elastic lower bound to grain boundary free energy")
print(f"  γ₀: paper mean {g0s.mean():.3f} vs OSTZ 0.10–0.12; within {abs(g0s.mean()-0.10)/0.10*100:.0f}% of GB value")"""

# ── Insert cells ──────────────────────────────────────────────────────────────
with open('OSET_derivations.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

def make_md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src}

def make_code(src):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src}

# Insert after cell 35 (§12 code) → at index 36
nb['cells'].insert(36, make_md(sec13_md))
nb['cells'].insert(37, make_code(sec13_code))

print(f"Total cells: {len(nb['cells'])}")
for i in [34, 35, 36, 37, 38]:
    src_p = ''.join(nb['cells'][i]['source'])[:60].encode('ascii', 'replace').decode('ascii')
    print(f"  Cell {i} [{nb['cells'][i]['cell_type'][:4]}]: {src_p}")

with open('OSET_derivations.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
print("Notebook saved.")
