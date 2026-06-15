import json

with open('OSET_derivations.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

def code_cell(src):
    return {"cell_type":"code","metadata":{},"source":src,"outputs":[],"execution_count":None}

def md_cell(src):
    return {"cell_type":"markdown","metadata":{},"source":src}

# ── Section 11 markdown ────────────────────────────────────────────────────────
cell_md = r"""---
## 11. Multi-Material Dislocation Property Comparison

OSET predictions ($W = b$, $\gamma_0 = 0.12$, $\varepsilon_0 = 0.05$, $\alpha = 0.5$)
for crystal-interior dislocations vs. literature:

| Quantity | OSET formula (crystal interior, $W=b$) |
|---|---|
| $\beta_1$ | $1 - 2S_{1313}(\nu, \alpha{=}0.5)$ |
| $\tau_P/G$ | $\tfrac{1}{2}\beta_1\gamma_0\exp(-2\pi) \times 2/(1-\nu)$ |
| $\gamma_{SF}$ | $\gamma_0^2 G b/4$ &nbsp; (with $\beta_1^{\rm eff}{=}1$, reproduces §10.2 table) |
| $E_{\rm core}/(2b)$ | $\pi\beta_1 b^2 \gamma_0 G/6$ per unit length (eV/Å) |
| $\Delta F_0$ | $\tfrac{1}{2}(\beta_1\gamma_0^2 + \beta_2\varepsilon_0^2)G V_0$, &nbsp; $V_0 = \tfrac{2\pi}{3}b^3$ |
| $N_c$ | $b/(\gamma_0 b) = 1/\gamma_0 = 8.33$ (FCC, $W{=}b$) |
"""

# ── Section 11 code cell ───────────────────────────────────────────────────────
cell_code = r"""# ── §10.2 Multi-Metal OSET Table with Literature Comparison ─────────────────
# Crystal interior: W = b, gamma0 = 0.12, eps0 = 0.05, alpha = 0.5
# beta1^eff = 1 for SFE and Peierls (per manuscript §10.2 header)

# ── Material parameters (room-temperature values) ────────────────────────────
# [struct, b(m), G(Pa), nu]
metals_db = {
    'Cu': ('FCC', 0.2556e-9, 48.3e9, 0.343),
    'Al': ('FCC', 0.2863e-9, 26.2e9, 0.347),
    'Ni': ('FCC', 0.2492e-9, 76.0e9, 0.276),
    'Ag': ('FCC', 0.2889e-9, 30.3e9, 0.367),
    'Au': ('FCC', 0.2884e-9, 27.0e9, 0.440),
    'Fe': ('BCC', 0.2482e-9, 82.0e9, 0.291),
    'W':  ('BCC', 0.2741e-9,161.0e9, 0.280),
}

# ── Literature reference values ───────────────────────────────────────────────
# SFE (mJ/m²): experimental values (Hirth & Lothe 1982; Murr 1975)
sfe_exp = {'Cu':45, 'Al':166, 'Ni':125, 'Ag':16, 'Au':32, 'Fe':None, 'W':None}

# Peierls stress τ_P/G: experimental estimates (Caillard & Martin 2003; Suzuki 1968)
# FCC values from low-T CRSS; BCC from kink-pair analysis
tau_P_exp_G = {
    'Cu': (5e-5, 2e-4),    # FCC, small but measurable
    'Al': (1e-6, 5e-5),    # FCC, anomalously low
    'Ni': (1e-4, 5e-4),    # FCC, intermediate
    'Ag': (2e-5, 1e-4),    # FCC, low (wide stacking fault)
    'Au': (2e-5, 1e-4),    # FCC, low
    'Fe': (4e-3, 2e-2),    # BCC, strong Peierls (kink-pair)
    'W':  (2e-2, 5e-2),    # BCC, very strong
}

# Core energy E_core/(2b) (eV/Å): DFT literature
# (Clouet 2019 for Al; Woodward 2002 for Cu/Ni; Frederiksen 2003 for Fe/W)
Ecore_dft = {
    'Cu': (0.05, 0.15), 'Al': (0.03, 0.08), 'Ni': (0.10, 0.25),
    'Ag': (0.03, 0.10), 'Au': (0.04, 0.12),
    'Fe': (0.20, 0.50), 'W':  (0.50, 1.20),
}

# ── OSET computation ──────────────────────────────────────────────────────────
g0 = 0.12;   e0 = 0.05

def S1313v(nu_v, alpha=0.5):
    e = np.sqrt(1-alpha**2)
    I1 = (2*np.pi*alpha/e**3)*(np.arccos(alpha)-alpha*e)
    I3 = 4*np.pi-2*I1
    I13 = (I3-I1)/(1-alpha**2)
    return ((1+alpha**2)*I13+(1-2*nu_v)*(I1+I3))/(16*np.pi*(1-nu_v))

res = {}
for name, (struct, b, G, nu) in metals_db.items():
    b1_esh = 1 - 2*S1313v(nu)           # Eshelby β₁
    b2     = 4*(1+nu)/(9*(1-nu))        # β₂ sphere approx
    W = b
    V0 = (2/3)*np.pi*b**3

    # Peierls stress (with β₁^Esh and HS factor 2/(1-ν))
    tau_th   = b1_esh*g0*G/2            # threshold stress
    tau_P    = tau_th*np.exp(-2*np.pi)*(2/(1-nu))

    # SFE: gamma0^2 * G * b / 4  (β₁^eff = 1, per manuscript §10.2)
    sfe_oset = g0**2*G*b/4 * 1e3        # mJ/m²

    # Core energy per 2b
    Ecore = np.pi*b1_esh*b*g0*G*W**2/3  # J
    Ecore_per2b = (Ecore/eV) / (2*b*1e10)  # eV/Å

    # Activation energy (crystal interior, W=b)
    dF0 = 0.5*(b1_esh*g0**2 + b2*e0**2)*G*V0/eV

    res[name] = dict(struct=struct, b_nm=b*1e9, G=G/1e9, nu=nu,
                     b1_esh=b1_esh, Nc=1/g0,
                     tau_P_G=tau_P/G, sfe=sfe_oset,
                     Ecore_per2b=Ecore_per2b, dF0=dF0)

# ── Print table 1: reproduce manuscript §10.2 ─────────────────────────────────
print("Table 1: OSET predictions reproducing manuscript §10.2")
print("="*100)
hdr = f"{'Metal':5s} {'struct':5s} {'G(GPa)':7s} {'ν':6s} {'b(nm)':7s} " \
      f"{'β₁(Esh)':8s} {'Nc':5s} {'τ_P/G':10s} {'γ_SF(OSET)':12s} {'γ_SF(exp)':10s} {'Match':6s}"
print(hdr)
print('-'*100)
for n, r in res.items():
    sfe_e = sfe_exp.get(n)
    match = '?' if sfe_e is None else ('✓' if abs(r['sfe']-sfe_e)/sfe_e < 0.20 else '~' if abs(r['sfe']-sfe_e)/sfe_e < 1.0 else '✗')
    sfe_str = f"{sfe_e:5d}" if sfe_e else "  N/A"
    print(f"{n:5s} {r['struct']:5s} {r['G']:7.1f} {r['nu']:6.3f} {r['b_nm']:7.4f} "
          f"{r['b1_esh']:8.4f} {r['Nc']:5.1f} {r['tau_P_G']:10.2e} "
          f"{r['sfe']:12.1f} {sfe_str:>10s} {match:6s}")

print()
print("  Notes: τ_P/G uses β₁^Esh×γ₀/2×exp(−2π)×2/(1−ν); γ_SF uses β₁^eff=1 (manuscript convention).")
print()

# ── Print table 2: additional dislocation quantities with literature ───────────
print("Table 2: Extended comparison — core energy, ΔF₀, Peierls stress")
print("="*110)
hdr2 = f"{'Metal':5s} {'E_c/(2b)(eV/Å)':16s} {'DFT range':14s} {'ΔF₀(eV)':10s} " \
       f"{'τ_P/G(OSET)':14s} {'τ_P/G lit.range':20s} {'HS factor':10s}"
print(hdr2)
print('-'*110)
for n, r in res.items():
    dft_lo, dft_hi = Ecore_dft[n]
    dft_str = f"{dft_lo:.2f}–{dft_hi:.2f}"
    ec_match = '✓' if dft_lo<=r['Ecore_per2b']<=dft_hi else '~' if dft_lo/3<=r['Ecore_per2b']<=dft_hi*3 else '✗'
    tp_lo, tp_hi = tau_P_exp_G[n]
    _, _, nu = metals_db[n][1], metals_db[n][2], metals_db[n][3]
    hs = 2/(1-nu)
    tp_str = f"{tp_lo:.0e}–{tp_hi:.0e}"
    tp_match = '✓' if tp_lo<=r['tau_P_G']<=tp_hi else '~' if tp_lo/5<=r['tau_P_G']<=tp_hi*5 else '✗'
    print(f"{n:5s} {r['Ecore_per2b']:10.4f} {ec_match:2s} {dft_str:14s} {r['dF0']:10.4f} "
          f"{r['tau_P_G']:14.2e} {tp_str:20s} {hs:.3f} {tp_match:2s}")

print()
print("  E_c/(2b) formula: π β₁^Esh γ₀ G b² / 6  (per unit OSTZ diameter 2W=2b)")
print("  ΔF₀ formula:      ½(β₁γ₀² + β₂ε₀²) G V₀  with V₀=(2/3)πb³")
print("  DFT references:   Clouet(2019)=Al; Woodward(2002)=Cu,Ni; Frederiksen(2003)=Fe,W")
print()

# ── Spot-checks against manuscript §10.2 ─────────────────────────────────────
check('γ_SF^OSET Cu  (mJ/m²)',  res['Cu']['sfe'],  43,  tol=0.03)
check('γ_SF^OSET Al  (mJ/m²)',  res['Al']['sfe'],  26,  tol=0.03)
check('γ_SF^OSET Ni  (mJ/m²)',  res['Ni']['sfe'],  67,  tol=0.03)
check('γ_SF^OSET Ag  (mJ/m²)',  res['Ag']['sfe'],  31,  tol=0.03)
check('γ_SF^OSET Au  (mJ/m²)',  res['Au']['sfe'],  27,  tol=0.03)
check('γ_SF^OSET Fe  (mJ/m²)',  res['Fe']['sfe'],  72,  tol=0.03)
check('γ_SF^OSET W   (mJ/m²)',  res['W']['sfe'],  155,  tol=0.03)
check('β₁(Cu, ν=0.343)',   res['Cu']['b1_esh'], 0.4487, tol=2e-3)
check('β₁(Ni, ν=0.276)',   res['Ni']['b1_esh'], 0.4290, tol=2e-3)
check('E_core/(2b) Cu (eV/Å)', res['Cu']['Ecore_per2b'], 0.124, tol=0.02)

# ── Comparison plot ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('OSET vs. Literature: Crystal-Interior Dislocation Properties ($W=b$, $\\gamma_0=0.12$)',
             fontsize=12, fontweight='bold')

names_list = list(metals_db.keys())
structs    = [res[n]['struct'] for n in names_list]
colors     = ['royalblue' if s=='FCC' else 'tomato' if s=='BCC' else 'seagreen' for s in structs]

# ── Plot 1: SFE ──────────────────────────────────────────────────────────────
ax = axes[0]
x   = np.arange(len(names_list))
oset_sfe = [res[n]['sfe'] for n in names_list]
exp_sfe  = [sfe_exp[n] if sfe_exp[n] else 0 for n in names_list]
bars = ax.bar(x-0.2, oset_sfe, 0.35, label='OSET', color=colors, alpha=0.85)
ax.bar(x+0.2, exp_sfe, 0.35, label='Experiment', color='silver', edgecolor='black', alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(names_list)
ax.set_ylabel('SFE (mJ/m²)'); ax.set_title('Stacking Fault Energy')
ax.legend(); ax.set_ylim(0, 200)
for xi, (o, e) in zip(x, zip(oset_sfe, exp_sfe)):
    if e > 0:
        ratio = o/e
        ax.text(xi, max(o,e)+5, f'{ratio:.1f}×', ha='center', fontsize=7.5, color='darkblue')

# ── Plot 2: Peierls stress τ_P/G ─────────────────────────────────────────────
ax = axes[1]
oset_tP = [res[n]['tau_P_G'] for n in names_list]
exp_tP_mid = [(tau_P_exp_G[n][0]*tau_P_exp_G[n][1])**0.5 for n in names_list]
exp_tP_lo  = [tau_P_exp_G[n][0] for n in names_list]
exp_tP_hi  = [tau_P_exp_G[n][1] for n in names_list]
ax.bar(x-0.2, oset_tP, 0.35, label='OSET', color=colors, alpha=0.85)
ax.bar(x+0.2, exp_tP_mid, 0.35, label='Lit. geometric mean', color='silver', edgecolor='black', alpha=0.85)
err_lo = [m-lo for m,lo in zip(exp_tP_mid,exp_tP_lo)]
err_hi = [hi-m for m,hi in zip(exp_tP_mid,exp_tP_hi)]
ax.errorbar(x+0.2, exp_tP_mid, yerr=[err_lo,err_hi], fmt='none', color='black', capsize=4)
ax.set_xticks(x); ax.set_xticklabels(names_list)
ax.set_yscale('log'); ax.set_ylabel('$\\tau_P / G$'); ax.set_title('Peierls Stress / Shear Modulus')
ax.legend()

# ── Plot 3: Core energy E_core/(2b) ─────────────────────────────────────────
ax = axes[2]
oset_Ec  = [res[n]['Ecore_per2b'] for n in names_list]
dft_mid  = [(Ecore_dft[n][0]*Ecore_dft[n][1])**0.5 for n in names_list]
dft_lo   = [Ecore_dft[n][0] for n in names_list]
dft_hi   = [Ecore_dft[n][1] for n in names_list]
ax.bar(x-0.2, oset_Ec, 0.35, label='OSET', color=colors, alpha=0.85)
ax.bar(x+0.2, dft_mid, 0.35, label='DFT', color='silver', edgecolor='black', alpha=0.85)
err_lo2 = [m-lo for m,lo in zip(dft_mid,dft_lo)]
err_hi2 = [hi-m for m,hi in zip(dft_mid,dft_hi)]
ax.errorbar(x+0.2, dft_mid, yerr=[err_lo2,err_hi2], fmt='none', color='black', capsize=4)
ax.set_xticks(x); ax.set_xticklabels(names_list)
ax.set_ylabel('$E_{\\rm core}/(2b)$ (eV/Å)'); ax.set_title('Core Energy per Unit Length')
ax.legend()

plt.tight_layout()
plt.savefig('fig_materials_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: fig_materials_comparison.png")
"""

# ── Insert before summary markdown (cell index 32) ────────────────────────────
cells = nb['cells']
# Find the summary markdown cell (contains "## 10. Summary Table")
idx_summary_md = None
for i, c in enumerate(cells):
    if c['cell_type'] == 'markdown' and '## 10. Summary' in ''.join(c['source']):
        idx_summary_md = i
        break
print(f"Summary markdown found at index {idx_summary_md}")

# Insert two cells before the summary markdown
cells.insert(idx_summary_md, code_cell(cell_code))
cells.insert(idx_summary_md, md_cell(cell_md))

with open('OSET_derivations.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Notebook now has {len(cells)} cells.")
# Print the new cell positions
for i, c in enumerate(cells):
    src = ''.join(c['source'])[:60].replace('\n', ' ')
    print(f"  [{i:2d}] {c['cell_type'][:4]:4s} : {src[:58]}")
