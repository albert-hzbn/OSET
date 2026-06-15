import json

sec12_md = r"""---
## 12. Grain-Boundary Free Energy from OSET

The OSTZ is an oblate spheroid with semi-major axis $R_0 = W$ and $\alpha = c/R_0 = 0.5$, giving:

$$V_0 = \tfrac{4}{3}\pi R_0^2 c = \tfrac{4}{3}\pi\alpha W^3 = \tfrac{2}{3}\pi W^3 \qquad A_0 = \pi R_0^2 = \pi W^2$$

$$\boxed{\gamma_\text{GB}^\text{OSET} = \frac{U}{A_0} = \frac{1}{3}\!\left(\beta_1\gamma_0^2 + \beta_2\varepsilon_0^2\right)GW}$$

Three boundary types are compared with literature:

| Type | $W$ | $\gamma_0$ | Formula |
|------|-----|-----------|---------|
| Coherent twin (CTB) | $b$ | 0.12 | $\gamma_0^2 Gb/4$ ($\beta_1^\mathrm{eff}=1$, SFE convention) |
| Low-angle tilt $\theta$ | $b$ | 0.12 | $(\theta/\theta_\mathrm{sat})\cdot\tfrac{1}{3}(\beta_1\gamma_0^2+\beta_2\varepsilon_0^2)Gb$ |
| Random HAGB | $2.5b$ | 0.10 | $\tfrac{1}{3}(\beta_1\gamma_0^2+\beta_2\varepsilon_0^2)G\cdot 2.5b$ |"""

sec12_code = r"""# ─────────────────────────────────────────────────────────────────────────────
# SS12  Grain-Boundary Free Energy from OSET
# ─────────────────────────────────────────────────────────────────────────────
# Derivation (alpha=0.5, R0=W):
#   V0 = (4/3)*pi*R0^2*c = (4/3)*pi*alpha*W^3 = (2/3)*pi*W^3
#   A0 = pi*R0^2 = pi*W^2
#   U  = 0.5*(b1*g0^2 + b2*e0^2)*G*V0
#   gamma_GB = U/A0 = (1/3)*(b1*g0^2 + b2*e0^2)*G*W

# SS12.1  Derivation spot-checks (pure geometry, no material data) -------------
alpha_ostz = 0.5
W_unit     = 1.0   # normalised to 1 for algebra checks

V0_spheroid = (4/3)*np.pi * alpha_ostz * W_unit**3   # oblate spheroid formula
V0_shortcut = (2/3)*np.pi * W_unit**3                # simplified for alpha=0.5
A0_val      = np.pi * W_unit**2
V0_A0_ratio = V0_shortcut / A0_val                   # should be (2/3)*W

b1_ref = 0.4456;  g0_ref = 0.12;  G_ref = 1.0
U_over_A0   = 0.5 * b1_ref * g0_ref**2 * G_ref * V0_shortcut / A0_val
formula_val = (1/3) * b1_ref * g0_ref**2 * G_ref * W_unit

check('V0 = (4/3)*pi*alpha*W^3 = (2/3)*pi*W^3  [alpha=0.5]',
      V0_spheroid, V0_shortcut, tol=1e-12)
check('V0/A0 = (2/3)*W  [geometric ratio, alpha=0.5]',
      V0_A0_ratio, (2/3)*W_unit, tol=1e-12)
check('U/A0 = (1/3)*b1*g0^2*G*W  [energy-per-area derivation]',
      U_over_A0, formula_val, tol=1e-12)

# SS12.2  OSET GB energy function ----------------------------------------------
def gamma_GB_oset(b1_v, b2_v, g0_v, e0_v, G_v, W_v):
    # gamma_GB = (1/3)*(b1*g0^2 + b2*e0^2)*G*W  [J/m^2]
    return (1/3) * (b1_v*g0_v**2 + b2_v*e0_v**2) * G_v * W_v

g0_cry  = 0.12;  g0_gb = 0.10;  e0_v = 0.05
theta_LA  = np.radians(10.0)
theta_sat = np.radians(17.0)    # Read-Shockley saturation angle ~17 deg

# SS12.3  Literature GB energies (J/m^2): (lo, hi) ranges ---------------------
# CTB = coherent Sigma3 twin ~ ISF: Hirth & Lothe 1982; Murr 1975
gb_ctb_lit = {'Cu':(0.041,0.050),'Al':(0.130,0.175),'Ni':(0.110,0.140),
              'Ag':(0.013,0.020),'Au':(0.028,0.036),'Fe':None,'W':None}
# Low-angle tilt at ~10 deg: Olmsted et al. 2009 MD; Sutton & Balluffi 1995
gb_la_lit  = {'Cu':(0.200,0.400),'Al':(0.150,0.300),'Ni':(0.250,0.500),
              'Ag':(0.150,0.300),'Au':(0.150,0.300),
              'Fe':(0.350,0.600),'W':(0.500,0.900)}
# Random HAGB: Rohrer 2011; Olmsted et al. 2009; Sutton & Balluffi 1995
gb_ha_lit  = {'Cu':(0.500,0.900),'Al':(0.300,0.600),'Ni':(0.700,1.100),
              'Ag':(0.350,0.650),'Au':(0.350,0.650),
              'Fe':(0.800,1.200),'W':(1.000,2.000)}

# SS12.4  Compute predictions --------------------------------------------------
gb12 = {}
for nm_v, (struct_v, b_v, G_v, nu_v) in metals_db.items():
    b1_v = 1 - 2*S1313v(nu_v)
    b2_v = 4*(1+nu_v)/(9*(1-nu_v))

    # CTB: beta1_eff = 1, W = b  (SFE formula, same as SS10.2)
    g_ctb = g0_cry**2 * G_v * b_v / 4

    # LA: partial OSTZ coverage proportional to theta/theta_sat
    g_la_full = gamma_GB_oset(b1_v, b2_v, g0_cry, e0_v, G_v, b_v)
    g_la      = g_la_full * (theta_LA / theta_sat)

    # HA: GB regime W = 2.5b, gamma0 = 0.10
    g_ha = gamma_GB_oset(b1_v, b2_v, g0_gb, e0_v, G_v, 2.5*b_v)

    gb12[nm_v] = dict(ctb=g_ctb, la=g_la, ha=g_ha, b1=b1_v, b2=b2_v,
                      G=G_v, b=b_v)

# SS12.5  Print table ----------------------------------------------------------
def rng_s(d):
    return f"{d[0]*1e3:.0f}-{d[1]*1e3:.0f}" if d else " N/A"
def flg(v, d):
    if d is None: return '-'
    return 'V' if d[0]<=v<=d[1] else '~' if d[0]/6<=v<=d[1]*6 else 'X'

print("SS12  OSET Grain-Boundary Free Energy vs. Literature  (mJ/m^2)")
print("="*108)
print(f"{'Metal':5s} | {'CTB OSET':8s} {'CTB lit':10s}  | {'LA(10deg) OSET':14s} {'LA lit':10s}  | {'HA OSET':8s} {'HA lit':10s}")
print('-'*108)
for nm_v in metals_db:
    r = gb12[nm_v]
    print(f"{nm_v:5s} | {r['ctb']*1e3:8.1f} {rng_s(gb_ctb_lit[nm_v]):10s} {flg(r['ctb'],gb_ctb_lit[nm_v])}"
          f"  | {r['la']*1e3:14.1f} {rng_s(gb_la_lit[nm_v]):10s} {flg(r['la'],gb_la_lit[nm_v])}"
          f"  | {r['ha']*1e3:8.1f} {rng_s(gb_ha_lit[nm_v]):10s} {flg(r['ha'],gb_ha_lit[nm_v])}")
print()
print("  OSET formula captures the correct G*b scaling (shear elastic component only).")
print("  Systematic underprediction by 5-15x for LA/HA GBs reflects missing")
print("  core/misfit contributions. CTB ~ SFE: within 2-3x for FCC noble metals.")
print()

# SS12.6  Spot-checks ----------------------------------------------------------
# Geometry derivation
check('V0 formula (alpha=0.5): (4/3)*pi*alpha*W^3 == (2/3)*pi*W^3',
      (4/3)*np.pi*0.5*2.0**3, (2/3)*np.pi*2.0**3, tol=1e-12)
check('V0/A0 == (2/3)*W for any W',
      (2/3)*np.pi*3.0**3 / (np.pi*3.0**2), (2/3)*3.0, tol=1e-12)

# CTB must equal SFE from SS11 (same formula, same parameters)
check('gamma_CTB Cu == gamma_SF Cu  (mJ/m^2)', gb12['Cu']['ctb']*1e3, res['Cu']['sfe'], tol=1e-5)
check('gamma_CTB Al == gamma_SF Al  (mJ/m^2)', gb12['Al']['ctb']*1e3, res['Al']['sfe'], tol=1e-5)
check('gamma_CTB Ni == gamma_SF Ni  (mJ/m^2)', gb12['Ni']['ctb']*1e3, res['Ni']['sfe'], tol=1e-5)

# HA hand-computed values (verified externally)
check('gamma_HA Cu  (mJ/m^2)',  gb12['Cu']['ha']*1e3,  69.5,  tol=0.01)
check('gamma_HA Fe  (mJ/m^2)',  gb12['Fe']['ha']*1e3, 107.6,  tol=0.01)
check('gamma_HA W   (mJ/m^2)',  gb12['W']['ha']*1e3,  230.7,  tol=0.01)

# G*b scaling: gamma_HA / (G*b) ~ constant across metals (ν variation < 2%)
Cu_ratio = gb12['Cu']['ha'] / (gb12['Cu']['G'] * gb12['Cu']['b'])
Al_ratio = gb12['Al']['ha'] / (gb12['Al']['G'] * gb12['Al']['b'])
check('gamma_HA/(G*b) ratio Cu ~ Al  [Gb scaling]', Cu_ratio, Al_ratio, tol=0.02)

# Monotonicity and positivity
for nm_v in metals_db:
    assert gb12[nm_v]['ha'] > gb12[nm_v]['la'] > 0, \
        f"Expected gamma_HA > gamma_LA > 0 for {nm_v}"
print("V PASS  Monotonicity: gamma_HA > gamma_LA > 0  for all 7 metals")

# W -> 0 limit
assert abs(gamma_GB_oset(0.45, 0.89, 0.12, 0.05, 1e10, 0.0)) < 1e-20
print("V PASS  gamma_GB -> 0 as W -> 0  [correct limiting behaviour]")

# SS12.7  Plot -----------------------------------------------------------------
nm_lst  = list(metals_db.keys())
str_lst = [metals_db[n][0] for n in nm_lst]
clr_lst = ['royalblue' if s=='FCC' else 'tomato' for s in str_lst]
xp_12   = np.arange(len(nm_lst));  bw_12 = 0.35

fig12, axes12 = plt.subplots(1, 3, figsize=(15, 5))
fig12.suptitle('OSET Grain-Boundary Free Energy vs. Literature',
               fontweight='bold', fontsize=12)
panels_12 = [
    ('CTB / SFE',           'ctb', gb_ctb_lit),
    ('Low-Angle Tilt 10deg', 'la',  gb_la_lit),
    ('Random HAGB',           'ha',  gb_ha_lit),
]
for ax12, (ttl, key_v, lit_v) in zip(axes12, panels_12):
    oset_v  = [gb12[n][key_v]*1e3 for n in nm_lst]
    lit_mid  = [(l[0]+l[1])/2*1e3 if l else 0 for l in [lit_v[n] for n in nm_lst]]
    half_rng = [(l[1]-l[0])/2*1e3 if l else 0 for l in [lit_v[n] for n in nm_lst]]
    ax12.bar(xp_12-bw_12/2, oset_v,  bw_12, label='OSET',       color=clr_lst, alpha=0.85)
    ax12.bar(xp_12+bw_12/2, lit_mid, bw_12, label='Literature',
             color='silver', edgecolor='k', alpha=0.85)
    ax12.errorbar(xp_12+bw_12/2, lit_mid, yerr=[half_rng, half_rng],
                  fmt='none', color='k', capsize=4)
    ax12.set_xticks(xp_12); ax12.set_xticklabels(nm_lst)
    ax12.set_ylabel(r'$\gamma_{GB}$ (mJ/m$^2$)'); ax12.set_title(ttl)
    ax12.legend(fontsize=8)
plt.tight_layout()
plt.savefig('fig_GB_energy.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: fig_GB_energy.png")"""

# ── Insert cells at index 34 ─────────────────────────────────────────────────
with open('OSET_derivations.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

def make_md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src}

def make_code(src):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src}

nb['cells'].insert(34, make_md(sec12_md))
nb['cells'].insert(35, make_code(sec12_code))

print(f"Total cells: {len(nb['cells'])}")
for i in [33, 34, 35, 36]:
    src_p = ''.join(nb['cells'][i]['source'])[:55].encode('ascii','replace').decode('ascii')
    print(f"  Cell {i} [{nb['cells'][i]['cell_type'][:4]}]: {src_p}")

with open('OSET_derivations.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
print("Notebook saved.")
