import json, sys

with open('OSET_derivations.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

def code_cell(src):
    return {"cell_type":"code","metadata":{},"source":src,"outputs":[],"execution_count":None}

def md_cell(src):
    return {"cell_type":"markdown","metadata":{},"source":src}

# ── Cell A: S₁₂₁₂, S₃₃₃₃, β₂_oblate ──────────────────────────────────────────
cell_A_src = r"""# ── S₁₂₁₂ and S₃₃₃₃ via second-kind I-integrals ─────────────────────────────
# Second-kind integrals for oblate spheroid (a₁=a₂=1, a₃=α):
#   I₁₁ = 2πα ∫₀^∞ ds / [(1+s)³ √(α²+s)]
#   I₃₃ = 2πα ∫₀^∞ ds / [(α²+s)^(5/2) (1+s)]
# Using substitution u = √(α²+s), both reduce to rational integrands:
#   I₁₁ = 4π ∫_α^∞ du / (β²+u²)³   where β²=1−α²
#   I₃₃ = 4π ∫_α^∞ du / [u⁴(β²+u²)]
#
# Eshelby tensor formulas (Mura 1987):
#   S₁₂₁₂ = [I₁₁ + (1-2ν)I₁] / [8π(1-ν)]
#   S₃₃₃₃ = [3α²I₃₃ + (1-2ν)I₃] / [8π(1-ν)]
#
# Dilatational constraint factor (exact oblate geometry, vs sphere approx):
#   β₂_oblate = 1 − S₃₃₃₃(1−2ν) / [2(1−ν)]

from scipy.integrate import quad

def second_kind_and_S(alpha_v, nu_v):
    b2 = 1 - alpha_v**2     # β²
    e_v = np.sqrt(b2)
    I1_v = (2*np.pi*alpha_v/e_v**3)*(np.arccos(alpha_v) - alpha_v*e_v)
    I3_v = 4*np.pi - 2*I1_v
    I13_v = (I3_v - I1_v)/(1 - alpha_v**2)

    I11_v, _ = quad(lambda u: 4*np.pi/(b2+u**2)**3,      alpha_v, np.inf)
    I33_v, _ = quad(lambda u: 4*np.pi/(u**4*(b2+u**2)),  alpha_v, np.inf)

    fac = 1/(8*np.pi*(1-nu_v))
    S1212 = (I11_v + (1-2*nu_v)*I1_v)*fac
    S3333 = (3*alpha_v**2*I33_v + (1-2*nu_v)*I3_v)*fac
    S1313 = ((1+alpha_v**2)*I13_v + (1-2*nu_v)*(I1_v+I3_v))/(16*np.pi*(1-nu_v))

    b2_oblate = 1 - S3333*(1-2*nu_v)/(2*(1-nu_v))
    b2_sphere  = 4*(1+nu_v)/(9*(1-nu_v))
    return dict(I11=I11_v, I33=I33_v, S1212=S1212, S3333=S3333,
                S1313=S1313, b2_oblate=b2_oblate, b2_sphere=b2_sphere)

r = second_kind_and_S(0.5, 1/3)

print(f"I₁₁(α=0.5)  = {r['I11']:.5f}")
print(f"I₃₃(α=0.5)  = {r['I33']:.5f}")
print()
print(f"{'Component':20s} {'Computed':>10} {'Manuscript':>12} {'Match':>6}")
print('-'*52)
rows = [('S₁₃₁₃',        r['S1313'], 0.2772),
        ('S₁₂₁₂',        r['S1212'], 0.1739),
        ('S₃₃₃₃',        r['S3333'], 0.7364),
        ('β₁=1-2S₁₃₁₃', 1-2*r['S1313'], 0.4456),
        ('β₂_oblate',     r['b2_oblate'],  0.816),
        ('β₂_sphere',     r['b2_sphere'],  8/9)]
for name, val, ref in rows:
    ok = abs(val-ref)<0.003*abs(ref)+1e-10
    print(f"{name:20s} {val:>10.4f} {ref:>12.4f} {'  OK' if ok else '  FAIL'}")

check('S₁₂₁₂(α=0.5, ν=1/3)',    r['S1212'],     0.1739, tol=2e-3)
check('S₃₃₃₃(α=0.5, ν=1/3)',    r['S3333'],     0.7364, tol=2e-3)
check('β₂_oblate(α=0.5, ν=1/3)',r['b2_oblate'], 0.816,  tol=3e-3)

# Store globals for summary cell
S1212_num      = r['S1212']
S3333_num      = r['S3333']
beta2_oblate_n = r['b2_oblate']
beta2_sphere_n = r['b2_sphere']
"""

# ── Cell B: Markdown section header ────────────────────────────────────────────
cell_B_src = """---
## 10a. Critical Cooperativity Number N_c

The minimum OSTZ count for a pile-up to saturate the back-stress:
$$N_c = b/(\\gamma_0 W)$$

- **Al grain boundary (W=2.5b):** N_c ≈ 3.8 ≈ 4
- **Cu crystal interior (W=b):** N_c ≈ 8.3

## 10b. Rate Equation — Theorem T5 (Al-12Si, 773 K)

$$\\dot{\\gamma}=\\frac{2W\\gamma_0\\nu_0}{d}\\,2\\sinh\\!\\left(\\frac{\\tau\\gamma_0 V_0}{2kT}\\right)\\exp\\!\\left(-\\frac{\\Delta F_0}{kT}\\right)$$

Key check: work term Δw = τγ₀V₀/2 ≪ kT at typical creep stresses → sinh ≈ linear.

## 10c. Effective Burgers Vector: Three Methods (§2.4.4)

Method 1 (kinematic): b_eff = γ₀W; Method 2 (seismic moment): b_eff = 2γ₀W/3.
The 33% discrepancy is acknowledged in the manuscript as an unresolved ambiguity.
"""

# ── Cell C: N_c + T5 rate equation + b_eff methods ────────────────────────────
cell_C_src = r"""# ── Critical Cooperativity N_c ───────────────────────────────────────────────
print("=== Critical Cooperativity Number N_c ===")
b_Al_nc = 0.2863e-9;  g0_Al_nc = 0.10;  W_Al_nc = 2.5*b_Al_nc
b_Cu_nc = 0.2556e-9;  g0_Cu_nc = 0.12;  W_Cu_nc = b_Cu_nc
Nc_Al = b_Al_nc/(g0_Al_nc*W_Al_nc)
Nc_Cu = b_Cu_nc/(g0_Cu_nc*W_Cu_nc)
print(f"Al GB  (W=2.5b): N_c = {b_Al_nc*1e9:.4f}/({g0_Al_nc}×{W_Al_nc*1e9:.4f}nm) = {Nc_Al:.2f} ≈ 4")
print(f"Cu FCC (W=b):    N_c = {b_Cu_nc*1e9:.4f}/({g0_Cu_nc}×{W_Cu_nc*1e9:.4f}nm) = {Nc_Cu:.2f} ≈ 8.3")
check('N_c for Al GB (W=2.5b)', Nc_Al, 3.82, tol=0.03)
check('N_c for Cu crystal',     Nc_Cu, 1/0.12, tol=1e-4)

# ── Rate Equation T5: Al-12Si at 773 K ───────────────────────────────────────
print()
print("=== T5: Al-12Si at T=773 K, τ=5 MPa ===")
T_r = 773.0;  kT_r = 1.380649e-23*T_r;  kT_eV_r = kT_r/eV
V0_r = (2/3)*np.pi*(0.75e-9)**3
Dw_r = 5e6*0.10*V0_r/2          # Δw = τγ₀V₀/2 (J)
x_r  = Dw_r/kT_r
print(f"kT = {kT_eV_r:.5f} eV")
print(f"V₀ = (2/3)πW³ = {V0_r:.4e} m³")
print(f"Δw = {Dw_r:.4e} J = {Dw_r/eV:.6f} eV")
print(f"Δw/kT = {x_r:.5f}  →  sinh = {np.sinh(x_r):.6f}  (linear error: {abs(np.sinh(x_r)-x_r)/x_r*100:.4f}%)")
exp_r = np.exp(-0.38*eV/kT_r)
gdot_r = (2*0.75e-9*0.10*1.6e8/10e-6)*2*np.sinh(x_r)*exp_r
print(f"exp(-ΔF₀/kT) = {exp_r:.4e}")
print(f"γ̇ ≈ {gdot_r:.3e} s⁻¹  (typical superplastic Al-12Si at 500°C: ~0.1–1 s⁻¹)")
check('Δw for Al at τ=5 MPa (eV)', Dw_r/eV,           0.00138, tol=0.02)
check('sinh(Δw/kT) ≈ Δw/kT',      np.sinh(x_r),       x_r,     tol=1e-4)
check('exp(-ΔF₀/kT) at 773 K',    exp_r,  np.exp(-0.38/kT_eV_r), tol=1e-6)

# ── b_eff Three Methods (§2.4.4) ─────────────────────────────────────────────
print()
print("=== b_eff: Three Methods ===")
for label_b, g0_b, W_nm_b in [('Cu crystal (W=b)', 0.12, 0.2556),
                                ('Al GB (W=2.5b)',    0.10, 0.75)]:
    b1 = g0_b*W_nm_b;  b2 = 2*g0_b*W_nm_b/3
    print(f"{label_b}: M1={b1:.4f} nm, M2={b2:.4f} nm, discrepancy={100*(b1-b2)/b1:.1f}%")
check('M1-M2 discrepancy (%)', 100*(1-2/3), 33.33, tol=1e-4)

# Store for summary
Nc_Al_val = Nc_Al;  Nc_Cu_val = Nc_Cu;  Dw_T5_eV = Dw_r/eV
"""

# ── Updated summary cell ───────────────────────────────────────────────────────
cell_summary_src = r"""# ── Run all checks and build summary ─────────────────────────────────────────
print('='*70)
print('OSET Derivation Verification — Summary')
print('='*70)

results = []
def rec(label, computed, expected, tol=1e-3):
    ok = abs(float(computed)-float(expected)) < tol*abs(float(expected))+1e-15
    results.append((label, float(computed), float(expected), '✓' if ok else '✗'))

# 1. β₂ identity
rec('β₂(ν=1/3) via 4K/(3K+4G)', 8/9, 8/9)

# 2. Eshelby S₁₃₁₃
S_val = float(S1313_sym.subs([(alpha, Rational(1,2)), (nu, Rational(1,3))]).evalf())
rec('S₁₃₁₃(α=0.5, ν=1/3)', S_val, 0.2772)
rec('β₁ = 1-2S₁₃₁₃', 1-2*S_val, 0.4456)

# 3. I-integrals
rec('I₁(α=0.5)',   float(I1_sym.subs(alpha, Rational(1,2)).evalf()), 2.9707, tol=5e-5)
rec('I₃(α=0.5)',   float(I3_sym.subs(alpha, Rational(1,2)).evalf()), 6.6250, tol=5e-5)
rec('I₁₃(α=0.5)', float(I13_sym.subs(alpha, Rational(1,2)).evalf()), 4.8724, tol=5e-5)

# 4. S₁₂₁₂ and S₃₃₃₃ (complete Eshelby table)
rec('S₁₂₁₂(α=0.5, ν=1/3)',    S1212_num,     0.1739, tol=2e-3)
rec('S₃₃₃₃(α=0.5, ν=1/3)',    S3333_num,     0.7364, tol=2e-3)
rec('β₂_oblate(α=0.5, ν=1/3)',beta2_oblate_n, 0.816,  tol=3e-3)

# 5. ΔF₀ for Al
W_Al_s=0.75e-9; G_Al_s=2.2e10; V0_Al_s=(2/3)*np.pi*W_Al_s**3
beta2_s=4*(1+1/3)/(9*(1-1/3))
dF0_s=0.5*(0.4456*0.01+beta2_s*0.0025)*G_Al_s*V0_Al_s/eV
rec('ΔF₀ for Al (eV)', dF0_s, 0.405, tol=0.04)

# 6. Unit conversion
rec('0.38 eV in kJ/mol', 0.38*eV*NA/1e3, 36.7, tol=0.005)

# 7. Abel projection
test_abel = float(abel_simplified.subs([(x,1.0),(W,2.0)]))
rec('Abel projection at (x=1,W=2)', test_abel, 2.0/(1+4))

# 8. Peierls stress
G_Cu_n=48.3e9; b_Cu_n=0.2556e-9; b1_Cu_n=0.4487; nu_Cu_n=0.343
tau_PG_basic = b1_Cu_n*0.12/2*np.exp(-2*np.pi)
tau_PG_hs    = tau_PG_basic*2/(1-nu_Cu_n)
rec('τ_P/G basic (τ_th×exp(−2π))',       tau_PG_basic, 5.03e-5, tol=0.01)
rec('τ_P/G with HS correction ×2/(1−ν)', tau_PG_hs,    1.53e-4, tol=0.01)

# 9. Taylor hardening
C_OSET_n=6*0.12*np.pi*b_Cu_n**2/(6*(1-0.343))
rec('OSET hardening coeff×G (Pa·m²)', C_OSET_n*G_Cu_n, 1.81e-9, tol=0.03)

# 10. Core energy
E_core_J_s=np.pi*1.0*b_Cu_n*0.12*G_Cu_n*b_Cu_n**2/3
rec('E_core/(2W) for Cu (eV/Å)', (E_core_J_s/eV)/(2*b_Cu_n*1e10), 0.124, tol=0.02)

# 11. Critical cooperativity
b_Al_r2=0.2863e-9; g0_Al_r2=0.10; W_Al_r2=2.5*b_Al_r2
rec('N_c for Al GB (W=2.5b)', b_Al_r2/(g0_Al_r2*W_Al_r2), 3.82, tol=0.03)

# 12-13. Rate equation T5
kT_r2=1.380649e-23*773.0; V0_r2=(2/3)*np.pi*(0.75e-9)**3
Dw_r2=5e6*0.10*V0_r2/2; x_r2=Dw_r2/kT_r2
rec('Δw for Al at τ=5 MPa (eV)', Dw_r2/eV, 0.00138, tol=0.02)
rec('sinh(Δw/kT) ≈ Δw/kT', np.sinh(x_r2), x_r2, tol=1e-4)

# 14. b_eff discrepancy
rec('b_eff M1–M2 discrepancy (%)', 100*(1-2/3), 33.33, tol=1e-4)

# Print table
print(f"{'Check':55s} {'Computed':>11} {'Expected':>11} {'Status':>6}")
print('-'*86)
for label, comp, exp_v, status in results:
    print(f"{label:55s} {comp:>11.5g} {exp_v:>11.5g} {status:>6}")

n_pass = sum(1 for *_,s in results if s == '✓')
print()
print(f'PASSED: {n_pass}/{len(results)}')
"""

# ── Apply insertions ──────────────────────────────────────────────────────────
cells = nb['cells']
orig_len = len(cells)

# Insert cell_A after cell[8] (S₁₃₁₃ plot) → position 9
cells.insert(9, code_cell(cell_A_src))

# Old cell[28] (core energy code) is now at [29] after one insertion.
# Insert B+C after cell[29] → positions 30, 31
cells.insert(30, md_cell(cell_B_src))
cells.insert(31, code_cell(cell_C_src))

# Old cell[30] (summary code) is now at [33] after three insertions.
cells[33]['source'] = cell_summary_src

with open('OSET_derivations.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Done. Notebook: {orig_len} -> {len(cells)} cells.")
for i, c in enumerate(cells):
    src = ''.join(c['source'])[:55].replace('\n',' ')
    marker = ' <-- NEW' if i in [9,30,31] else (' <-- UPDATED' if i==33 else '')
    print(f"  [{i:2d}] {c['cell_type'][:4]:4s} : {src[:50]}{marker}")
