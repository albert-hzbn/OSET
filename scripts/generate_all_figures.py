"""
generate_all_figures.py
========================
Regenerates every figure used in the OSET manuscript / Supporting Information
at consistent, publication-quality settings, and extends the grain-boundary
verification to additional (HCP) crystal systems plus a consolidated
application of the Unified Rate Equation across material classes.

All literature numbers used below are taken from the references already
cited in the manuscript/SI, or (for the two new HCP metals, Mg and Ti) from
published DFT/experimental values cited in comments next to each entry.
No numbers are invented; ranges/medians are reported as found.

Run with:  python generate_all_figures.py
Outputs go to ../manuscript/figures/
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "manuscript", "figures")
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Global publication-quality style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 400,
    "font.size": 13,
    "axes.titlesize": 14.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 13,
    "legend.fontsize": 10.5,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.7",
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "axes.linewidth": 1.0,
    "axes.edgecolor": "0.3",
    "lines.linewidth": 2.0,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.06,
    "figure.facecolor": "white",
})

# Crystal-structure palette is used ONLY to colour x-axis tick labels, never the
# bars, so the bars unambiguously encode OSET (one colour) vs literature (grey).
FCC_COLOR, BCC_COLOR, HCP_COLOR = "#2f6fc0", "#d6502a", "#1f8a4c"
OSET_COLOR = "#34495e"   # single dark slate-blue for every OSET bar
LIT_COLOR = "#c7c7c7"    # neutral grey for literature/reference bars


def struct_color(s):
    return {"FCC": FCC_COLOR, "BCC": BCC_COLOR, "HCP": HCP_COLOR}.get(s, "gray")


def _color_xticklabels(ax, structs):
    """Colour each x-tick label by crystal structure (bars stay OSET vs lit)."""
    for lbl, s in zip(ax.get_xticklabels(), structs):
        lbl.set_color(struct_color(s))
        lbl.set_fontweight("bold")


def _struct_legend(fig, structs=("FCC", "BCC", "HCP"), y=0.5, loc="center left"):
    """Shared crystal-structure colour legend (explains the x-label colours)."""
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="s", ls="", markersize=9,
                      markerfacecolor=struct_color(s), markeredgecolor="0.3",
                      label=f"{s} (axis label)") for s in structs]
    return fig.legend(handles=handles, frameon=False, fontsize=11,
                      loc=loc, bbox_to_anchor=(1.0, y), title="crystal\nstructure")


# ===========================================================================
# Shared Eshelby helper
# ===========================================================================
def S1313v(nu_v, alpha=0.5):
    e = np.sqrt(1 - alpha ** 2)
    I1 = (2 * np.pi * alpha / e ** 3) * (np.arccos(alpha) - alpha * e)
    I3 = 4 * np.pi - 2 * I1
    I13 = (I3 - I1) / (1 - alpha ** 2)
    return ((1 + alpha ** 2) * I13 + (1 - 2 * nu_v) * (I1 + I3)) / (16 * np.pi * (1 - nu_v))


eV = 1.602176634e-19

# ===========================================================================
# FIGURE 1: S_1313 vs aspect ratio
# ===========================================================================
def fig_S1313():
    alpha_arr = np.linspace(0.05, 0.999, 400)
    nu_list = [0.30, 1 / 3, 0.35, 0.40]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for nu_v, col in zip(nu_list, colors):
        S_vals = [S1313v(nu_v, a) for a in alpha_arr]
        label = f"$\\nu={nu_v:.3g}$" if abs(nu_v - 1 / 3) > 0.001 else "$\\nu=1/3$"
        ax.plot(alpha_arr, S_vals, color=col, lw=1.8, label=label)
    ax.axhline(0.5, ls="--", color="gray", lw=2.2, label="thin-disk limit (1/2)")
    ax.axvline(0.5, ls=":", color="black", lw=2.2, label=r"OSTZ $\alpha=0.5$")
    ax.set_xlabel(r"Aspect ratio $\alpha=c/a$")
    ax.set_ylabel(r"$S_{1313}$")
    ax.set_title(r"Eshelby tensor component $S_{1313}$ for an oblate spheroid")
    ax.set_ylim(top=0.565)
    ax.legend(frameon=True, loc="upper right", facecolor="white",
              framealpha=0.92, edgecolor="0.7")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "fig_S1313.pdf"))
    plt.close(fig)
    print("Saved fig_S1313.pdf")


# ===========================================================================
# FIGURE 2: OSET vs Taylor hardening
# ===========================================================================
def fig_Taylor():
    rho_arr = np.logspace(10, 17, 400)
    G_cu, b_cu, z = 48.3e9, 0.2556e-9, 6
    tau_OSET_coeff = z * G_cu * b_cu * 0.12 / 8
    tau_Taylor_coeff = 0.3 * G_cu * b_cu
    tau_OSET = tau_OSET_coeff * rho_arr / 1e6
    tau_Taylor = tau_Taylor_coeff * np.sqrt(rho_arr) / 1e6

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    ax = axes[0]
    ax.loglog(rho_arr, tau_OSET, lw=2, label=r"OSET: $\tau\propto\rho$ (back-stress)")
    ax.loglog(rho_arr, tau_Taylor, lw=2, ls="--", label=r"Taylor: $\tau\propto\sqrt{\rho}$ (forest)")
    ax.axvspan(1e10, 1e14, alpha=0.08, color="green", label="Typical range (annealed-CW)")
    ax.set_xlabel(r"Dislocation density $\rho$ (m$^{-2}$)")
    ax.set_ylabel("Flow-stress increment (MPa)")
    ax.set_title("OSET vs. Taylor hardening (Cu, $z=6$)")
    ax.legend(frameon=True, fontsize=10.5, loc="lower right",
              facecolor="white", framealpha=0.92, edgecolor="0.7")
    ax.set_ylim(1e-3, 2e8)

    ax2 = axes[1]
    rho_hi = np.logspace(14, 17, 300)
    ax2.plot(rho_hi * 1e-15, tau_OSET_coeff * rho_hi / 1e6, lw=2, label="OSET (linear)")
    ax2.plot(rho_hi * 1e-15, tau_Taylor_coeff * np.sqrt(rho_hi) / 1e6, lw=2, ls="--", label=r"Taylor ($\sqrt{\rho}$)")
    ax2.set_xlabel(r"$\rho$ ($\times10^{15}\,\mathrm{m}^{-2}$)")
    ax2.set_ylabel("Flow-stress increment (MPa)")
    ax2.set_title("High-density regime")
    ax2.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "fig_Taylor_comparison.pdf"))
    plt.close(fig)
    print("Saved fig_Taylor_comparison.pdf")


# ===========================================================================
# FIGURE 3: partition function vs incorrect sinh occupation
# ===========================================================================
def fig_partition():
    dw_arr = np.linspace(0, 2, 300)
    n_corr = np.exp(2 * dw_arr)
    n_wrong = 2 * np.sinh(dw_arr)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    ax1, ax2 = axes
    ax1.plot(dw_arr, n_corr, lw=2, label=r"correct: $e^{2\Delta w/kT}$")
    ax1.plot(dw_arr, n_wrong, lw=2, ls="--", label=r"incorrect: $2\sinh(\Delta w/kT)$")
    ax1.set_xlabel(r"$\Delta w/kT$")
    ax1.set_ylabel("relative OSTZ occupation")
    ax1.set_title("Partition function vs. sinh occupation")
    ax1.legend(frameon=False)

    ax2.semilogy(dw_arr, np.abs(n_corr - n_wrong) / n_corr, color="firebrick", lw=2)
    ax2.set_xlabel(r"$\Delta w/kT$")
    ax2.set_ylabel("relative error")
    ax2.set_title("Error of the incorrect sinh expression")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "fig_partition_vs_sinh.pdf"))
    plt.close(fig)
    print("Saved fig_partition_vs_sinh.pdf")


# ===========================================================================
# Material parameter databases (crystal-interior dislocation properties)
# ===========================================================================
metals_db = {
    "Cu": ("FCC", 0.2556e-9, 48.3e9, 0.343),
    "Al": ("FCC", 0.2863e-9, 26.2e9, 0.347),
    "Ni": ("FCC", 0.2492e-9, 76.0e9, 0.276),
    "Ag": ("FCC", 0.2889e-9, 30.3e9, 0.367),
    "Au": ("FCC", 0.2884e-9, 27.0e9, 0.440),
    "Fe": ("BCC", 0.2482e-9, 82.0e9, 0.291),
    "W": ("BCC", 0.2741e-9, 161.0e9, 0.280),
    # HCP additions for basal-plane SFE comparison (same G, b, nu as the
    # GB-energy verification): Mg G=17.0 GPa, nu=0.290; Ti G=41.4 GPa,
    # nu=0.340 (handbook polycrystalline values).
    "Mg": ("HCP", 0.320e-9, 17.0e9, 0.290),
    "Ti": ("HCP", 0.295e-9, 41.4e9, 0.340),
}

# Basal-plane stacking-fault energy (I1/I2-type), literature range, mJ/m^2:
#  Mg: DFT I1 = 14.4 (Wang et al., Mater. Res. Lett. 2 (2014) 29-36);
#      experimental estimates 50-78 reported in the same source -- the wide
#      DFT-vs-experiment spread is itself noted in the literature as a
#      known difficulty for Mg basal SFE. Range given as (DFT, experiment).
#  Ti: DFT I2 = 259-327, cross-validated across 5 independent first-principles
#      studies compiled in Yu, SCIREA J. Phys. 1(1) (2016) Table 2 (this
#      work's DFT=327; Wu et al., Appl. Surf. Sci. 256 (2010) 3409, I2=287;
#      Kwasniak et al., Phys. Rev. B 89 (2014) 144105, I2=307; Benoit &
#      Tarrat, Modelling Simul. Mater. Sci. Eng. 21 (2013) 015009, I2=259;
#      Domain & Legris (2004), I2=291).
# No fabricated/estimated values are used for BCC metals beyond Fe, W: BCC
# metals lack a single well-defined low-energy stacking fault analogous to
# FCC/HCP, and the literature "generalized planar fault energy" is highly
# method- and plane-dependent (see main text/SI discussion); Fe, W are kept
# in the comparison with no experimental SFE bar, exactly as before.
sfe_exp = {"Cu": (45, 45), "Al": (166, 166), "Ni": (125, 125), "Ag": (16, 16),
           "Au": (32, 32), "Fe": None, "W": None,
           "Mg": (14.4, 78), "Ti": (259, 327)}
tau_P_exp_G = {
    "Cu": (5e-5, 2e-4), "Al": (1e-6, 5e-5), "Ni": (1e-4, 5e-4),
    "Ag": (2e-5, 1e-4), "Au": (2e-5, 1e-4), "Fe": (4e-3, 2e-2), "W": (2e-2, 5e-2),
}
Ecore_dft = {
    "Cu": (0.05, 0.15), "Al": (0.03, 0.08), "Ni": (0.10, 0.25),
    "Ag": (0.03, 0.10), "Au": (0.04, 0.12), "Fe": (0.20, 0.50), "W": (0.50, 1.20),
}

e0 = 0.05

# Material-specific OSTZ parameters (root-cause fix): gamma0 is the ideal shear
# strain (gamma-surface maximum, Ogata et al., Science 298 (2002) 807) and W/b
# is the gamma-surface-controlled core-width factor (>1 wide planar FCC cores,
# <1 compact non-planar BCC cores). The legacy universal (0.12, W=b) is the Cu
# value misapplied to all metals -- the cause of the Al/Ag/BCC discrepancies.
# gamma0 = LITERATURE relaxed ideal shear strain s_m on the primary slip system
# (Ogata, Li, Hirosaki, Shibutani & Yip, Phys. Rev. B 70 (2004) 104104, Tab. II),
# used directly with no fit. W (structural) = b for the crystal interior; W_P is
# the Peierls glide-misfit width (FCC planar wide, BCC non-planar narrow).
ostz_params = {  # name: (gamma0_Ogata2004, Wb_struct, WbP_peierls)
    "Cu": (0.137, 1.0, 1.08), "Al": (0.200, 1.0, 1.55), "Ni": (0.140, 1.0, 1.18),
    "Ag": (0.145, 1.0, 1.11), "Au": (0.105, 1.0, 1.14),
    "Fe": (0.178, 1.0, 0.64), "W": (0.179, 1.0, 0.52),
    # HCP: Mg basal s_m=0.152; Ti only prismatic s_m=0.099 in Ogata Tab. II
    # (Ti's primary slip is prismatic, so basal SFE is underpredicted -- see SI).
    "Mg": (0.152, 1.0, 1.0), "Ti": (0.099, 1.0, 1.0),
}


def materials_results():
    res = {}
    for name, (struct, b, G, nu) in metals_db.items():
        b1_esh = 1 - 2 * S1313v(nu)
        b2 = 4 * (1 + nu) / (9 * (1 - nu))
        g0, Wb, WbP = ostz_params[name]
        W = Wb * b            # structural OSTZ width (energy, SFE, volume)
        WP = WbP * b          # Peierls glide-misfit width
        V0 = (2 / 3) * np.pi * W ** 3
        tau_P = (2 / (1 - nu)) * np.exp(-2 * np.pi * WP / (b * (1 - nu)))
        # SFE uses the (sqrt(3)-1) cluster prefactor with beta1_eff = 1
        sfe_oset = (np.sqrt(3) - 1) * g0 ** 2 * G * W / 3 * 1e3
        Ecore = np.pi * b1_esh * b * g0 * G * W ** 2 / 3
        Ecore_per2b = (Ecore / eV) / (2 * W * 1e10)
        dF0 = 0.5 * (b1_esh * g0 ** 2 + b2 * e0 ** 2) * G * V0 / eV
        res[name] = dict(struct=struct, b_nm=b * 1e9, G=G / 1e9, nu=nu,
                          b1_esh=b1_esh, Nc=1 / (g0 * Wb), tau_P_G=tau_P,
                          sfe=sfe_oset, Ecore_per2b=Ecore_per2b, dF0=dF0)
    return res


# ===========================================================================
# FIGURE 4: multi-metal comparison (SFE, Peierls, core energy)
# ===========================================================================
def fig_materials_comparison():
    res = materials_results()
    # SFE panel: only metals with a comparable literature SFE (5 FCC + 2 HCP).
    # BCC metals (Fe, W) are excluded: a BCC crystal has no single well-defined
    # low-energy stacking fault, so there is no unambiguous experimental value.
    # Peierls-stress and core-energy panels: restricted to the 7 FCC+BCC metals.
    names_sfe = [n for n in metals_db if sfe_exp.get(n) is not None]
    names_other = [n for n in metals_db if n not in ("Mg", "Ti")]
    struct_sfe = [res[n]["struct"] for n in names_sfe]
    struct_other = [res[n]["struct"] for n in names_other]
    x_sfe = np.arange(len(names_sfe))
    x_other = np.arange(len(names_other))
    bw = 0.38

    # Per-metal gamma0 values are NOT written on the plot (they crowd the axis);
    # they are listed with their source (Ogata et al. 2004) in SI Table S2.
    def _bars(ax, x, oset, mid, lo, hi, litlabel):
        ax.bar(x - bw / 2, oset, bw, color=OSET_COLOR, alpha=0.95,
               edgecolor="0.15", linewidth=0.6, label="OSET", zorder=3)
        ax.bar(x + bw / 2, mid, bw, color=LIT_COLOR, edgecolor="0.25",
               linewidth=0.6, hatch="//", label=litlabel, zorder=3)
        ax.errorbar(x + bw / 2, mid,
                    yerr=[[m - l for m, l in zip(mid, lo)], [h - m for h, m in zip(hi, mid)]],
                    fmt="none", ecolor="0.15", capsize=4, lw=1.2, zorder=4)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0))
    fig.suptitle(r"OSET vs. literature: crystal-interior dislocation properties "
                 r"(literature $\gamma_0$ = ideal shear strain, Ogata et al. 2004)",
                 fontweight="bold", fontsize=14.5)

    # --- Panel 1: SFE ---
    ax = axes[0]
    oset_sfe = [res[n]["sfe"] for n in names_sfe]
    mid = [(sfe_exp[n][0] * sfe_exp[n][1]) ** 0.5 for n in names_sfe]
    lo = [sfe_exp[n][0] for n in names_sfe]
    hi = [sfe_exp[n][1] for n in names_sfe]
    _bars(ax, x_sfe, oset_sfe, mid, lo, hi, "Literature (exp./DFT)")
    ax.set_xticks(x_sfe); ax.set_xticklabels(names_sfe)
    _color_xticklabels(ax, struct_sfe)
    ax.set_ylabel(r"$\gamma_\mathrm{SF}$ (mJ/m$^2$)"); ax.set_title("Stacking-fault energy")
    ax.legend(loc="upper left"); ax.set_ylim(0, max(max(oset_sfe), max(hi)) * 1.28)
    for xi, o, m in zip(x_sfe, oset_sfe, mid):
        ax.annotate(f"{o/m:.1f}$\\times$", (xi - bw / 2, o), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=10, color=OSET_COLOR, fontweight="bold")

    # --- Panel 2: Peierls stress ---
    ax = axes[1]
    oset_tP = [res[n]["tau_P_G"] for n in names_other]
    mid = [(tau_P_exp_G[n][0] * tau_P_exp_G[n][1]) ** 0.5 for n in names_other]
    lo = [tau_P_exp_G[n][0] for n in names_other]
    hi = [tau_P_exp_G[n][1] for n in names_other]
    _bars(ax, x_other, oset_tP, mid, lo, hi, "Literature range")
    ax.set_xticks(x_other); ax.set_xticklabels(names_other)
    _color_xticklabels(ax, struct_other)
    ax.set_yscale("log"); ax.set_ylabel(r"$\tau_P/G$"); ax.set_title("Peierls stress / shear modulus")
    ax.set_ylim(3e-7, 2e-1); ax.legend(loc="upper left")

    # --- Panel 3: core energy ---
    ax = axes[2]
    oset_Ec = [res[n]["Ecore_per2b"] for n in names_other]
    dft_mid = [(Ecore_dft[n][0] * Ecore_dft[n][1]) ** 0.5 for n in names_other]
    dft_lo = [Ecore_dft[n][0] for n in names_other]
    dft_hi = [Ecore_dft[n][1] for n in names_other]
    _bars(ax, x_other, oset_Ec, dft_mid, dft_lo, dft_hi, "DFT range")
    ax.set_xticks(x_other); ax.set_xticklabels(names_other)
    _color_xticklabels(ax, struct_other)
    ax.set_ylabel(r"$E_\mathrm{core}/(2W)$ (eV/Å)"); ax.set_title("Dislocation core energy")
    ax.set_ylim(0, max(max(oset_Ec), max(dft_hi)) * 1.15); ax.legend(loc="upper left")

    _struct_legend(fig, y=0.5)
    fig.tight_layout(rect=[0, 0, 0.94, 0.92])
    fig.savefig(os.path.join(OUTDIR, "fig_materials_comparison.pdf"))
    plt.close(fig)
    print("Saved fig_materials_comparison.pdf")


# ===========================================================================
# FIGURE 5: Grain-boundary energy, extended to HCP metals (Mg, Ti)
# ===========================================================================
# HCP additions -- literature sourced as noted; isotropic-elasticity inputs.
#  Mg: G=17.0 GPa, nu=0.290 (handbook polycrystalline values)
#  Ti: G=41.4 GPa, nu=0.340 (handbook polycrystalline values, alpha-Ti)
#  Mg {10-12} twin boundary energy: DFT 118 mJ/m^2, MD 125 mJ/m^2
#      (used as the coherent-boundary/CTB analogue; range 110-125 mJ/m^2)
#  Ti {10-12}/{11-2} twin boundary energy: DFT range ~100-300 mJ/m^2
#      (mode- and method-dependent; reported range used directly)
#  Mg average computed GB energy: DFT range 205-584 mJ/m^2 across distinct
#      boundaries (used as the HA analogue, since dedicated low-/high-angle
#      compilations as extensive as Olmsted (2009) for cubic metals are not
#      available for HCP metals)
#  Ti average GB energy: ~0.79 J/m^2 at 500 K (MD); range 0.70-1.00 J/m^2
#      adopted to allow for the lower-temperature value being somewhat higher
metals_db_gb = dict(metals_db)
metals_db_gb["Mg"] = ("HCP", 0.320e-9, 17.0e9, 0.290)
metals_db_gb["Ti"] = ("HCP", 0.295e-9, 41.4e9, 0.340)

gb_ctb_lit = {"Cu": (0.041, 0.050), "Al": (0.130, 0.175), "Ni": (0.110, 0.140),
              "Ag": (0.013, 0.020), "Au": (0.028, 0.036), "Fe": None, "W": None,
              "Mg": (0.110, 0.125), "Ti": (0.100, 0.300)}
gb_la_lit = {"Cu": (0.200, 0.400), "Al": (0.150, 0.300), "Ni": (0.250, 0.500),
             "Ag": (0.150, 0.300), "Au": (0.150, 0.300), "Fe": (0.350, 0.600), "W": (0.500, 0.900),
             "Mg": None, "Ti": None}
gb_ha_lit = {"Cu": (0.500, 0.900), "Al": (0.300, 0.600), "Ni": (0.700, 1.100),
             "Ag": (0.350, 0.650), "Au": (0.350, 0.650), "Fe": (0.800, 1.200), "W": (1.000, 2.000),
             "Mg": (0.205, 0.584), "Ti": (0.700, 1.000)}

theta_m = np.radians(15.0)
theta_LA = np.radians(10.0)
A_rs = 1.0 + np.log(theta_m)


def E0_rs(G_v, b_v, nu_v):
    return G_v * b_v / (4 * np.pi * (1 - nu_v))


def gamma_RS(theta, G_v, b_v, nu_v):
    return E0_rs(G_v, b_v, nu_v) * theta * (A_rs - np.log(theta))


def gb_results():
    out = {}
    for nm, (struct, b, G, nu) in metals_db_gb.items():
        E0 = E0_rs(G, b, nu)
        g0_m = ostz_params.get(nm, (0.12, 1.0))[0]   # material-specific ideal shear strain
        g_ctb = g0_m ** 2 * G * b / 4
        g_la = gamma_RS(theta_LA, G, b, nu) if gb_la_lit[nm] is not None else None
        g_ha = E0 * theta_m
        out[nm] = dict(struct=struct, ctb=g_ctb, la=g_la, ha=g_ha, E0=E0, G=G, b=b, nu=nu)
    return out


def fig_GB_energy():
    gb = gb_results()
    names = list(metals_db_gb.keys())
    colors = [struct_color(gb[n]["struct"]) for n in names]
    x = np.arange(len(names))
    bw = 0.36

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    fig.suptitle("OSET grain-boundary free energy vs. literature, FCC + BCC + HCP (dislocation-array model)",
                 fontweight="bold", fontsize=14.5)
    panels = [("Coherent twin / SF boundary", "ctb", gb_ctb_lit),
              (r"Low-angle tilt $10^\circ$", "la", gb_la_lit),
              ("Random high-angle boundary", "ha", gb_ha_lit)]
    for ax, (title, key, lit) in zip(axes, panels):
        # Only include metals that actually have a literature value for THIS
        # boundary type -- metals without an estimate are dropped from the panel
        # (no empty x-axis slots).
        present = [n for n in names if lit[n] is not None]
        oset_v = [gb[n][key] * 1e3 for n in present]
        lit_mid = [(lit[n][0] + lit[n][1]) / 2 * 1e3 for n in present]
        half_rng = [(lit[n][1] - lit[n][0]) / 2 * 1e3 for n in present]
        xt = np.arange(len(present))
        ax.bar(xt - bw / 2, oset_v, bw, label="OSET", color=OSET_COLOR,
               alpha=0.95, edgecolor="0.15", linewidth=0.6, zorder=3)
        ax.bar(xt + bw / 2, lit_mid, bw, label="Literature", color=LIT_COLOR,
               edgecolor="0.25", linewidth=0.6, hatch="//", zorder=3)
        ax.errorbar(xt + bw / 2, lit_mid, yerr=[half_rng, half_rng], fmt="none",
                    ecolor="0.15", capsize=4, lw=1.2, zorder=4)
        ax.set_xticks(xt); ax.set_xticklabels(present)
        _color_xticklabels(ax, [gb[n]["struct"] for n in present])
        ax.set_ylabel(r"$\gamma_\mathrm{GB}$ (mJ/m$^2$)"); ax.set_title(title)
        ax.legend(loc="upper left", fontsize=10.5)

    _struct_legend(fig, y=0.5)
    fig.tight_layout(rect=[0, 0, 0.965, 0.93])
    fig.savefig(os.path.join(OUTDIR, "fig_GB_energy.pdf"))
    plt.close(fig)
    print("Saved fig_GB_energy.pdf (now FCC+BCC+HCP, 9 systems)")
    return gb


# ===========================================================================
# Harisankar & Padmanabhan (2025) 41-system data (unchanged from SI Table S5)
# ===========================================================================
HP25 = [
    ("Zn22Al-2.5", "Zn", 2.5, 423, 0.558, 0.059, 0.034, 77.49, 51.38, 2902.33, 1.29, 0.262, 14.21),
    ("Zn22Al-0.9", "Zn", 0.9, 453, 0.597, 0.061, 0.035, 103.18, 49.39, 3614.76, 1.31, 0.248, 14.05),
    ("Al33Cu0.4Zr", "Al", 7.6, 713, 0.869, 0.066, 0.038, 170.35, 41.65, 146.16, 1.61, 0.332, 15.76),
    ("Al-MgScMn-3", "Al", 3.0, 723, 0.775, 0.089, 0.051, 135.44, 20.29, 2154.32, 1.57, 0.293, 15.57),
    ("Al-MgScMn-1", "Al", 1.0, 523, 0.561, 0.083, 0.048, 131.86, 22.78, 13365.50, 0.92, 0.255, 14.53),
    ("Al3Mg0.2Sc", "Al", 0.2, 573, 0.614, 0.086, 0.050, 118.56, 21.87, 5798.52, 1.46, 0.250, 14.77),
    ("Al-ZnMgSc", "Al", 0.7, 493, 0.528, 0.080, 0.046, 134.19, 25.27, 26323.70, 1.47, 0.249, 14.29),
    ("Al5Mg-24", "Al", 24.0, 748, 0.802, 0.090, 0.052, 141.75, 19.87, 1350.44, 0.87, 0.347, 19.57),
    ("Al17Si", "Al", 1.4, 763, 0.818, 0.089, 0.051, 116.48, 23.53, 9024.68, 0.88, 0.278, 15.78),
    ("Mg6Zn0.8Zr", "Mg", 0.7, 448, 0.485, 0.080, 0.046, 128.70, 18.04, 9296.55, 1.54, 0.245, 13.98),
    ("Mg4Y0.7Zr", "Mg", 2.0, 597, 0.647, 0.087, 0.050, 132.28, 15.21, 5021.66, 1.35, 0.273, 15.04),
    ("Mg5.8Zn1Y", "Mg", 17.5, 673, 0.728, 0.088, 0.051, 113.57, 15.38, 3528.13, 1.37, 0.358, 18.11),
    ("Ti6Al4V", "Ti", 0.9, 1023, 0.529, 0.116, 0.067, 124.29, 11.56, 4066.21, 1.14, 0.286, 16.73),
    ("Cu2.8Al-7", "Cu", 7.0, 723, 0.533, 0.096, 0.056, 106.59, 24.61, 8134.25, 1.48, 0.328, 15.70),
    ("Cu2.8Al-3", "Cu", 3.0, 673, 0.497, 0.095, 0.055, 107.17, 25.19, 6939.53, 0.99, 0.289, 15.41),
    ("IN836", "Ni", 2.5, 735, 0.652, 0.081, 0.047, 108.91, 33.01, 19591.86, 1.34, 0.289, 15.61),
    ("Ti43Al", "TiAl", 5.0, 1273, 0.695, 0.144, 0.088, 154.81, 7.23, 7394.89, 1.45, 0.341, 15.07),
    ("Ti48Al", "TiAl", 0.9, 1163, 0.654, 0.138, 0.080, 160.71, 7.73, 23813.29, 1.20, 0.291, 17.00),
    ("Ti46Al2Cr", "TiAl", 0.8, 1073, 0.599, 0.125, 0.072, 146.66, 9.39, 22138.44, 1.07, 0.287, 16.88),
    ("Co3Ti", "Co", 24.0, 1173, 0.776, 0.096, 0.055, 180.25, 23.78, 26514.51, 0.86, 0.389, 12.15),
    ("Ni3Si", "Ni", 15.0, 1323, 0.882, 0.092, 0.053, 147.16, 31.01, 6179.96, 0.87, 0.406, 10.97),
    ("ZrO2", "ZrO2", 0.07, 1273, 0.426, 0.082, 0.047, 112.41, 17.79, 13330.00, 1.65, 0.282, 17.56),
    ("ZrO2-3Y", "ZrO2", 0.51, 1523, 0.510, 0.088, 0.051, 122.47, 15.38, 1572.93, 0.85, 0.288, 17.30),
    ("ZrO2-4Y", "ZrO2", 0.75, 1573, 0.526, 0.089, 0.051, 126.12, 15.05, 592.37, 1.57, 0.291, 17.05),
    ("Al2O3-ZrO2-Si", "Ox", 0.4, 1673, 0.723, 0.067, 0.039, 144.83, 33.26, 1919.62, 0.86, 0.284, 17.16),
    ("Al2O3-NiAl-Zr", "Ox", 1.3, 1623, 0.702, 0.064, 0.037, 153.98, 40.41, 4690.84, 1.35, 0.298, 16.49),
    ("6061-20SiC", "Al", 0.8, 773, 0.590, 0.071, 0.041, 114.57, 37.50, 3954.22, 1.05, 0.272, 15.87),
    ("7075-20SiC", "Al", 5.0, 753, 0.596, 0.069, 0.040, 117.39, 39.33, 7220.88, 0.97, 0.314, 15.62),
    ("Zr65BMG", "Zr", 10., 653, 0.587, 0.074, 0.043, 124.71, 26.27, 13821.94, 1.35, 0.340, 16.22),
    ("Zr52BMG", "Zr", 10., 683, 0.621, 0.075, 0.043, 133.19, 27.60, 19572.74, 0.78, 0.343, 16.14),
    ("La55BMG", "La", 10., 483, 0.687, 0.059, 0.034, 139.26, 30.18, 14096.72, 1.60, 0.322, 16.57),
    ("La60BMG", "La", 10., 460, 0.713, 0.073, 0.042, 144.00, 30.35, 2170.54, 0.82, 0.320, 16.60),
    ("Fe72BMG", "Fe", 10., 863, 0.694, 0.079, 0.046, 112.65, 38.33, 29436.88, 0.78, 0.358, 15.48),
    ("Limestone", "Geo", 4.2, 973, 0.604, 0.073, 0.042, 152.98, 12.89, 13761.67, 0.57, 0.321, 15.80),
    ("AnDi-dry", "Geo", 3.1, 1323, 0.726, 0.062, 0.036, 165.03, 24.03, 10762.47, 1.12, 0.321, 15.89),
    ("AnDi-wet", "Geo", 3.1, 1223, 0.671, 0.062, 0.036, 142.75, 24.36, 22018.63, 1.27, 0.319, 16.03),
    ("Ice-10", "Ice", 10., 220, 0.806, 0.115, 0.067, 117.39, 6.65, 6227.10, 0.79, 0.289, 16.67),
    ("Ice-1700", "Ice", 1700., 241, 0.883, 0.117, 0.068, 144.41, 6.57, 1391.93, 1.44, 0.414, 2.30),
    ("Si3N4", "Ox", 0.07, 1723, 0.793, 0.085, 0.049, 140.67, 54.29, 5945.67, 1.05, 0.278, 17.36),
    ("ZrO2-Al2O3", "ZrO2", 0.06, 1623, 0.584, 0.069, 0.040, 109.91, 27.94, 2321.02, 0.94, 0.280, 17.55),
    ("ZrO2-Spinel", "ZrO2", 0.05, 1573, 0.574, 0.075, 0.043, 120.89, 21.78, 5912.42, 1.62, 0.281, 17.62),
]

b_nm = {"Zn": 0.267, "Al": 0.286, "Mg": 0.320, "Ti": 0.295, "TiAl": 0.285, "Cu": 0.256,
        "Ni": 0.250, "Co": 0.255, "Fe": 0.248, "Zr": 0.275, "La": 0.300, "ZrO2": 0.360,
        "Ox": 0.350, "Geo": 0.500, "Ice": 0.452}
nu_cls = {"Zn": 0.25, "Al": 0.345, "Mg": 0.29, "Ti": 0.32, "TiAl": 0.23, "Cu": 0.34,
          "Ni": 0.28, "Co": 0.31, "Fe": 0.29, "Zr": 0.36, "La": 0.36, "ZrO2": 0.30,
          "Ox": 0.25, "Geo": 0.25, "Ice": 0.33}
gB_exp = {"Zn": 0.34, "Al": 0.40, "Mg": 0.35, "Ti": 0.55, "TiAl": 0.50, "Cu": 0.60,
          "Ni": 0.87, "Co": 0.65, "Fe": 0.80, "Zr": 0.45, "La": 0.30, "ZrO2": 0.90,
          "Ox": 1.00, "Geo": 0.50, "Ice": 0.065}

NAv = 6.0221e23
beta1_13 = 1 - 2 * S1313v(1 / 3)
beta2_13 = 4 * (1 + 1 / 3) / (9 * (1 - 1 / 3))
theta_m13 = np.radians(15.0)


def hp25_results():
    res = []
    for lbl, cls, d_um, T, Thom, g0v, e0v, QoTm, GoTm, t0oTm, gB_lit, a_fit, Na_lit in HP25:
        b_v = b_nm.get(cls, 0.28) * 1e-9
        nu_v = nu_cls.get(cls, 0.30)
        Tm = T / Thom
        G_v = GoTm * Tm * 1e6
        V0 = (2 / 3) * np.pi * (2.5 * b_v) ** 3
        EF = beta1_13 * g0v ** 2 + beta2_13 * e0v ** 2
        dF0_eV = 0.5 * EF * G_v * V0 / eV
        Q_eV = QoTm * Tm / NAv / eV
        E0 = G_v * b_v / (4 * np.pi * (1 - nu_v))
        gB_oset = E0 * theta_m13
        res.append(dict(lbl=lbl, cls=cls, d=d_um, T=T, Tm=Tm, g0=g0v, e0=e0v, EF=EF,
                         G_GPa=G_v / 1e9, dF0=dF0_eV, Q_eV=Q_eV, gB_lit=gB_lit,
                         gB_oset=gB_oset, gB_exp=gB_exp.get(cls, 0.5), Na_lit=Na_lit))
    return res


def fig_HP2025_comparison():
    res13 = hp25_results()
    g0s = np.array([r["g0"] for r in res13]); e0s = np.array([r["e0"] for r in res13])
    dF0s = np.array([r["dF0"] for r in res13]); Qs = np.array([r["Q_eV"] for r in res13])
    gBl = np.array([r["gB_lit"] for r in res13]); gBo = np.array([r["gB_oset"] for r in res13])
    gBe = np.array([r["gB_exp"] for r in res13])
    ratG = gBl / gBo
    Th = np.array([r["T"] / r["Tm"] for r in res13])

    fig, ax = plt.subplots(2, 3, figsize=(14, 9.6))
    fig.suptitle("OSTZ theory vs. Harisankar & Padmanabhan (2025) -- 41 systems", fontsize=14.5, fontweight="bold")

    a0 = ax[0, 0]
    a0.hist(g0s, bins=14, color="steelblue", alpha=0.78, edgecolor="k")
    a0.axvline(0.12, color="red", ls="--", lw=2, label="OSTZ crystal 0.12")
    a0.axvline(0.10, color="darkorange", ls="--", lw=2, label="OSTZ GB 0.10")
    a0.axvline(g0s.mean(), color="navy", lw=2, label=f"mean {g0s.mean():.3f}")
    a0.set_xlabel(r"$\gamma_0$"); a0.set_ylabel("count"); a0.set_title(r"$\gamma_0$ distribution"); a0.legend(fontsize=9.5, frameon=False)

    a1 = ax[0, 1]
    a1.hist(e0s, bins=12, color="teal", alpha=0.78, edgecolor="k")
    a1.axvline(0.05, color="red", ls="--", lw=2, label="OSTZ 0.05")
    a1.axvline(e0s.mean(), color="navy", lw=2, label=f"mean {e0s.mean():.3f}")
    a1.set_xlabel(r"$\varepsilon_0$"); a1.set_ylabel("count"); a1.set_title(r"$\varepsilon_0$ distribution"); a1.legend(fontsize=9.5, frameon=False)

    a2 = ax[0, 2]
    a2.scatter(dF0s, Qs, c="purple", alpha=0.75, s=48, edgecolors="k", lw=0.5)
    lim = max(dF0s.max(), Qs.max()) * 1.1
    a2.plot([0, lim], [0, lim], "k--", lw=1, label="1:1")
    a2.plot([0, lim], [0, 4 * lim], "gray", ls=":", lw=2, label=r"$Q=4\Delta F_0$")
    a2.set_xlabel(r"$\Delta F_0$ OSTZ (eV)"); a2.set_ylabel(r"$Q$ paper (eV)")
    a2.set_title(r"Activation energy ($Q/\Delta F_0\approx4$)"); a2.legend(fontsize=9.5, frameon=False)

    a3 = ax[1, 0]
    a3.scatter(gBo, gBl, c="crimson", alpha=0.75, s=48, edgecolors="k", lw=0.5)
    lim = max(gBl.max(), gBo.max()) * 1.1
    a3.plot([0, lim], [0, lim], "k--", lw=1, label="1:1")
    for fac, ls in [(2, "--"), (5, ":")]:
        a3.plot([0, lim], [0, fac * lim], "gray", lw=1.8, ls=ls, label=f"{fac}$\\times$")
    a3.set_xlabel(r"$\gamma_B$ OSET Read--Shockley (J/m$^2$)"); a3.set_ylabel(r"$\gamma_B$ paper (J/m$^2$)")
    a3.set_title("GB energy: dislocation-array model"); a3.legend(fontsize=9.5, frameon=False)

    a4 = ax[1, 1]
    a4.scatter(gBe, gBl, c="green", alpha=0.75, s=48, edgecolors="k", lw=0.5)
    lim = max(gBl.max(), gBe.max()) * 1.1
    a4.plot([0, lim], [0, lim], "k--", lw=1, label="1:1")
    a4.set_xlabel(r"measured HAGB energy (J/m$^2$)"); a4.set_ylabel(r"$\gamma_B$ paper fit (J/m$^2$)")
    a4.set_title("Paper's $\\gamma_B$ vs. experiment"); a4.legend(fontsize=10.5, frameon=False)

    a5 = ax[1, 2]
    m20 = ratG < 20
    a5.scatter(Th[m20], ratG[m20], c="darkorange", alpha=0.78, s=48, edgecolors="k", lw=0.5, label="metals/ceramics")
    a5.scatter(Th[~m20], ratG[~m20], c="cornflowerblue", alpha=0.85, s=58, marker="X", edgecolors="k", lw=0.5, label="ice (anomalous)")
    a5.axhline(np.median(ratG[m20]), color="red", ls="--", lw=1.5, label=f"median {np.median(ratG[m20]):.1f}$\\times$")
    a5.axhline(1.0, color="green", ls=":", lw=2, label="perfect agreement")
    a5.set_xlabel(r"homologous temperature $T/T_m$"); a5.set_ylabel(r"$\gamma_B$ ratio (paper / OSET)")
    a5.set_title("Ratio is a constant offset (correct $Gb$ scaling)"); a5.legend(fontsize=9, frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUTDIR, "fig_HP2025_comparison.pdf"))
    plt.close(fig)
    print("Saved fig_HP2025_comparison.pdf")
    return res13


def fig_HP2025_universal(res13):
    g0s = np.array([r["g0"] for r in res13]); e0s = np.array([r["e0"] for r in res13])
    gBl = np.array([r["gB_lit"] for r in res13]); gBo = np.array([r["gB_oset"] for r in res13])
    Td = np.array([r["T"] for r in res13])
    T_rng = np.linspace(300, 2000, 300)
    g0_fit = 0.0827 + 1.3422 / T_rng
    e0_fit = 0.0408 + 0.0117 / T_rng
    gB_fit = 1.075 - 19.962 / T_rng

    fig, axb = plt.subplots(1, 3, figsize=(15, 6.8))
    fig.suptitle("Universal temperature dependences (Harisankar & Padmanabhan 2025, Eqs. 9, 10, 12)",
                 fontsize=14.5, fontweight="bold")

    # Legends are placed BELOW each panel (under the x-label) so they can never
    # overlap the scattered 41-system data inside the axes.
    def _below_legend(ax, ncol):
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=ncol,
                  fontsize=10, frameon=False, handlelength=1.6,
                  columnspacing=1.2, borderaxespad=0.0)

    a = axb[0]
    a.plot(T_rng, g0_fit, "-", color="#c0392b", lw=2.4, label=r"fit $\gamma_0=0.0827+1.342/T$")
    a.scatter(Td, g0s, c="0.45", s=30, zorder=5, edgecolors="white", lw=0.4, label="41 systems")
    a.axhline(0.12, color="#2f6fc0", ls="--", lw=1.6, label="OSTZ crystal 0.12")
    a.axhline(0.10, color="navy", ls=":", lw=1.6, label="OSTZ GB 0.10")
    a.set_xlabel("temperature $T$ (K)"); a.set_ylabel(r"$\gamma_0$")
    a.set_title(r"Shear eigenstrain $\gamma_0(T)$")
    a.set_ylim(0.04, max(g0_fit.max(), g0s.max()) * 1.12)
    _below_legend(a, ncol=2)

    a = axb[1]
    a.plot(T_rng, e0_fit, "-", color="#1f8a70", lw=2.4, label=r"fit $\varepsilon_0=0.0408+0.0117/T$")
    a.scatter(Td, e0s, c="0.45", s=30, zorder=5, edgecolors="white", lw=0.4, label="41 systems")
    a.axhline(0.05, color="#2f6fc0", ls="--", lw=1.6, label="OSTZ 0.05")
    a.set_xlabel("temperature $T$ (K)"); a.set_ylabel(r"$\varepsilon_0$")
    a.set_title(r"Dilatational eigenstrain $\varepsilon_0(T)$")
    a.set_ylim(0.03, max(e0_fit.max(), e0s.max()) * 1.12)
    _below_legend(a, ncol=2)

    a = axb[2]
    a.plot(T_rng, gB_fit, "-", color="#8e44ad", lw=2.4, label=r"fit $\gamma_B=1.075-19.96/T$")
    a.scatter(Td, gBl, c="0.45", s=30, zorder=5, edgecolors="white", lw=0.4, label="paper fit (41 sys.)")
    a.scatter(Td, gBo, c="#c0392b", s=34, marker="^", zorder=6, edgecolors="white", lw=0.4,
              label="OSET Read--Shockley")
    a.set_xlabel("temperature $T$ (K)"); a.set_ylabel(r"$\gamma_B$ (J/m$^2$)")
    a.set_title(r"Boundary energy $\gamma_B(T)$")
    a.set_ylim(0, max(gB_fit.max(), gBl.max(), gBo.max()) * 1.18)
    _below_legend(a, ncol=2)

    fig.tight_layout(rect=[0, 0.10, 1, 0.93])
    fig.savefig(os.path.join(OUTDIR, "fig_HP2025_universal.pdf"))
    plt.close(fig)
    print("Saved fig_HP2025_universal.pdf")


# ===========================================================================
# FIGURE: Unified Rate Equation applied across material classes
# (consolidates already-cited, independently published comparisons into one
#  forest-style figure; no new numbers invented)
# ===========================================================================
def fig_rate_equation():
    # (label, ratio_lo, ratio_hi, n_systems, source)
    entries = [
        (r"Bulk metallic glasses (8 systems)", 1.1, 2.5, 8, "buenz2015"),
        (r"Nanocrystalline Ni", 0.4, 1.6, 1, "pad2018"),
        (r"High-rate Al/composite alloys (8 systems)", 0.5, 2.0, 8, "padbasariya2009"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.4),
                              gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.32})

    ax = axes[0]
    y = np.arange(len(entries))[::-1]
    for yi, (lbl, lo, hi, n, src) in zip(y, entries):
        mid = (lo * hi) ** 0.5
        ax.plot([lo, hi], [yi, yi], color="steelblue", lw=3, solid_capstyle="round")
        ax.plot(mid, yi, "o", color="navy", ms=7, zorder=5)
        ax.text(hi * 1.08, yi, f"$n={n}$", va="center", fontsize=10.5, color="0.25")
    ax.axvline(1.0, color="green", ls=":", lw=1.3, label="perfect agreement")
    ax.set_xscale("log")
    ax.set_xlim(0.2, 6)
    # explicit, non-overlapping log ticks (default LogLocator crowds this range)
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_yticks(y); ax.set_yticklabels([e[0] for e in entries], fontsize=10.5)
    ax.set_ylim(y.min() - 0.6, y.max() + 0.6)
    ax.set_xlabel(r"strain-rate ratio  $\dot\gamma_\mathrm{pred}/\dot\gamma_\mathrm{exp}$")
    ax.set_title("Predicted-to-measured strain-rate ratio", fontsize=13.5)

    ax2 = axes[1]
    win_lo, win_hi = 1e-3, 1.0e-1
    ax2.axhspan(win_lo, win_hi, color="green", alpha=0.15,
                label="typical superplastic\nwindow $10^{-3}$--$10^{-1}\\,\\mathrm{s}^{-1}$")
    ax2.plot([0], [0.16], marker="D", color="crimson", ms=10, zorder=5,
             label=r"OSET prediction: $0.16\,\mathrm{s}^{-1}$")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-4, 1.0)
    ax2.set_xlim(-1, 1)
    ax2.set_xticks([0]); ax2.set_xticklabels(["Al-12Si"])
    ax2.set_ylabel(r"$\dot\gamma$ (s$^{-1}$)")
    ax2.set_title("Worked example (773 K, $d{=}10\\,\\mu$m)", fontsize=13.5)
    ax2.legend(frameon=False, fontsize=9.5, loc="lower left")

    fig.suptitle("Unified Rate Equation applied across material classes", fontweight="bold", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(os.path.join(OUTDIR, "fig_rate_equation.pdf"))
    plt.close(fig)
    print("Saved fig_rate_equation.pdf (new)")


# ===========================================================================
# Run everything
# ===========================================================================
if __name__ == "__main__":
    fig_S1313()
    fig_Taylor()
    fig_partition()
    fig_materials_comparison()
    gb = fig_GB_energy()
    res13 = fig_HP2025_comparison()
    fig_HP2025_universal(res13)
    fig_rate_equation()

    res_mat = materials_results()
    print("\nSFE table (mJ/m^2):")
    print(f"{'Metal':5s} {'struct':5s} {'OSET':8s} {'Lit.lo':8s} {'Lit.hi':8s}")
    for n in metals_db:
        lit = sfe_exp[n]
        lit_s = f"{lit[0]:8.1f} {lit[1]:8.1f}" if lit else "     N/A      N/A"
        print(f"{n:5s} {res_mat[n]['struct']:5s} {res_mat[n]['sfe']:8.1f} {lit_s}")

    print("\nGB energy table (mJ/m^2):")
    print(f"{'Metal':5s} {'struct':5s} {'CTB':8s} {'LA(10)':8s} {'HA':8s}")
    for n, r in gb.items():
        la_s = f"{r['la']*1e3:8.1f}" if r["la"] is not None else "     N/A"
        print(f"{n:5s} {r['struct']:5s} {r['ctb']*1e3:8.1f} {la_s} {r['ha']*1e3:8.1f}")

    print("\nAll figures regenerated in", os.path.abspath(OUTDIR))
