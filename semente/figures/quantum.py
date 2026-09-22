"""Figuras da Fase 4: GR vs Einstein-Cartan vs LQC com condicoes iniciais equivalentes."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..core.units import M_NEUTRON, M_PLANCK, M_SUN  # noqa: E402
from ..cosmology.seed_universe import SeedUniverse  # noqa: E402
from ..quantum.comparison import compare_models  # noqa: E402
from ..quantum.einstein_cartan import bounce_density_planck  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402

COL = {"GR": C_SING, "Einstein-Cartan": C_WH, "LQC efetiva": C_OTHER}


def fig_gr_ec_lqc():
    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    rows = {}
    for row, (w, lbl) in enumerate([(0.0, "poeira (w = 0)"), (1 / 3, "radiacao (w = 1/3)")]):
        comps = compare_models(rho0=1.0, w=w, rho_star=20.0, t_span=(0, 1.5))
        rows[lbl] = [dict(name=c.name, a_min=c.a_min, rho_max=c.rho_max, ricci_max=c.ricci_max, t_bounce=c.t_bounce,
                          Sigma_b=c.anisotropy_Sigma_b, efolds=c.efolds_after_bounce) for c in comps]
        ax = axes[row, 0]
        for c in comps:
            ls = "--" if c.name == "LQC efetiva" else "-"
            ax.plot(c.solution.t, c.solution.a, color=COL[c.name], lw=2.2, ls=ls, label=c.name)
        ax.set_xlabel("t (unidades de Planck)")
        ax.set_ylabel("a(t)")
        ax.set_title(f"{lbl}: fator de escala, mesmas condicoes iniciais (a0 = 1, contracao)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax = axes[row, 1]
        for c in comps:
            ls = "--" if c.name == "LQC efetiva" else "-"
            ax.semilogy(c.solution.t, c.solution.rho, color=COL[c.name], lw=2, ls=ls, label=f"rho, {c.name}")
        ax.axhline(20.0, color="white", ls=":", lw=1, label="rho_* = 20")
        ax.set_ylim(1e-1, 1e4)
        ax.set_xlabel("t")
        ax.set_title("densidade de materia: GR diverge; EC/LQC identicas para poeira,\ndiferentes para radiacao (EC escala com n^2 ~ a^-6)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
        ax = axes[row, 2]
        for c in comps:
            ls = "--" if c.name == "LQC efetiva" else "-"
            ax.semilogy(c.solution.t, np.abs(c.solution.ricci_scalar) + 1e-12, color=COL[c.name], lw=2, ls=ls, label=c.name)
        ax.set_ylim(1e-1, 1e6)
        ax.set_xlabel("t")
        ax.set_title("|escalar de Ricci|: curvatura maxima finita com correcao")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("GR vs Einstein-Cartan (sigma = rho0^2/rho_*) vs LQC efetiva (rho_c = rho_*), rho0 = 1, rho_* = 20, k = 0", fontsize=11)
    fig.tight_layout()
    rec = RunRecord(model="quantum/comparison", parameters=dict(rho0=1.0, rho_star=20.0, k=0))
    stamp_figure(fig, rec)
    path = _save(fig, "27_gr_vs_ec_vs_lqc.png")
    save_json(os.path.join(OUT, "27_gr_vs_ec_vs_lqc.json"), rows, rec)
    return path


def fig_progenitor_to_child():
    masses = np.logspace(0, 10, 60) * M_SUN
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    ax = axes[0]
    for mech, col, lbl in [("einstein_cartan", C_WH, "Einstein-Cartan (neutrons)"), ("lqc", C_OTHER, "LQC efetiva (rho_c = 0.41 rho_Pl)")]:
        Rb = [SeedUniverse(m, mechanism=mech).table()["R_areal no ricochete (m)"] for m in masses]
        ax.loglog(masses / M_SUN, Rb, color=col, lw=2.2, label=lbl)
    ax.axhline(1.616e-35, color="white", ls=":", lw=1)
    ax.text(1.5, 2.5e-35, "comprimento de Planck", color="white", fontsize=8)
    ax.set_xlabel("massa do buraco negro pai (M_sol)")
    ax.set_ylabel("raio areal do ricochete R_b (m)")
    ax.set_title("Progenitor -> filho: raio do ricochete R_b = (3M / 4 pi rho_*)^(1/3)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[1]
    mf = np.logspace(-3, 2, 100) * M_NEUTRON
    ax.loglog(mf / M_NEUTRON, [bounce_density_planck(m) for m in mf], color=C_WH, lw=2.2, label="rho_b(EC) = 4 m^2 / pi")
    ax.axhline(0.41, color=C_OTHER, ls="--", label="rho_c (LQC)")
    ax.set_xlabel("massa do fermion / massa do neutron")
    ax.set_ylabel("densidade de ricochete (rho_Pl)")
    ax.set_title("Densidade critica: EC depende da massa do fermion; LQC e universal")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    rec = RunRecord(model="quantum/progenitor_child", parameters=dict(m_fermion=M_NEUTRON))
    stamp_figure(fig, rec)
    return _save(fig, "28_progenitor_para_filho.png")


def make_all():
    fig_gr_ec_lqc()
    fig_progenitor_to_child()
