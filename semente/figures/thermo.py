"""Figuras da Fase 8: termodinamica de horizontes e curvas de Page (exploratorio)."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..collapse.dynamic import compare_classical_vs_bounce  # noqa: E402
from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..geometry.static import StaticSphericalMetric  # noqa: E402
from ..information.page_curve import page_curves, timescales  # noqa: E402
from ..thermodynamics.horizon import apparent_horizon_entropy, event_horizon_entropy, static_thermo  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402


def fig_thermo():
    ls = np.linspace(0.0, 1.98, 60)
    T, S, res = [], [], []
    for l in ls:
        m = StaticSphericalMetric("1 - 2*M/sqrt(r**2 + l**2)", "sqrt(r**2 + l**2)", ("M", "l"), dict(M=1.0, l=max(l, 1e-9)), name="SV")
        th = static_thermo(m)
        T.append(th.temperature)
        S.append(th.entropy)
        res.append(th.first_law_residual)
    c, b = compare_classical_vs_bounce(M=1.0, R0=8.0, Rb_over_R0=0.05, n_shells=40)
    sc, sb = apparent_horizon_entropy(c), apparent_horizon_entropy(b)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    ax.plot(ls, np.array(T) * 8 * np.pi, color=C_WH, lw=2, label="T / T_Schwarzschild")
    ax.plot(ls, np.array(S) / (4 * np.pi), color=C_OUR, lw=2, label="S / S_Schwarzschild (= 1: R(r_h) = 2M)")
    ax.plot(ls, res, color=C_SING, lw=2, ls="--", label="residuo da 1a lei: 1 - T dS/dM")
    ax.plot(ls, 1 - np.sqrt(1 - ls**2 / 4), color="white", lw=1, ls=":", label="1 - sqrt(1 - l^2/4M^2) (analitico)")
    ax.set_xlabel("l / M")
    ax.set_title("Simpson-Visser: S = A/4 nao muda com l, T cai,\ne dM = T dS falha (termo de trabalho em l ou S != A/4)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    se = event_horizon_entropy(c)
    ax.plot(se["tau"], se["S"], color=C_OUR, lw=2.5, label="classico: S do horizonte de EVENTOS cresce (teorema da area)")
    ax.plot(sc["tau"], sc["S"], color=C_SING, lw=1.5, ls="--", label="classico: S do horizonte APARENTE (propaga para dentro)")
    ax.plot(sb["tau"], sb["S"], color=C_WH, lw=2, label="ricochete: S_AH transiente (volta a zero)")
    ax.axhline(4 * np.pi, color="white", ls=":", lw=1, label="4 pi M^2")
    ax.set_xlabel("tau [M]")
    ax.set_ylabel("S = pi R^2  [M^2]")
    ax.set_title("Horizontes no colapso dinamico: o teorema da area vale para o de eventos;\nno ricochete a regiao presa e transiente (NEC violada)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    Ms = np.logspace(0, 10, 50)
    ts = [timescales(m) for m in Ms]
    ax.loglog(Ms, [t["t_bounce_s"] for t in ts], color=C_WH, lw=2, label="ricochete interior ~ pi M")
    ax.loglog(Ms, [t["t_page_s"] for t in ts], color=C_OTHER, lw=2, label="tempo de Page ~ 0.54 t_ev")
    ax.loglog(Ms, [t["t_evap_s"] for t in ts], color=C_BH, lw=2, label="evaporacao t_ev")
    ax.axhline(4.35e17, color="white", ls=":", lw=1)
    ax.text(1.5, 8e17, "idade do universo", color="white", fontsize=8)
    ax.set_xlabel("massa (M_sol)")
    ax.set_ylabel("tempo (s)")
    ax.set_title("Escalas de tempo: o destino do interior e decidido\n~1e78 vezes antes do tempo de Page")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    rec = RunRecord(model="thermodynamics/horizon", parameters=dict(M=1.0))
    stamp_figure(fig, rec)
    return _save(fig, "35_termodinamica.png")


def fig_page():
    pc = page_curves(M0_planck=1e3, beta=1.5, t_star_over_tev=0.3)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    x = pc.t_over_tev
    ax.plot(x, pc.S_BH / pc.S0, color="white", lw=2, label="S_BH(t) = A/4")
    ax.plot(x, pc.S_hawking / pc.S0, color=C_SING, lw=2, ls="--", label="Hawking (radiacao termica): cresce sempre")
    ax.plot(x, pc.S_page / pc.S0, color=C_OUR, lw=2.5, label="curva de Page (evaporacao unitaria)")
    ax.plot(x, pc.S_baby / pc.S0, color=C_OTHER, lw=2.5, ls="-.", label="universo-filho: congela em t_* (parceiros inacessiveis)")
    ax.axvline(pc.t_page_over_tev, color=C_OUR, ls=":", lw=1)
    ax.text(pc.t_page_over_tev + 0.01, 0.9, f"t_Page = {pc.t_page_over_tev:.2f} t_ev", color=C_OUR, fontsize=8)
    ax.axvline(pc.t_star_over_tev, color=C_OTHER, ls=":", lw=1)
    ax.text(pc.t_star_over_tev + 0.01, 0.75, "t_* (ilustrativo)", color=C_OTHER, fontsize=8)
    ax.set_xlabel("t / t_evaporacao")
    ax.set_ylabel("entropia / S_BH(0)")
    ax.set_title("EXPLORATORIO: entropia de emaranhamento da radiacao em tres cenarios\n(beta = 1.5; nenhum e derivado de gravidade quantica)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    rec = RunRecord(model="information/page_curve", parameters=dict(M0_planck=1e3, beta=1.5, t_star=0.3),
                    assumptions=list(pc.assumptions), notes="EXPLORATORIO")
    stamp_figure(fig, rec, extra="EXPLORATORIO")
    save_json(os.path.join(OUT, "36_page_curves.json"), dict(t_page_over_tev=pc.t_page_over_tev, timescales_10Msun=timescales(10.0)), rec)
    return _save(fig, "36_curvas_de_page.png")


def make_all():
    fig_thermo()
    fig_page()
