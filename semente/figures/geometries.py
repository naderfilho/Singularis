"""Figuras da Fase 7: geometrias com carga e rotacao."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..geometry.charged import charged_black_bounce, charged_bounce_phase, reissner_nordstrom  # noqa: E402
from ..geometry.kerr import Kerr  # noqa: E402
from ..stability.energy_conditions import static_metric_report  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402

PHASE_COLORS = {"regular BH com horizonte interno": 0, "black-bounce": 1, "buraco de minhoca": 2, "sem horizontes": 2,
                "Reissner-Nordstrom": 3, "singularidade nua": 4}


def _phase_index(s):
    for k, v in PHASE_COLORS.items():
        if s.startswith(k):
            return v
    return 5


def fig_charged():
    Qs = np.linspace(0.0, 1.3, 66)
    ls = np.linspace(0.0, 2.6, 66)
    grid = np.zeros((ls.size, Qs.size))
    for i, l in enumerate(ls):
        for j, Q in enumerate(Qs):
            grid[i, j] = _phase_index(charged_bounce_phase(1.0, Q, l))
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    im = ax.imshow(grid, origin="lower", aspect="auto", extent=[Qs[0], Qs[-1], ls[0], ls[-1]], cmap="viridis", vmin=0, vmax=5)
    Qc = np.linspace(0, 1, 200)
    ax.plot(Qc, 1 - np.sqrt(1 - Qc**2), color="white", lw=1.2, label="l = R_- (horizonte interno some)")
    ax.plot(Qc, 1 + np.sqrt(1 - Qc**2), color=C_WH, lw=1.2, label="l = R_+ (horizonte externo some)")
    ax.axvline(1.0, color=C_SING, ls="--", lw=1, label="Q = M")
    ax.set_xlabel("Q / M")
    ax.set_ylabel("l / M")
    ax.set_title("Black-bounce carregado (Franzin et al. 2021): fases no plano (Q, l)\n0: BH regular c/ horizonte interno, 1: black-bounce, 2: minhoca, 3: RN, 4: singularidade nua")
    ax.legend(fontsize=7, loc="upper right")
    ax = axes[1]
    r = np.linspace(-6, 6, 600)
    for l, col in [(1e-6, "white"), (0.3, C_OUR), (0.6, C_WH), (1.9, C_OTHER)]:
        cb = charged_black_bounce(1.0, 0.5, l)
        ax.semilogy(r, np.abs(cb.kretschmann(r)) + 1e-12, color=col, lw=1.6, label=f"l = {l:g} M (Q = 0.5 M)")
    ax.set_ylim(1e-3, 1e8)
    ax.set_xlabel("r / M")
    ax.set_ylabel("|Kretschmann| [M^-4]")
    ax.set_title("Curvatura: RN diverge em r = 0; com l > 0 e finita na garganta")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    rep = static_metric_report("charged_black_bounce", r, M=1.0, l=0.6, Q=0.5)
    ax.plot(r, rep.rho, color=C_OUR, lw=1.6, label="rho")
    ax.plot(r, rep.p_r, color=C_WH, lw=1.6, label="p_r")
    ax.plot(r, rep.p_t, color=C_OTHER, lw=1.6, label="p_t")
    for a, b in rep.violation_intervals("NEC"):
        ax.axvspan(a, b, color=C_SING, alpha=0.15)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_ylim(-0.02, 0.02)
    ax.set_xlabel("r / M")
    ax.set_title("Fluido efetivo (Q = 0.5, l = 0.6): NEC violada nas faixas vermelhas")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    rec = RunRecord(model="geometry/charged", parameters=dict(M=1.0))
    stamp_figure(fig, rec)
    return _save(fig, "33_geometrias_carregadas.png")


def fig_rotating():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    chis = np.linspace(0, 0.999, 200)
    ks = [Kerr(1.0, c) for c in chis]
    ax.plot(chis, [k.kerr_horizons[1] for k in ks], color="white", lw=2, label="r_+ (horizonte)")
    ax.plot(chis, [k.kerr_horizons[0] for k in ks], color=C_BH, lw=2, label="r_- (Cauchy)")
    ax.plot(chis, [k.ergosurface(np.pi / 2) for k in ks], color=C_WH, lw=2, ls="--", label="ergosuperficie (equador)")
    ax.plot(chis, [k.isco() for k in ks], color=C_OUR, lw=2, label="ISCO prograda")
    ax.plot(chis, [k.isco(False) for k in ks], color=C_OUR, lw=1, ls=":", label="ISCO retrograda")
    ax.set_xlabel("chi = a / M")
    ax.set_ylabel("raio / M")
    ax.set_title("Kerr: horizontes, ergosfera e ISCO (Bardeen-Press-Teukolsky 1972)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    r = np.linspace(1.85, 12, 400)
    for a, col in [(0.3, C_OUR), (0.6, C_WH), (0.9, C_OTHER)]:
        k = Kerr(1.0, a)
        rr = r[r > k.kerr_horizons[1]]
        ax.semilogy(rr, k.frame_dragging(rr), color=col, lw=2, label=f"omega(r), a = {a} M")
        ax.axhline(k.Omega_H, color=col, ls=":", lw=1)
    ax.set_xlabel("r / M (equador)")
    ax.set_ylabel("arrasto do referencial omega [1/M]")
    ax.set_title("Frame dragging (ZAMO) -> Omega_H no horizonte (pontilhado)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[2]
    ls = np.linspace(0.0, 2.2, 45)
    rows = []
    K0 = []
    hs = []
    for l in ls:
        kb = Kerr(1.0, 0.6, l)
        h = kb.horizons
        hs.append(h[-1] if h else np.nan)
        K0.append(float(kb.kretschmann(0.0)) if l > 0 else np.nan)
        rows.append(dict(l=float(l), horizons=list(h), kind=kb.classify()["kind"], K_throat=K0[-1]))
    ax.plot(ls, hs, color="white", lw=2, label="horizonte externo r_h(l)")
    ax.axvline(1.8, color=C_WH, ls="--", lw=1, label="l = r_+(Kerr): horizontes desaparecem")
    ax.set_xlabel("l / M")
    ax.set_ylabel("r_h / M")
    ax2 = ax.twinx()
    ax2.semilogy(ls, K0, color=C_OTHER, lw=2, label="Kretschmann na garganta (r = 0, equador)")
    ax2.set_ylabel("K(r=0) [M^-4]", color=C_OTHER)
    ax.set_title("Black-bounce rotativo (a = 0.6 M, Mazza-Franzin-Liberati 2021):\nhorizonte encolhe com l; curvatura finita na garganta")
    ax.legend(fontsize=8, loc="upper right")
    ax2.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    rec = RunRecord(model="geometry/kerr", parameters=dict(M=1.0, a=0.6))
    stamp_figure(fig, rec)
    path = _save(fig, "34_geometrias_rotativas.png")
    save_json(os.path.join(OUT, "34_black_bounce_rotativo.json"), rows, rec)
    return path


def make_all():
    fig_charged()
    fig_rotating()
