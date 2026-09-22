"""Figuras da Fase 5: espectros atraves do ricochete (fundo generico) e expansao emergente."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..cosmology.expansion import analyze_expansion, scan_parameters  # noqa: E402
from ..cosmology.friedmann import einstein_cartan_dust, lqc_fluid  # noqa: E402
from ..cosmology.modes import ModeSolver, bounce_background  # noqa: E402
from ..observations.bridge import compare, table  # noqa: E402
from ..observations.constraints import PLANCK  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402


def dust_background(a_b=0.05, a0=2.0e4):
    return bounce_background(einstein_cartan_dust(1.0, a_b**3), a0=a0)


def radiation_background(rho_star=20.0, a0=300.0):
    return bounce_background(lqc_fluid(1.0, 1 / 3, rho_star), a0=a0)


def fig_spectra():
    ks = np.logspace(-0.7, 0.4, 14)
    sd = dust_background()
    sr = radiation_background()
    md, mr = ModeSolver(sd), ModeSolver(sr)
    rd = md.spectrum(ks)
    rr = mr.spectrum(np.logspace(-0.5, 0.3, 10))
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    ax.loglog(rd.k, rd.P_test / rd.P_test[0], color=C_OUR, lw=2.2, marker="o", ms=4, label=f"poeira + torcao: n_s = {rd.n_s:.3f}")
    ax.loglog(rr.k, rr.P_test / rr.P_test[0], color=C_OTHER, lw=2.2, marker="s", ms=4, label=f"radiacao + LQC: n_s = {rr.n_s:.2f}")
    kk = rd.k
    ax.loglog(kk, (kk / kk[0]) ** (PLANCK["n_s"] - 1), color="white", ls=":", lw=1.5, label=f"Planck: n_s = {PLANCK['n_s']}")
    ax.axvline(0.1 * rd.k_bounce, color=C_SING, ls="--", lw=1, label="0.1 k_ricochete (fim do plato)")
    ax.set_xlabel("k (unidades do fundo)")
    ax.set_ylabel("P(k) / P(k_min)")
    ax.set_title("Espectro do campo de teste (= tensorial) atraves do ricochete\nfundo generico (ModeSolver), vacuo de Bunch-Davies na contracao")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[1]
    ok = np.isfinite(rd.r) & (rd.r > 0)
    ax.semilogx(rd.k[ok], rd.r[ok], color=C_WH, lw=2.2, marker="o", ms=4, label="r(k) = 16 eps(t_exit): poeira")
    ax.axhline(PLANCK["r_max_95"], color=C_SING, ls="--", label=f"BICEP/Keck 2021: r < {PLANCK['r_max_95']} (95%)")
    ax.set_ylim(1e-2, 50)
    ax.set_yscale("log")
    ax.set_xlabel("k")
    ax.set_ylabel("razao tensor/escalar r")
    ax.set_title("r = 16 epsilon na saida do horizonte: 24 para contracao de poeira\n(hipotese: zeta transportada pelo ricochete como o campo de teste)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[2]
    for sol, col, lbl in [(sd, C_OUR, "poeira + torcao"), (sr, C_OTHER, "radiacao + LQC")]:
        rep = analyze_expansion(sol)
        m = (sol.t > sol.t_bounce - 0.6) & (sol.t < sol.t_bounce + 0.6)
        ax.plot(sol.t[m] - sol.t_bounce, sol.accel[m] * sol.a[m] ** 0 , color=col, lw=2, label=f"{lbl}: N_acel = {rep.N_accelerated:.2f}")
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel("t - t_ricochete")
    ax.set_ylabel("a'' / a")
    ax.set_title("Expansao emergente: a fase acelerada dura ~1/k_b e rende N < 1 e-fold\n(o codigo determina: NAO ha inflacao)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    rec = RunRecord(model="cosmology/modes", parameters=dict(a_b=0.05, rho_star_rad=20.0))
    stamp_figure(fig, rec)
    path = _save(fig, "29_perturbacoes_fundo_generico.png")
    rows = [compare("poeira+torcao (SEMENTE)", "n_s", rd.n_s, 0.02, rd.assumptions, note="incerteza do modelo: sistematico de extracao validado contra a referencia exata (~2% em n_s)"),
            compare("poeira+torcao (SEMENTE)", "alpha_s", rd.alpha_s, 0.05, rd.assumptions),
            compare("poeira+torcao (SEMENTE)", "r", float(np.nanmedian(rd.r[ok])), 0.0, rd.assumptions),
            compare("matter bounce (literatura)", "f_NL_local", -35 / 8, None, ["Cai, Xue, Brandenberger & Zhang 2009"], literature=True),
            compare("radiacao+LQC", "n_s", rr.n_s, 0.1, rr.assumptions)]
    save_json(os.path.join(OUT, "29_perturbacoes_ponte.json"), [r.as_row() for r in rows], rec)
    with open(os.path.join(OUT, "29_perturbacoes_ponte.md"), "w", encoding="utf-8") as fh:
        fh.write(table(rows))
    return path


def fig_expansion_scan():
    rows = scan_parameters(rho_stars=(2, 5, 20, 100, 1000, 1e4), ws=(0.0, 1 / 3))
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for w, col, lbl in [(0.0, C_OUR, "poeira"), (1 / 3, C_OTHER, "radiacao")]:
        rs = [r for r in rows if r["w"] == w]
        ax[0].semilogx([r["rho_star"] for r in rs], [r["N_accelerated"] for r in rs], color=col, lw=2, marker="o", label=lbl)
        ax[1].loglog([r["rho_star"] for r in rs], [r["duration"] for r in rs], color=col, lw=2, marker="o", label=lbl)
    ax[0].axhline(60, color=C_SING, ls="--", label="N = 60 (inflacao)")
    ax[0].set_ylim(0, 2)
    ax[0].set_xlabel("rho_* (unidades de rho_0)")
    ax[0].set_ylabel("e-folds da fase acelerada pos-ricochete")
    ax[0].set_title("N acelerado e independente de rho_* e fica < 1:\nnenhuma inflacao emergente (LQC efetiva, k = 0)")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("rho_*")
    ax[1].set_ylabel("duracao da fase acelerada (t)")
    ax[1].set_title("duracao ~ rho_*^(-1/2): a fase acelerada e a escala do ricochete")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    rec = RunRecord(model="cosmology/expansion.scan", parameters=dict(mechanism="lqc"))
    stamp_figure(fig, rec)
    path = _save(fig, "30_expansao_emergente.png")
    save_json(os.path.join(OUT, "30_expansao_emergente.json"), rows, rec)
    return path


def make_all():
    fig_spectra()
    fig_expansion_scan()
