"""Figuras do modulo nascimento (previsoes do modelo vs o nosso universo)."""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..bounce import M_SUN  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402
from ..nascimento import PLANCK, BounceSpectrum, curvature_constraints, verdict  # noqa: E402


def fig_espectro(out=None):
    out = verdict(verbose=False) if out is None else out
    bs = BounceSpectrum(a_b=0.05)
    ks, P = out["ks"], out["P"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    m = np.abs(bs.eta) < 3
    ax.plot(bs.eta[m], bs.app_over_a[m], color=C_WH, lw=2, label="a''/a  (poeira + torcao)")
    eta = bs.eta[m]
    with np.errstate(divide="ignore"):
        ax.plot(eta, 2 / eta**2, color="white", ls=":", lw=1, label="2/eta^2 (poeira pura, singular)")
    ax.set_yscale("log")
    ax.set_ylim(1, 1e3)
    ax.set_xlabel("tempo conforme eta (ricochete em 0)")
    ax.set_title("O 'potencial' das perturbacoes atraves do ricochete:\nfinito na garganta, 2/eta^2 longe dela")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.loglog(ks, P / P[0], color=C_OUR, lw=2.5, marker="o", ms=4, label="P_phi(k) numerico (Bunch-Davies -> ricochete -> expansao)")
    ax.axvline(bs.k_bounce, color=C_SING, ls="--", label=f"k_ricochete = {bs.k_bounce:.1f}")
    kk = np.logspace(np.log10(ks.min()), np.log10(0.3 * bs.k_bounce), 50)
    ax.loglog(kk, (kk / ks.min()) ** (PLANCK["n_s"] - 1), color=C_OTHER, ls=":", lw=2, label=f"Planck: n_s = {PLANCK['n_s']}")
    ax.set_xlabel("k (unidades do fundo)")
    ax.set_ylabel("P(k) / P(k_min)")
    ax.set_title(f"Espectro primordial gerado pela contracao de poeira:\nn_s = {out['n_s_numerico']:.3f} no plato (invariante de escala SEM inflacao)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[2]
    for k, col in [(0.3, C_OUR), (3.0, C_WH), (30.0, C_SING)]:
        # solucao completa de um modo
        from scipy.integrate import solve_ivp
        eta0 = bs.eta_min
        v0 = np.exp(-1j * k * eta0) / np.sqrt(2 * k)
        y0 = [v0.real, v0.imag, (-1j * k * v0).real, (-1j * k * v0).imag]
        sol = solve_ivp(lambda e, y: [y[2], y[3], -(k * k - bs._app(e)) * y[0], -(k * k - bs._app(e)) * y[1]],
                        (eta0, 6.0), y0, method="DOP853", rtol=1e-8, atol=1e-11, max_step=0.02, dense_output=True)
        e = np.linspace(-6, 6, 3000)
        y = sol.sol(e)
        phi = np.abs(y[0] + 1j * y[1]) / bs._a(e)
        ax.semilogy(e, phi * k**1.5, color=col, lw=1.6, label=f"k = {k:g}")
    ax.set_xlabel("eta")
    ax.set_ylabel("k^{3/2} |phi_k|")
    ax.set_title("Modos atraves do ricochete: os de k pequeno congelam\n(amplitude constante = o que vira estrutura no filho)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, "20_nascimento_espectro.png")


def fig_curvatura_veredito(out=None):
    out = verdict(verbose=False) if out is None else out
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw=dict(width_ratios=[1, 1.5]))
    ax = axes[0]
    Oks = -np.logspace(-4, -1, 80)
    Mmin = []
    for ok in Oks:
        cc = curvature_constraints(ok, 0.0, n_sigma=0.0)
        Mmin.append(cc["M_parent_min_Msol"])
    ax.loglog(-Oks, Mmin, color=C_OTHER, lw=2.5)
    ax.axvspan(1e-4, 0.0031, color=C_OUR, alpha=0.15, label="permitido por Planck+BAO (2 sigma)")
    ax.axvline(0.011, color=C_WH, ls="--", label="Planck sozinho (central, fechado)")
    ax.axhline(1e23, color="white", ls=":", lw=1)
    ax.text(1.2e-4, 1.3e23, "massa do universo observavel", color="white", fontsize=8)
    ax.set_xlabel("-Omega_k  (universo fechado)")
    ax.set_ylabel("massa minima do buraco negro pai (M_sol)")
    ax.set_title("Curvatura observada -> tamanho da esfera-3 -> massa do pai\n(so conservacao de massa: consistente, nao evidencia)")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3, which="both")
    ax = axes[1]
    ax.axis("off")
    cells = [[r["item"], r["previsto"], r["observado"]] for r in out["rows"]]
    status = []
    for r in out["rows"]:
        v = r["veredito"].split(" ")[0].strip(":;,")
        status.append(v)
    tbl = ax.table(cellText=[[c[0], c[1][:48], c[2][:48], s] for c, s in zip(cells, status)],
                   colLabels=["previsao do modelo", "valor previsto", "observado", "veredito"],
                   loc="center", cellLoc="left", colWidths=[0.34, 0.26, 0.26, 0.14])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 2.2)
    colors = {"COMPATIVEL": C_OUR, "CONSISTENTE": C_OUR, "FALSIFICADO": C_SING}
    for (row, col), c in tbl.get_celld().items():
        c.set_facecolor("#151925")
        c.set_edgecolor("#3a4152")
        c.get_text().set_color("#e6e9f0")
        if row > 0 and col == 3:
            c.get_text().set_color(colors.get(status[row - 1], C_WH))
            c.get_text().set_fontweight("bold")
    ax.set_title("Veredito: o que o modelo preve para um universo nascido do ricochete\nvs o que medimos no nosso", fontsize=11, fontweight="bold")
    fig.tight_layout()
    return _save(fig, "21_nascimento_veredito.png")


def make_all():
    out = verdict(verbose=False)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "nascimento_veredito.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(n_s=out["n_s_numerico"], k_bounce=out["k_bounce"], curvatura=out["curvatura"],
                       curvatura_planck_only=out["curvatura_planck_only"], r_min=out["r_min"], rows=out["rows"]),
                  fh, indent=2, default=float, ensure_ascii=False)
    fig_espectro(out)
    fig_curvatura_veredito(out)
    return out


if __name__ == "__main__":
    make_all()
