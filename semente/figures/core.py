"""
Todas as figuras do projeto.  `python -m semente.figures` gera tudo em ./output.

Cada figura e calculada do zero pelos modulos fisicos; nenhuma e desenhada "a mao".
"""
from __future__ import annotations

import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

from ..bounce import (LQCCosmology, SeedUniverse, TorsionCosmology, classical_collapse_reference,  # noqa: E402
                     observable_universe_as_black_hole, torsion_dust_exact, M_SUN)
from ..collapse import OppenheimerSnyder  # noqa: E402
from ..geodesics import GeodesicSolver, deflection_angle  # noqa: E402
from ..geometry import BlackBounce, Schwarzschild, penrose_from_kruskal  # noqa: E402
from ..interior import (kantowski_sachs_black_bounce, kantowski_sachs_schwarzschild,  # noqa: E402
                       radial_infall_black_bounce)
from ..selection import CosmicSelection, fecundity_landscape  # noqa: E402

OUT = "output"
DARK = "#0b0d12"
plt.rcParams.update({
    "figure.facecolor": DARK, "axes.facecolor": DARK, "savefig.facecolor": DARK,
    "axes.edgecolor": "#8a93a6", "axes.labelcolor": "#e6e9f0", "xtick.color": "#c7cbd6",
    "ytick.color": "#c7cbd6", "text.color": "#e6e9f0", "grid.color": "#2a2f3d", "font.size": 10,
    "axes.titleweight": "bold", "legend.facecolor": "#151925", "legend.edgecolor": "#3a4152",
})
C_BH, C_WH, C_OUR, C_OTHER, C_SING = "#4b5bff", "#ffd166", "#7bd389", "#c77dff", "#ff4d6d"


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("  ->", path)
    return path


# ----------------------------------------------------------------------------
def fig_kruskal():
    """Extensao maximal: as quatro regioes, o buraco branco que 'aparece sem querer',
    a estrela de Oppenheimer-Snyder (que cobre so parte do diagrama) e um observador em queda."""
    from scipy.integrate import cumulative_trapezoid
    s = Schwarzschild(1.0)
    fig, ax = plt.subplots(figsize=(9, 9))
    lim = 3.2
    X = np.linspace(-lim, lim, 600)
    ax.fill_between(X, np.abs(X), np.sqrt(1 + X**2), color=C_BH, alpha=0.18, label="II: buraco NEGRO (r<2M)")
    ax.fill_between(X, -np.sqrt(1 + X**2), -np.abs(X), color=C_WH, alpha=0.18, label="III: buraco BRANCO")
    ax.fill_between(X[X > 0], -X[X > 0], X[X > 0], color=C_OUR, alpha=0.10, label="I: nosso universo")
    ax.fill_between(X[X < 0], X[X < 0], -X[X < 0], color=C_OTHER, alpha=0.10, label="IV: outro universo")
    ax.plot(X, np.sqrt(1 + X**2), color=C_SING, lw=2.5, label="singularidade r = 0")
    ax.plot(X, -np.sqrt(1 + X**2), color=C_SING, lw=2.5)
    ax.plot([-lim, lim], [-lim, lim], color="white", lw=1.2, alpha=0.8)
    ax.plot([-lim, lim], [lim, -lim], color="white", lw=1.2, alpha=0.8, label="horizontes r = 2M")
    for r in [2.5, 3.0, 4.0, 5.0]:
        T, Xc = s.const_r_curve(r, tmax=14, n=600)
        ax.plot(Xc, T, color="#9aa4bd", lw=0.7, alpha=0.7)
        ax.plot(-Xc, T, color="#9aa4bd", lw=0.7, alpha=0.7)
    for r in [0.5, 1.0, 1.5]:
        T, Xc = s.const_r_curve(r, tmax=14, n=600)
        ax.plot(Xc, T, color="#9aa4bd", lw=0.7, alpha=0.7, ls="--")
        ax.plot(Xc, -T, color="#9aa4bd", lw=0.7, alpha=0.7, ls="--")
    for tt in [-4, -2, 0, 2, 4]:
        ax.plot([-lim, lim], np.tanh(tt / 4) * np.array([-lim, lim]), color="#5e6a86", lw=0.5, alpha=0.7)

    # --- estrela de Oppenheimer-Snyder: escolhemos a origem de t para que a superficie cruze
    # o horizonte em V = e^{v/4M} = 1.25 (translacao temporal = simetria do exterior)
    R0 = 4.5
    os_ = OppenheimerSnyder(1.0, R0)
    sk0 = os_.surface_kruskal(n=4000, eta_max=np.pi * 0.9995)
    v_h = np.interp(2.0, sk0["R"][::-1], sk0["v"][::-1])
    shift = 4 * np.log(1.25) - v_h
    sk = os_.surface_kruskal(n=4000, eta_max=np.pi * 0.9995, v_shift=shift)
    t_static = np.linspace(-60, 0, 1500) + shift
    Ts, Xs = s.kruskal_from_tr(t_static, np.full_like(t_static, R0))
    Xw = np.concatenate([Xs, sk["X"]])
    Tw = np.concatenate([Ts, sk["T"]])
    ax.plot(Xw, Tw, color="#ff9f1c", lw=2.5, label="superficie da estrela (Oppenheimer-Snyder)")
    ax.fill_betweenx(Tw, -lim, Xw, color="#ff9f1c", alpha=0.14, hatch="//", edgecolor="#ff9f1c", linewidth=0)

    # --- observador em queda livre a partir de r = 4M (geodesica radial via v de Eddington-Finkelstein)
    Rq = 4.0
    eta = np.linspace(0, np.pi * 0.999, 3000)
    Rg = Rq / 2 * (1 + np.cos(eta))
    tau = np.sqrt(Rq**3 / 8) * (eta + np.sin(eta))
    E = np.sqrt(1 - 2 / Rq)
    f = 1 - 2 / Rg
    v = cumulative_trapezoid(1 / (E + np.sqrt(np.maximum(E**2 - f, 0))), tau, initial=0)
    v_h = np.interp(2.0, Rg[::-1], v[::-1])
    v = v - v_h + 4 * np.log(1.9)  # cruza o horizonte em V = 1.9
    V = np.exp(v / 4)
    U = -(Rg / 2 - 1) * np.exp(Rg / 2) / V
    ax.plot((V - U) / 2, (V + U) / 2, color="#40e0d0", lw=2.2, label="observador em queda livre (r = 4M -> 0)")
    for c in [0.4, 0.9]:
        ax.plot([c - lim, c + lim], [-lim, lim], color=C_WH, lw=0.8, ls=":")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("X (Kruskal)")
    ax.set_ylabel("T (Kruskal)")
    ax.set_title("Extensao maximal de Kruskal-Szekeres: o buraco branco vem junto com o buraco negro")
    ax.text(1.9, 0.0, "I", fontsize=18, ha="center", color=C_OUR)
    ax.text(-1.9, 0.0, "IV", fontsize=18, ha="center", color=C_OTHER)
    ax.text(0.0, 1.6, "II", fontsize=18, ha="center", color=C_BH)
    ax.text(0.0, -1.7, "III", fontsize=18, ha="center", color=C_WH)
    nota = chr(10).join([
        "Hachurado = interior da estrela: la a metrica nao e Schwarzschild,",
        "e FRW (um universo fechado em contracao). Num colapso real tudo",
        "a esquerda da superficie - inclusive III (buraco branco) e IV - e",
        "substituido pela estrela. O buraco branco so existe na solucao",
        "de vacuo ETERNA... ou se a singularidade for substituida por um ricochete."])
    ax.text(-3.1, 3.05, nota,
            fontsize=8.3, va="top", color="#ffd9a8")
    ax.legend(loc="lower right", fontsize=8)
    return _save(fig, "01_kruskal.png")


def fig_penrose():
    """Diagramas de Penrose: (a) Schwarzschild eterno; (b) a 'escada' do black-bounce."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    ax = axes[0]
    lim = 40.0
    X = np.linspace(-lim, lim, 4000)
    for sign, col in [(1, C_SING), (-1, C_SING)]:
        x, y = penrose_from_kruskal(sign * np.sqrt(1 + X**2), X)
        ax.plot(x, y, color=col, lw=2.5)
    for tt in np.linspace(-lim, lim, 9):
        x, y = penrose_from_kruskal(np.full_like(X, tt), X)
    s = Schwarzschild(1.0)
    for r in [2.2, 3, 5, 10]:
        T, Xc = s.const_r_curve(r, tmax=60, n=3000)
        for sx in (1, -1):
            x, y = penrose_from_kruskal(T, sx * Xc)
            ax.plot(x, y, color="#9aa4bd", lw=0.6, alpha=0.8)
    for r in [0.3, 1.0, 1.7]:
        T, Xc = s.const_r_curve(r, tmax=60, n=3000)
        for st in (1, -1):
            x, y = penrose_from_kruskal(st * T, Xc)
            ax.plot(x, y, color="#9aa4bd", lw=0.6, alpha=0.8, ls="--")
    q = np.pi / 4
    # horizontes (do ponto de bifurcacao ate a singularidade) e infinitos nulos scri+-
    ax.plot([-q, q], [-q, q], color="white", lw=1.2)
    ax.plot([-q, q], [q, -q], color="white", lw=1.2)
    for sx in (1, -1):
        ax.plot([sx * q, sx * 2 * q, sx * q], [q, 0, -q], color="white", lw=1.2)
    ax.add_patch(Polygon([[0, 0], [q, q], [2 * q, 0], [q, -q]], closed=True, color=C_OUR, alpha=0.10))
    ax.add_patch(Polygon([[0, 0], [-q, q], [-2 * q, 0], [-q, -q]], closed=True, color=C_OTHER, alpha=0.10))
    x, y = penrose_from_kruskal(np.sqrt(1 + X**2), X)
    ax.fill_between(x, np.abs(x) - 0 * x, y, where=y > np.abs(x) - 1e-9, color=C_BH, alpha=0.18)
    x, y = penrose_from_kruskal(-np.sqrt(1 + X**2), X)
    ax.fill_between(x, y, -np.abs(x), color=C_WH, alpha=0.18)
    for lbl, xy, col in [("i+", (q, q + 0.07), "w"), ("i-", (q, -q - 0.07), "w"), ("i0", (2 * q + 0.09, 0), "w"),
                         ("scri+", (1.5 * q + 0.12, 0.5 * q + 0.08), "w"), ("scri-", (1.5 * q + 0.12, -0.5 * q - 0.08), "w"),
                         ("BURACO NEGRO", (0, 0.45), C_BH), ("BURACO BRANCO", (0, -0.5), C_WH),
                         ("nosso universo", (q, 0), C_OUR), ("outro universo", (-q, 0), C_OTHER),
                         ("singularidade futura r=0", (0, 0.86), C_SING), ("singularidade passada r=0", (0, -0.9), C_SING)]:
        ax.text(*xy, lbl, ha="center", va="center", fontsize=9, color=col, fontweight="bold")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.0, 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(a) Penrose: Schwarzschild eterno")

    # (b) escada do black-bounce (Simpson-Visser, 0 < l < 2M)
    ax = axes[1]
    ax.set_title("(b) Penrose: black-bounce (0 < l < 2M) - universos empilhados")
    for k in range(-1, 2):
        y0 = 2 * k
        # dois losangos exteriores (universo k, lado esquerdo e direito) e uma faixa interna
        for sx in (1, -1):
            ax.add_patch(Polygon([[0, y0 - 1], [sx * 1, y0], [0, y0 + 1]], closed=True,
                                 facecolor=(C_OUR if sx > 0 else C_OTHER), alpha=0.12, edgecolor="white", lw=1))
        ax.add_patch(Polygon([[-1, y0], [0, y0 + 1], [1, y0]], closed=True, facecolor=C_BH, alpha=0.20, edgecolor="none"))
        ax.add_patch(Polygon([[-1, y0], [0, y0 - 1], [1, y0]], closed=True, facecolor=C_WH, alpha=0.20, edgecolor="none"))
        ax.plot([-1, 1], [y0 + 1, y0 + 1], color="#40e0d0", lw=2.5)
        ax.text(0, y0 + 1.08, "garganta r = 0 (superficie ESPACIAL regular: o ricochete)", ha="center", fontsize=8, color="#40e0d0")
        ax.text(0, y0 + 0.55, "BN", ha="center", color=C_BH, fontsize=11, fontweight="bold")
        ax.text(0, y0 - 0.6, "BB", ha="center", color=C_WH, fontsize=11, fontweight="bold")
        ax.text(0.62, y0, f"universo {k + 1}", ha="center", fontsize=8, color=C_OUR)
        ax.text(-0.62, y0, f"universo {k + 1}'", ha="center", fontsize=8, color=C_OTHER)
    ax.annotate("", xy=(0.15, 2.3), xytext=(0.15, -1.6), arrowprops=dict(arrowstyle="->", color="#ff9f1c", lw=2))
    ax.text(0.22, 0.2, "materia que cai no BN\ndo universo 1 emerge\ndo BB no universo 2:\no 'Big Bang' dele", fontsize=8, color="#ff9f1c")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-2.6, 3.2)
    ax.set_aspect("equal")
    ax.axis("off")
    return _save(fig, "02_penrose.png")


def fig_flamm():
    """Ponte de Einstein-Rosen (Flamm) e garganta lisa do black-bounce (l > 2M)."""
    fig = plt.figure(figsize=(14, 6))
    s = Schwarzschild(1.0)
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    r = np.linspace(2.0, 12.0, 120)
    ph = np.linspace(0, 2 * np.pi, 120)
    R, PH = np.meshgrid(r, ph)
    Z = s.flamm_paraboloid(R)
    for sign, cmap in [(1, "Blues_r"), (-1, "Purples_r")]:
        ax.plot_surface(R * np.cos(PH), R * np.sin(PH), sign * Z, cmap=cmap, alpha=0.85, linewidth=0, antialiased=True)
    ax.set_title("Paraboloide de Flamm: a fatia t=const de Schwarzschild\n(garganta r=2M = ponte de Einstein-Rosen para a regiao IV)")
    ax.set_facecolor(DARK)
    ax.set_axis_off()
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    bb = BlackBounce(1.0, 2.5)
    rr = np.linspace(-12, 12, 240)
    Rr = bb.R(rr)
    f = bb.f(rr)
    dzdr = np.sqrt(np.maximum(1 / f - (rr / Rr) ** 2, 0))
    from scipy.integrate import cumulative_trapezoid
    z = cumulative_trapezoid(dzdr, rr, initial=0)
    z -= z[len(z) // 2]
    RR, PH = np.meshgrid(Rr, ph)
    ZZ = np.broadcast_to(z, RR.shape)
    col = np.where(np.broadcast_to(rr, RR.shape) > 0, 0.85, 0.15)
    ax.plot_surface(RR * np.cos(PH), RR * np.sin(PH), ZZ, facecolors=plt.cm.cool(col), alpha=0.9, linewidth=0)
    ax.set_title("Black-bounce com l = 2.5M (> 2M): buraco de minhoca\nliso ligando dois universos pela garganta R = l")
    ax.set_facecolor(DARK)
    ax.set_axis_off()
    return _save(fig, "03_flamm_wormhole.png")


def fig_interior():
    ks = kantowski_sachs_schwarzschild(1.0)
    kb = kantowski_sachs_black_bounce(1.0, 0.5)
    inf = radial_infall_black_bounce(1.0, 0.5, r0=8.0, tau_max=60.0)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    ax = axes[0, 0]
    ax.plot(ks["tau"], ks["a_perp"], color=C_BH, lw=2, label="a_perp = r (esferas)")
    ax.plot(ks["tau"], ks["a_par"], color=C_SING, lw=2, label="a_par = sqrt(2M/r - 1) (direcao t)")
    ax.set_yscale("log")
    ax.set_xlabel("tempo proprio desde o horizonte  (tau / M)")
    ax.set_title("Dentro de Schwarzschild: cosmologia de Kantowski-Sachs\n(uma direcao estica ao infinito, as esferas colapsam em tau = pi M)")
    ax.axvline(np.pi, color="white", ls=":", lw=1)
    ax.text(np.pi, 1e-2, " tau = pi M", color="white")
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[0, 1]
    ax.plot(kb["tau"], kb["a_perp"], color=C_BH, lw=2, label="a_perp = R = sqrt(r^2 + l^2)")
    ax.plot(kb["tau"], kb["a_par"], color=C_SING, lw=2, label="a_par = sqrt(2M/R - 1)")
    ax.axhline(0.5, color=C_WH, ls="--", lw=1, label="a_perp,min = l")
    ax.axvline(kb["tau_bounce"], color="#40e0d0", ls=":", lw=1.2, label="ricochete (r = 0)")
    ax.set_xlabel("tempo proprio desde o horizonte do BN  (tau / M)")
    ax.set_title("Dentro de um black-bounce (l = 0.5M): o mesmo universo\nCONTRAI, RICOCHETEIA em R = l e RE-EXPANDE no lado do buraco branco")
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[1, 0]
    ax.plot(ks["tau"], ks["kretschmann"], color=C_SING, lw=2, label="Schwarzschild: K = 48M^2/r^6 -> infinito")
    ax.plot(kb["tau"], kb["kretschmann"], color="#40e0d0", lw=2, label="black-bounce l=0.5M: K finito (max. na garganta)")
    ax.set_yscale("log")
    ax.set_ylim(1e-1, 1e9)
    ax.set_xlabel("tau / M")
    ax.set_ylabel("escalar de Kretschmann K")
    ax.set_title("A singularidade (curvatura infinita) e substituida por curvatura finita")
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[1, 1]
    ax.plot(inf["tau"], inf["r"], color="#ff9f1c", lw=2)
    for h in BlackBounce(1.0, 0.5).horizons:
        ax.axhline(h, color="white", ls="--", lw=0.8)
    ax.axhline(0, color="#40e0d0", ls=":", lw=1.2)
    ax.text(1, 0.15, "garganta r = 0", color="#40e0d0")
    ax.text(1, 2.1, "horizonte do BN (nosso lado)", color="white", fontsize=8)
    ax.text(1, -2.5, "horizonte do BB (outro universo)", color="white", fontsize=8)
    ax.set_xlabel("tempo proprio do observador em queda  (tau / M)")
    ax.set_ylabel("coordenada r  (r < 0 = outro universo)")
    ax.set_title("Queda livre desde r = 8M atravessando a garganta:\nemerge do BURACO BRANCO no outro universo (e volta a cair no BN de la)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, "04_interior_universo.png")


def fig_collapse():
    os_ = OppenheimerSnyder(1.0, 8.0)
    eta = np.linspace(0, np.pi * 0.9999, 2000)
    tau = os_.tau(eta)
    R = os_.R_surface(eta)
    eh = os_.event_horizon_interior()
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    ax = axes[0]
    ax.plot(tau, R, color="#ff9f1c", lw=2.5, label="superficie da estrela R(tau)")
    ax.axhline(2, color="white", ls="--", lw=1, label="r = 2M")
    ax.plot(eh["tau"], eh["R"], color=C_BH, lw=2.5, label="horizonte de EVENTOS (nasce no centro!)")
    ax.fill_between(tau, 0, R, color="#ff9f1c", alpha=0.12)
    ax.axvline(os_.tau_s, color=C_SING, ls=":", lw=1.5)
    ax.text(os_.tau_s, 4, " singularidade", color=C_SING, rotation=90, va="bottom")
    ax.set_xlabel("tempo proprio da superficie (tau / M)")
    ax.set_ylabel("raio areal R / M")
    ax.set_title("Colapso de Oppenheimer-Snyder (R0 = 8M)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    a = os_.a(eta)
    ax.plot(tau, a / a.max(), color=C_OUR, lw=2.5, label="a(tau) dentro da estrela: FRW fechado em CONTRACAO")
    ax.plot(2 * os_.tau_s - tau, a / a.max(), color=C_WH, lw=2.5, ls="--", label="o mesmo a(tau) com o tempo invertido: um BIG BANG")
    ax.set_xlabel("tau / M")
    ax.set_ylabel("fator de escala a / a_max")
    ax.set_title("O interior da estrela E um universo de Friedmann\n(a mesma equacao do Big Bang, ao contrario)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    ax.semilogy(tau, os_.density(eta), color=C_SING, lw=2)
    ax.set_xlabel("tau / M")
    ax.set_ylabel("densidade rho (unidades G=c=1, M=1)")
    ax.set_title("Densidade -> infinito na singularidade classica.\nSe algo (torcao, LQC) a limita, a contracao vira expansao.")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, "05_colapso_oppenheimer_snyder.png")


def fig_bounce():
    tc = TorsionCosmology(rho_m0=1.0, sigma=0.05)
    st = tc.solve(a0=1.0, t_span=(0, 1.2))
    ex, tb = torsion_dust_exact(st["t"], 1.0, 0.05)
    lq = LQCCosmology(rho0=0.2, rho_c=0.41)
    sl = lq.solve(a0=1.0, t_span=(0, 1.2))
    cl = classical_collapse_reference(rho_m0=1.0)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    ax.plot(cl["t"], cl["a"], color=C_SING, lw=2, label="GR pura: a -> 0 (singularidade)")
    ax.plot(st["t"], st["a"], color=C_WH, lw=2.5, label="Einstein-Cartan (torcao/spin): ricochete")
    ax.plot(st["t"], ex, color="white", lw=1, ls=":", label="solucao exata  a^3 = sigma/rho + 6 pi rho (t - t_b)^2")
    ax.plot(sl["t"], sl["a"], color=C_OTHER, lw=2.5, label="LQC (rho_c = 0.41 rho_Pl): ricochete")
    ax.set_xlabel("tempo t (unidades de Planck)")
    ax.set_ylabel("fator de escala a(t)")
    ax.set_title("Contracao -> ricochete -> expansao (um Big Bang)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.plot(st["t"], st["H"], color=C_WH, lw=2, label="H(t) torcao")
    ax.plot(sl["t"], sl["H"], color=C_OTHER, lw=2, label="H(t) LQC")
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel("H = a'/a")
    ax.set_title("H cruza zero suavemente: o 'Big Bang' do filho")
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[2]
    ax.semilogy(st["t"], st["rho"], color=C_SING, lw=2, label="rho (materia)")
    ax.semilogy(st["t"], np.abs(st["rho_eff"]), color=C_WH, lw=2, label="|rho_eff| = |rho - sigma a^-6|  (torcao)")
    ax.axhline(0.41, color=C_OTHER, ls="--", label="rho_c (LQC)")
    ax.set_xlabel("t")
    ax.set_title("A densidade satura: torcao = repulsao ~ n^2")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, "06_ricochete.png")


def fig_seed_table():
    masses = [3, 10, 30, 1e3, 4e6, 6.5e9]
    rows = []
    for m in masses:
        rows.append(SeedUniverse(m * M_SUN).table())
    keys = list(rows[0].keys())
    fig, ax = plt.subplots(figsize=(15, 4.6))
    ax.axis("off")
    cell = [[f"{r[k]:.3g}" for r in rows] for k in keys]
    tbl = ax.table(cellText=cell, rowLabels=keys, colLabels=[f"{m:g} M_sol" for m in masses], loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.4)
    for (_, _), c in tbl.get_celld().items():
        c.set_facecolor("#151925")
        c.set_edgecolor("#3a4152")
        c.get_text().set_color("#e6e9f0")
    ax.set_title("Modelo SEMENTE: o universo-filho em funcao da massa do buraco negro pai\n(interior de Oppenheimer-Snyder feito de neutrons + ricochete de torcao de Einstein-Cartan)",
                 fontsize=11, fontweight="bold")
    path = _save(fig, "07_tabela_universo_filho.png")
    with open(os.path.join(OUT, "07_tabela_universo_filho.json"), "w", encoding="utf-8") as fh:
        json.dump({f"{m:g} M_sol": {k: float(v) for k, v in r.items()} for m, r in zip(masses, rows)}, fh, indent=2)
    return path


def fig_selection():
    res = CosmicSelection(n_pop=4000).run(60)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    ax.plot(res["fecundity"], color=C_WH, lw=2.5)
    ax.set_xlabel("geracao")
    ax.set_ylabel("<N_bh> por universo")
    ax.set_title("Selecao natural cosmologica (Smolin):\na fecundidade media cresce geracao a geracao")
    ax.grid(alpha=0.3)
    ax = axes[1]
    for i, lbl, col in [(0, "escala de massa dos fermions", C_OUR), (1, "constante cosmologica", C_OTHER), (2, "amplitude Q", C_BH)]:
        ax.plot(res["mean"][:, i], color=col, lw=2, label=lbl)
        ax.fill_between(range(len(res["mean"])), res["mean"][:, i] - res["std"][:, i], res["mean"][:, i] + res["std"][:, i], color=col, alpha=0.15)
    ax.set_xlabel("geracao")
    ax.set_title("Os parametros da populacao convergem\npara o maximo de producao de buracos negros")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    g = np.linspace(-3, 3, 200)
    lam = np.linspace(-3, 3, 200)
    G, LM = np.meshgrid(g, lam)
    Q = np.full(G.size, 2.0)
    Z = fecundity_landscape(np.stack([G.ravel(), LM.ravel(), Q], 1)).reshape(G.shape)
    ax.contourf(G, LM, Z, levels=30, cmap="magma")
    for gen, th in res["snapshots"]:
        ax.scatter(th[:300, 0], th[:300, 1], s=4, alpha=0.6, label=f"geracao {gen}")
    ax.set_xlabel("escala de massa dos fermions")
    ax.set_ylabel("constante cosmologica")
    ax.set_title("Paisagem de fecundidade N_bh(theta) e a populacao migrando")
    ax.legend(fontsize=7, markerscale=3)
    fig.tight_layout()
    return _save(fig, "08_selecao_cosmologica.png")


def fig_deflection():
    b = np.linspace(5.25, 40, 160)
    alpha = np.array([deflection_angle(bb) for bb in b])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(b, np.degrees(alpha), color=C_WH, lw=2.5, label="integracao de u'' + u = 3Mu^2")
    ax.plot(b, np.degrees(4 / b), color=C_OUR, ls="--", label="Einstein (campo fraco): 4M/b")
    ax.axvline(3 * np.sqrt(3), color=C_SING, ls=":", label="b = 3 sqrt(3) M (esfera de fotons): deflexao infinita")
    ax.set_yscale("log")
    ax.set_xlabel("parametro de impacto b / M")
    ax.set_ylabel("angulo de deflexao (graus)")
    ax.set_title("Lente gravitacional: a mesma equacao que gera as imagens")
    ax.legend()
    ax.grid(alpha=0.3)
    return _save(fig, "09_deflexao.png")


def fig_geodesic_gallery():
    """Orbitas de fotons e particulas em Schwarzschild e no black-bounce, calculadas com o integrador
    geral (Christoffel via sympy)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    ax = axes[0]
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(2 * np.cos(th), 2 * np.sin(th), color="black")
    ax.plot(3 * np.cos(th), 3 * np.sin(th), color=C_SING, ls=":", lw=1, label="esfera de fotons r=3M")
    ax.plot(6 * np.cos(th), 6 * np.sin(th), color=C_OUR, ls=":", lw=1, label="ISCO r=6M")
    from ..geometry.geodesics import orbit_rhs_null
    for b, col in [(5.3, "#ff9f1c"), (5.6, "#ffd166"), (7.0, "#7bd389"), (10.0, "#40e0d0"), (15.0, "#c77dff")]:
        u, up, phi = 0.0, 1 / b, 0.0
        us, ps = [], []
        dphi = 2e-3
        while phi < 6 * np.pi and 0 <= u < 0.5:
            k1u, k1p = up, orbit_rhs_null(u, 1)
            k2u, k2p = up + .5 * dphi * k1p, orbit_rhs_null(u + .5 * dphi * k1u, 1)
            k3u, k3p = up + .5 * dphi * k2p, orbit_rhs_null(u + .5 * dphi * k2u, 1)
            k4u, k4p = up + dphi * k3p, orbit_rhs_null(u + dphi * k3u, 1)
            u += dphi * (k1u + 2 * k2u + 2 * k3u + k4u) / 6
            up += dphi * (k1p + 2 * k2p + 2 * k3p + k4p) / 6
            phi += dphi
            if u > 1 / 60:
                us.append(u)
                ps.append(phi)
        r = 1 / np.array(us)
        ax.plot(r * np.cos(ps), r * np.sin(ps), color=col, lw=1.3, label=f"foton, b = {b}M")
    ax.set_aspect("equal")
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    ax.set_title("Raios de luz em Schwarzschild (b -> 3 sqrt(3) M: orbitam varias vezes)")
    ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    g = GeodesicSolver("black_bounce", 1.0, 0.6)
    bb = BlackBounce(1.0, 0.6)
    for rh in bb.horizons:
        ax.plot(abs(rh) * np.cos(th), abs(rh) * np.sin(th), color="white", ls="--", lw=0.8)
    for L, col in [(3.9, "#ff9f1c"), (4.3, "#40e0d0"), (6.0, "#c77dff")]:
        # particula massiva vinda de r=15 com momento angular L (usa E fixada pela normalizacao com dr<0)
        r0 = 15.0
        x0 = [0.0, r0, np.pi / 2, 0.0]
        v0 = g.initial_velocity(x0, dr=-0.35, dphi=L / r0**2)
        sol = g.integrate(x0, v0, 400.0, stop_r=bb.horizons[1] * 1.001 if bb.horizons else None, max_step=0.1)
        r, ph = sol.y[1], sol.y[3]
        ax.plot(r * np.cos(ph), r * np.sin(ph), color=col, lw=1.3, label=f"particula, L = {L}M (ate o horizonte)")
    ax.fill(bb.horizons[1] * np.cos(th), bb.horizons[1] * np.sin(th), color="#1a1030")
    ax.set_aspect("equal")
    ax.set_xlim(-16, 16)
    ax.set_ylim(-16, 16)
    ax.set_title("Black-bounce l = 0.6M: geodesicas do integrador geral\n(simbolos de Christoffel gerados por sympy)")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    return _save(fig, "10_geodesicas.png")


def fig_pinn():
    try:
        from ..ml.pinn import BouncePINN, PhotonOrbitPINN, reference_photon_orbit
    except ImportError:
        print("  (torch ausente: pulando figura da PINN)")
        return None
    p = BouncePINN(rho_m0=1.0, sigma=0.05, t_max=1.2)
    hist = p.train(epochs=3000, verbose=False, lbfgs_steps=400)
    ref = TorsionCosmology(1.0, 0.0, 0.05).solve(a0=1.0, t_span=(0, 1.2))
    pred = p.predict(ref["t"])
    q = PhotonOrbitPINN(b=5.5)
    hq = q.train(epochs=1500)
    phi, u = reference_photon_orbit(5.5)
    uq = q.predict(phi)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    ax.plot(ref["t"], ref["a"], color="white", lw=4, alpha=0.35, label="integrador numerico (DOP853)")
    ax.plot(ref["t"], pred, color=C_WH, lw=1.8, label="PINN (so viu a equacao)")
    ax.set_xlabel("t")
    ax.set_ylabel("a(t)")
    ax.set_title(f"PINN aprende o ricochete de torcao\nerro maximo |a_PINN - a_num| = {np.abs(pred - ref['a']).max():.1e}")
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.semilogy(hist, color=C_OTHER, lw=1)
    ax.set_xlabel("iteracao (Adam + L-BFGS)")
    ax.set_ylabel("perda = <residuo^2> + <vinculo^2>")
    ax.set_title("Treino: o residuo da equacao de Friedmann cai ~10 ordens")
    ax.grid(alpha=0.3)
    ax = axes[2]
    r = 1 / np.maximum(u, 1e-3)
    rq = 1 / np.maximum(uq, 1e-3)
    ax.plot(r * np.cos(phi), r * np.sin(phi), color="white", lw=4, alpha=0.35, label="RK45")
    ax.plot(rq * np.cos(phi), rq * np.sin(phi), color=C_WH, lw=1.5, label="PINN")
    th = np.linspace(0, 2 * np.pi, 100)
    ax.fill(2 * np.cos(th), 2 * np.sin(th), color="black")
    ax.set_aspect("equal")
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)
    ax.set_title(f"PINN da orbita de um foton (b = 5.5M): u'' + u = 3Mu^2\nerro maximo = {np.abs(uq - u).max():.1e}")
    ax.legend()
    fig.tight_layout()
    return _save(fig, "11_pinn.png")


def renders(width=1280, height=720):
    from ..raytracing.raytracer import Camera, Scene, render, save_png
    os.makedirs(OUT, exist_ok=True)
    paths = []
    for l, name, title in [(0.0, "12_render_schwarzschild.png", "Schwarzschild"),
                           (1.0, "13_render_black_bounce.png", "black-bounce l=1M"),
                           (2.6, "14_render_wormhole.png", "buraco de minhoca l=2.6M")]:
        t0 = time.time()
        img = render(Scene(l=l), Camera(width=width, height=height), verbose=False)
        p = os.path.join(OUT, name)
        save_png(img, p)
        print(f"  -> {p}  ({title}, {time.time() - t0:.0f}s)")
        paths.append(p)
    # vista quase de cima, sem disco, para ver so a lente + o outro universo
    img = render(Scene(l=1.2, disk=False, exposure=0.6), Camera(width=width, height=height, inclination_deg=25), verbose=False)
    p = os.path.join(OUT, "15_render_lente_outro_universo.png")
    save_png(img, p)
    paths.append(p)
    print("  ->", p)
    return paths


def make_all(with_renders=True, render_size=(1280, 720)):
    t0 = time.time()
    print("Gerando figuras...")
    fig_kruskal()
    fig_penrose()
    fig_flamm()
    fig_interior()
    fig_collapse()
    fig_bounce()
    fig_seed_table()
    fig_selection()
    fig_deflection()
    fig_geodesic_gallery()
    fig_pinn()
    if with_renders:
        renders(*render_size)
    from .fronteira import make_all as fronteira
    fronteira()
    from .nascimento import make_all as nascimento
    nascimento()
    from .collapse import make_all as collapse
    collapse()
    from .stability import make_all as stability
    stability()
    from .quantum import make_all as quantum
    quantum()
    from .perturbations import make_all as perturbations
    perturbations()
    from .gw import make_all as gw
    gw()
    from .geometries import make_all as geometries
    geometries()
    info = observable_universe_as_black_hole()
    with open(os.path.join(OUT, "universo_observavel_como_buraco_negro.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2)
    print(f"Pronto em {time.time() - t0:.0f}s. Raio de Schwarzschild do universo observavel / raio: {info['razao_rs_R']:.2f}")


if __name__ == "__main__":
    import sys
    make_all(with_renders="--sem-render" not in sys.argv)
