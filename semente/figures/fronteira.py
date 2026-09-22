"""Figuras dos experimentos de fronteira (semente.fronteira)."""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..bounce import M_SUN  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402
from ..fronteira import (K_B_ENTROPIES, BudgetedSelection, TorsionCollapse, genealogy_depth,  # noqa: E402
                        l_from_black_hole_mass, parent_mass_lower_bound, smolin_neutron_star_test)
from ..selection import CosmicSelection  # noqa: E402


def fig_juncao():
    tc = TorsionCollapse(M=1.0, R0=8.0, Rb_over_R0=0.05)
    sol = tc.solve()
    js = tc.junction_schwarzschild(sol)
    jl = tc.junction_black_bounce(sol, tc.R_b)
    ls, frac = tc.scan_l(sol)
    tau, R = sol["tau"], sol["R"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    ax = axes[0, 0]
    ax.plot(tau, R, color="#ff9f1c", lw=2.2, label="raio areal da superficie R(tau)")
    ax.axhline(2, color="white", ls="--", lw=0.8, label="r = 2M")
    ax.axhline(tc.R_b, color="#40e0d0", ls=":", lw=1.2, label=f"R_b = {tc.R_b} M (ricochete)")
    inv = ~js["valid"]
    ax.fill_between(tau, 0, R.max(), where=inv, color=C_SING, alpha=0.25, label="juncao com Schwarzschild IMPOSSIVEL")
    ax.set_xlabel("tempo proprio da superficie (tau / M)")
    ax.set_ylabel("R / M")
    ax.set_yscale("log")
    ax.set_title("Colapso com torcao: ricochete dentro do horizonte.\nFaixa vermelha: nenhuma camada cola o interior a Schwarzschild(M).")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    ax = axes[0, 1]
    ax.plot(tau, js["Es2"], color=C_SING, lw=2, label="1 + Rdot^2 - 2M/R  (exterior Schwarzschild)")
    ax.plot(tau, jl["F"], color="#40e0d0", lw=2, label="f + rdot^2  (exterior black-bounce, l = R_b)")
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlim(sol["tau_b"] - 1.2, sol["tau_b"] + 1.2)
    ax.set_ylim(-3, 3)
    ax.set_xlabel("tau / M")
    ax.set_title("Termo sob a raiz da curvatura extrinseca exterior:\nnegativo = impossivel; com l = R_b fica >= 0 sempre.")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1, 0]
    ax.plot(ls / tc.R_b, 100 * frac, color=C_WH, lw=2.5, marker="o", ms=3)
    ax.set_xlabel("l / R_b  (parametro do black-bounce exterior)")
    ax.set_ylabel("% do tempo em que a juncao e impossivel")
    ax.set_title("Varredura em l: so l = R_b funciona (garganta = raio do ricochete).\nO parametro 'livre' de Simpson-Visser fica fixado pelo interior.")
    ax.grid(alpha=0.3)
    ax = axes[1, 1]
    ms = np.logspace(0, 10.5, 200) * M_SUN
    ax.loglog(ms / M_SUN, l_from_black_hole_mass(ms), color=C_OTHER, lw=2.5)
    for name, m in [("10 M_sol", 10), ("Sgr A*", 4e6), ("M87*", 6.5e9)]:
        lm = l_from_black_hole_mass(m * M_SUN)
        ax.scatter([m], [lm], color="white", zorder=5)
        ax.annotate(f"{name}: l = {lm:.1e} m", (m, lm), textcoords="offset points", xytext=(8, -12), fontsize=8)
    ax.set_xlabel("massa do buraco negro (M_sol)")
    ax.set_ylabel("l previsto (m)")
    ax.set_title("Previsao:  l^3 = 3M / (4 pi rho_b),  rho_b = 4 m_n^2 / pi\n(garganta sub-nanometrica; camada de fronteira pesa ~ R_b c^2/G)")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    return _save(fig, "16_fronteira_juncao.png")


def fig_entropia():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    ax = axes[0]
    names = list(K_B_ENTROPIES)
    S = np.array([K_B_ENTROPIES[k] for k in names])
    Mmin = np.array([parent_mass_lower_bound(s)[1] for s in S])
    ax.barh(names, np.log10(Mmin), color=[C_WH, C_WH, C_BH, C_BH])
    for i, (m, s) in enumerate(zip(Mmin, S)):
        ax.text(np.log10(m) + 0.2, i, f"M_pai >= {m:.1e} M_sol   (S = {s:.1e} k_B)", va="center", fontsize=8.5)
    ax.set_xlim(0, 20)
    ax.set_xlabel("log10( massa minima do buraco negro pai / M_sol )")
    ax.set_title("Se o horizonte do pai limita a entropia do filho (4 pi M^2 >= S_filho):\npiso para a massa do buraco negro em que o NOSSO universo nasceu")
    ax.grid(alpha=0.3, axis="x")
    ax = axes[1]
    N = np.logspace(1, 25, 200)
    for M0, col in [(10, C_OUR), (1e6, C_WH), (1e10, C_OTHER)]:
        ax.semilogx(N, genealogy_depth(M0 * M_SUN, N), color=col, lw=2.2, label=f"buraco negro raiz de {M0:g} M_sol")
    ax.axvline(1e20, color="white", ls=":", lw=1)
    ax.text(1.3e20, 30, "N_bh do nosso universo ~ 1e20", color="white", fontsize=8, rotation=90, va="center")
    ax.set_xlabel("N = buracos negros por universo (por geracao)")
    ax.set_ylabel("k_max = geracoes ate a massa de Planck")
    ax.set_title("Arvore de universos sob o limite holografico: profundidade finita\nk_max = ln(M_0 / m_Pl) / ln(sqrt N)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    return _save(fig, "17_fronteira_entropia.png")


def fig_selecao_orcamento():
    free = CosmicSelection(n_pop=3000).run(120)
    hist = BudgetedSelection().run(120)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    ax.semilogy(free["fecundity"], color=C_WH, lw=2.2, label="Smolin classico: <N_bh> (sem limite)")
    ax.semilogy([h["N_fit_mean"] for h in hist], color=C_SING, lw=2.2, label="com orcamento holografico: <filhos que cabem>")
    ax.set_xlabel("geracao")
    ax.set_title("A selecao perde a alavanca: com orcamento cada universo\ncabe ~1 filho, e sem diferenca de fecundidade nao ha otimizacao")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.semilogy([h["m_median"] for h in hist], color=C_BH, lw=2.2)
    ax.set_xlabel("geracao")
    ax.set_ylabel("massa mediana dos buracos negros (m_Pl)")
    ax.set_title("Deriva neutra: a massa tipica dos buracos negros cai\nlentamente geracao a geracao (orcamento herdado)")
    ax.grid(alpha=0.3)
    ax = axes[2]
    ax.plot([h["alive"] for h in hist], color=C_OUR, lw=2.2)
    ax.set_xlabel("geracao")
    ax.set_ylabel("universos com pelo menos 1 filho")
    ax.set_ylim(0, 3200)
    ax.set_title("Populacao sobrevive, mas sem amplificacao:\nreproducao ~1:1, o oposto da selecao de Smolin")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, "18_fronteira_selecao_orcamento.png")


def fig_estrelas_neutrons():
    rows = smolin_neutron_star_test()
    fig, ax = plt.subplots(figsize=(9, 4.8))
    names = list(rows)
    m = [rows[n]["massa"] for n in names]
    e = [rows[n]["erro"] for n in names]
    ax.errorbar(m, range(len(names)), xerr=e, fmt="o", color=C_WH, capsize=4, lw=2)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.axvline(1.6, color=C_SING, ls="--", lw=1.5, label="Smolin 1992: M_max ~ 1.6 M_sol")
    ax.axvline(2.0, color=C_OTHER, ls=":", lw=1.5, label="Smolin (revisao): ~2.0 M_sol")
    ax.set_xlabel("massa medida (M_sol)")
    ax.set_title("Teste de dados: a previsao original da selecao natural cosmologica\nesta excluida por >4 sigma em todos os pulsares; a revisada esta no limite")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, axis="x")
    return _save(fig, "19_fronteira_estrelas_neutrons.png")


def make_all():
    from ..fronteira import run_all
    out = run_all(verbose=False)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "fronteira_resultados.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, default=float)
    fig_juncao()
    fig_entropia()
    fig_selecao_orcamento()
    fig_estrelas_neutrons()
    return out


if __name__ == "__main__":
    make_all()
