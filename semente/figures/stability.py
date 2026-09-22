"""Figuras da Fase 3: condicoes de energia e estabilidade."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..cosmology.friedmann import einstein_cartan_dust, gr_dust  # noqa: E402
from ..stability.bounce import anisotropy_robustness  # noqa: E402
from ..stability.energy_conditions import friedmann_report, static_metric_report  # noqa: E402
from ..stability.perturbations import analyze, master_potential, simpson_visser_metric, tortoise_grid  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402


def fig_energy_conditions_map():
    ls = np.linspace(0.05, 3.0, 60)
    r = np.linspace(-8, 8, 400)
    maps = {c: np.zeros((ls.size, r.size)) for c in ("NEC", "WEC", "SEC", "DEC")}
    for i, l in enumerate(ls):
        rep = static_metric_report("simpson_visser", r, M=1.0, l=l)
        for c in maps:
            maps[c][i] = getattr(rep, c)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.8), sharey=True)
    for ax, c in zip(axes, maps):
        ax.imshow(maps[c], origin="lower", aspect="auto", extent=[r[0], r[-1], ls[0], ls[-1]], cmap="RdYlGn", vmin=0, vmax=1)
        rh = np.sqrt(np.maximum(4 - ls**2, 0))
        ax.plot(rh, ls, color="white", lw=1, ls="--")
        ax.plot(-rh, ls, color="white", lw=1, ls="--")
        ax.axhline(2.0, color=C_WH, ls=":", lw=1)
        ax.set_title(f"{c}: verde = satisfeita, vermelho = violada")
        ax.set_xlabel("r [M]  (r < 0: outro lado da garganta)")
    axes[0].set_ylabel("l / M  (parametro de Simpson-Visser)")
    axes[0].text(-7.8, 2.1, "l = 2M: horizontes desaparecem", color=C_WH, fontsize=8)
    axes[0].text(-7.8, 0.15, "tracejado: horizontes", color="white", fontsize=8)
    fig.suptitle("Condicoes de energia do fluido efetivo de Simpson-Visser (M = 1): NEC violada em todo o exterior para l > 0", fontsize=11)
    fig.tight_layout()
    rec = RunRecord(model="energy_conditions/simpson_visser", parameters=dict(M=1.0, l_min=0.05, l_max=3.0))
    stamp_figure(fig, rec)
    return _save(fig, "24_condicoes_energia_black_bounce.png")


def fig_master_potentials():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    for l, col in [(0.0, "white"), (0.8, C_OUR), (1.6, C_WH), (2.5, C_OTHER)]:
        m = simpson_visser_metric(1.0, l) if l > 0 else simpson_visser_metric(1.0, 1e-9)
        if m.horizons:
            r, rs = tortoise_grid(m, m.horizons[1] * 1.001, 40, 3000)
        else:
            r, rs = tortoise_grid(m, -40, 40, 3000)
        V = master_potential(m, r, spin=0, ell=2)
        rs = rs - rs[int(np.argmax(V))]   # r_* e definida a menos de constante: pico do potencial em r_* = 0
        ax.plot(rs, V, color=col, lw=1.8, label=f"l = {l:g} M")
    ax.set_xlim(-40, 40)
    ax.set_xlabel("coordenada tartaruga r_* [M]")
    ax.set_ylabel("V_0 (escalar, l_mult = 2)")
    ax.set_title("Potencial mestre escalar (M = 1), pico em r_* = 0: sem regioes negativas")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ls = np.linspace(0.1, 3.5, 25)
    rows = []
    for spin, col in [(0, C_OUR), (1, C_WH)]:
        eigs, vmins = [], []
        for l in ls:
            res = analyze(simpson_visser_metric(1.0, l), spin=spin, ell=2, n=(1500, 3000))
            eigs.append(res.lowest_eigenvalue)
            vmins.append(res.potential_min)
            rows.append(dict(l=float(l), spin=spin, classification=res.classification, lowest_eigenvalue=res.lowest_eigenvalue,
                             potential_min=res.potential_min, converged=res.convergence["converged"]))
        ax.plot(ls, eigs, color=col, lw=2, marker="o", ms=3, label=f"menor autovalor, spin {spin}")
        ax.plot(ls, vmins, color=col, lw=1, ls="--", label=f"min V, spin {spin}")
    ax.axhline(0, color=C_SING, lw=1)
    ax.axvline(2.0, color="white", ls=":", lw=1)
    ax.set_xlabel("l / M")
    ax.set_title("Criterio de estado ligado: autovalor negativo <=> modo crescente.\nTodos positivos: campos de teste linearmente estaveis")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    sol = einstein_cartan_dust(1.0, 0.05).solve(1.0, (0, 1.2))
    rep = friedmann_report(sol)
    ax.plot(sol.t, rep.rho + rep.p_r, color=C_WH, lw=2, label="rho_eff + p_eff (NEC)")
    ax.plot(sol.t, rep.rho, color=C_OUR, lw=1.5, label="rho_eff (WEC)")
    ax.plot(sol.t, rep.rho + 3 * rep.p_r, color=C_OTHER, lw=1.5, label="rho_eff + 3 p_eff (SEC)")
    for a, b in rep.violation_intervals("NEC"):
        ax.axvspan(a, b, color=C_SING, alpha=0.2)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel("t (unidades de Planck)")
    ax.set_title("Fluido efetivo total do ricochete de torcao (k = 0):\nNEC e SEC violadas em torno do ricochete (faixa vermelha)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    rec = RunRecord(model="stability/master_potentials", parameters=dict(M=1.0, ell=2))
    stamp_figure(fig, rec)
    path = _save(fig, "25_estabilidade_potenciais.png")
    save_json(os.path.join(OUT, "25_estabilidade_black_bounce.json"), rows, rec)
    return path


def fig_bounce_anisotropy():
    sol = einstein_cartan_dust(1.0, 0.05).solve(1.0, (0, 1.2))
    sig = np.logspace(-8, 0, 60)
    reps = [anisotropy_robustness(sol, rho_star=20.0, sigma0_sq=s) for s in sig]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(sig, [r.Sigma_bounce for r in reps], color=C_WH, lw=2.5)
    ax.axhline(1, color=C_SING, ls="--", label="Sigma_b = 1: anisotropia domina (BKL)")
    ax.axhline(1e-2, color=C_OUR, ls=":", label="Sigma_b = 0.01: ricochete isotropico robusto")
    ax.axvline(reps[0].sigma0_sq_max, color=C_OTHER, ls="-.", label=f"sigma_0^2 max = {reps[0].sigma0_sq_max:.2e}")
    ax.set_xlabel("anisotropia inicial sigma_0^2 (unidades de rho_0)")
    ax.set_ylabel("Sigma_b = sigma^2(a_b) / rho_*")
    ax.set_title(f"Robustez do ricochete isotropico a cisalhamento (a_0/a_b = {reps[0].a_0 / reps[0].a_b:.2f}, rho_* = 20 rho_0)\ncriterio conservador: cisalhamento ~ a^-6 nao regularizado")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    rec = RunRecord(model="stability/bounce_anisotropy", parameters=dict(sigma=0.05, rho_star=20.0))
    stamp_figure(fig, rec)
    return _save(fig, "26_ricochete_anisotropia.png")


def make_all():
    fig_energy_conditions_map()
    fig_master_potentials()
    fig_bounce_anisotropy()
