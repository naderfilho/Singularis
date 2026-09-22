"""Figuras da Fase 9: mapas de fase do espaco de parametros e benchmark das PINNs."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..parameter_space.explorer import (CLASSES, CLASS_INDEX, adaptive_refine, classify_bounce, classify_static,  # noqa: E402
                                        monte_carlo, scan_2d)
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402

CMAP = ListedColormap(["#ff4d6d", "#9aa4bd", "#4b5bff", "#ffd166", "#c77dff", "#ff9f1c", "#7f1d1d", "#7bd389", "#40e0d0", "#333333"])


def fig_phase_maps():
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))
    # (a) familia estatica (Q, l), varredura 2D
    s2 = scan_2d(classify_static, ("Q", np.linspace(0, 1.3, 40)), ("l", np.linspace(0, 2.6, 40)), fixed=dict(M=1.0))
    G = s2.grid_labels()
    ax = axes[0]
    ax.imshow(G, origin="lower", aspect="auto", extent=[0, 1.3, 0, 2.6], cmap=CMAP, vmin=0, vmax=9, interpolation="nearest")
    ax.set_xlabel("Q / M")
    ax.set_ylabel("l / M")
    ax.set_title("Familia estatica (M, Q, l): classes de solucao\n(varredura 2D, 1600 pontos)")
    # (b) refinamento adaptativo: pontos avaliados
    ad = adaptive_refine(classify_static, ("Q", np.linspace(0, 1.3, 9)), ("l", np.linspace(0, 2.6, 9)), fixed=dict(M=1.0), levels=2)
    ax = axes[1]
    ax.scatter(ad.points[:, 0], ad.points[:, 1], c=[CLASS_INDEX[l] for l in ad.labels], cmap=CMAP, vmin=0, vmax=9, s=14)
    ax.set_xlabel("Q / M")
    ax.set_ylabel("l / M")
    ax.set_title(f"Refinamento adaptativo nas fronteiras de fase\n({ad.extra['n_evaluations']} avaliacoes vs 1600 da grade)")
    # (c) ricochete cosmologico (rho_*, sigma_0^2): Monte Carlo log-uniforme
    mc = monte_carlo(classify_bounce, dict(rho_star=(1.5, 1e4), sigma0_sq=(1e-9, 1e2)), n=400, seed=11, log_scale=("rho_star", "sigma0_sq"))
    ax = axes[2]
    ax.scatter(mc.points[:, 0], mc.points[:, 1], c=[CLASS_INDEX[l] for l in mc.labels], cmap=CMAP, vmin=0, vmax=9, s=12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rho_* / rho_0")
    ax.set_ylabel("anisotropia inicial sigma_0^2 / rho_0")
    ax.set_title("Ricochete LQC efetivo (poeira): Monte Carlo com semente 11 (400 pontos)\nverde: universo-filho em expansao; vermelho escuro: ricochete instavel (BKL)")
    handles = [plt.Line2D([], [], marker="s", ls="", color=CMAP(i), label=c) for i, c in enumerate(CLASSES)]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    rec = RunRecord(model="parameter_space/phase_maps", parameters=dict(M=1.0), seed=11)
    stamp_figure(fig, rec)
    path = _save(fig, "37_mapas_de_fase.png")
    s2.save(os.path.join(OUT, "37_scan_estatico.json"))
    mc.save(os.path.join(OUT, "37_monte_carlo_ricochete.json"))
    return path


def fig_pinn_benchmark():
    from ..ml.benchmarks import benchmark_bounce_pinn
    seeds = (0, 1, 2, 3)
    b0 = benchmark_bounce_pinn(seeds=seeds, epochs=2000, lbfgs_steps=200, data_weight=0.0)
    b1 = benchmark_bounce_pinn(seeds=seeds, epochs=2000, lbfgs_steps=200, data_weight=1.0)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    x = np.arange(len(seeds))
    ax = axes[0]
    ax.bar(x - 0.2, b0.abs_error, 0.4, color=C_OUR, label="so fisica (L_boundary = 0 por construcao)")
    ax.bar(x + 0.2, b1.abs_error, 0.4, color=C_OTHER, label="hibrido: + 1.0 L_data (16 pontos do solver)")
    ax.axhline(b0.numeric_vs_exact_abs, color="white", ls=":", label=f"solver numerico vs exato: {b0.numeric_vs_exact_abs:.1e}")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"semente {s}" for s in seeds])
    ax.set_ylabel("max |a_PINN - a_exato|")
    ax.set_title("Erro absoluto vs solucao exata, por semente\n(estabilidade do treino)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    ax.bar(x - 0.2, b0.conservation, 0.4, color=C_OUR, label="vinculo de Friedmann (normalizado)")
    ax.bar(x + 0.2, b0.constraint, 0.4, color=C_WH, label="residuo da EDO (normalizado)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"semente {s}" for s in seeds])
    ax.set_title("Erro de conservacao e violacao do vinculo da PINN (so fisica)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax = axes[2]
    ax.bar(x - 0.2, b0.time_s, 0.4, color=C_OUR, label="so fisica")
    ax.bar(x + 0.2, b1.time_s, 0.4, color=C_OTHER, label="hibrido")
    ax.set_xticks(x)
    ax.set_xticklabels([f"semente {s}" for s in seeds])
    ax.set_ylabel("tempo de treino (s, CPU)")
    ax.set_title(f"Custo: PINN ~ {np.mean(b0.time_s):.0f} s vs solver DOP853 < 0.1 s\n(ML como ferramenta de verificacao, nao substituto)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    rec = RunRecord(model="ml/benchmark_bounce_pinn", parameters=b0.parameters, seed=list(seeds)[0], notes=f"seeds={list(seeds)}")
    stamp_figure(fig, rec)
    save_json(os.path.join(OUT, "38_pinn_benchmark.json"), dict(physics_only=b0.summary(), hybrid=b1.summary()), rec)
    return _save(fig, "38_pinn_benchmark.png")


def make_all():
    fig_phase_maps()
    fig_pinn_benchmark()
