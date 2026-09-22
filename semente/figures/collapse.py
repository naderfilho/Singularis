"""Figuras do colapso dinamico (Fase 2): classico vs ricochete e estrutura causal."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..collapse.dynamic import CollapseResult, compare_classical_vs_bounce  # noqa: E402
from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402


def _record(res: CollapseResult) -> RunRecord:
    return RunRecord(model=f"SphericalCollapse/{res.model}", parameters=res.parameters, solver=res.solver_metadata["solver"],
                     assumptions=["poeira sem pressao", "simetria esferica", "dinamica de camada com densidade media (modelo efetivo)"])


def fig_collapse_comparison(M=1.0, R0=8.0, Rb_over_R0=0.05, n_shells=40):
    c, b = compare_classical_vs_bounce(M=M, R0=R0, Rb_over_R0=Rb_over_R0, n_shells=n_shells)
    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))
    for col, (res, title, colr) in enumerate([(c, "RG classica (Oppenheimer-Snyder)", C_SING), (b, f"ricochete efetivo, rho_* = {b.parameters['rho_star']:.3g}", C_WH)]):
        ax = axes[0, col]
        for j in range(0, res.chi.size, max(1, res.chi.size // 10)):
            ax.plot(res.tau, res.R[:, j], color=colr, lw=0.9, alpha=0.75)
        ax.plot(res.tau, res.surface_R, color="white", lw=2, label="superficie")
        ax.plot(res.tau, res.apparent_horizon_R, color=C_BH, lw=2.2, label="horizonte aparente (externo)")
        ax.axhline(2 * M, color="#9aa4bd", ls="--", lw=0.8, label="R = 2M")
        ax.set_yscale("log")
        ax.set_xlabel("tempo proprio comovel tau [M]")
        ax.set_ylabel("raio areal R [M]")
        ax.set_title(f"{title}\nR(tau) das camadas, horizonte aparente")
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(alpha=0.3)
    ax = axes[0, 2]
    ax.semilogy(c.tau, np.nanmax(c.rho, axis=1), color=C_SING, lw=2, label="rho_max classico")
    ax.semilogy(b.tau, np.nanmax(b.rho, axis=1), color=C_WH, lw=2, label="rho_max ricochete")
    ax.axhline(b.parameters["rho_star"], color=C_OTHER, ls=":", label="rho_*")
    ax.set_xlabel("tau [M]")
    ax.set_ylabel("densidade [M^-2]")
    ax.set_title("Densidade maxima: diverge (classico) vs satura em rho_* (ricochete)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1, 0]
    ax.plot(c.tau, c.m_MS[:, -1], color=C_SING, lw=2, label="classico")
    ax.plot(b.tau, b.m_MS[:, -1], color=C_WH, lw=2, label="ricochete")
    ax.set_xlabel("tau [M]")
    ax.set_ylabel("massa de Misner-Sharp na superficie [M]")
    ax.set_title("Massa quasi-local: conservada em RG; no ricochete cai para ~R_b/R0 do valor\n(o deficit e a energia da correcao, cf. juncao de Israel)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1, 1]
    if c.kretschmann_homog is not None:
        ax.semilogy(c.tau, np.nanmax(c.kretschmann_homog, axis=1), color=C_SING, lw=2, label="K classico")
        ax.semilogy(b.tau, np.nanmax(b.kretschmann_homog, axis=1), color=C_WH, lw=2, label="K ricochete")
    ax.set_xlabel("tau [M]")
    ax.set_ylabel("escalar de Kretschmann [M^-4]")
    ax.set_title("Curvatura: singular vs limitada")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1, 2]
    ax.plot(b.tau, b.trapped.sum(axis=1) / b.chi.size, color=C_BH, lw=2, label="fracao de camadas presas (ricochete)")
    ax.plot(c.tau, c.trapped.sum(axis=1) / c.chi.size, color=C_SING, lw=2, ls="--", label="fracao presa (classico)")
    ax.axvline(b.summary()["bounce_tau_surface"], color=C_WH, ls=":", label="ricochete")
    ax.set_xlabel("tau [M]")
    ax.set_title("Regiao presa: permanente (classico) vs transiente (ricochete)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.suptitle(f"Colapso dinamico de poeira, M = {M}, R0 = {R0} M, {n_shells} camadas  |  unidades G = c = 1", fontsize=11)
    fig.tight_layout()
    stamp_figure(fig, _record(b))
    path = _save(fig, "22_colapso_dinamico.png")
    save_json(os.path.join(OUT, "22_colapso_dinamico.json"), dict(classico=c.summary(), ricochete=b.summary()), _record(b))
    return path


def fig_causal_structure(M=1.0, R0=8.0, Rb_over_R0=0.05, n_shells=40):
    c, b = compare_classical_vs_bounce(M=M, R0=R0, Rb_over_R0=Rb_over_R0, n_shells=n_shells)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    for ax, res, title in [(axes[0], c, "RG classica"), (axes[1], b, "ricochete efetivo")]:
        T, X = np.meshgrid(res.tau, res.chi, indexing="ij")
        ax.contourf(X, T, res.trapped.astype(float), levels=[0.5, 1.5], colors=[C_BH], alpha=0.35)
        ax.contour(X, T, res.trapped.astype(float), levels=[0.5], colors=[C_BH], linewidths=1.5)
        # raios nulos de saida (dchi/dtau = sqrt(1+2E)/R') e de entrada (sinal oposto) a partir de varios pontos
        prof_E = res.R_prime  # placeholder for shape
        E = -res.m_MS[0] / res.R[0]  # E aproximada da grade inicial (perfil em repouso)
        for t0 in np.linspace(res.tau[0], res.tau[-1] * 0.9, 9):
            for sign, colr in [(1, C_WH), (-1, C_OTHER)]:
                cc = res.chi[0] if sign > 0 else res.chi[-1]
                i0 = int(np.searchsorted(res.tau, t0))
                xs, ts = [cc], [res.tau[i0]]
                for i in range(i0, res.tau.size - 1):
                    Rp = np.interp(cc, res.chi, res.R_prime[i])
                    if Rp <= 0:
                        break
                    cc = cc + sign * np.sqrt(max(1 + 2 * np.interp(cc, res.chi, E), 0)) / Rp * (res.tau[i + 1] - res.tau[i])
                    if not (res.chi[0] <= cc <= res.chi[-1]):
                        break
                    xs.append(cc)
                    ts.append(res.tau[i + 1])
                ax.plot(xs, ts, color=colr, lw=0.8, alpha=0.8)
        ax.axvline(res.chi[-1], color="white", lw=2)
        if np.isfinite(res.bounce_tau).any():
            ax.axhline(np.nanmean(res.bounce_tau), color=C_WH, ls=":", lw=1.2)
            ax.text(0.05, np.nanmean(res.bounce_tau) + 0.4, "ricochete (R minimo)", color=C_WH, fontsize=8)
        ax.set_xlabel("coordenada comovel chi (superficie em chi = 1)")
        ax.set_ylabel("tempo proprio tau [M]")
        ax.set_title(f"{title}: raios nulos de saida (amarelo) e de entrada (violeta),\nregiao presa (azul)")
        ax.grid(alpha=0.2)
    fig.suptitle(f"Estrutura causal do interior, M = {M}, R0 = {R0} M  |  G = c = 1", fontsize=11)
    fig.tight_layout()
    stamp_figure(fig, _record(b))
    return _save(fig, "23_estrutura_causal_colapso.png")


def make_all():
    fig_collapse_comparison()
    fig_causal_structure()
