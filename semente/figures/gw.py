"""Figuras da Fase 6: QNMs, ringdown no dominio do tempo, ecos e ponte observacional de GW."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.provenance import RunRecord, save_json, stamp_figure  # noqa: E402
from ..core.units import Scale  # noqa: E402
from ..gravitational_waves.qnm import LEAVER_SCHWARZSCHILD, kerr_bcw, wkb3  # noqa: E402
from ..gravitational_waves.ringdown import classical_template, evolve  # noqa: E402
from ..observations.bridge import compare, table  # noqa: E402
from ..observations.constraints import CONSTRAINTS  # noqa: E402
from ..stability.perturbations import schwarzschild_metric, simpson_visser_metric  # noqa: E402
from .core import C_BH, C_OTHER, C_OUR, C_SING, C_WH, OUT, _save  # noqa: E402


def fig_qnm():
    ls = np.concatenate([[1e-4], np.linspace(0.1, 1.95, 30)])
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    rows = []
    q0 = {s: wkb3(schwarzschild_metric(1.0), spin=s, ell=2) for s in (0, 1)}
    for spin, col in [(0, C_OUR), (1, C_WH)]:
        fr, dr = [], []
        for l in ls:
            q = wkb3(simpson_visser_metric(1.0, l), spin=spin, ell=2)
            fr.append(q.frequency / q0[spin].frequency - 1)
            dr.append(q.damping_rate / q0[spin].damping_rate - 1)
            rows.append(dict(l=float(l), spin=spin, omega_re=q.frequency, omega_im=-q.damping_rate))
        axes[0].plot(ls, np.array(fr) * 100, color=col, lw=2, label=f"spin {spin}: delta f / f")
        axes[0].plot(ls, np.array(dr) * 100, color=col, lw=2, ls="--", label=f"spin {spin}: delta (1/tau) / (1/tau)")
    axes[0].axhline(0, color="white", lw=0.8)
    axes[0].set_xlabel("l / M (black-bounce; l = 2M: sem horizonte)")
    axes[0].set_ylabel("desvio relativo em relacao a Schwarzschild (%)")
    axes[0].set_title("QNM fundamental (l_mult = 2, WKB 3a ordem): frequencia quase insensivel a l,\namortecimento cai ate ~35% perto de l = 2M")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    ax = axes[1]
    chis = np.linspace(0, 0.98, 50)
    ax.plot(chis, [kerr_bcw(c).frequency for c in chis], color=C_OTHER, lw=2, label="M omega_R (Kerr, ajuste BCW 2006)")
    ax.plot(chis, [kerr_bcw(c).damping_rate for c in chis], color=C_OTHER, lw=2, ls="--", label="M omega_I")
    ax.axhline(LEAVER_SCHWARZSCHILD[(2, 2, 0)].real, color="white", ls=":", lw=1, label="Schwarzschild (Leaver)")
    ax.set_xlabel("spin adimensional chi = a/M")
    ax.set_title("Kerr l = m = 2 (literatura): a rotacao muda a frequencia em ate 2x.\nSem rotacao, nenhum template reproduz GW150914")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    Ms = np.logspace(0.5, 3, 100)
    for chi, col in [(0.0, "white"), (0.67, C_OTHER)]:
        ax.loglog(Ms, [classical_template(m, chi)["f_Hz"] for m in Ms], color=col, lw=2, label=f"f_ringdown, chi = {chi}")
    c = CONSTRAINTS["GW150914_f_ringdown_Hz"]
    ax.errorbar([62.0], [c.value], yerr=[c.sigma], fmt="o", color=C_SING, capsize=4, label="GW150914 (LIGO/Virgo 2016)")
    ax.set_xlabel("massa final (M_sol)")
    ax.set_ylabel("frequencia de ringdown (Hz)")
    ax.set_title("Templates classicos vs GW150914")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    rec = RunRecord(model="gravitational_waves/qnm", parameters=dict(M=1.0, ell=2))
    stamp_figure(fig, rec)
    path = _save(fig, "31_qnm_black_bounce_kerr.png")
    save_json(os.path.join(OUT, "31_qnm_black_bounce.json"), rows, rec)
    return path


def fig_ringdown_echoes():
    wf_bh = evolve(schwarzschild_metric(1.0), spin=0, ell=2, n_grid=4000, t_max=300.0)
    wf_bb = evolve(simpson_visser_metric(1.0, 1.0), spin=0, ell=2, n_grid=4000, t_max=300.0)
    wf_wh = evolve(simpson_visser_metric(1.0, 2.5), spin=0, ell=2, n_grid=4000, t_max=500.0)
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    ax = axes[0, 0]
    for wf, col, lbl in [(wf_bh, "white", "Schwarzschild"), (wf_bb, C_WH, "black-bounce l = 1M"), (wf_wh, C_OTHER, "buraco de minhoca l = 2.5M")]:
        ax.semilogy(wf.t, np.abs(wf.psi) / np.max(np.abs(wf.psi)) + 1e-9, color=col, lw=1.2, label=lbl)
    for e in wf_wh.echoes[:5]:
        ax.axvline(e[0], color=C_OTHER, ls=":", lw=0.8)
    ax.set_ylim(1e-6, 2)
    ax.set_xlabel("t [M]")
    ax.set_ylabel("|Psi| / max")
    ax.set_title(f"Ringdown (campo escalar, l_mult = 2): ecos so no buraco de minhoca\n(atraso previsto 2 dr_* = {wf_wh.echo_delay_predicted:.1f} M; medido ~ {np.mean(np.diff([e[0] for e in wf_wh.echoes[:4]])):.1f} M)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[0, 1]
    ax.plot(wf_bh.t, wf_bh.psi / np.max(np.abs(wf_bh.psi)), color="white", lw=1, label="Schwarzschild")
    ax.plot(wf_bb.t, wf_bb.psi / np.max(np.abs(wf_bb.psi)), color=C_WH, lw=1, ls="--", label="black-bounce l = 1M")
    ax.set_xlim(wf_bh.t[np.argmax(np.abs(wf_bh.psi))] - 10, wf_bh.t[np.argmax(np.abs(wf_bh.psi))] + 80)
    q_bh, q_bb = wf_bh.fitted_omega, wf_bb.fitted_omega
    ax.set_title(f"Forma de onda: omega_fit = {q_bh.real:.3f} - {-q_bh.imag:.3f}i (Schw), {q_bb.real:.3f} - {-q_bb.imag:.3f}i (l = 1M)\nWKB Schw: {wkb3(schwarzschild_metric(1.0), 0, 2).omega:.3f}")
    ax.set_xlabel("t [M]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    for ax, wf, lbl in [(axes[1, 0], wf_bh, "Schwarzschild"), (axes[1, 1], wf_wh, "buraco de minhoca l = 2.5M")]:
        f, tt, Z = wf.spectrogram(nperseg=512)
        s = Scale.from_solar_masses(62.0)
        ax.pcolormesh(tt * s.time_unit * 1e3, f / s.time_unit, np.log10(Z + 1e-12), shading="auto", cmap="magma")
        ax.set_ylim(0, 800)
        ax.set_xlabel("tempo (ms) para M = 62 M_sol")
        ax.set_ylabel("frequencia (Hz)")
        ax.set_title(f"Espectrograma, {lbl}: escala de GW150914")
    fig.tight_layout()
    rec = RunRecord(model="gravitational_waves/ringdown", parameters=dict(M=1.0, ell=2, spin=0))
    stamp_figure(fig, rec)
    path = _save(fig, "32_ringdown_ecos.png")
    # ponte observacional (massa no referencial do detector: M_fonte (1 + z))
    z = CONSTRAINTS["GW150914_z"].value
    M_det = CONSTRAINTS["GW150914_Mf_Msun"].value * (1 + z)
    dM = CONSTRAINTS["GW150914_Mf_Msun"].sigma / CONSTRAINTS["GW150914_Mf_Msun"].value
    t0 = classical_template(M_det, 0.0)
    t1 = classical_template(M_det, CONSTRAINTS["GW150914_chi_f"].value)
    rows = [compare("Schwarzschild (sem rotacao)", "GW150914_f_ringdown_Hz", t0["f_Hz"], t0["f_Hz"] * dM,
                    note=f"M_det = {M_det:.1f} M_sol; modelo sem rotacao: a discrepancia mede o spin, nao o ricochete"),
            compare("Schwarzschild (sem rotacao)", "GW150914_tau_ms", t0["tau_ms"], t0["tau_ms"] * dM),
            compare("Kerr BCW (literatura)", "GW150914_f_ringdown_Hz", t1["f_Hz"], t1["f_Hz"] * dM, literature=True),
            compare("Kerr BCW (literatura)", "GW150914_tau_ms", t1["tau_ms"], t1["tau_ms"] * dM, literature=True),
            compare("black-bounce l = R_b (juncao), sem rotacao", "GW150914_f_ringdown_Hz", t0["f_Hz"], t0["f_Hz"] * dM,
                    note="l/M ~ 1e-14: desvio de QNM < 1e-6, indistinguivel do classico; falta rotacao (fase 7)")]
    save_json(os.path.join(OUT, "32_gw_ponte.json"), [r.as_row() for r in rows], rec)
    with open(os.path.join(OUT, "32_gw_ponte.md"), "w", encoding="utf-8") as fh:
        fh.write(table(rows))
    return path


def make_all():
    fig_qnm()
    fig_ringdown_echoes()
