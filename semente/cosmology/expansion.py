"""
Expansao emergente apos o ricochete: o codigo DETERMINA se ha uma fase acelerada e quanto ela dura.
Nada e assumido.

Para um BackgroundSolution:
    H(t),  epsilon(t) = -Hdot/H^2,  N = int H dt,
    fases aceleradas = intervalos com addot > 0 (equivale a epsilon < 1 quando H > 0),
    entrada/saida: instantes em que epsilon cruza 1 (para H > 0).

Criterio de "inflacao" (documentado): fase acelerada pos-ricochete com N >= 60 e-folds.  Um ricochete
generico produz uma fase de super-aceleracao curta (Hdot > 0) em torno de t_b com N ~ O(1): isso NAO e
inflacao e e reportado como tal.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .friedmann import BackgroundSolution, Fluid, FriedmannModel
from ..quantum.einstein_cartan import EinsteinCartanCorrection
from ..quantum.lqc import LQCCorrection


@dataclass
class ExpansionReport:
    background: str
    t_bounce: Optional[float]
    phases: list                      # [(t_ini, t_fim, N, eps_min)]
    N_total_after_bounce: float
    N_accelerated: float
    inflation_like: bool
    eps_min_after_bounce: float
    duration_accelerated: float
    entry_exit: list                  # instantes em que epsilon cruza 1 (H>0)
    note: str = ""


def analyze_expansion(sol: BackgroundSolution, N_inflation=60.0) -> ExpansionReport:
    ib = sol.i_bounce
    if ib is None:
        return ExpansionReport(sol.model, None, [], float("nan"), 0.0, False, float("nan"), 0.0, [],
                               note="sem ricochete: nao ha fase de expansao emergente")
    t, eps, H, acc = sol.t[ib:], sol.epsilon[ib:], sol.H[ib:], sol.accel[ib:]
    phases = []
    N_acc = 0.0
    dur = 0.0
    i0 = None
    for i in range(t.size):
        on = acc[i] > 0
        if on and i0 is None:
            i0 = i
        if (not on or i == t.size - 1) and i0 is not None:
            N = float(np.trapezoid(H[i0:i + 1], t[i0:i + 1]))
            em = float(np.nanmin(eps[i0:i + 1])) if np.isfinite(eps[i0:i + 1]).any() else float("nan")
            phases.append((float(t[i0]), float(t[i]), N, em))
            N_acc += N
            dur += float(t[i] - t[i0])
            i0 = None
    N_tot = float(np.trapezoid(H, t))
    with np.errstate(invalid="ignore"):
        pos = (H > 0) & np.isfinite(eps)
    crossings = []
    e = np.where(pos, eps, np.nan)
    for i in range(1, t.size):
        if np.isfinite(e[i - 1]) and np.isfinite(e[i]) and (e[i - 1] - 1) * (e[i] - 1) < 0:
            crossings.append(float(t[i]))
    infl = N_acc >= N_inflation
    note = ("fase acelerada satisfaz N >= 60: inflacao emergente" if infl else
            f"super-aceleracao curta em torno do ricochete com N = {N_acc:.2f} e-folds: NAO e inflacao")
    return ExpansionReport(sol.model, sol.t_bounce, phases, N_tot, N_acc, infl,
                           float(np.nanmin(eps)) if np.isfinite(eps).any() else float("nan"), dur, crossings, note)


def scan_parameters(rho_stars=(2, 5, 20, 100, 1000), ws=(0.0, 1 / 3), mechanism="lqc", rho0=1.0, t_span=(0, 3.0)):
    """Dependencia dos e-folds acelerados com rho_* e w.  Devolve lista de dicts."""
    out = []
    for w in ws:
        for rs in rho_stars:
            fl = (Fluid(rho0, w),)
            corr = LQCCorrection(rho_c=rs) if mechanism == "lqc" else EinsteinCartanCorrection(sigma=rho0**2 / rs)
            m = FriedmannModel(fl, 0.0, corr)
            sol = m.solve(a0=1.0, t_span=t_span, contracting=True, n_samples=3000)
            rep = analyze_expansion(sol)
            out.append(dict(mechanism=mechanism, w=w, rho_star=rs, a_min=sol.a_min, N_accelerated=rep.N_accelerated,
                            duration=rep.duration_accelerated, eps_min=rep.eps_min_after_bounce, inflation_like=rep.inflation_like))
    return out
