"""
Condicoes de energia (NEC, WEC, SEC, DEC) avaliadas automaticamente.

1) Metricas estaticas esfericamente simetricas  ds^2 = -f dt^2 + dr^2/f + R(r)^2 dOmega^2.
   O tensor de Einstein e derivado SIMBOLICAMENTE (sympy) para f(r), R(r) genericos e depois
   avaliado para cada metrica; nao ha formulas transcritas de memoria.  O "fluido efetivo" e
       rho = -G^t_t / 8 pi,   p_r = G^r_r / 8 pi,   p_t = G^theta_theta / 8 pi   (base ortonormal).
   Classificacao: RESULTADO MATEMATICO/NUMERICO sobre a metrica dada.  Para Schwarzschild tudo e zero;
   para Simpson-Visser (l > 0) a NEC e violada perto da garganta (Simpson & Visser 2019, sec. 4).

2) Fundos de Friedmann (BackgroundSolution): fluido efetivo TOTAL (materia + correcao)
       rho_eff = 3 (H^2 + k/a^2) / 8 pi,     p_eff = -(2 addot/a + H^2 + k/a^2) / 8 pi,
   avaliado ao longo de t.  Num ricochete com k = 0, H = 0 e Hdot > 0 => rho_eff + p_eff = -Hdot/4pi < 0:
   a NEC e necessariamente violada (consistente com os teoremas de singularidade).

Definicoes (fluido anisotropico, base ortonormal):
   NEC: rho + p_r >= 0  e  rho + p_t >= 0
   WEC: NEC  e  rho >= 0
   SEC: NEC  e  rho + p_r + 2 p_t >= 0
   DEC: rho >= |p_r|  e  rho >= |p_t|
Referencia: Hawking & Ellis (1973), cap. 4; Visser, "Lorentzian Wormholes" (1995), cap. 12.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------------
# Tensor de Einstein simbolico para -f dt^2 + dr^2/f + R^2 dOmega^2
# ----------------------------------------------------------------------------
_r, _th = sp.symbols("r theta", real=True)
_f = sp.Function("f")(_r)
_R = sp.Function("R")(_r)


@lru_cache(maxsize=None)
def _einstein_mixed_symbolic():
    """Devolve (G^t_t, G^r_r, G^theta_theta) como expressoes em f, f', f'', R, R', R''."""
    t, ph = sp.symbols("t phi", real=True)
    x = (t, _r, _th, ph)
    g = sp.diag(-_f, 1 / _f, _R**2, _R**2 * sp.sin(_th) ** 2)
    ginv = g.inv()
    n = 4
    Gam = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                Gam[a][b][c] = sum(ginv[a, d] * (sp.diff(g[d, c], x[b]) + sp.diff(g[d, b], x[c]) - sp.diff(g[b, c], x[d]))
                                   for d in range(n)) / 2
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            expr = 0
            for a in range(n):
                expr += sp.diff(Gam[a][b][d], x[a]) - sp.diff(Gam[a][b][a], x[d])
                for e in range(n):
                    expr += Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a]
            Ric[b, d] = expr
    Rs = sum(ginv[a, a] * Ric[a, a] for a in range(n))
    G = Ric - g * Rs / 2
    Gmixed = [sp.simplify(ginv[i, i] * G[i, i]) for i in range(3)]  # G^t_t, G^r_r, G^th_th
    return tuple(Gmixed)


@lru_cache(maxsize=None)
def effective_fluid_functions(f_expr_str: str, R_expr_str: str, param_names: tuple):
    """Compila (rho, p_r, p_t)(r, *params) para f(r) e R(r) dados como strings sympy."""
    params = sp.symbols(param_names, positive=True) if param_names else ()
    local = {"r": _r, **{str(p): p for p in params}}
    f_expr = sp.sympify(f_expr_str, locals=local)
    R_expr = sp.sympify(R_expr_str, locals=local)
    Gtt, Grr, Gthth = _einstein_mixed_symbolic()
    subs = {_f: f_expr, _R: R_expr}
    out = []
    for expr in (Gtt, Grr, Gthth):
        e = expr.subs(subs).doit()
        out.append(sp.simplify(e))
    rho = -out[0] / (8 * sp.pi)
    pr = out[1] / (8 * sp.pi)
    pt = out[2] / (8 * sp.pi)
    fn = sp.lambdify((_r, *params), (rho, pr, pt), modules="numpy")
    return fn, (rho, pr, pt)


METRIC_LIBRARY = {
    # nome: (f(r), R(r), parametros)
    "schwarzschild": ("1 - 2*M/r", "r", ("M",)),
    "simpson_visser": ("1 - 2*M/sqrt(r**2 + l**2)", "sqrt(r**2 + l**2)", ("M", "l")),
    "reissner_nordstrom": ("1 - 2*M/r + Q**2/r**2", "r", ("M", "Q")),
    "charged_black_bounce": ("1 - 2*M/sqrt(r**2 + l**2) + Q**2/(r**2 + l**2)", "sqrt(r**2 + l**2)", ("M", "l", "Q")),
}


@dataclass
class EnergyConditionReport:
    x: np.ndarray            # coordenada (r ou t)
    rho: np.ndarray
    p_r: np.ndarray
    p_t: np.ndarray
    NEC: np.ndarray          # bool
    WEC: np.ndarray
    SEC: np.ndarray
    DEC: np.ndarray
    label: str

    def violation_intervals(self, cond: str):
        """Lista de intervalos [x_ini, x_fim] onde a condicao `cond` e violada."""
        v = ~getattr(self, cond)
        out = []
        i = 0
        while i < v.size:
            if v[i]:
                j = i
                while j + 1 < v.size and v[j + 1]:
                    j += 1
                out.append((float(self.x[i]), float(self.x[j])))
                i = j + 1
            else:
                i += 1
        return out

    def summary(self) -> dict:
        return {c: dict(satisfied_everywhere=bool(getattr(self, c).all()), violated_fraction=float((~getattr(self, c)).mean()),
                        intervals=self.violation_intervals(c)[:6]) for c in ("NEC", "WEC", "SEC", "DEC")}


def _conditions(rho, pr, pt, tol=1e-12):
    nec = (rho + pr >= -tol) & (rho + pt >= -tol)
    wec = nec & (rho >= -tol)
    sec = nec & (rho + pr + 2 * pt >= -tol)
    dec = (rho >= np.abs(pr) - tol) & (rho >= np.abs(pt) - tol)
    return nec, wec, sec, dec


def static_metric_report(metric: str, r: np.ndarray, **params) -> EnergyConditionReport:
    f_s, R_s, names = METRIC_LIBRARY[metric]
    fn, _ = effective_fluid_functions(f_s, R_s, names)
    args = [params[n] for n in names]
    r = np.asarray(r, float)
    with np.errstate(all="ignore"):
        rho, pr, pt = fn(r, *args)
    rho, pr, pt = (np.broadcast_to(np.asarray(v, float), r.shape).copy() for v in (rho, pr, pt))
    nec, wec, sec, dec = _conditions(rho, pr, pt)
    label = f"{metric}(" + ", ".join(f"{n}={params[n]:g}" for n in names) + ")"
    return EnergyConditionReport(r, rho, pr, pt, nec, wec, sec, dec, label)


def friedmann_report(sol) -> EnergyConditionReport:
    """Fluido efetivo total de um BackgroundSolution (materia + correcao), isotropico (p_r = p_t)."""
    H2k = sol.H**2 + sol.k / sol.a**2
    rho = 3 * H2k / (8 * np.pi)
    p = -(2 * sol.accel + H2k) / (8 * np.pi)
    nec, wec, sec, dec = _conditions(rho, p, p)
    return EnergyConditionReport(sol.t, rho, p, p, nec, wec, sec, dec, f"FRW efetivo: {sol.model}")
