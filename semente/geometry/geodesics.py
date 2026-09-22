"""
Geodesicas: (1) integrador GERAL para qualquer metrica diagonal, com os
simbolos de Christoffel derivados simbolicamente (sympy) e compilados;
(2) a forma orbital u(phi) para Schwarzschild, usada em massa pelo ray-tracer.

Equacao da geodesica:   d^2x^mu/dlambda^2 + Gamma^mu_{ab} dx^a/dlambda dx^b/dlambda = 0.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


# ----------------------------------------------------------------------------
# Metricas simbolicas
# ----------------------------------------------------------------------------
t, r, th, ph = sp.symbols("t r theta phi", real=True)
M_s, l_s = sp.symbols("M ell", positive=True)
COORDS = (t, r, th, ph)


def metric_schwarzschild():
    f = 1 - 2 * M_s / r
    return sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)


def metric_black_bounce():
    R2 = r**2 + l_s**2
    f = 1 - 2 * M_s / sp.sqrt(R2)
    return sp.diag(-f, 1 / f, R2, R2 * sp.sin(th) ** 2)


def christoffel(g: sp.Matrix):
    """Gamma^a_{bc} = 1/2 g^{ad} (d_b g_{dc} + d_c g_{db} - d_d g_{bc})."""
    ginv = g.inv()
    n = 4
    Gam = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                expr = sum(
                    ginv[a, d]
                    * (sp.diff(g[d, c], COORDS[b]) + sp.diff(g[d, b], COORDS[c]) - sp.diff(g[b, c], COORDS[d]))
                    for d in range(n)
                ) / 2
                expr = sp.simplify(expr)
                Gam[a][b][c] = expr
                Gam[a][c][b] = expr
    return Gam


def kretschmann_symbolic(g: sp.Matrix):
    """Escalar de Kretschmann R_{abcd} R^{abcd} calculado do zero (lento, mas exato)."""
    n = 4
    Gam = christoffel(g)
    ginv = g.inv()
    # R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}
    Riem = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    expr = sp.diff(Gam[a][d][b], COORDS[c]) - sp.diff(Gam[a][c][b], COORDS[d])
                    expr += sum(Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b] for e in range(n))
                    Riem[(a, b, c, d)] = expr
    # baixar / levantar indices (metrica diagonal simplifica)
    K = sp.S(0)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    low = g[a, a] * Riem[(a, b, c, d)]  # R_{abcd}
                    up = ginv[b, b] * ginv[c, c] * ginv[d, d] * ginv[a, a] * low  # R^{abcd}
                    K += low * up
    return sp.simplify(K)


@lru_cache(maxsize=None)
def compiled_geodesic_rhs(name: str):
    """Compila o lado direito da equacao da geodesica para a metrica `name`."""
    g = {"schwarzschild": metric_schwarzschild, "black_bounce": metric_black_bounce}[name]()
    Gam = christoffel(g)
    v = sp.symbols("v0:4", real=True)
    acc = []
    for a in range(4):
        expr = -sum(Gam[a][b][c] * v[b] * v[c] for b in range(4) for c in range(4))
        acc.append(sp.simplify(expr))
    params = (M_s, l_s)
    fn = sp.lambdify((COORDS, v, params), acc, modules="numpy")
    return fn


@dataclass
class GeodesicSolver:
    """Integra geodesicas em qualquer metrica registrada (Schwarzschild ou black-bounce)."""

    metric: str = "schwarzschild"
    M: float = 1.0
    l: float = 0.0

    def rhs(self, lam, y):
        x, v = y[:4], y[4:]
        a = compiled_geodesic_rhs(self.metric)(x, v, (self.M, self.l))
        return np.concatenate([v, np.asarray(a, float)])

    def initial_velocity(self, x0, dr, dphi, kind="timelike", dth=0.0):
        """Fixa dt/dlambda pela normalizacao g_{ab} v^a v^b = -1 (tipo tempo) ou 0 (nula)."""
        if self.metric == "schwarzschild":
            rr = x0[1]
            f = 1 - 2 * self.M / rr
            R2 = rr**2
        else:
            R2 = x0[1] ** 2 + self.l**2
            f = 1 - 2 * self.M / np.sqrt(R2)
        eps = -1.0 if kind == "timelike" else 0.0
        spatial = dr**2 / f + R2 * (dth**2 + np.sin(x0[2]) ** 2 * dphi**2)
        dt = np.sqrt((spatial - eps) / f)
        return np.array([dt, dr, dth, dphi])

    def integrate(self, x0, v0, lam_max, stop_r=None, max_step=0.05, rtol=1e-9, atol=1e-11):
        y0 = np.concatenate([np.asarray(x0, float), np.asarray(v0, float)])
        events = []
        if stop_r is not None:
            def ev(lam, y):
                return y[1] - stop_r
            ev.terminal = True
            events.append(ev)
        sol = solve_ivp(self.rhs, (0, lam_max), y0, method="DOP853", max_step=max_step,
                        rtol=rtol, atol=atol, events=events, dense_output=True)
        return sol


# ----------------------------------------------------------------------------
# Forma orbital de Schwarzschild: u = 1/r em funcao de phi
# ----------------------------------------------------------------------------
def orbit_rhs_null(u, M):
    """d^2u/dphi^2 = -u + 3 M u^2   (fotons)."""
    return -u + 3 * M * u**2


def orbit_rhs_timelike(u, M, L):
    """d^2u/dphi^2 = -u + M/L^2 + 3 M u^2   (particulas massivas)."""
    return -u + M / L**2 + 3 * M * u**2


def deflection_angle(b, M=1.0, dphi=1e-3, max_turns=3.0):
    """Angulo de deflexao de um raio de luz com parametro de impacto b (integracao numerica de u(phi)).

    Para b >> M: alpha ~ 4M/b (Einstein).  Para b -> 3*sqrt(3) M (esfera de fotons) diverge (aneis).
    """
    u = 0.0
    up = 1.0 / b  # du/dphi em u=0:  (du/dphi)^2 = 1/b^2 - u^2 + 2Mu^3
    phi = 0.0
    turned = False
    while phi < max_turns * 2 * np.pi:
        u_prev = u
        # RK4
        k1u, k1p = up, orbit_rhs_null(u, M)
        k2u, k2p = up + 0.5 * dphi * k1p, orbit_rhs_null(u + 0.5 * dphi * k1u, M)
        k3u, k3p = up + 0.5 * dphi * k2p, orbit_rhs_null(u + 0.5 * dphi * k2u, M)
        k4u, k4p = up + dphi * k3p, orbit_rhs_null(u + dphi * k3u, M)
        u += dphi * (k1u + 2 * k2u + 2 * k3u + k4u) / 6
        up += dphi * (k1p + 2 * k2p + 2 * k3p + k4p) / 6
        phi += dphi
        if u > 0.5 / M:  # cruzou o horizonte
            return np.nan
        if u < 0 and turned:
            # interpola o cruzamento u=0 dentro do ultimo passo
            frac = u_prev / (u_prev - u)
            return (phi - dphi + frac * dphi) - np.pi
        if up < 0:
            turned = True
    return np.nan
