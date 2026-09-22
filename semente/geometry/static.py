"""
Metricas estaticas esfericamente simetricas GENERICAS:  ds^2 = -f(r) dt^2 + dr^2/f(r) + R(r)^2 dOmega^2,
com f e R dados como expressoes sympy.  Fornece, para qualquer membro da familia:
   horizontes (raizes de f), classificacao causal, esfera de fotons (extremos de f/R^2),
   escalar de Kretschmann (base ortonormal, derivadas simbolicas), condicao da garganta.

Kretschmann (RESULTADO MATEMATICO padrao para esta forma de metrica, ja testado em geometry/spherical.py):
   A = f''/2,  B = f' R'/(2R),  C = -(f R'' + f' R'/2)/R,  D = (1 - f R'^2)/R^2,   K = 4A^2 + 8B^2 + 8C^2 + 4D^2.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Sequence

import numpy as np
import sympy as sp
from scipy.optimize import brentq

_r = sp.symbols("r", real=True)


@lru_cache(maxsize=None)
def _compile(f_str: str, R_str: str, names: tuple):
    params = sp.symbols(names, positive=True) if names else ()
    loc = {"r": _r, **{str(p): p for p in params}}
    f = sp.sympify(f_str, locals=loc)
    R = sp.sympify(R_str, locals=loc)
    fp, fpp = sp.diff(f, _r), sp.diff(f, _r, 2)
    Rp, Rpp = sp.diff(R, _r), sp.diff(R, _r, 2)
    A = fpp / 2
    B = fp * Rp / (2 * R)
    C = -(f * Rpp + fp * Rp / 2) / R
    D = (1 - f * Rp**2) / R**2
    K = 4 * A**2 + 8 * B**2 + 8 * C**2 + 4 * D**2
    g = f / R**2
    args = (_r, *params)
    return dict(f=sp.lambdify(args, f, "numpy"), R=sp.lambdify(args, R, "numpy"), fp=sp.lambdify(args, fp, "numpy"),
                K=sp.lambdify(args, K, "numpy"), dg=sp.lambdify(args, sp.diff(g, _r), "numpy"))


@dataclass
class StaticSphericalMetric:
    f_expr: str
    R_expr: str
    param_names: tuple
    params: dict
    name: str = "metrica"
    r_domain: tuple = (-np.inf, np.inf)   # (-inf, inf) para black-bounces; (0, inf) para r = R

    def __post_init__(self):
        self._c = _compile(self.f_expr, self.R_expr, tuple(self.param_names))
        self._args = [self.params[n] for n in self.param_names]

    def f(self, r):
        return np.asarray(self._c["f"](np.asarray(r, float), *self._args), float)

    def R(self, r):
        return np.asarray(self._c["R"](np.asarray(r, float), *self._args), float)

    def fprime(self, r):
        return np.asarray(self._c["fp"](np.asarray(r, float), *self._args), float)

    def kretschmann(self, r):
        with np.errstate(all="ignore"):
            return np.asarray(self._c["K"](np.asarray(r, float), *self._args), float)

    def horizons(self, r_max=None, n=20000) -> list:
        """Raizes de f(r) no dominio (varredura + brentq)."""
        M = self.params.get("M", 1.0)
        r_max = 20 * M if r_max is None else r_max
        lo = 1e-6 * M if self.r_domain[0] == 0 else -r_max
        grid = np.linspace(lo, r_max, n)
        with np.errstate(all="ignore"):
            fv = self.f(grid)
        roots = []
        for i in range(n - 1):
            if np.isfinite(fv[i]) and np.isfinite(fv[i + 1]) and fv[i] * fv[i + 1] < 0:
                roots.append(float(brentq(lambda x: float(self.f(x)), grid[i], grid[i + 1])))
        return sorted(roots)

    def photon_spheres(self, r_max=None, n=20000) -> list:
        """Extremos de f/R^2 (orbitas circulares de fotons) fora dos horizontes."""
        M = self.params.get("M", 1.0)
        r_max = 20 * M if r_max is None else r_max
        lo = 1e-6 * M if self.r_domain[0] == 0 else -r_max
        grid = np.linspace(lo, r_max, n)
        with np.errstate(all="ignore"):
            dg = np.asarray(self._c["dg"](grid, *self._args), float)
        out = []
        for i in range(n - 1):
            if np.isfinite(dg[i]) and np.isfinite(dg[i + 1]) and dg[i] * dg[i + 1] < 0:
                out.append(float(grid[i]))
        return out

    def surface_gravity(self, r_h):
        return float(abs(self.fprime(r_h)) / 2)

    def classify(self) -> dict:
        """Estrutura causal a partir dos horizontes e do sinal de f na garganta (r = 0) quando aplicavel."""
        hs = self.horizons()
        M = self.params.get("M", 1.0)
        info = dict(horizons=hs, n_horizons=len(hs))
        if self.r_domain[0] == 0:
            # familia r = R (Schwarzschild, Reissner-Nordstrom): singularidade em r = 0
            K0 = float(self.kretschmann(1e-3 * M))
            info.update(regular_center=bool(np.isfinite(K0) and K0 < 1e12), throat=None)
            if len(hs) == 0:
                info["kind"] = "singularidade nua" if not info["regular_center"] else "solucao regular sem horizonte"
            elif len(hs) == 1:
                info["kind"] = "buraco negro (um horizonte)"
            else:
                info["kind"] = "buraco negro com horizonte de Cauchy interno"
        else:
            f0 = float(self.f(0.0))
            K0 = float(self.kretschmann(0.0))
            info.update(regular_center=bool(np.isfinite(K0)), throat_f=f0, kretschmann_throat=K0)
            pos = [h for h in hs if h > 0]
            if len(hs) == 0:
                info["kind"] = "buraco de minhoca atravessavel" if f0 > 0 else "sem horizonte (f<0 na garganta: nao fisico)"
            elif f0 < 0:
                info["kind"] = ("black-bounce: garganta espacial (ricochete BN -> BB)" if len(pos) == 1
                                else "black-bounce com horizonte interno: garganta espacial entre horizontes de Cauchy")
            else:
                info["kind"] = "buraco negro regular com garganta tipo-tempo (buraco de minhoca escondido por horizontes)"
        return info
