"""
Estabilidade linear de metricas estaticas esfericamente simetricas por potenciais mestres.

Para ds^2 = -f dt^2 + dr^2/f + R(r)^2 dOmega^2, com coordenada tartaruga dr_* = dr/f, perturbacoes
de campos de teste obedecem a equacao mestre

    -d^2 Psi/dt^2 + d^2 Psi/dr_*^2 - V(r) Psi = 0 .

Potenciais implementados (RESULTADOS MATEMATICOS padrao para esta forma de metrica):
   spin 0 (escalar sem massa, Psi = R phi):  V_0 = f l(l+1)/R^2 + (f/R) d/dr( f dR/dr )
   spin 1 (eletromagnetico):                 V_1 = f l(l+1)/R^2
   spin 2 axial, SO para regioes de VACUO (Schwarzschild, R = r):
                                             V_2 = f [ l(l+1)/r^2 - 6M/r^3 ]          (Regge-Wheeler 1957)
   Para metricas com materia efetiva (black-bounce) o potencial gravitacional axial depende de como
   a materia responde a perturbacao; NAO esta implementado aqui (limitacao documentada).  Usa-se
   spin 0 e 1 como campos de teste.

Criterio de estabilidade (documentado, nao arbitrario):
   Psi = e^{-i omega t} psi(r_*)  =>  -psi'' + V psi = omega^2 psi.  Um modo com omega^2 < 0 cresce como
   e^{|omega| t}.  Logo existe instabilidade linear (para esse campo) se e somente se o operador de
   Schrodinger  H = -d^2/dr_*^2 + V  tem um autovalor negativo (estado ligado).  Se V >= 0 em todo
   lugar, H >= 0 e nao ha modos crescentes.  O menor autovalor e obtido por diferencas finitas numa grade
   uniforme em r_*; a classificacao exige CONVERGENCIA sob refinamento da grade.

Classes de saida:
   stable                  V >= 0 em todo lugar, ou menor autovalor >= -tol com convergencia
   weakly_unstable         autovalor negativo convergido com tempo de crescimento 1/|omega| > 100 M
   strongly_unstable       autovalor negativo convergido com 1/|omega| <= 100 M
   numerically_unresolved  resultados nao convergem entre resolucoes

Referencias: Regge & Wheeler (1957) Phys. Rev. 108, 1063; Chandrasekhar, "The Mathematical Theory of Black
Holes" (1983); Bronnikov, Konoplya & Zhidenko (2012) Phys. Rev. D 86, 024028 (criterio de estado ligado
para buracos de minhoca); Simpson & Visser (2019).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import eigh_tridiagonal


@dataclass
class StaticMetric:
    """f(r), R(r) e derivadas numericas (ou analiticas quando dadas)."""

    f: Callable
    R: Callable
    name: str = "metrica"
    M: float = 1.0
    horizons: tuple = ()

    def dR(self, r, h=1e-5):
        return (self.R(r + h) - self.R(r - h)) / (2 * h)

    def d(self, g, r, h=1e-5):
        return (g(r + h) - g(r - h)) / (2 * h)


def simpson_visser_metric(M=1.0, l=0.5) -> StaticMetric:
    from ..geometry.spherical import BlackBounce
    bb = BlackBounce(M, l)
    return StaticMetric(f=lambda r: bb.f(r), R=lambda r: bb.R(r), name=f"Simpson-Visser(M={M:g}, l={l:g})", M=M, horizons=bb.horizons)


def schwarzschild_metric(M=1.0) -> StaticMetric:
    return StaticMetric(f=lambda r: 1 - 2 * M / np.asarray(r, float), R=lambda r: np.asarray(r, float),
                        name=f"Schwarzschild(M={M:g})", M=M, horizons=(2 * M,))


def master_potential(metric: StaticMetric, r, spin: int = 0, ell: int = 2):
    r = np.asarray(r, float)
    f = metric.f(r)
    R = metric.R(r)
    if spin == 0:
        g = lambda x: metric.f(x) * metric.dR(x)      # noqa: E731
        return f * ell * (ell + 1) / R**2 + f / R * metric.d(g, r)
    if spin == 1:
        return f * ell * (ell + 1) / R**2
    if spin == 2:
        if not np.allclose(R, r):
            raise NotImplementedError("Potencial gravitacional axial so implementado para vacuo (R = r).")
        return f * (ell * (ell + 1) / r**2 - 6 * metric.M / r**3)
    raise ValueError("spin deve ser 0, 1 ou 2")


def tortoise_grid(metric: StaticMetric, r_min: float, r_max: float, n: int = 4000):
    """Grade em r e a coordenada tartaruga r_*(r) = int dr/f (regiao sem horizontes no intervalo)."""
    r = np.linspace(r_min, r_max, n)
    f = metric.f(r)
    if np.any(f <= 0):
        raise ValueError("intervalo contem horizonte (f <= 0); escolha uma regiao exterior ou o dominio completo do buraco de minhoca")
    rs = cumulative_trapezoid(1 / f, r, initial=0.0)
    return r, rs


@dataclass
class StabilityResult:
    metric: str
    spin: int
    ell: int
    classification: str
    lowest_eigenvalue: float
    growth_rate: Optional[float]
    growth_time_M: Optional[float]
    potential_min: float
    potential_negative_region: Optional[tuple]
    convergence: dict
    criterion: str = "estado ligado de -d2/dr*2 + V (autovalor negativo <=> modo crescente)"


def lowest_eigenvalue(V_rs: np.ndarray, rs: np.ndarray) -> float:
    """Menor autovalor de -d^2/dr_*^2 + V com Dirichlet nas bordas (grade uniforme em r_*)."""
    h = rs[1] - rs[0]
    diag = 2.0 / h**2 + V_rs
    off = -np.ones(V_rs.size - 1) / h**2
    w = eigh_tridiagonal(diag, off, select="i", select_range=(0, 0), eigvals_only=True)
    return float(w[0])


def analyze(metric: StaticMetric, spin=0, ell=2, r_range=None, n=(2000, 4000, 8000), tol=1e-6,
            weak_threshold_M=100.0) -> StabilityResult:
    """Analise completa com verificacao de convergencia em tres resolucoes."""
    if r_range is None:
        if metric.horizons:
            rh = max(metric.horizons)
            r_range = (rh * (1 + 1e-3), rh + 60 * metric.M)
        else:
            r_range = (-60 * metric.M, 60 * metric.M)   # buraco de minhoca: dominio completo em r
    eigs = []
    Vmin = None
    neg = None
    for nn in n:
        r, rs = tortoise_grid(metric, r_range[0], r_range[1], nn)
        V = master_potential(metric, r, spin, ell)
        # reamostra V numa grade UNIFORME em r_*
        rs_u = np.linspace(rs[0], rs[-1], nn)
        V_u = np.interp(rs_u, rs, V)
        eigs.append(lowest_eigenvalue(V_u, rs_u))
        if Vmin is None:
            Vmin = float(np.min(V))
            idx = np.where(V < 0)[0]
            neg = (float(r[idx[0]]), float(r[idx[-1]])) if idx.size else None
    e_last = eigs[-1]
    converged = abs(eigs[-1] - eigs[-2]) < max(tol, 0.05 * abs(eigs[-1])) if len(eigs) > 1 else True
    if Vmin >= -1e-14:
        cls = "stable"
        rate = None
    elif not converged and e_last < -tol:
        cls = "numerically_unresolved"
        rate = None
    elif e_last >= -tol:
        cls = "stable"
        rate = None
    else:
        rate = float(np.sqrt(-e_last))
        cls = "weakly_unstable" if 1 / rate > weak_threshold_M * metric.M else "strongly_unstable"
    return StabilityResult(metric=metric.name, spin=spin, ell=ell, classification=cls, lowest_eigenvalue=e_last,
                           growth_rate=rate, growth_time_M=(None if rate is None else 1 / rate / metric.M),
                           potential_min=Vmin, potential_negative_region=neg,
                           convergence=dict(n=list(n), eigenvalues=eigs, converged=bool(converged)))
