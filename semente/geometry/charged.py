"""
Geometrias com carga: Reissner-Nordstrom e o black-bounce carregado de Franzin et al. (2021).

Reissner-Nordstrom (RESULTADO MATEMATICO):  f = 1 - 2M/r + Q^2/r^2,  R = r.
   horizontes r_+- = M +- sqrt(M^2 - Q^2)  (Q < M: dois; Q = M: extremal; Q > M: singularidade nua)
   K = 8 (6 M^2 r^2 - 12 M Q^2 r + 7 Q^4) / r^8
   fluido efetivo eletromagnetico rho = -p_r = p_t = Q^2/(8 pi r^4)

Black-bounce carregado (MODELO DA LITERATURA, Franzin, Liberati, Mazza, Simpson & Visser 2021, JCAP 07, 036):
   f = 1 - 2M/sqrt(r^2 + l^2) + Q^2/(r^2 + l^2),   R = sqrt(r^2 + l^2)
   horizontes onde R(r) = R_+- = M +- sqrt(M^2 - Q^2):  r_h = +- sqrt(R_+-^2 - l^2)  (se R_+- > l).
   Classificacao em (Q, l):  l < R_-: dois pares de horizontes (garganta entre os horizontes de Cauchy);
   R_- < l < R_+: um par (garganta espacial, ricochete BN -> BB);  l > R_+: buraco de minhoca.
"""
from __future__ import annotations

import numpy as np

from .static import StaticSphericalMetric


def reissner_nordstrom(M=1.0, Q=0.5) -> StaticSphericalMetric:
    return StaticSphericalMetric("1 - 2*M/r + Q**2/r**2", "r", ("M", "Q"), dict(M=M, Q=Q),
                                 name=f"Reissner-Nordstrom(M={M:g}, Q={Q:g})", r_domain=(0, np.inf))


def charged_black_bounce(M=1.0, Q=0.5, l=0.5) -> StaticSphericalMetric:
    return StaticSphericalMetric("1 - 2*M/sqrt(r**2 + l**2) + Q**2/(r**2 + l**2)", "sqrt(r**2 + l**2)", ("M", "Q", "l"),
                                 dict(M=M, Q=Q, l=l), name=f"charged black-bounce(M={M:g}, Q={Q:g}, l={l:g})")


def rn_horizons_analytic(M=1.0, Q=0.5):
    if Q > M:
        return ()
    s = np.sqrt(M**2 - Q**2)
    return (M - s, M + s)


def rn_kretschmann_analytic(r, M=1.0, Q=0.5):
    r = np.asarray(r, float)
    return 8 * (6 * M**2 * r**2 - 12 * M * Q**2 * r + 7 * Q**4) / r**8


def charged_bounce_phase(M=1.0, Q=0.5, l=0.5) -> str:
    """Fase do black-bounce carregado no plano (Q, l)."""
    if Q > M:
        return "sem horizontes (Q > M): buraco de minhoca regular" if l > 0 else "singularidade nua (RN, Q > M)"
    Rm, Rp = rn_horizons_analytic(M, Q)
    if l <= 0:
        return "Reissner-Nordstrom"
    if l < Rm:
        return "regular BH com horizonte interno (garganta entre horizontes de Cauchy)"
    if l < Rp:
        return "black-bounce (garganta espacial, ricochete BN -> BB)"
    return "buraco de minhoca atravessavel"
