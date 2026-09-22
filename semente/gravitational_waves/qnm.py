"""
Modos quasi-normais (QNM) por WKB de 3a ordem (Iyer & Will 1987) para potenciais mestres de metricas
estaticas esfericamente simetricas, e ajustes publicados para Kerr (Berti, Cardoso & Will 2006).

WKB (Iyer & Will 1987, eq. 1.5), com derivadas de V em relacao a coordenada tartaruga no pico:
    omega^2 = [V0 + sqrt(-2 V0'') Lambda] - i nu sqrt(-2 V0'') (1 + Omega),   nu = n + 1/2,  alpha = nu^2
    Lambda = (1/sqrt(-2V0'')) [ (1/8)(V0''''/V0'')(1/4 + alpha) - (1/288)(V0'''/V0'')^2 (7 + 60 alpha) ]
    Omega  = (1/(-2V0'')) [ (5/6912)(V0'''/V0'')^4 (77 + 188 alpha) - (1/384)(V0'''^2 V0''''/V0''^3)(51 + 100 alpha)
             + (1/2304)(V0''''/V0'')^2 (67 + 68 alpha) + (1/288)(V0''' V0'''''/V0''^2)(19 + 28 alpha)
             - (1/288)(V0''''''/V0'')(5 + 4 alpha) ]
Classificacao: RESULTADO NUMERICO (aproximacao WKB, precisao ~1% para n = 0, l >= 2).

Kerr, ajuste de Berti-Cardoso-Will (2006, Tab. VIII), modo l = m = 2, n = 0:
    M omega_R = f1 + f2 (1 - chi)^f3,   Q = q1 + q2 (1 - chi)^q3,   omega_I = omega_R / (2 Q)
    f = (1.5251, -1.1568, 0.1292),  q = (0.7000, 1.4187, -0.4990).
Classificacao: MODELO/AJUSTE DA LITERATURA (nao calculado aqui).

Valores de referencia (Leaver 1985) para validacao, M = 1, n = 0:
    gravitacional l=2: 0.37367 - 0.08896 i;  escalar l=2: 0.48364 - 0.09676 i;  EM l=2: 0.45760 - 0.09500 i.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..stability.perturbations import StaticMetric, master_potential, tortoise_grid

LEAVER_SCHWARZSCHILD = {(2, 2, 0): 0.37367 - 0.08896j, (0, 2, 0): 0.48364 - 0.09676j, (1, 2, 0): 0.45760 - 0.09500j,
                        (2, 3, 0): 0.59944 - 0.09270j}
BCW_FIT = dict(f=(1.5251, -1.1568, 0.1292), q=(0.7000, 1.4187, -0.4990))


@dataclass
class QNM:
    metric: str
    spin: int
    ell: int
    n: int
    omega: complex          # unidades 1/M
    method: str
    note: str = ""

    @property
    def frequency(self) -> float:
        return float(self.omega.real)

    @property
    def damping_rate(self) -> float:
        return float(-self.omega.imag)

    @property
    def damping_time_M(self) -> float:
        return 1.0 / self.damping_rate

    @property
    def quality(self) -> float:
        return self.frequency / (2 * self.damping_rate)


def wkb3(metric: StaticMetric, spin=0, ell=2, n=0, r_range=None, n_grid=6000) -> QNM:
    if r_range is None:
        if metric.horizons:
            rh = max(metric.horizons)
            r_range = (rh * (1 + 1e-4), rh + 40 * metric.M)
        else:
            r_range = (-40 * metric.M, 40 * metric.M)
    r, rs = tortoise_grid(metric, r_range[0], r_range[1], n_grid)
    V = master_potential(metric, r, spin, ell)
    # derivadas ate 6a ordem no pico: ajuste polinomial (grau 10) de V(r_*) numa janela uniforme em r_*
    # em torno do maximo (splines dao derivadas altas ruidosas; o polinomio e suave e analitico)
    i0 = int(np.argmax(V))
    half = 5.0 * metric.M
    xs = np.linspace(rs[i0] - half, rs[i0] + half, 801)
    Vx = np.interp(xs, rs, V)
    P = np.polynomial.Polynomial.fit(xs, Vx, 10)
    dP = P.deriv(1)
    roots = dP.roots()
    real = roots[np.abs(roots.imag) < 1e-8].real
    real = real[(real > xs[0]) & (real < xs[-1])]
    x0 = float(real[np.argmax(P(real))]) if real.size else float(xs[int(np.argmax(Vx))])
    V0 = float(P(x0))
    V2, V3, V4, V5, V6 = (float(P.deriv(k)(x0)) for k in range(2, 7))
    nu = n + 0.5
    al = nu**2
    s2 = np.sqrt(-2 * V2)
    Lam = (1 / s2) * ((1 / 8) * (V4 / V2) * (1 / 4 + al) - (1 / 288) * (V3 / V2) ** 2 * (7 + 60 * al))
    Om = (1 / (-2 * V2)) * ((5 / 6912) * (V3 / V2) ** 4 * (77 + 188 * al)
                            - (1 / 384) * (V3**2 * V4 / V2**3) * (51 + 100 * al)
                            + (1 / 2304) * (V4 / V2) ** 2 * (67 + 68 * al)
                            + (1 / 288) * (V3 * V5 / V2**2) * (19 + 28 * al)
                            - (1 / 288) * (V6 / V2) * (5 + 4 * al))
    om2 = (V0 + s2 * Lam) - 1j * nu * s2 * (1 + Om)
    om = np.sqrt(om2)
    if om.imag > 0:
        om = -om
    return QNM(metric.name, spin, ell, n, complex(om), "WKB 3a ordem (Iyer-Will 1987)")


def kerr_bcw(chi: float, M: float = 1.0) -> QNM:
    """Modo l = m = 2, n = 0 de Kerr pelo ajuste de Berti-Cardoso-Will (2006)."""
    f1, f2, f3 = BCW_FIT["f"]
    q1, q2, q3 = BCW_FIT["q"]
    wr = (f1 + f2 * (1 - chi) ** f3) / M
    Q = q1 + q2 * (1 - chi) ** q3
    wi = wr / (2 * Q)
    return QNM(f"Kerr(chi={chi:g})", 2, 2, 0, complex(wr, -wi), "ajuste BCW 2006 (literatura)",
               note="ajuste publicado a QNMs de Kerr; nao calculado neste codigo")
