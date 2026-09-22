"""
Estabilidade dinamica do ricochete homogeneo.

1) Perturbacao homogenea e isotropica do fator de escala, a -> a + delta a.
   RESULTADO MATEMATICO: para um fundo FRW com fluido barotropico, a perturbacao homogenea de a(t)
   satisfazendo o vinculo e uma translacao temporal (delta a = adot * delta t): e GAUGE, nao fisica.
   Por isso NAO usamos "crescimento de delta a" como criterio.  A funcao `homogeneous_mode_is_gauge`
   verifica numericamente que delta a(t) = adot(t) delta t resolve a equacao linearizada.

2) Anisotropia (cisalhamento de Bianchi I).  Numa contracao, o cisalhamento cresce como
   sigma^2 ~ a^-6, mais rapido que qualquer fluido com w < 1: e a instabilidade BKL do ricochete
   isotropico.  Criterio (documentado): o ricochete isotropico e robusto se, no instante do ricochete,
       Sigma_b = sigma_b^2 / (8 pi rho_eff-scale)  =  sigma_0^2 a_b^-6 / rho_*   << 1 .
   Aqui rho_* e a densidade que limita o ricochete (rho_c em LQC; rho_b em EC).  Em LQC efetiva de
   Bianchi I (Gupt & Singh 2012) o cisalhamento tambem e limitado, sigma^2 <= sigma_max^2 ~ 10 rho_c/...
   mas no modelo isotropico usado aqui a anisotropia entra como fluido rigido (w = 1) NAO regularizado:
   esse e um LIMITE do modelo, e por isso o criterio e conservador.  Reportamos:
     * Sigma_b para uma anisotropia inicial dada;
     * a anisotropia inicial maxima  sigma_0^2 < rho_* a_b^6  para que o ricochete isotropico sobreviva;
     * a razao entre o ganho de anisotropia e o de densidade entre o inicio e o ricochete: (a_0/a_b)^3.

3) Perturbacoes tensoriais (ondas gravitacionais) atravessam o ricochete com amplitude finita — ver
   cosmology/perturbations.py; a ausencia de crescimento exponencial ali e o teste de estabilidade
   linear inomogenea.

Referencias: Belinskii, Khalatnikov & Lifshitz (1970) Adv. Phys. 19, 525; Gupt & Singh (2012) PRD 85,
044011; Erickson, Wesley, Steinhardt & Turok (2004) PRD 69, 063514 (crescimento de anisotropias em
contracao).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..cosmology.friedmann import BackgroundSolution, FriedmannModel


def homogeneous_mode_is_gauge(model: FriedmannModel, sol: BackgroundSolution, delta_t: float = 1e-3) -> float:
    """Residuo maximo (relativo) da equacao linearizada  d2(delta a)/dt2 = d(A a)/da * delta a
    para delta a = adot * delta_t.  Deve ser ~0 (a translacao temporal e solucao)."""
    a = sol.a
    da = sol.adot * delta_t
    h = 1e-6 * a
    dAa = (model.accel(a + h) * (a + h) - model.accel(a - h) * (a - h)) / (2 * h)
    lhs = np.gradient(np.gradient(da, sol.t), sol.t)
    rhs = dAa * da
    scale = np.max(np.abs(rhs)) + 1e-300
    inner = slice(5, -5)
    return float(np.max(np.abs(lhs[inner] - rhs[inner])) / scale)


@dataclass
class AnisotropyReport:
    a_0: float
    a_b: float
    rho_star: float
    sigma0_sq: float
    Sigma_bounce: float          # sigma^2(a_b) / rho_*
    sigma0_sq_max: float         # anisotropia inicial maxima para Sigma_bounce < 1
    amplification: float         # (a_0/a_b)^3 : quanto sigma^2/rho cresce ate o ricochete
    verdict: str


def anisotropy_robustness(sol: BackgroundSolution, rho_star: float, sigma0_sq: float = 1e-6) -> AnisotropyReport:
    """Criterio BKL conservador para um ricochete isotropico limitado por rho_*."""
    a0 = float(sol.a[0])
    ab = sol.a_min
    sig_b = sigma0_sq * (a0 / ab) ** 6
    Sigma = sig_b / rho_star
    sig_max = rho_star * ab**6 / a0**6
    amp = (a0 / ab) ** 3
    if Sigma < 1e-2:
        v = "ricochete isotropico robusto para esta anisotropia (Sigma_b << 1)"
    elif Sigma < 1:
        v = "anisotropia comparavel a densidade no ricochete: resultado isotropico nao confiavel (requer Bianchi I)"
    else:
        v = "anisotropia domina no ricochete: o modelo isotropico e invalido aqui (instabilidade BKL)"
    return AnisotropyReport(a0, ab, rho_star, sigma0_sq, Sigma, sig_max, amp, v)
