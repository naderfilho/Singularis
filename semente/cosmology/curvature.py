"""
Curvatura observada -> tamanho da esfera-3 -> massa minima do buraco negro pai; razao
tensor/escalar do ricochete de materia minimo.  RESULTADO COM CONEXAO OBSERVACIONAL
(depende dos dados em observations/constraints.py).
"""
from __future__ import annotations

import numpy as np

from ..core.units import C_SI, G_SI, GYR, M_SUN, MPC
from ..observations.constraints import PLANCK


# ----------------------------------------------------------------------------
# 2-4. Curvatura, massa do pai, idade
# ----------------------------------------------------------------------------
def curvature_constraints(Omega_k=None, Omega_k_err=None, n_sigma=2.0):
    """Do limite observado em Omega_k (fechado: Omega_k < 0) ao piso da massa do pai.

    Universo fechado de poeira:  a_m = (8 pi G / 3 c^2) rho_m0 a0^3,  a0 = c / (H0 sqrt(-Omega_k)).
    O universo observavel (raio comovel D) ocupa chi_obs = D / a0 da esfera-3; a estrela-pai
    cobria chi0 >= chi_obs, e  M_pai = (a_m / 2) sin^3(chi0) c^2/G  >=  (a_m/2) sin^3(chi_obs) c^2/G.
    """
    Ok = PLANCK["Omega_k"] if Omega_k is None else Omega_k
    Oke = PLANCK["Omega_k_err"] if Omega_k_err is None else Omega_k_err
    H0 = PLANCK["H0"] * 1e3 / MPC
    rho_c = 3 * H0**2 / (8 * np.pi * G_SI)
    rho_m0 = PLANCK["Omega_m"] * rho_c
    D_obs = 4.4e26  # m, raio comovel do universo observavel (~46.5 Gly)
    # o valor mais negativo permitido a n_sigma:
    Ok_min = Ok - n_sigma * Oke
    if Ok_min >= 0:
        return dict(compatible_closed=False)
    a0 = C_SI / (H0 * np.sqrt(-Ok_min))  # raio de curvatura MINIMO permitido
    a_m = (8 * np.pi * G_SI / (3 * C_SI**2)) * rho_m0 * a0**3  # em metros (raio maximo do universo fechado)
    chi_obs = D_obs / a0
    M_min_kg = (a_m / 2) * np.sin(chi_obs) ** 3 * C_SI**2 / G_SI
    t_recollapse = np.pi * a_m / 2 / C_SI  # s (poeira pura, sem Lambda)
    # massa da esfera-3 inteira (chi0 = pi/2, se a estrela cobrisse meio universo) para referencia
    M_half = (a_m / 2) * C_SI**2 / G_SI
    return dict(compatible_closed=True, Omega_k_used=Ok_min, a0_min_m=a0, a_m_m=a_m, chi_obs=chi_obs,
                M_parent_min_kg=M_min_kg, M_parent_min_Msol=M_min_kg / M_SUN,
                M_parent_halfsphere_Msol=M_half / M_SUN,
                t_recollapse_Gyr=t_recollapse / GYR, age_Gyr=PLANCK["age_Gyr"])


def matter_bounce_tensor_ratio():
    """Ricochete de materia minimo: r = 16 epsilon, epsilon = 3/2 (poeira)  =>  r = 24."""
    return 16 * 1.5


