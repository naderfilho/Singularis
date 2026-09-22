"""
Limite entropico e genealogia de universos.

HIPOTESE ESPECULATIVA explicitamente marcada: que a entropia de Bekenstein-Hawking do pai
(4 pi M^2) limite a entropia total do universo-filho.  Os numeros sao consequencias dessa
hipotese, nao evidencia dela.  Ver docs/02_fronteira.md (secao B).
"""
from __future__ import annotations

import numpy as np

from ..core.units import M_PLANCK, M_SUN


# ----------------------------------------------------------------------------
# B. Limite entropico e profundidade da arvore
# ----------------------------------------------------------------------------
K_B_ENTROPIES = {
    # estimativas de Egan & Lineweaver (2010), ApJ 710, 1825, em unidades de k_B
    "CMB (fotons)": 2.03e89,
    "neutrinos cosmicos": 5.2e89,
    "buracos negros supermassivos": 3.1e104,
    "buracos negros estelares": 5.9e97,
}


def parent_mass_lower_bound(S_child_kB):
    """S_BH = 4 pi M^2 / m_Pl^2  >=  S_child   =>   M >= m_Pl sqrt(S_child / 4 pi).  Retorna kg e M_sol."""
    M = M_PLANCK * np.sqrt(S_child_kB / (4 * np.pi))
    return M, M / M_SUN


def genealogy_depth(M0_kg, N_per_generation, m_min_kg=M_PLANCK):
    """Sob o limite holografico, um filho so cabe se  sum_i m_i^2 <= M_pai^2.
    Com N buracos negros iguais por geracao: m_{k+1} = m_k / sqrt(N).
    A linhagem acaba quando m_k < m_min:  k_max = ln(M0/m_min) / ln(sqrt(N))."""
    return np.log(M0_kg / m_min_kg) / np.log(np.sqrt(N_per_generation))


