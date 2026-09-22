"""
Einstein–Cartan–Sciama–Kibble (ECSK) com fluido de spin: a fisica que entra na
equacao de Friedmann modificada.

Classificacao: MODELO TEORICO EXISTENTE NA LITERATURA (nao e resultado deste
codigo).  Equacoes usadas (unidades de Planck, G = c = hbar = 1):

    H^2 + k/a^2   = (8 pi/3) (rho - alpha n^2)              [Poplawski 2010, eq. (5)-(6)]
    a''/a         = -(4 pi/3) (rho + 3 p - 4 alpha n^2)
    rho_eff = rho - alpha n^2,  p_eff = p - alpha n^2,  alpha = kappa / 32 = pi / 4  (spin 1/2)

onde n e a densidade numerica dos fermions (n ~ a^-3) e o termo -alpha n^2 vem da
media do quadrado da densidade de spin, s^2 = n^2 / 8, com kappa = 8 pi.  A
conservacao vale separadamente: d(alpha n^2)/dt = -6 H alpha n^2 e
3 H (rho_eff + p_eff) reproduz esse termo.

Referencias
-----------
Poplawski, N. J. (2010). Phys. Lett. B 694, 181.  Poplawski, N. J. (2012). Phys. Rev. D 85, 107502.
Hehl, F. W., von der Heyde, P., Kerlick, G. D. & Nester, J. M. (1976). Rev. Mod. Phys. 48, 393.
Gasperini, M. (1986). Phys. Rev. Lett. 56, 2873 (spin-dominated bounce).

Hipoteses/limitacoes (documentadas):
* Fluido de spin de Weyssenhoff com spins nao polarizados (media <s^2> = n^2/8).
* Torcao tratada como contribuicao efetiva a materia em RG; condicoes de juncao
  proprias de EC podem adicionar termos de spin na superficie (nao incluidos).
* Poeira nao relativistica ate o ricochete: a densidade de ricochete rho_b = 4 m^2/pi
  para neutrons (~4e58 kg/m^3) esta muito acima da densidade nuclear, onde a equacao
  de estado real nao e a de poeira.  O valor deve ser lido como ordem de grandeza.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.units import M_NEUTRON, M_PLANCK, RHO_PLANCK

ALPHA_SPIN_HALF = np.pi / 4.0  # kappa/32 em unidades de Planck


def spin_torsion_coefficient(kappa: float = 8 * np.pi) -> float:
    """alpha = kappa / 32 (spin 1/2, spins nao polarizados)."""
    return kappa / 32.0


def bounce_density_planck(m_fermion_kg: float = M_NEUTRON) -> float:
    """Densidade de ricochete para poeira de fermions de massa m: rho_b = 4 m^2 / pi (Planck).

    Derivacao: rho_eff = 0  <=>  rho = alpha n^2 = alpha (rho/m)^2  <=>  rho = m^2/alpha = 4 m^2/pi."""
    m = m_fermion_kg / M_PLANCK
    return 4 * m**2 / np.pi


def bounce_density_si(m_fermion_kg: float = M_NEUTRON) -> float:
    return bounce_density_planck(m_fermion_kg) * RHO_PLANCK


def sigma_from_number_density(n0: float, alpha: float = ALPHA_SPIN_HALF) -> float:
    """Termo de torcao em a = 1:  sigma = alpha n0^2  (entra como -sigma a^-6 em rho_eff)."""
    return alpha * n0**2


def sigma_from_matter(rho_m0: float, m_fermion_planck: float, alpha: float = ALPHA_SPIN_HALF) -> float:
    """sigma a partir da densidade de massa rho_m0 (a=1) e da massa por particula m (Planck)."""
    return alpha * (rho_m0 / m_fermion_planck) ** 2


@dataclass(frozen=True)
class EinsteinCartanCorrection:
    """Estrategia de densidade efetiva para o fundo de Friedmann.

    rho_eff = rho - sigma a^-6,  p_eff = p - sigma a^-6.
    """

    sigma: float
    name: str = "Einstein-Cartan (spin fluid)"

    def rho_eff(self, a, rho, p):
        return rho - self.sigma * a**-6

    def p_eff(self, a, rho, p):
        return p - self.sigma * a**-6

    def H2_contribution(self, a, rho, p):
        return (8 * np.pi / 3) * self.rho_eff(a, rho, p)

    def accel(self, a, rho, p):
        """a''/a."""
        return -(4 * np.pi / 3) * (self.rho_eff(a, rho, p) + 3 * self.p_eff(a, rho, p))

    def classical_limit(self) -> "EinsteinCartanCorrection":
        return EinsteinCartanCorrection(0.0, name="GR (sigma = 0)")
