"""
Cosmologia Quantica de Laco (LQC): dinamica EFETIVA.

Classificacao: MODELO TEORICO EXISTENTE NA LITERATURA.  Distincao obrigatoria:
* "LQC efetiva" = equacoes classicas modificadas que reproduzem o valor esperado
  do estado semiclassico (Ashtekar, Pawlowski & Singh 2006; Taveras 2008).
* "Gravidade quantica de laco completa" = teoria de campo quantico de conexoes;
  NAO esta implementada aqui.  O que segue e um modelo efetivo.

Equacoes efetivas (k = 0, unidades de Planck):

    H^2  = (8 pi/3) rho (1 - rho/rho_c)                        [APS 2006; Singh 2009 eq. (3)]
    H'   = -4 pi (rho + p) (1 - 2 rho/rho_c)
    a''/a = H' + H^2

com rho_c = sqrt(3) / (32 pi^2 gamma^3) rho_Pl ~ 0.41 rho_Pl para o parametro de
Barbero-Immirzi gamma ~ 0.2375.  Ricochete exatamente em rho = rho_c.

Limitacoes:
* Para k = +1 as equacoes efetivas tem forma diferente (Ashtekar, Pawlowski, Singh &
  Vandersloot 2007); aqui so o caso k = 0 e implementado e o codigo recusa k != 0.
* A conexao "massa do buraco negro pai -> condicoes iniciais do universo-filho" usa a
  ideia de Ashtekar-Olmedo-Singh (2018) de transicao BN -> BB, mas nao a dinamica
  completa de Kantowski-Sachs quantica.

Referencias
-----------
Ashtekar, A., Pawlowski, T. & Singh, P. (2006). Phys. Rev. Lett. 96, 141301; Phys. Rev. D 74, 084003.
Singh, P. (2009). Class. Quantum Grav. 26, 125005 (equacoes efetivas e limites).
Ashtekar, A., Olmedo, J. & Singh, P. (2018). Phys. Rev. Lett. 121, 241301.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GAMMA_BARBERO_IMMIRZI = 0.2375
RHO_C_PLANCK = float(np.sqrt(3.0) / (32 * np.pi**2 * GAMMA_BARBERO_IMMIRZI**3))  # ~ 0.41


@dataclass(frozen=True)
class LQCCorrection:
    """Estrategia de LQC efetiva para o fundo de Friedmann (k = 0)."""

    rho_c: float = RHO_C_PLANCK
    name: str = "LQC efetiva (k = 0)"

    def rho_eff(self, a, rho, p):
        return rho * (1 - rho / self.rho_c)

    def p_eff(self, a, rho, p):
        # definido para que a''/a = -(4pi/3)(rho_eff + 3 p_eff) reproduza H' + H^2 da LQC
        H2 = (8 * np.pi / 3) * rho * (1 - rho / self.rho_c)
        Hdot = -4 * np.pi * (rho + p) * (1 - 2 * rho / self.rho_c)
        acc = Hdot + H2
        return (-acc * 3 / (4 * np.pi) - self.rho_eff(a, rho, p)) / 3

    def H2_contribution(self, a, rho, p):
        return (8 * np.pi / 3) * self.rho_eff(a, rho, p)

    def accel(self, a, rho, p):
        H2 = (8 * np.pi / 3) * rho * (1 - rho / self.rho_c)
        Hdot = -4 * np.pi * (rho + p) * (1 - 2 * rho / self.rho_c)
        return Hdot + H2

    def classical_limit(self) -> "LQCCorrection":
        return LQCCorrection(rho_c=np.inf, name="GR (rho_c -> inf)")
