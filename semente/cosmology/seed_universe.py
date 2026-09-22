"""
Universo-filho (modelo SEMENTE) e o argumento de Pathria.

O ricochete em si vive em cosmology/friedmann.py (GR, Einstein-Cartan, LQC).

Duas fisicas concretas (ambas publicadas, ambas resolvidas aqui numericamente):

1) Einstein-Cartan-Sciama-Kibble (ECSK) com fluido de spin  [Poplawski 2010, 2012]
   A torcao do espaco-tempo, gerada pelo spin dos fermions, adiciona um termo
   repulsivo proporcional ao quadrado da densidade de particulas:

       H^2 + k/a^2 = (8 pi/3) [ rho - alpha n^2 ],       n ~ a^{-3}
       a''/a       = -(4 pi/3) [ rho + 3p - 4 alpha n^2 ]

   Em unidades de Planck (G = c = hbar = 1), para fermions de spin 1/2:
   alpha = kappa/32 = pi/4 (s^2 = n^2/8, kappa = 8 pi).  O termo ~ a^{-6}
   domina em a pequeno e forca um ricochete em rho_b = 4 m^2/pi (poeira de massa m).

2) Cosmologia Quantica de Laco (LQC)  [Ashtekar-Pawlowski-Singh 2006]
       H^2 = (8 pi/3) rho (1 - rho/rho_c),   H' = -4 pi (rho + p)(1 - 2 rho/rho_c)
   com rho_c ~ 0.41 rho_Planck.  Ricochete quando rho = rho_c.

3) Modelo SEMENTE (sintese deste projeto): pega o interior FRW EXATO do colapso
   de Oppenheimer-Snyder (semente.collapse), poe torcao nele, e devolve
   as propriedades do universo-filho em funcao da massa do buraco negro pai.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.units import (C_SI, G_SI, HBAR_SI, L_PLANCK, M_NEUTRON, M_PLANCK, M_SUN, RHO_PLANCK, T_PLANCK)  # noqa: F401
from ..quantum.einstein_cartan import bounce_density_planck


# ----------------------------------------------------------------------------
# Modelo SEMENTE: universo-filho em funcao da massa do buraco negro pai
# ----------------------------------------------------------------------------
@dataclass
class SeedUniverse:
    """Propriedades do universo que nasce dentro de um buraco negro de massa M (kg),
    supondo o interior de Oppenheimer-Snyder feito de fermions de massa m e o
    ricochete de torcao de Einstein-Cartan.

    Em unidades de Planck:  rho_b = 4 m^2 / pi.
    """

    M_kg: float
    m_fermion_kg: float = M_NEUTRON
    R0_over_rs: float = 100.0  # raio inicial da estrela em unidades do raio de Schwarzschild
    mechanism: str = "einstein_cartan"   # "einstein_cartan" (rho_b = 4 m^2/pi) ou "lqc" (rho_c ~ 0.41 rho_Pl)

    def __post_init__(self):
        from ..quantum.lqc import RHO_C_PLANCK
        self.M_pl = self.M_kg / M_PLANCK  # massa em unidades de Planck (G=c=hbar=1 => comprimento = massa)
        self.m_pl = self.m_fermion_kg / M_PLANCK
        self.rs_pl = 2 * self.M_pl
        if self.mechanism == "einstein_cartan":
            self.rho_bounce_pl = bounce_density_planck(self.m_fermion_kg)
        elif self.mechanism == "lqc":
            self.rho_bounce_pl = RHO_C_PLANCK
        else:
            raise ValueError("mechanism deve ser 'einstein_cartan' ou 'lqc'")
        # interior OS: rho a^3 = const  =>  rho(a) = 3 M / (4 pi R^3) com R = a sin(chi0)
        # bounce quando rho = rho_b  =>  R_b^3 = 3 M / (4 pi rho_b)
        self.R_bounce_pl = (3 * self.M_pl / (4 * np.pi * self.rho_bounce_pl)) ** (1 / 3)
        R0 = self.R0_over_rs * self.rs_pl
        self.chi0 = float(np.arcsin(np.sqrt(2 * self.M_pl / R0)))
        self.a_m_pl = R0 / np.sin(self.chi0)
        self.a_bounce_pl = self.R_bounce_pl / np.sin(self.chi0)
        # numero de particulas dentro do buraco negro
        self.N_particles = self.M_kg / self.m_fermion_kg
        # tempo proprio EXATO ao longo da superficie de Oppenheimer-Snyder, do horizonte ate R_b:
        #   R = R0 (1 + cos eta)/2,  tau = sqrt(R0^3 / 8M) (eta + sin eta)
        eta_h = np.arccos(4 * self.M_pl / R0 - 1.0)
        eta_b = np.arccos(2 * self.R_bounce_pl / R0 - 1.0)
        self.tau_horizon_to_bounce_s = np.sqrt(R0**3 / (8 * self.M_pl)) * ((eta_b + np.sin(eta_b)) - (eta_h + np.sin(eta_h))) * T_PLANCK
        # densidade do ricochete em SI e comparacao com a nuclear
        self.rho_bounce_SI = self.rho_bounce_pl * RHO_PLANCK
        # "quao maior" o universo-filho pode ficar depois: fator de escala maximo a_m (universo fechado) / a_bounce
        self.expansion_factor = self.a_m_pl / self.a_bounce_pl
        # raio proprio maximo do universo-filho (em metros): a_m * chi0
        self.R_max_child_m = self.a_m_pl * self.chi0 * L_PLANCK
        self.R_bounce_child_m = self.a_bounce_pl * self.chi0 * L_PLANCK

    def table(self) -> dict:
        return {
            "M_pai (M_sol)": self.M_kg / M_SUN,
            "raio de Schwarzschild (m)": self.rs_pl * L_PLANCK,
            "N fermions": self.N_particles,
            "rho_ricochete (kg/m^3)": self.rho_bounce_SI,
            "rho_ricochete / rho_Planck": self.rho_bounce_pl,
            "R_areal no ricochete (m)": self.R_bounce_pl * L_PLANCK,
            "R_areal no ricochete / r_s": self.R_bounce_pl / self.rs_pl,
            "tempo horizonte->ricochete (s)": self.tau_horizon_to_bounce_s,
            "raio proprio do filho no Big Bang (m)": self.R_bounce_child_m,
            "raio proprio maximo do filho (m)": self.R_max_child_m,
            "fator de expansao a_max/a_bounce": self.expansion_factor,
        }


def observable_universe_as_black_hole():
    """Argumento de Pathria (1972) / Good (1972): o raio de Schwarzschild da massa contida na
    esfera de Hubble e igual ao proprio raio de Hubble.

    ATENCAO (honestidade): para um universo plano isso e uma IDENTIDADE da equacao de Friedmann,
    nao uma coincidencia:  rho_c = 3H^2/(8 pi G)  =>  r_s = 2G M_H / c^2 = c/H = R_H  exatamente.
    Por isso a ideia "o universo e o interior de um buraco negro" e sugestiva, mas essa conta
    sozinha nao a demonstra.  O que a sustenta e o restante deste projeto (interior FRW do
    colapso, ricochete, black-bounce)."""
    H0 = 67.4e3 / 3.0857e22  # s^-1
    rho_c = 3 * H0**2 / (8 * np.pi * G_SI)
    R_H = C_SI / H0
    M_H = rho_c * (4 / 3) * np.pi * R_H**3
    rs = 2 * G_SI * M_H / C_SI**2
    return dict(rho_critica_kg_m3=rho_c, R_Hubble_m=R_H, M_Hubble_kg=M_H, M_Hubble_Msol=M_H / M_SUN,
                rs_m=rs, razao_rs_R=rs / R_H, nota="razao = 1 e uma identidade de Friedmann para k=0")
