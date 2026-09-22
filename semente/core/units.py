"""
Unidades e constantes.  Fonte unica para todo o pacote.

Convencoes
----------
* Modulos geometricos trabalham em unidades geometricas G = c = 1, com a massa M
  do buraco negro como escala de comprimento e tempo.
* Modulos cosmologicos/quanticos trabalham em unidades de Planck G = c = hbar = 1.
* A classe `Scale` converte entre geometricas (em unidades de M), Planck e SI.

Valores CODATA 2018 / IAU 2015.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --- SI ---------------------------------------------------------------------------
G_SI = 6.67430e-11          # m^3 kg^-1 s^-2
C_SI = 2.99792458e8         # m/s
HBAR_SI = 1.054571817e-34   # J s
K_B_SI = 1.380649e-23       # J/K
M_SUN = 1.98892e30          # kg
M_NEUTRON = 1.67492749804e-27  # kg
MPC = 3.0857e22             # m
GYR = 3.15576e16            # s
PC = 3.0857e16              # m

# --- Planck -----------------------------------------------------------------------
L_PLANCK = float(np.sqrt(HBAR_SI * G_SI / C_SI**3))   # 1.616e-35 m
T_PLANCK = L_PLANCK / C_SI                             # 5.391e-44 s
M_PLANCK = float(np.sqrt(HBAR_SI * C_SI / G_SI))       # 2.176e-8 kg
RHO_PLANCK = M_PLANCK / L_PLANCK**3                    # 5.155e96 kg/m^3
E_PLANCK = M_PLANCK * C_SI**2                          # J
TEMP_PLANCK = E_PLANCK / K_B_SI                        # K


@dataclass(frozen=True)
class Scale:
    """Escala geometrica definida por uma massa M (em kg).

    length_unit  = G M / c^2  (metros por unidade de comprimento geometrica)
    time_unit    = G M / c^3  (segundos)
    density_unit = c^6 / (G^3 M^2)  (kg/m^3 por unidade geometrica de densidade)
    """

    M_kg: float

    @classmethod
    def from_solar_masses(cls, m: float) -> "Scale":
        return cls(m * M_SUN)

    @property
    def length_unit(self) -> float:
        return G_SI * self.M_kg / C_SI**2

    @property
    def time_unit(self) -> float:
        return G_SI * self.M_kg / C_SI**3

    @property
    def density_unit(self) -> float:
        return C_SI**6 / (G_SI**3 * self.M_kg**2)

    @property
    def M_planck(self) -> float:
        """Massa em unidades de Planck (= comprimento em l_Pl, ja que G=c=hbar=1)."""
        return self.M_kg / M_PLANCK

    def length_to_si(self, x_geo: float) -> float:
        return x_geo * self.length_unit

    def time_to_si(self, t_geo: float) -> float:
        return t_geo * self.time_unit

    def density_to_si(self, rho_geo: float) -> float:
        return rho_geo * self.density_unit

    def frequency_to_hz(self, omega_geo: float) -> float:
        """omega em unidades 1/M -> frequencia em Hz: f = omega / (2 pi t_unit)."""
        return omega_geo / (2 * np.pi * self.time_unit)

    def schwarzschild_radius_m(self) -> float:
        return 2 * self.length_unit


def planck_density_to_si(rho_pl: float) -> float:
    return rho_pl * RHO_PLANCK


def si_mass_to_planck(m_kg: float) -> float:
    return m_kg / M_PLANCK


def kelvin_to_planck(T_K: float) -> float:
    return T_K / TEMP_PLANCK
