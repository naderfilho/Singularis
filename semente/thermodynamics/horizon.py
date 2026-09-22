"""
Termodinamica de horizontes.

RESULTADOS MATEMATICOS (padrao):
   area A = 4 pi R(r_h)^2,  entropia S = A/4 (Bekenstein-Hawking, G = hbar = c = k_B = 1),
   temperatura T = kappa / 2 pi,  kappa = |f'(r_h)| / 2  para  -f dt^2 + dr^2/f + R^2 dOmega^2.
   Schwarzschild: T = 1/(8 pi M), S = 4 pi M^2.   Kerr: S = pi (r_+^2 + a^2), primeira lei dM = T dS + Omega_H dJ.

RESULTADOS DESTE CODIGO:
   * Simpson-Visser: R(r_h) = 2M para todo l < 2M  =>  S = 4 pi M^2 (igual a Schwarzschild), mas
     T = sqrt(4M^2 - l^2)/(16 pi M^2) < T_Schw.  Com S = A/4 a primeira lei dM = T dS FALHA por um fator
     sqrt(1 - l^2/4M^2): a diferenca e um termo de trabalho associado a l (ou a entropia nao e A/4).
     Reportamos o residuo, sem escolher uma das duas interpretacoes.
   * Horizontes no colapso dinamico.  O teorema da area (Hawking 1971) vale para o horizonte de EVENTOS:
     S_EH(tau) = pi R_EH^2 cresce monotonicamente no colapso classico (verificado).  O horizonte APARENTE
     de Oppenheimer-Snyder nasce na superficie quando R_s = 2M e se propaga para DENTRO (R_AH = 2 m(chi_AH)
     decresce): sua area nao obedece ao teorema, e isso e classico e conhecido.  No ricochete a regiao presa
     e transiente e S_AH volta a zero, e o horizonte de eventos classico deixa de estar definido (docs/03).

Referencias: Bekenstein (1973) PRD 7, 2333; Hawking (1975) CMP 43, 199; Bardeen, Carter & Hawking (1973)
CMP 31, 161; Hayward (1994) PRD 49, 6467 (horizontes aprisionados dinamicos).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.units import HBAR_SI, C_SI, G_SI, K_B_SI, M_SUN, Scale
from ..geometry.kerr import Kerr
from ..geometry.static import StaticSphericalMetric


@dataclass
class HorizonThermo:
    metric: str
    r_h: float
    area: float
    entropy: float
    temperature: float
    surface_gravity: float
    first_law_residual: float | None = None   # 1 - T dS/dM (a l fixo) quando calculavel
    note: str = ""

    def temperature_kelvin(self, M_solar: float) -> float:
        """T em kelvin para uma massa M_solar (T_geo em unidades 1/M)."""
        # T = hbar kappa_SI / (2 pi k_B), com kappa_SI = kappa_geo / time_unit; self.temperature = kappa_geo / 2 pi
        s = Scale.from_solar_masses(M_solar)
        return self.temperature * HBAR_SI / (K_B_SI * s.time_unit)


def static_thermo(metric: StaticSphericalMetric, dM=1e-5) -> HorizonThermo:
    hs = [h for h in metric.horizons() if h >= 0]
    if not hs:
        return HorizonThermo(metric.name, float("nan"), 0.0, 0.0, 0.0, 0.0, None, "sem horizonte")
    rh = hs[-1]
    Rh = float(metric.R(rh))
    A = 4 * np.pi * Rh**2
    S = A / 4
    kappa = metric.surface_gravity(rh)
    T = kappa / (2 * np.pi)
    res = None
    if "M" in metric.params:
        # dS/dM numerico a parametros restantes fixos
        p2 = dict(metric.params)
        p2["M"] = metric.params["M"] + dM
        m2 = StaticSphericalMetric(metric.f_expr, metric.R_expr, metric.param_names, p2, metric.name, metric.r_domain)
        h2 = [h for h in m2.horizons() if h >= 0]
        if h2:
            S2 = np.pi * float(m2.R(h2[-1])) ** 2
            res = float(1 - T * (S2 - S) / dM)
    return HorizonThermo(metric.name, rh, A, S, T, kappa, res)


def kerr_thermo(k: Kerr) -> dict:
    if not k.kerr_horizons:
        return dict(note="sem horizonte")
    S = k.horizon_area / 4
    T = k.surface_gravity / (2 * np.pi)
    J = k.a * k.M
    # primeira lei: dM = T dS + Omega_H dJ (verificacao numerica a J fixo e a M fixo)
    dM = 1e-5
    k2 = Kerr(k.M + dM, J / (k.M + dM))     # J fixo
    dS = k2.horizon_area / 4 - S
    res = 1 - T * dS / dM
    return dict(entropy=S, temperature=T, Omega_H=k.Omega_H, J=J, first_law_residual_fixed_J=float(res))


def hawking_temperature_kelvin(M_solar: float) -> float:
    """T_H = hbar c^3 / (8 pi G M k_B)."""
    return HBAR_SI * C_SI**3 / (8 * np.pi * G_SI * M_solar * M_SUN * K_B_SI)


def event_horizon_entropy(res) -> dict:
    """S_EH(tau) = pi R_EH^2 ao longo do gerador do horizonte de eventos (raio nulo de saida lancado do
    centro em tau_birth_center, do CollapseResult).  Verifica monotonicidade (teorema da area)."""
    info = res.event_horizon
    tb = info.get("tau_birth_center")
    if tb is None:
        return dict(defined=False, note=info.get("note", ""))
    tau, chi = res.tau, res.chi
    E = -res.m_MS[0] / res.R[0]   # perfil em repouso: E = -m/R0
    c = chi[0]
    i0 = int(np.searchsorted(tau, tb))
    taus, Rs = [tau[i0]], [float(np.interp(c, chi, res.R[i0]))]
    for i in range(i0, tau.size - 1):
        Rp = np.interp(c, chi, res.R_prime[i])
        if Rp <= 0:
            break
        c = c + np.sqrt(max(1 + 2 * np.interp(c, chi, E), 0.0)) / Rp * (tau[i + 1] - tau[i])
        if c >= chi[-1]:
            break
        taus.append(tau[i + 1])
        Rs.append(float(np.interp(c, chi, res.R[i + 1])))
    taus, Rs = np.array(taus), np.array(Rs)
    S = np.pi * Rs**2
    dS = np.diff(S)
    return dict(defined=True, tau=taus, R=Rs, S=S, monotonic=bool(np.all(dS >= -1e-9 * S.max())),
                S_final=float(S[-1]), S_final_over_4piM2=float(S[-1] / (4 * np.pi * res.M**2)))


def apparent_horizon_entropy(res) -> dict:
    """S_AH(tau) = pi R_AH^2 de um CollapseResult (camada presa mais externa).  NAO se espera monotonicidade:
    o horizonte aparente de OS se propaga para dentro.  Serve para ver a regiao presa transiente do ricochete."""
    R = res.apparent_horizon_R
    S = np.where(np.isfinite(R), np.pi * R**2, 0.0)
    dS = np.diff(S)
    m = np.isfinite(dS)
    tol = 1e-6 * max(S.max(), 1e-300)
    decreases = np.where(dS[m] < -tol)[0]
    return dict(tau=res.tau, S=S, S_max=float(S.max()), monotonic=bool(decreases.size == 0),
                first_decrease_tau=float(res.tau[1:][m][decreases[0]]) if decreases.size else None,
                final_S=float(S[-1]))


def baby_universe_budget(M_solar: float, m_particle_kg: float) -> dict:
    """Graus de liberdade: entropia do horizonte do pai vs entropia (ordem N) da materia que vira o filho."""
    s = Scale.from_solar_masses(M_solar)
    S_BH = 4 * np.pi * s.M_planck**2
    N = s.M_kg / m_particle_kg
    return dict(S_BH_kB=S_BH, N_particles=N, S_matter_order_N=N, ratio_matter_over_horizon=N / S_BH)
