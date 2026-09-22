"""
O ricochete: como o universo em contracao dentro do buraco negro vira um Big Bang.

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
from scipy.integrate import solve_ivp

# --- constantes fisicas (SI) para converter ------------------------------------------------
G_SI = 6.67430e-11
C_SI = 2.99792458e8
HBAR_SI = 1.054571817e-34
M_SUN = 1.98892e30
L_PLANCK = np.sqrt(HBAR_SI * G_SI / C_SI**3)  # 1.616e-35 m
T_PLANCK = L_PLANCK / C_SI
M_PLANCK = np.sqrt(HBAR_SI * C_SI / G_SI)  # 2.176e-8 kg
RHO_PLANCK = M_PLANCK / L_PLANCK**3  # 5.16e96 kg/m^3
M_NEUTRON = 1.67492749804e-27  # kg


@dataclass
class TorsionCosmology:
    """Friedmann com torcao (ECSK) para poeira + radiacao + termo de spin sigma a^{-6}."""

    rho_m0: float = 1.0  # densidade de materia em a=1
    rho_r0: float = 0.0  # radiacao em a=1
    sigma: float = 1e-3  # alpha n0^2: termo de torcao em a=1
    k: float = 0.0  # curvatura (+1 fechado, como o interior de Oppenheimer-Snyder)

    def rho(self, a):
        return self.rho_m0 * a**-3 + self.rho_r0 * a**-4

    def p(self, a):
        return self.rho_r0 * a**-4 / 3

    def rho_eff(self, a):
        return self.rho(a) - self.sigma * a**-6

    def H2(self, a):
        return (8 * np.pi / 3) * self.rho_eff(a) - self.k / a**2

    def a_bounce(self):
        """Raiz de rho_eff(a) a^2 = 3k/(8pi): resolvida numericamente (para k=0 e so poeira: a_b = (sigma/rho_m0)^(1/3))."""
        from scipy.optimize import brentq
        f = lambda a: self.H2(a)  # noqa: E731
        lo = 1e-9
        hi = 1.0
        # H2 e negativo em a muito pequeno (torcao) e positivo depois
        while f(hi) < 0 and hi < 1e6:
            hi *= 2
        return brentq(f, lo, hi)

    def rhs(self, t, y):
        a, adot = y
        acc = -(4 * np.pi / 3) * (self.rho(a) + 3 * self.p(a) - 4 * self.sigma * a**-6) * a
        return [adot, acc]

    def solve(self, a0=1.0, t_span=(0.0, 6.0), contracting=True, n=4000, **kw):
        """Integra a(t) comecando em a0 na fase de contracao (como o interior do colapso)."""
        H2 = self.H2(a0)
        if H2 < 0:
            raise ValueError("a0 esta abaixo do ricochete.")
        adot0 = -np.sqrt(H2) * a0 if contracting else np.sqrt(H2) * a0
        sol = solve_ivp(self.rhs, t_span, [a0, adot0], method="DOP853", rtol=1e-11, atol=1e-13,
                        dense_output=True, **kw)
        t = np.linspace(t_span[0], sol.t[-1], n)
        a, adot = sol.sol(t)
        return dict(t=t, a=a, adot=adot, H=adot / a, rho=self.rho(a), rho_eff=self.rho_eff(a),
                    a_bounce=self.a_bounce())


@dataclass
class LQCCosmology:
    """Dinamica efetiva da Cosmologia Quantica de Laco (k=0, fluido w = p/rho constante)."""

    rho0: float = 1.0
    w: float = 0.0
    rho_c: float = 0.41  # em densidades de Planck

    def rho(self, a):
        return self.rho0 * a ** (-3 * (1 + self.w))

    def H2(self, a):
        r = self.rho(a)
        return (8 * np.pi / 3) * r * (1 - r / self.rho_c)

    def a_bounce(self):
        return (self.rho0 / self.rho_c) ** (1.0 / (3 * (1 + self.w)))

    def rhs(self, t, y):
        a, H = y
        r = self.rho(a)
        Hdot = -4 * np.pi * (1 + self.w) * r * (1 - 2 * r / self.rho_c)
        return [a * H, Hdot]

    def solve(self, a0=1.0, t_span=(0.0, 6.0), n=4000):
        if self.H2(a0) < 0:
            raise ValueError("rho(a0) > rho_c: a0 ja esta alem do ricochete; use rho0 menor.")
        H0 = -np.sqrt(self.H2(a0))
        sol = solve_ivp(self.rhs, t_span, [a0, H0], method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
        t = np.linspace(t_span[0], sol.t[-1], n)
        a, H = sol.sol(t)
        return dict(t=t, a=a, H=H, rho=self.rho(a), a_bounce=self.a_bounce())


def torsion_dust_exact(t, rho_m0=1.0, sigma=1e-3, a0=1.0):
    """Solucao EXATA de poeira + torcao (k=0): na variavel de volume w = a^3,
        w'^2 = 24 pi (rho_m0 w - sigma)   =>   w'' = 12 pi rho_m0   (uma parabola!)
        a(t)^3 = sigma/rho_m0 + 6 pi rho_m0 (t - t_b)^2,
    com t_b fixado por a(0) = a0 na fase de contracao.  Vale para qualquer sigma > 0."""
    t = np.asarray(t, float)
    w0 = a0**3
    tb = np.sqrt((w0 - sigma / rho_m0) / (6 * np.pi * rho_m0))
    w = sigma / rho_m0 + 6 * np.pi * rho_m0 * (t - tb) ** 2
    return w ** (1 / 3), tb


def classical_collapse_reference(rho_m0=1.0, a0=1.0, n=4000):
    """Poeira em GR pura (sem torcao), k=0: a(t) = a0 (1 - t/t_s)^{2/3} atinge a=0 em t_s = 2/(3 H0)."""
    H0 = np.sqrt(8 * np.pi * rho_m0 / 3) * a0 ** (-1.5)
    ts = 2 / (3 * H0)
    t = np.linspace(0, ts * 0.9999, n)
    return dict(t=t, a=a0 * (1 - t / ts) ** (2 / 3), t_singularity=ts)


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

    def __post_init__(self):
        self.M_pl = self.M_kg / M_PLANCK  # massa em unidades de Planck (G=c=hbar=1 => comprimento = massa)
        self.m_pl = self.m_fermion_kg / M_PLANCK
        self.rs_pl = 2 * self.M_pl
        self.rho_bounce_pl = 4 * self.m_pl**2 / np.pi
        # interior OS: rho a^3 = const  =>  rho(a) = 3 M / (4 pi R^3) com R = a sin(chi0)
        # bounce quando rho = rho_b  =>  R_b^3 = 3 M / (4 pi rho_b)
        self.R_bounce_pl = (3 * self.M_pl / (4 * np.pi * self.rho_bounce_pl)) ** (1 / 3)
        R0 = self.R0_over_rs * self.rs_pl
        self.chi0 = float(np.arcsin(np.sqrt(2 * self.M_pl / R0)))
        self.a_m_pl = R0 / np.sin(self.chi0)
        self.a_bounce_pl = self.R_bounce_pl / np.sin(self.chi0)
        # numero de particulas dentro do buraco negro
        self.N_particles = self.M_kg / self.m_fermion_kg
        # tempo proprio do horizonte ate o ricochete ~ pi M (classico ate r~0; o ricochete ocorre em R_b << 2M)
        self.tau_horizon_to_bounce_s = np.pi * self.M_pl * T_PLANCK * (1 - (self.R_bounce_pl / self.rs_pl) ** 1.5)
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
