"""
Juncao de Israel entre um interior FRW (com ou sem ricochete) e um exterior estatico
(Schwarzschild ou black-bounce de Simpson-Visser).

RESULTADO DESTE CODIGO (tratamento efetivo): a juncao com Schwarzschild falha num intervalo
finito em torno do ricochete; na familia black-bounce so l = R_b e compativel.
Ver docs/02_fronteira.md (secao A) para a derivacao e as ressalvas.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from ..core.units import L_PLANCK, M_NEUTRON, M_PLANCK, M_SUN


# ----------------------------------------------------------------------------
# A. Colapso com torcao + juncao de Israel
# ----------------------------------------------------------------------------
@dataclass
class TorsionCollapse:
    """Interior FRW fechado de poeira + termo de torcao s a^-6, ajustado para ricochetear em R_b.

    OS:  sin^2(chi0) = 2M/R0,  a_m = R0/sin(chi0),  rho a^3 = 3 a_m/(8 pi).
    Torcao:  H^2 + 1/a^2 = (8 pi/3)(rho - s a^-6).  Bounce em a_b = R_b/sin(chi0):
             s = (3/8pi) a_b^3 (a_m - a_b).
    """

    M: float = 1.0
    R0: float = 8.0
    Rb_over_R0: float = 0.05

    def __post_init__(self):
        self.chi0 = float(np.arcsin(np.sqrt(2 * self.M / self.R0)))
        self.sin, self.cos = np.sin(self.chi0), np.cos(self.chi0)
        self.a_m = self.R0 / self.sin
        self.R_b = self.Rb_over_R0 * self.R0
        self.a_b = self.R_b / self.sin
        self.rho_c = 3 * self.a_m / (8 * np.pi)  # rho a^3
        self.s = 3 / (8 * np.pi) * self.a_b**3 * (self.a_m - self.a_b)

    def rho(self, a):
        return self.rho_c / a**3

    def rho_eff(self, a):
        return self.rho(a) - self.s / a**6

    def H2(self, a):
        return (8 * np.pi / 3) * self.rho_eff(a) - 1 / a**2

    def rhs(self, tau, y):
        a, adot = y
        acc = -(4 * np.pi / 3) * (self.rho(a) - 4 * self.s / a**6) * a
        return [adot, acc]

    def solve(self, n=6000):
        """Da estrela em repouso (a = a_m, adot = 0-) atraves do ricochete ate voltar a a_m."""
        from scipy.optimize import brentq
        # com torcao o ponto de retorno maximo fica ligeiramente abaixo de a_m: achamos a raiz de H^2
        self.a_max = brentq(self.H2, 1.5 * self.a_b, self.a_m)
        a0 = self.a_max * (1 - 1e-7)
        adot0 = -np.sqrt(max(self.H2(a0), 0.0)) * a0

        def back(tau, y):  # para no maximo de expansao depois do ricochete (adot cruza zero descendo)
            return y[1] if y[0] > 2 * self.a_b else 1.0
        back.terminal = True
        back.direction = -1
        sol = solve_ivp(self.rhs, (0, 1e5), [a0, adot0], method="DOP853", rtol=1e-11, atol=1e-13,
                        dense_output=True, events=back, max_step=self.a_m / 50)
        tau = np.linspace(0, sol.t[-1], n)
        a, adot = sol.sol(tau)
        R = a * self.sin
        Rdot = adot * self.sin
        m_in = (4 * np.pi / 3) * self.rho_eff(a) * R**3  # massa de Misner-Sharp interior na superficie
        ib = int(np.argmin(a))
        return dict(tau=tau, a=a, adot=adot, R=R, Rdot=Rdot, m_in=m_in, tau_b=tau[ib], R_min=R[ib],
                    tau_h=float(np.interp(2 * self.M, R[:ib][::-1], tau[:ib][::-1])) if R[0] > 2 * self.M else 0.0)

    # --- juncoes ---------------------------------------------------------------
    def junction_schwarzschild(self, sol):
        """Camada fina necessaria para colar o interior a um exterior de Schwarzschild(M).

        K^theta_theta interior = cos(chi0)/R  (exato para FRW, usa a eq. de Friedmann);
        K^theta_theta exterior = sqrt(1 + Rdot^2 - 2M/R)/R.
        sigma = [cos(chi0) - sqrt(1 + Rdot^2 - 2M/R)] / (4 pi R).
        Quando 1 + Rdot^2 - 2M/R < 0 NAO existe camada (de nenhuma energia) que faca a juncao.
        """
        R, Rdot = sol["R"], sol["Rdot"]
        Es2 = 1 + Rdot**2 - 2 * self.M / R
        ok = Es2 >= 0
        sigma = np.full_like(R, np.nan)
        sigma[ok] = (self.cos - np.sqrt(Es2[ok])) / (4 * np.pi * R[ok])
        m_shell = 4 * np.pi * R**2 * sigma
        return dict(Es2=Es2, valid=ok, sigma=sigma, m_shell=m_shell, frac_time_invalid=float((~ok).mean()),
                    tau_fail=(sol["tau"][~ok].min(), sol["tau"][~ok].max()) if (~ok).any() else None)

    def junction_black_bounce(self, sol, l):
        """Camada fina necessaria para colar o interior a um exterior de Simpson-Visser(M, l).

        Exterior: R = sqrt(r^2 + l^2), f = 1 - 2M/R.  K^theta_theta = (r/R^2) sqrt(f + rdot^2),
        com rdot = (R/r) Rdot.  Requer l <= R(tau) sempre e f + rdot^2 >= 0.
        """
        R, Rdot = sol["R"], sol["Rdot"]
        if l > R.min() * (1 + 1e-9):
            return dict(valid=np.zeros_like(R, bool), frac_time_invalid=1.0, sigma=np.full_like(R, np.nan),
                        m_shell=np.full_like(R, np.nan), reason="l maior que o raio minimo do interior")
        r2 = np.maximum(R**2 - l**2, 0.0)
        r = np.sqrt(r2)
        f = 1 - 2 * self.M / R
        with np.errstate(divide="ignore", invalid="ignore"):
            rdot2 = np.where(r2 > 0, R**2 * Rdot**2 / r2, np.nan)
        F = f + rdot2
        # limite regular na garganta (r -> 0): R^2 Rdot^2/(R^2-l^2) -> R_b * Rddot_b
        ib = int(np.argmin(R))
        if r2[ib] <= 1e-12 * R[ib] ** 2:
            Rdd = np.gradient(np.gradient(R, sol["tau"]), sol["tau"])[ib]
            F[ib] = f[ib] + R[ib] * Rdd
            # pontos vizinhos com r pequeno: usar o mesmo limite se numericamente instaveis
            near = (r2 < 1e-8 * R**2)
            F[near] = f[near] + R[near] * Rdd
        ok = np.isfinite(F) & (F >= -1e-9)
        Kext = np.full_like(R, np.nan)
        Kext[ok] = (r[ok] / R[ok] ** 2) * np.sqrt(np.maximum(F[ok], 0))
        sigma = np.full_like(R, np.nan)
        sigma[ok] = (self.cos / R[ok] - Kext[ok]) / (4 * np.pi)
        return dict(valid=ok, F=F, sigma=sigma, m_shell=4 * np.pi * R**2 * sigma,
                    frac_time_invalid=float((~ok).mean()), reason="")

    def scan_l(self, sol, n=60):
        """Varre l em [0, R_b]: fracao do tempo em que a juncao e impossivel."""
        ls = np.linspace(0, sol["R_min"], n)
        frac = np.array([self.junction_black_bounce(sol, l)["frac_time_invalid"] for l in ls])
        return ls, frac

    def l_predicted(self):
        return self.R_b


def l_from_black_hole_mass(M_kg, m_fermion_kg=M_NEUTRON):
    """Previsao do modelo: l = R_b = (3M / (4 pi rho_b))^{1/3},  rho_b = 4 m^2/pi (Planck).  Retorna metros."""
    M_pl = M_kg / M_PLANCK
    m_pl = m_fermion_kg / M_PLANCK
    rho_b = 4 * m_pl**2 / np.pi
    return (3 * M_pl / (4 * np.pi * rho_b)) ** (1 / 3) * L_PLANCK


