"""
O interior de um buraco negro E um universo.

Dentro do horizonte de Schwarzschild (r < 2M) a coordenada r e TEMPORAL e t e
ESPACIAL.  Renomeando T := r, a metrica fica

    ds^2 = -dT^2 / (2M/T - 1) + (2M/T - 1) dt^2 + T^2 dOmega^2,

que e uma cosmologia homogenea e anisotropica de Kantowski-Sachs com dois
fatores de escala:

    a_par(T)  = sqrt(2M/T - 1)   (direcao t, que "estica")
    a_perp(T) = T                (esferas, que "esmagam")

e tempo proprio  dtau = dT / sqrt(2M/T - 1).  O tempo proprio total do
horizonte ate a singularidade e exatamente pi*M.

Na metrica black-bounce de Simpson-Visser a mesma construcao vale com
R = sqrt(r^2 + l^2) no lugar de T: as esferas encolhem ate a_perp = l e
voltam a crescer -> o universo interno RICOCHETEIA e emerge no lado do
buraco branco.  Este modulo calcula essas duas cosmologias.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad

from .geometry import BlackBounce, Schwarzschild


def proper_time_to_singularity(M: float = 1.0) -> float:
    """tau = int_0^{2M} dT / sqrt(2M/T - 1) = pi M (calculado numericamente para conferir)."""
    val, _ = quad(lambda T: 1.0 / np.sqrt(2 * M / T - 1.0), 0.0, 2 * M, limit=200)
    return val


def kantowski_sachs_schwarzschild(M: float = 1.0, n: int = 4000):
    """Fatores de escala e taxas de Hubble direcionais dentro de um buraco negro de Schwarzschild.

    Retorna dict com T (=r), tau (tempo proprio desde o horizonte), a_par, a_perp,
    H_par, H_perp, volume (a_par * a_perp^2) e o escalar de Kretschmann.
    """
    # parametrizacao exata: T = M (1 + cos eta), tau = M (eta + sin eta), eta in (0, pi)
    # (resolve dtau = dT / sqrt(2M/T - 1) em forma fechada; tau_total = pi M exatamente)
    eta = np.linspace(1e-4, np.pi - 1e-4, n)
    T = M * (1 + np.cos(eta))
    tau = M * (eta + np.sin(eta))
    a_par = np.sqrt(2 * M / T - 1.0)
    a_perp = T.copy()
    dT_dtau = -np.sqrt(2 * M / T - 1.0)
    da_par_dT = -(M / T**2) / np.sqrt(2 * M / T - 1.0)
    H_par = (da_par_dT * dT_dtau) / a_par
    H_perp = dT_dtau / a_perp
    vol = a_par * a_perp**2
    K = Schwarzschild(M).kretschmann(T)
    return dict(T=T, tau=tau, a_par=a_par, a_perp=a_perp, H_par=H_par, H_perp=H_perp, volume=vol, kretschmann=K)


def kantowski_sachs_black_bounce(M: float = 1.0, l: float = 0.5, n: int = 6000):
    """Cosmologia interna de um black-bounce: das esferas que encolhem ate a_perp = l e voltam a crescer.

    A coordenada r percorre o interior de +r_h ate -r_h (r_h = sqrt(4M^2 - l^2)); r e temporal la dentro.
    """
    bb = BlackBounce(M, l)
    if not bb.horizons:
        raise ValueError("l >= 2M: nao ha horizonte; e um buraco de minhoca atravessavel.")
    rh = bb.horizons[1]
    # grade r = rh cos(theta): adensa perto dos horizontes, onde 1/sqrt(g) diverge (integravel)
    theta = np.linspace(1e-4, np.pi - 1e-4, n)
    r = rh * np.cos(theta)
    R = bb.R(r)
    g = 2 * M / R - 1.0  # > 0 dentro
    a_par = np.sqrt(g)
    a_perp = R
    dtau_dr = 1.0 / np.sqrt(g)
    tau = -cumulative_trapezoid(dtau_dr, r, initial=0.0)
    dr_dtau = -np.sqrt(g)
    dR_dr = r / R
    H_perp = (dR_dr * dr_dtau) / a_perp
    dg_dr = -2 * M * r / R**3
    da_par_dr = 0.5 * dg_dr / np.sqrt(g)
    H_par = (da_par_dr * dr_dtau) / a_par
    K = bb.kretschmann(r)
    return dict(r=r, R=R, tau=tau, a_par=a_par, a_perp=a_perp, H_par=H_par, H_perp=H_perp,
                volume=a_par * a_perp**2, kretschmann=K, r_h=rh, tau_bounce=float(np.interp(0.0, -r, tau)))


def radial_infall_black_bounce(M=1.0, l=0.5, r0=8.0, tau_max=60.0, n=6000):
    """Queda radial livre a partir do repouso em r0 atravessando a garganta.

    (dr/dtau)^2 = E^2 - f(r),  E^2 = f(r0).  Integramos a forma de 2a ordem
    r'' = -f'(r)/2 (valida em todo lado, inclusive nos horizontes e na garganta).
    """
    bb = BlackBounce(M, l)
    from scipy.integrate import solve_ivp

    def rhs(tau, y):
        r, v = y
        return [v, bb.timelike_accel(r, 0.0)]

    sol = solve_ivp(rhs, (0, tau_max), [r0, 0.0], max_step=tau_max / n, rtol=1e-10, atol=1e-12, dense_output=True)
    tau = np.linspace(0, sol.t[-1], n)
    y = sol.sol(tau)
    return dict(tau=tau, r=y[0], v=y[1], R=bb.R(y[0]), f=bb.f(y[0]), kretschmann=bb.kretschmann(y[0]))
