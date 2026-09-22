"""
Perturbacoes cosmologicas atraves do ricochete: equacao de Mukhanov-Sasaki para um campo
de teste sem massa (mesma equacao dos modos tensoriais), com extracao exata do modo constante.

RESULTADO DESTE CODIGO: n_s ~ 1 para a contracao de poeira + ricochete de torcao
(reproduz a dualidade de Wands 1999 num fundo nao singular).  Ver docs/03_nascimento.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp


# ----------------------------------------------------------------------------
# 1. Espectro primordial atraves do ricochete de torcao
# ----------------------------------------------------------------------------
@dataclass
class BounceSpectrum:
    """Fundo: poeira + torcao (k=0), solucao exata  a^3 = a_b^3 + 6 pi rho0 t^2  (rho0 = 1).
    Perturbacao: campo escalar de teste sem massa, v = a phi."""

    a_b: float = 0.05
    t_max: float = 1e6
    n_t: int = 400000

    def __post_init__(self):
        t = np.concatenate([-np.logspace(np.log10(self.t_max), -6, self.n_t // 2), [0.0],
                            np.logspace(-6, np.log10(self.t_max), self.n_t // 2)])
        a = (self.a_b**3 + 6 * np.pi * t**2) ** (1 / 3)
        eta = cumulative_trapezoid(1 / a, t, initial=0.0)
        eta -= np.interp(0.0, t, eta)  # eta = 0 no ricochete
        # a''/a em tempo conforme = a^2 (H^2 + addot/a);  para a^3 = A + B t^2:
        adot = (2 * np.pi * t) * 2 / a**2  # d/dt (A + 6 pi t^2)^{1/3} = (1/3)(12 pi t) a^{-2}
        addot = 4 * np.pi / a**2 - 2 * adot**2 / a
        self.t, self.a, self.eta = t, a, eta
        self.app_over_a = a**2 * ((adot / a) ** 2 + addot / a)
        self.k_bounce = float(np.sqrt(self.app_over_a.max()))  # escala do ricochete
        self.eta_min, self.eta_max = eta[0], eta[-1]

    def _app(self, eta):
        return np.interp(eta, self.eta, self.app_over_a)

    def _a(self, eta):
        return np.interp(eta, self.eta, self.a)

    def mode(self, k, eta_eval=None):
        """Integra v_k de Bunch-Davies em eta_min ate eta_eval (por padrao, o fim da grade, ja no
        regime de poeira da expansao, a = a1 eta^2 com a1 = 2 pi/3)."""
        eta0 = self.eta_min
        # vacuo de Bunch-Davies na contracao de poeira: a solucao exata e f+ = e^{-ix}(1 - i/x)/sqrt(2k)
        x0 = k * eta0
        v0 = np.exp(-1j * x0) * (1 - 1j / x0) / np.sqrt(2 * k)
        dv0 = k * np.exp(-1j * x0) * (-1j * (1 - 1j / x0) + 1j / x0**2) / np.sqrt(2 * k)
        eta_eval = self.eta_max if eta_eval is None else eta_eval
        y0 = [v0.real, v0.imag, dv0.real, dv0.imag]

        def rhs(eta, y):
            w2 = k * k - self._app(eta)
            return [y[2], y[3], -w2 * y[0], -w2 * y[1]]

        sol = solve_ivp(rhs, (eta0, eta_eval), y0, method="DOP853", rtol=1e-9, atol=1e-12,
                        max_step=min(0.2 / k, 0.02))
        v = sol.y[0, -1] + 1j * sol.y[1, -1]
        dv = sol.y[2, -1] + 1j * sol.y[3, -1]
        return v, dv, eta_eval

    def growing_amplitude(self, k):
        """Decompoe v no fim da grade nas solucoes exatas do regime de poeira,
            f+ = e^{-ix}(1 - i/x),  f- = e^{+ix}(1 + i/x),  x = k eta,
        e devolve a amplitude CONSTANTE de phi = v/a fora do horizonte:
            f+ + f- = -(2/3) x^2 + O(x^4)   =>   phi_const = -(c+ + c-) k^2 / (3 a1),  a1 = 2 pi/3.
        Isso remove exatamente a contaminacao do modo decrescente (~ eta^-3)."""
        v, dv, eta = self.mode(k)
        x = k * eta
        fp = np.exp(-1j * x) * (1 - 1j / x)
        fm = np.exp(1j * x) * (1 + 1j / x)
        dfp = k * np.exp(-1j * x) * (-1j * (1 - 1j / x) + 1j / x**2)
        dfm = k * np.exp(1j * x) * (1j * (1 + 1j / x) - 1j / x**2)
        A = np.array([[fp, fm], [dfp, dfm]])
        cp, cm = np.linalg.solve(A, np.array([v, dv]))
        a1 = 2 * np.pi / 3
        return -(cp + cm) * k**2 / (3 * a1), cp, cm

    def spectrum(self, ks):
        """P_phi(k) = k^3 |phi_const|^2 / (2 pi^2): espectro do modo que sobrevive fora do horizonte."""
        P = []
        for k in ks:
            phi, _, _ = self.growing_amplitude(k)
            P.append(k**3 * abs(phi) ** 2 / (2 * np.pi**2))
        return np.array(P)

    def spectral_index(self, ks, P, kmin=None, kmax=None):
        """n_s - 1 = d ln P / d ln k ajustado no plato (k << k_bounce)."""
        kmin = kmin if kmin is not None else ks.min()
        kmax = kmax if kmax is not None else 0.1 * self.k_bounce
        m = (ks >= kmin) & (ks <= kmax)
        slope, _ = np.polyfit(np.log(ks[m]), np.log(P[m]), 1)
        return 1.0 + slope


