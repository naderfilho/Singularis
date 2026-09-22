"""
Colapso de Oppenheimer-Snyder (1939): uma estrela de poeira homogenea colapsa
e forma um buraco negro.  O fato central deste projeto:

    O INTERIOR DA ESTRELA EM COLAPSO E EXATAMENTE UM UNIVERSO FRW FECHADO
    (a mesma metrica de Friedmann do Big Bang, com o tempo invertido).

    interior:  ds^2 = -dtau^2 + a(eta)^2 [ dchi^2 + sin^2(chi) dOmega^2 ],   0 <= chi <= chi0
               a(eta) = a_m (1 + cos eta)/2,     tau(eta) = a_m (eta + sin eta)/2
    exterior:  Schwarzschild com massa M
    juncao:    sin^2(chi0) = 2M/R0,   a_m = R0 / sin(chi0)

Portanto "um universo em contracao dentro de um buraco negro" nao e uma
metafora: e a solucao exata de Einstein para o colapso.  Se algo impedir a
singularidade (torcao, gravidade quantica), esse universo ricocheteia e
vira um Big Bang.  Ver semente.bounce.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid

from ..geometry.spherical import Schwarzschild


@dataclass
class OppenheimerSnyder:
    M: float = 1.0
    R0: float = 10.0  # raio inicial da estrela (em repouso)

    def __post_init__(self):
        if self.R0 <= 2 * self.M:
            raise ValueError("A estrela precisa comecar fora do proprio raio de Schwarzschild.")
        self.chi0 = float(np.arcsin(np.sqrt(2 * self.M / self.R0)))
        self.a_m = self.R0 / np.sin(self.chi0)
        self.eta_h = float(np.arccos(4 * self.M / self.R0 - 1.0))  # superficie cruza R = 2M
        self.tau_h = self.tau(self.eta_h)
        self.tau_s = self.tau(np.pi)  # singularidade

    # --- interior FRW ---------------------------------------------------------
    def a(self, eta):
        return self.a_m * (1 + np.cos(eta)) / 2

    def tau(self, eta):
        return self.a_m * (eta + np.sin(eta)) / 2

    def R_surface(self, eta):
        """Raio areal da superficie da estrela: R = a(eta) sin(chi0) = R0 (1+cos eta)/2."""
        return self.a(eta) * np.sin(self.chi0)

    def density(self, eta):
        """Densidade homogenea de poeira: rho = 3M / (4 pi R^3)."""
        return 3 * self.M / (4 * np.pi * self.R_surface(eta) ** 3)

    def eta_from_tau(self, tau, n=20000):
        eta = np.linspace(0, np.pi, n)
        return np.interp(tau, self.tau(eta), eta)

    # --- horizonte de eventos dentro da estrela -------------------------------
    def event_horizon_interior(self, n=400):
        """Raio nulo de saida dchi/deta = +1 que alcanca a superficie exatamente quando R = 2M.

        chi_EH(eta) = chi0 - (eta_h - eta).  Nasce em chi = 0 em eta_0 = eta_h - chi0 (se > 0):
        o horizonte de eventos se forma no CENTRO da estrela ANTES de a superficie cruzar 2M,
        e cresce ate coincidir com r = 2M.
        """
        eta0 = max(self.eta_h - self.chi0, 0.0)
        eta = np.linspace(eta0, self.eta_h, n)
        chi = self.chi0 - (self.eta_h - eta)
        R = self.a(eta) * np.sin(chi)
        return dict(eta=eta, tau=self.tau(eta), chi=chi, R=R, eta_birth=eta0, tau_birth=self.tau(eta0))

    # --- linha de mundo da superficie em Kruskal (usa v de Eddington-Finkelstein) ----
    def surface_kruskal(self, n=3000, eta_max=None, v_shift=0.0):
        """Coordenadas de Kruskal (T, X) da superficie da estrela ao longo do colapso.

        dv/dtau = 1 / (E + sqrt(E^2 - f)),  E = sqrt(1 - 2M/R0)  (finito atraves do horizonte).
        Depois U = -(r/2M - 1) e^{r/2M} / V,  V = e^{v/4M}.
        """
        M = self.M
        eta_max = np.pi * 0.999 if eta_max is None else eta_max
        eta = np.linspace(0, eta_max, n)
        R = self.R_surface(eta)
        tau = self.tau(eta)
        E = np.sqrt(1 - 2 * M / self.R0)
        f = 1 - 2 * M / R
        dv_dtau = 1.0 / (E + np.sqrt(np.maximum(E**2 - f, 0.0)))
        v = cumulative_trapezoid(dv_dtau, tau, initial=0.0)
        v0 = self.R0 + 2 * M * np.log(self.R0 / (2 * M) - 1)  # t=0 em tau=0 -> v = r*(R0)
        v = v + v0 + v_shift  # v_shift = translacao temporal (simetria do exterior)
        V = np.exp(v / (4 * M))
        U = -(R / (2 * M) - 1.0) * np.exp(R / (2 * M)) / V
        T = (V + U) / 2
        X = (V - U) / 2
        return dict(eta=eta, tau=tau, R=R, T=T, X=X, v=v)

    def summary(self) -> dict:
        return dict(
            M=self.M, R0=self.R0, chi0=self.chi0, a_max=self.a_m,
            tau_horizonte=self.tau_h, tau_singularidade=self.tau_s,
            rho_inicial=float(self.density(0.0)),
            comentario="interior = FRW fechado (k=+1) de poeira; universo em contracao.",
        )
