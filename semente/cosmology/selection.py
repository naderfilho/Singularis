"""
Selecao natural cosmologica (Smolin, 1992): se cada buraco negro gera um
universo-filho com constantes ligeiramente diferentes das do pai, entao,
depois de muitas geracoes, a maioria dos universos tera constantes que
MAXIMIZAM a producao de buracos negros.  Previsao testavel: nosso universo
deve estar num maximo local de "fecundidade" (numero de buracos negros).

Este modulo implementa a dinamica populacional: universos com parametros
theta (adimensionais, ex.: escala de massa dos fermions, constante
cosmologica, amplitude de perturbacoes), fecundidade N_bh(theta) e mutacao
gaussiana.  A paisagem N_bh e um modelo de brinquedo, mas a dinamica
(replicador + mutacao) e exata.

Referencia: L. Smolin, "Did the universe evolve?", Class. Quantum Grav. 9, 173 (1992);
"The Life of the Cosmos" (1997).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.units import M_PLANCK, M_SUN  # noqa: F401


def fecundity_landscape(theta: np.ndarray) -> np.ndarray:
    """Numero esperado de buracos negros por universo, N_bh(theta) >= 0.

    theta[:,0] ~ log10 da razao (massa do nucleon / massa de Planck) deslocada:  estrelas so existem
                 numa faixa estreita  -> gaussiana estreita.
    theta[:,1] ~ constante cosmologica (em unidades arbitrarias): Lambda grande impede a formacao
                 de galaxias -> decaimento exponencial.
    theta[:,2] ~ amplitude de perturbacoes Q: pequena -> poucas estruturas; grande -> tudo vira BN cedo
                 -> pico assimetrico.
    """
    x, lam, q = theta[:, 0], theta[:, 1], theta[:, 2]
    stars = np.exp(-0.5 * (x / 0.6) ** 2)
    galaxies = np.exp(-np.maximum(lam, 0.0) / 0.8) * (1 / (1 + np.exp(-(lam + 2.5) / 0.3)))
    structure = (q**2) * np.exp(-q / 1.0)
    return 1e4 * stars * galaxies * structure


@dataclass
class CosmicSelection:
    n_pop: int = 4000
    mutation: float = 0.08
    rng: np.random.Generator = field(default_factory=lambda: np.random.default_rng(7))
    theta0_spread: float = 2.0

    def run(self, generations: int = 60):
        """Cada geracao: cada universo tem N_bh filhos (amostragem multinomial com ponderacao),
        os filhos herdam theta + ruido.  Populacao mantida constante (renormalizacao)."""
        theta = self.rng.normal(0.0, self.theta0_spread, size=(self.n_pop, 3))
        theta[:, 2] = np.abs(theta[:, 2]) + 0.05
        hist_mean, hist_std, hist_fec, snapshots = [], [], [], []
        for g in range(generations):
            fec = fecundity_landscape(theta)
            hist_mean.append(theta.mean(axis=0))
            hist_std.append(theta.std(axis=0))
            hist_fec.append(fec.mean())
            if g in (0, 3, 10, 30, generations - 1):
                snapshots.append((g, theta.copy()))
            w = fec + 1e-12
            w /= w.sum()
            parents = self.rng.choice(self.n_pop, size=self.n_pop, p=w)
            theta = theta[parents] + self.rng.normal(0.0, self.mutation, size=(self.n_pop, 3))
            theta[:, 2] = np.abs(theta[:, 2])
        return dict(mean=np.array(hist_mean), std=np.array(hist_std), fecundity=np.array(hist_fec),
                    snapshots=snapshots, final=theta)


# ----------------------------------------------------------------------------
# C. Selecao natural cosmologica com orcamento herdado
# ----------------------------------------------------------------------------
@dataclass
class BudgetedSelection:
    """Cada universo nasce dentro de um buraco negro de massa m_parent (Planck) e recebe
    capacidade C = m_parent^2.  Ele forma N(theta) buracos negros de massa m(theta) cada,
    mas so os que cabem no orcamento (sum m^2 <= C) geram filhos.

    m(theta) ~ massa de Chandrasekhar-like: m = m_ref * exp(-2 theta_0)  (theta_0 = log da escala
    de massa dos fermions em relacao ao nosso universo).  N(theta) = paisagem de fecundidade.
    """

    n_pop: int = 3000
    mutation: float = 0.08
    m_ref: float = 1e39  # ~10 M_sol em massas de Planck (2e31 kg / 2.2e-8 kg ~ 1e39)
    M_root: float = 1e48  # ~ 10^10 M_sol: o buraco negro raiz da arvore
    rng: np.random.Generator = field(default_factory=lambda: np.random.default_rng(3))

    def run(self, generations=40):
        theta = self.rng.normal(0.0, 1.0, size=(self.n_pop, 3))
        theta[:, 2] = np.abs(theta[:, 2]) + 0.05
        cap = np.full(self.n_pop, self.M_root**2)
        alive = np.ones(self.n_pop, bool)
        hist = []
        for g in range(generations):
            N = fecundity_landscape(theta)
            m = self.m_ref * np.exp(-2 * theta[:, 0])
            N_fit = np.minimum(N, cap / m**2)  # filhos que cabem no orcamento
            N_fit[~alive] = 0.0
            N_fit[m < 1.0] = 0.0  # buracos negros abaixo da massa de Planck nao existem
            hist.append(dict(gen=g, alive=int(alive.sum()), N_mean=float(N[alive].mean()) if alive.any() else 0.0,
                             N_fit_mean=float(N_fit[alive].mean()) if alive.any() else 0.0,
                             cap_median=float(np.median(cap[alive])) if alive.any() else 0.0,
                             m_median=float(np.median(m[alive])) if alive.any() else 0.0))
            if N_fit.sum() <= 0:
                break
            w = N_fit / N_fit.sum()
            parents = self.rng.choice(self.n_pop, size=self.n_pop, p=w)
            theta = theta[parents] + self.rng.normal(0.0, self.mutation, size=(self.n_pop, 3))
            theta[:, 2] = np.abs(theta[:, 2])
            cap = m[parents] ** 2  # o filho herda a capacidade do buraco negro que o gerou
            alive = N_fit[parents] > 0
        return hist


