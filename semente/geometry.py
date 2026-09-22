"""
Geometria: Schwarzschild, extensao maximal de Kruskal-Szekeres (onde o buraco
branco aparece "sem querer"), compactificacao de Penrose e a metrica
black-bounce de Simpson-Visser (buraco negro cujo interior atravessa uma
garganta e emerge como buraco branco em OUTRO universo).

Tudo em unidades geometricas G = c = 1.

Referencias
-----------
Kruskal, M. (1960) Phys. Rev. 119, 1743.  Szekeres, G. (1960) Publ. Mat. Debrecen 7, 285.
Simpson, A. & Visser, M. (2019) JCAP 02, 042  ("Black-bounce to traversable wormhole").
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import lambertw


# ----------------------------------------------------------------------------
# Schwarzschild
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class Schwarzschild:
    """Metrica de Schwarzschild  ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2,  f = 1 - 2M/r."""

    M: float = 1.0

    @property
    def rs(self) -> float:
        return 2.0 * self.M

    def f(self, r):
        return 1.0 - 2.0 * self.M / np.asarray(r, float)

    def tortoise(self, r):
        """r* = r + 2M ln|r/2M - 1|  (coordenada tartaruga)."""
        r = np.asarray(r, float)
        return r + 2 * self.M * np.log(np.abs(r / (2 * self.M) - 1.0))

    # --- Kruskal-Szekeres --------------------------------------------------
    def kruskal_from_tr(self, t, r):
        """(t, r) -> (T, X) de Kruskal-Szekeres.

        Regiao I (exterior, r>2M):  X = sqrt(r/2M-1) e^{r/4M} cosh(t/4M), T = ... sinh
        Regiao II (buraco negro, r<2M): T = sqrt(1-r/2M) e^{r/4M} cosh(t/4M), X = ... sinh
        Vale a identidade  X^2 - T^2 = (r/2M - 1) e^{r/2M}  em todas as regioes.
        """
        t = np.asarray(t, float)
        r = np.asarray(r, float)
        M = self.M
        amp = np.sqrt(np.abs(r / (2 * M) - 1.0)) * np.exp(r / (4 * M))
        outside = r > 2 * M
        X = np.where(outside, amp * np.cosh(t / (4 * M)), amp * np.sinh(t / (4 * M)))
        T = np.where(outside, amp * np.sinh(t / (4 * M)), amp * np.cosh(t / (4 * M)))
        return T, X

    def r_from_kruskal(self, T, X):
        """Inverte X^2 - T^2 = (r/2M-1) e^{r/2M} via a funcao W de Lambert.

        r/2M - 1 = W( (X^2 - T^2)/e ).  A singularidade r=0 e a hiperbole T^2 - X^2 = 1,
        que corresponde ao ponto de ramificacao W(-1/e) = -1.
        """
        w = (np.asarray(X, float) ** 2 - np.asarray(T, float) ** 2) / np.e
        w = np.maximum(w, -1.0 / np.e * (1 - 1e-12))  # alem da singularidade nao ha espaco-tempo
        return 2 * self.M * np.maximum(1.0 + np.real(lambertw(w, 0)), 0.0)

    def t_from_kruskal(self, T, X):
        """t de Schwarzschild a partir de Kruskal: t = 4M artanh(T/X) (I, IV) ou 4M artanh(X/T) (II, III)."""
        T = np.asarray(T, float)
        X = np.asarray(X, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            inside = np.abs(T) > np.abs(X)
            ratio = np.where(inside, X / T, T / X)
            return 4 * self.M * np.arctanh(np.clip(ratio, -1 + 1e-15, 1 - 1e-15))

    @staticmethod
    def region(T, X) -> np.ndarray:
        """Rotula a regiao de Kruskal: 1 = nosso exterior, 2 = buraco NEGRO, 3 = buraco BRANCO, 4 = outro exterior."""
        T = np.asarray(T, float)
        X = np.asarray(X, float)
        reg = np.full(T.shape, 1, dtype=int)
        reg[(T > np.abs(X))] = 2
        reg[(T < -np.abs(X))] = 3
        reg[(X < -np.abs(T))] = 4
        return reg

    def const_r_curve(self, r, tmax=6.0, n=400):
        """Hiperbole de r constante no diagrama de Kruskal (ramo de t crescente)."""
        t = np.linspace(-tmax * self.M, tmax * self.M, n)
        T, X = self.kruskal_from_tr(t, np.full_like(t, r))
        return T, X

    def flamm_paraboloid(self, r):
        """Paraboloide de Flamm z(r) = sqrt(8M (r-2M)): a fatia t=const, theta=pi/2 mergulhada em R^3.

        A garganta em r = 2M e uma ponte de Einstein-Rosen para a regiao IV (o outro exterior).
        """
        r = np.asarray(r, float)
        return np.sqrt(np.maximum(8 * self.M * (r - 2 * self.M), 0.0))

    def kretschmann(self, r):
        """Escalar de Kretschmann K = R_{abcd}R^{abcd} = 48 M^2 / r^6 (diverge em r=0: singularidade fisica)."""
        return 48.0 * self.M**2 / np.asarray(r, float) ** 6


# ----------------------------------------------------------------------------
# Compactificacao de Penrose
# ----------------------------------------------------------------------------
def penrose_from_kruskal(T, X):
    """Compactifica Kruskal: U = T - X, V = T + X  ->  u = arctan U, v = arctan V.

    Devolve coordenadas de desenho (x, y) = ((v - u)/2, (v + u)/2), em que os
    infinitos nulos (scri+/-), i0, i+ e i- ficam em posicoes finitas.
    """
    T = np.asarray(T, float)
    X = np.asarray(X, float)
    u = np.arctan(T - X)
    v = np.arctan(T + X)
    return (v - u) / 2.0, (v + u) / 2.0


# ----------------------------------------------------------------------------
# Simpson-Visser black-bounce  (buraco negro -> garganta -> buraco branco em outro universo)
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class BlackBounce:
    """Metrica de Simpson-Visser (2019):

        ds^2 = -f dt^2 + dr^2/f + (r^2 + l^2) dOmega^2,   f = 1 - 2M / sqrt(r^2 + l^2),
        r in (-inf, +inf).

    * l = 0            : Schwarzschild exato.
    * 0 < l < 2M       : buraco negro REGULAR. Horizontes em r = +-sqrt(4M^2 - l^2).
                         A singularidade r=0 e substituida por uma "garganta" ESPACIAL
                         (um ricochete): tudo que cai atravessa r=0 e emerge de um
                         BURACO BRANCO no universo r<0.
    * l = 2M           : extremal (garganta nula, horizonte unico em r=0).
    * l > 2M           : buraco de minhoca atravessavel entre os dois universos.

    A metrica e uma solucao exata das equacoes de Einstein com um tensor de
    energia-momento efetivo que viola a condicao de energia nula perto da
    garganta - o tipo de "materia" que efeitos de gravidade quantica fornecem.
    O escalar de Kretschmann e finito em todo lugar para l>0.
    """

    M: float = 1.0
    l: float = 0.5

    def R(self, r):
        """Raio areal R = sqrt(r^2 + l^2) (nunca chega a zero se l>0)."""
        return np.sqrt(np.asarray(r, float) ** 2 + self.l**2)

    def f(self, r):
        return 1.0 - 2.0 * self.M / self.R(r)

    @property
    def horizons(self):
        """Posicoes dos horizontes (vazio se l > 2M: buraco de minhoca atravessavel)."""
        if self.l >= 2 * self.M:
            return ()
        rh = float(np.sqrt(4 * self.M**2 - self.l**2))
        return (-rh, rh)

    @property
    def kind(self) -> str:
        if self.l == 0:
            return "Schwarzschild (singular)"
        if self.l < 2 * self.M:
            return "buraco negro regular com ricochete para buraco branco"
        if self.l == 2 * self.M:
            return "buraco negro extremal / garganta nula"
        return "buraco de minhoca atravessavel"

    def kretschmann(self, r):
        """Escalar de Kretschmann K = R_{abcd}R^{abcd} da metrica -f dt^2 + dr^2/f + R(r)^2 dOmega^2.

        Em base ortonormal, os componentes independentes de Riemann sao
            A = R_{trtr}          = f''/2
            B = R_{t th t th}     = f' R'/(2R)
            C = R_{r th r th}     = -(f R'' + f' R'/2)/R
            D = R_{th ph th ph}   = (1 - f R'^2)/R^2
        e K = 4A^2 + 8B^2 + 8C^2 + 4D^2.  Para l=0 reproduz 48 M^2/r^6.
        """
        r = np.asarray(r, float)
        R = self.R(r)
        l, M = self.l, self.M
        f = 1 - 2 * M / R
        Rp = r / R
        Rpp = l**2 / R**3
        fp = 2 * M * r / R**3
        fpp = 2 * M * (l**2 - 2 * r**2) / R**5
        A = fpp / 2
        B = fp * Rp / (2 * R)
        C = -(f * Rpp + fp * Rp / 2) / R
        D = (1 - f * Rp**2) / R**2
        return 4 * A**2 + 8 * B**2 + 8 * C**2 + 4 * D**2

    def null_accel(self, r, L):
        """Geodesicas nulas: (dr/dlambda)^2 = E^2 - L^2 f/(r^2+l^2)  =>  d^2r/dlambda^2 = V'(r)/2.

        V'(r)/2 = L^2 r (R - 3M) / R^5.  Para l=0 recupera r'' = L^2 (r - 3M)/r^4.
        """
        R = self.R(r)
        return L**2 * r * (R - 3 * self.M) / R**5

    def timelike_accel(self, r, L):
        """Geodesicas tipo tempo: (dr/dtau)^2 = E^2 - f (1 + L^2/R^2)  =>  r'' = -(1/2) d/dr[f(1+L^2/R^2)]."""
        R = self.R(r)
        M, l = self.M, self.l
        # d/dr [ (1 - 2M/R)(1 + L^2/R^2) ] = (r/R) d/dR[...]
        dR = 2 * M / R**2 * (1 + L**2 / R**2) + (1 - 2 * M / R) * (-2 * L**2 / R**3)
        return -0.5 * (r / R) * dR

    def photon_sphere(self):
        """Esfera de fotons: R = 3M  ->  r = +-sqrt(9M^2 - l^2)."""
        if self.l >= 3 * self.M:
            return ()
        rp = float(np.sqrt(9 * self.M**2 - self.l**2))
        return (-rp, rp)

    def describe(self) -> str:
        return {
            "Schwarzschild (singular)": "Diagrama de Kruskal padrao: regioes I, II (BN), III (BB), IV.",
            "buraco negro regular com ricochete para buraco branco":
                "Diagrama de Penrose em escada infinita: BN no nosso universo -> garganta r=0 "
                "(hipersuperficie espacial regular) -> BB no universo seguinte -> ...",
            "buraco negro extremal / garganta nula": "Horizonte duplo em r=0.",
            "buraco de minhoca atravessavel":
                "Dois universos assintoticamente planos ligados por uma garganta em r=0, sem horizontes.",
        }[self.kind]
