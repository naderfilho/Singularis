"""
Geometrias com rotacao: Kerr e o black-bounce rotativo de Mazza, Franzin & Liberati (2021).

Kerr (RESULTADO MATEMATICO), Boyer-Lindquist, G = c = 1, a = J/M:
   Sigma = r^2 + a^2 cos^2 th,   Delta = r^2 - 2 M r + a^2
   ds^2 = -(1 - 2Mr/Sigma) dt^2 - (4 M a r sin^2 th / Sigma) dt dphi + (Sigma/Delta) dr^2 + Sigma dth^2
          + [ r^2 + a^2 + 2 M a^2 r sin^2 th / Sigma ] sin^2 th dphi^2
   horizontes  r_+- = M +- sqrt(M^2 - a^2)          (a <= M)
   ergosuperficie  r_E(th) = M + sqrt(M^2 - a^2 cos^2 th)   (g_tt = 0)
   arrasto  omega(r, th) = 2 M a r / [ (r^2 + a^2)^2 - a^2 Delta sin^2 th ]   (ZAMO)
   Omega_H = a / (r_+^2 + a^2),   kappa = (r_+ - r_-) / (2 (r_+^2 + a^2)),   A = 4 pi (r_+^2 + a^2)
   ISCO (Bardeen, Press & Teukolsky 1972):  Z1 = 1 + (1-a^2)^{1/3} [(1+a)^{1/3} + (1-a)^{1/3}],
        Z2 = sqrt(3 a^2 + Z1^2),  r_isco = 3 + Z2 -+ sqrt((3 - Z1)(3 + Z1 + 2 Z2))  (M = 1; -: prograda)
   Kretschmann: K = 48 M^2 (r^2 - a^2 c^2)(r^4 - 14 a^2 r^2 c^2 + a^4 c^4) / (r^2 + a^2 c^2)^6,  c = cos th

Black-bounce rotativo (MODELO DA LITERATURA, Mazza, Franzin & Liberati 2021, JCAP 04, 082):
   substituicao r -> sqrt(r^2 + l^2) em Sigma e Delta (e no fator 2 M r -> 2 M sqrt(r^2 + l^2)):
   Sigma = r^2 + l^2 + a^2 cos^2 th,   Delta = r^2 + l^2 - 2 M sqrt(r^2 + l^2) + a^2.
   horizontes: sqrt(r^2 + l^2) = r_+-(Kerr)  =>  r_h = +- sqrt(r_+-^2 - l^2)   (se r_+- > l)
   ergosuperficie: sqrt(r^2 + l^2) = M + sqrt(M^2 - a^2 cos^2 th)
   Os autores mostram que a geometria e regular (invariantes de curvatura finitos) para l > 0; aqui o escalar
   de Kretschmann e calculado SIMBOLICAMENTE (sympy) e avaliado numericamente para verificar isso, e
   reduz-se ao de Kerr para l -> 0.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class Kerr:
    M: float = 1.0
    a: float = 0.5
    l: float = 0.0   # l > 0: black-bounce rotativo de Mazza-Franzin-Liberati

    # --- funcoes basicas ------------------------------------------------------------------
    def rho_bar(self, r):
        """Raio 'regularizado' sqrt(r^2 + l^2) (= r para Kerr)."""
        return np.sqrt(np.asarray(r, float) ** 2 + self.l**2)

    def Sigma(self, r, th):
        return self.rho_bar(r) ** 2 + self.a**2 * np.cos(th) ** 2

    def Delta(self, r):
        rb = self.rho_bar(r)
        return rb**2 - 2 * self.M * rb + self.a**2

    @property
    def kerr_horizons(self):
        if self.a > self.M:
            return ()
        s = np.sqrt(self.M**2 - self.a**2)
        return (self.M - s, self.M + s)

    @property
    def horizons(self):
        """Posicoes em r (coordenada) dos horizontes; para l > 0 so existem se r_+-(Kerr) > l."""
        out = []
        for rk in self.kerr_horizons:
            if rk > self.l:
                out.append(float(np.sqrt(rk**2 - self.l**2)))
        return tuple(sorted(out))

    def ergosurface(self, th):
        th = np.asarray(th, float)
        rk = self.M + np.sqrt(np.maximum(self.M**2 - self.a**2 * np.cos(th) ** 2, 0.0))
        return np.sqrt(np.maximum(rk**2 - self.l**2, 0.0))

    def frame_dragging(self, r, th=np.pi / 2):
        """Velocidade angular do ZAMO omega = -g_tphi / g_phiphi."""
        rb = self.rho_bar(r)
        S = self.Sigma(r, th)
        D = self.Delta(r)
        s2 = np.sin(th) ** 2
        return 2 * self.M * self.a * rb / ((rb**2 + self.a**2) ** 2 - self.a**2 * D * s2)

    @property
    def Omega_H(self):
        if not self.kerr_horizons:
            return None
        rp = self.kerr_horizons[1]
        return self.a / (rp**2 + self.a**2)

    @property
    def surface_gravity(self):
        if not self.kerr_horizons:
            return None
        rm, rp = self.kerr_horizons
        return (rp - rm) / (2 * (rp**2 + self.a**2))

    @property
    def horizon_area(self):
        if not self.kerr_horizons:
            return None
        rp = self.kerr_horizons[1]
        return 4 * np.pi * (rp**2 + self.a**2)

    def isco(self, prograde=True):
        """ISCO de Kerr (BPT 1972) em unidades de M; para l > 0 e o valor de Kerr em rho_bar (aproximacao
        documentada: a substituicao r -> rho_bar preserva as orbitas equatoriais em termos de rho_bar)."""
        a = self.a / self.M
        Z1 = 1 + (1 - a**2) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
        Z2 = np.sqrt(3 * a**2 + Z1**2)
        s = -1 if prograde else 1
        return self.M * (3 + Z2 + s * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))

    def kretschmann_kerr_analytic(self, r, th=np.pi / 2):
        r = np.asarray(r, float)
        c = np.cos(th)
        a = self.a
        return 48 * self.M**2 * (r**2 - a**2 * c**2) * (r**4 - 14 * a**2 * r**2 * c**2 + a**4 * c**4) / (r**2 + a**2 * c**2) ** 6

    def kretschmann(self, r, th=np.pi / 2):
        """Kretschmann do black-bounce rotativo (simbolico via sympy, cacheado); para l = 0 e Kerr."""
        fn = _kretschmann_rotating_bb()
        r = np.asarray(r, float).copy()
        th = np.clip(np.asarray(th, float), 1e-6, np.pi - 1e-6)   # polos: singularidade de coordenadas
        # onde g_tt = 0 (ergosuperficie) ou Delta = 0 (horizonte) a expressao simbolica tem 0/0 removivel:
        # desloca r por 1e-7 M (a curvatura e continua ali)
        rb = self.rho_bar(r)
        gtt = 1 - 2 * self.M * rb / (rb**2 + self.a**2 * np.cos(th) ** 2)
        bad = (np.abs(gtt) < 1e-9) | (np.abs(self.Delta(r)) < 1e-9)
        r = np.where(bad, r + 1e-7 * self.M, r)
        with np.errstate(all="ignore"):
            return np.asarray(fn(r, th, self.M, self.a, self.l), float)

    def classify(self) -> dict:
        hs = self.horizons
        info = dict(horizons=hs, ergosurface_equator=float(self.ergosurface(np.pi / 2)), l=self.l, a=self.a)
        if self.a > self.M:
            info["kind"] = "super-extremal (a > M): sem horizontes" + (" - anel singular nu" if self.l == 0 else " - regular")
        elif self.l == 0:
            info["kind"] = "Kerr" if self.a < self.M else "Kerr extremal"
        elif len(hs) == 2:
            info["kind"] = "black-bounce rotativo com horizonte interno"
        elif len(hs) == 1:
            info["kind"] = "black-bounce rotativo (garganta espacial entre horizontes externos)"
        else:
            info["kind"] = "buraco de minhoca rotativo atravessavel"
        return info


@lru_cache(maxsize=None)
def _kretschmann_rotating_bb():
    """Riemann completo da metrica de Mazza-Franzin-Liberati em Boyer-Lindquist (simbolico), lambdificado."""
    t, r, th, ph = sp.symbols("t r theta phi", real=True)
    M, a, l = sp.symbols("M a ell", positive=True)
    rb = sp.sqrt(r**2 + l**2)
    Sig = rb**2 + a**2 * sp.cos(th) ** 2
    Del = rb**2 - 2 * M * rb + a**2
    s2 = sp.sin(th) ** 2
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M * rb / Sig)
    g[0, 3] = g[3, 0] = -2 * M * a * rb * s2 / Sig
    g[1, 1] = Sig / Del
    g[2, 2] = Sig
    g[3, 3] = (rb**2 + a**2 + 2 * M * a**2 * rb * s2 / Sig) * s2
    x = (t, r, th, ph)
    ginv = g.inv()
    n = 4
    Gam = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(j, n):
                e = sum(ginv[i, d] * (sp.diff(g[d, k], x[j]) + sp.diff(g[d, j], x[k]) - sp.diff(g[j, k], x[d])) for d in range(n)) / 2
                Gam[i][j][k] = Gam[i][k][j] = e
    Riem = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for m in range(k + 1, n):
                    e = sp.diff(Gam[i][m][j], x[k]) - sp.diff(Gam[i][k][j], x[m])
                    e += sum(Gam[i][k][p] * Gam[p][m][j] - Gam[i][m][p] * Gam[p][k][j] for p in range(n))
                    Riem[(i, j, k, m)] = e
    # K = R_{ijkm} R^{ijkm}; baixar o primeiro indice com g e levantar os outros tres com ginv
    # (usa antissimetria em (k,m): soma sobre k<m com fator 2 em cada par, e sobre todos i,j)
    Rlow = {}
    for (i, j, k, m), e in Riem.items():
        Rlow[(i, j, k, m)] = sum(g[i, p] * Riem[(p, j, k, m)] for p in range(n))
    K = sp.S(0)
    idx = list(Riem.keys())
    for (i, j, k, m) in idx:
        # R^{ijkm} = ginv[i,p] ginv[j,q] ginv[k,s] ginv[m,u] R_{pqsu}; a metrica so mistura t e phi
        up = sp.S(0)
        for p in range(n):
            if ginv[i, p] == 0:
                continue
            for q in range(n):
                if ginv[j, q] == 0:
                    continue
                for s in range(n):
                    if ginv[k, s] == 0:
                        continue
                    for u in range(n):
                        if ginv[m, u] == 0:
                            continue
                        if s < u:
                            val = Rlow[(p, q, s, u)]
                        elif s > u:
                            val = -Rlow[(p, q, u, s)]
                        else:
                            continue
                        up += ginv[i, p] * ginv[j, q] * ginv[k, s] * ginv[m, u] * val
        K += 2 * Rlow[(i, j, k, m)] * up
    return sp.lambdify((r, th, M, a, l), K, "numpy")
