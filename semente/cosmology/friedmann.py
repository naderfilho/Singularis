"""
Fundo de Friedmann unificado: um unico integrador para GR, Einstein–Cartan e LQC efetiva.

    H^2 = (8 pi/3) F(a) - k/a^2,     a''/a = A(a)

onde F e A vem de uma estrategia `correction` (GR, EinsteinCartanCorrection, LQCCorrection)
aplicada ao conteudo de fluidos (rho0_i, w_i).  Unidades de Planck (G = c = hbar = 1).

Resultados classificados:
* `BackgroundSolution` e um RESULTADO NUMERICO do modelo escolhido.
* `exact_torsion_dust` e um RESULTADO MATEMATICO (solucao fechada) usado para validar.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import brentq

from ..core.solvers import DEFAULT_ODE, ODESettings, integrate
from ..quantum.einstein_cartan import EinsteinCartanCorrection
from ..quantum.lqc import LQCCorrection


@dataclass(frozen=True)
class Fluid:
    """Fluido perfeito barotropico: rho = rho0 a^{-3(1+w)}, p = w rho."""

    rho0: float
    w: float = 0.0
    name: str = "dust"


@dataclass(frozen=True)
class GRCorrection:
    """Sem correcao: relatividade geral."""

    name: str = "GR"

    def rho_eff(self, a, rho, p):
        return rho

    def p_eff(self, a, rho, p):
        return p

    def H2_contribution(self, a, rho, p):
        return (8 * np.pi / 3) * rho

    def accel(self, a, rho, p):
        return -(4 * np.pi / 3) * (rho + 3 * p)

    def classical_limit(self):
        return self


@dataclass
class BackgroundSolution:
    t: np.ndarray
    a: np.ndarray
    adot: np.ndarray
    H: np.ndarray
    rho: np.ndarray
    p: np.ndarray
    rho_eff: np.ndarray
    accel: np.ndarray                  # a''/a
    epsilon: np.ndarray                # -Hdot/H^2 (nan onde H ~ 0)
    ricci_scalar: np.ndarray           # R = 6 (a''/a + H^2 + k/a^2)
    constraint_violation: np.ndarray   # |adot^2 - H2(a) a^2| / max_t(adot^2) : verificacao de conservacao (vinculo)
    k: float
    model: str
    parameters: dict = field(default_factory=dict)
    solver_metadata: dict = field(default_factory=dict)

    @property
    def i_bounce(self) -> Optional[int]:
        i = int(np.argmin(self.a))
        return i if 0 < i < self.a.size - 1 else None

    @property
    def t_bounce(self) -> Optional[float]:
        i = self.i_bounce
        return None if i is None else float(self.t[i])

    @property
    def a_min(self) -> float:
        return float(self.a.min())

    @property
    def rho_max(self) -> float:
        return float(self.rho.max())

    @property
    def ricci_max(self) -> float:
        return float(np.nanmax(np.abs(self.ricci_scalar)))

    def efolds(self, t0: Optional[float] = None, t1: Optional[float] = None) -> float:
        """N = int H dt entre t0 e t1 (padrao: do ricochete ate o fim)."""
        t0 = self.t_bounce if t0 is None else t0
        t1 = self.t[-1] if t1 is None else t1
        if t0 is None:
            return float("nan")
        m = (self.t >= t0) & (self.t <= t1)
        return float(np.trapezoid(self.H[m], self.t[m]))

    def accelerated_phases(self):
        """Intervalos (t_ini, t_fim, N) em que a'' > 0 (expansao acelerada), depois do ricochete."""
        acc = self.accel > 0
        out = []
        i0 = None
        for i in range(self.a.size):
            if acc[i] and i0 is None:
                i0 = i
            if (not acc[i] or i == self.a.size - 1) and i0 is not None:
                i1 = i
                if self.i_bounce is None or i1 > self.i_bounce:
                    N = float(np.trapezoid(self.H[i0:i1 + 1], self.t[i0:i1 + 1]))
                    out.append((float(self.t[i0]), float(self.t[i1]), N))
                i0 = None
        return out


@dataclass
class FriedmannModel:
    fluids: Sequence[Fluid] = (Fluid(1.0, 0.0),)
    k: float = 0.0
    correction: object = field(default_factory=GRCorrection)

    def __post_init__(self):
        if isinstance(self.correction, LQCCorrection) and self.k != 0:
            raise ValueError("LQC efetiva implementada so para k = 0 (ver quantum/lqc.py: limitacoes).")

    # --- conteudo -------------------------------------------------------------------------
    def rho(self, a):
        a = np.asarray(a, float)
        return sum(f.rho0 * a ** (-3 * (1 + f.w)) for f in self.fluids)

    def p(self, a):
        a = np.asarray(a, float)
        return sum(f.w * f.rho0 * a ** (-3 * (1 + f.w)) for f in self.fluids)

    def rho_eff(self, a):
        return self.correction.rho_eff(a, self.rho(a), self.p(a))

    def H2(self, a):
        a = np.asarray(a, float)
        return self.correction.H2_contribution(a, self.rho(a), self.p(a)) - self.k / a**2

    def accel(self, a):
        return self.correction.accel(a, self.rho(a), self.p(a))

    @property
    def name(self) -> str:
        return f"Friedmann[{self.correction.name}, k={self.k:g}]"

    def parameters(self) -> dict:
        d = {f"rho0_{f.name}": f.rho0 for f in self.fluids}
        d.update({f"w_{f.name}": f.w for f in self.fluids})
        d["k"] = self.k
        for key in ("sigma", "rho_c"):
            if hasattr(self.correction, key):
                d[key] = getattr(self.correction, key)
        return d

    # --- ricochete --------------------------------------------------------------------------
    def turning_points(self, a_lo=1e-6, a_hi=1e4, n=4000):
        """Raizes de H^2(a) = 0 (ricochete e/ou recolapso) por varredura + brentq."""
        grid = np.logspace(np.log10(a_lo), np.log10(a_hi), n)
        with np.errstate(all="ignore"):
            h = self.H2(grid)
        roots = []
        for i in range(n - 1):
            if np.isfinite(h[i]) and np.isfinite(h[i + 1]) and h[i] * h[i + 1] < 0:
                roots.append(brentq(self.H2, grid[i], grid[i + 1]))
        return roots

    def bounce_scale(self) -> Optional[float]:
        r = self.turning_points()
        return min(r) if r else None

    # --- evolucao ---------------------------------------------------------------------------
    def rhs(self, t, y):
        a, adot = y
        return [adot, self.accel(a) * a]

    def solve(self, a0=1.0, t_span=(0.0, 6.0), contracting=True, n_samples=4000,
              settings: ODESettings = DEFAULT_ODE, stop_at_turnaround=False) -> BackgroundSolution:
        H2 = float(self.H2(a0))
        if H2 < 0:
            raise ValueError(f"H^2(a0={a0}) < 0: a0 esta numa regiao proibida (alem de um ponto de retorno).")
        adot0 = (-1 if contracting else 1) * np.sqrt(H2) * a0
        events = None
        if stop_at_turnaround:
            a_b = self.bounce_scale() or 0.0

            def turnaround(t, y):
                return y[1] if y[0] > 2 * a_b else 1.0
            turnaround.terminal = True
            turnaround.direction = -1
            events = [turnaround]
        sol = integrate(self.rhs, t_span, [a0, adot0], n_samples=n_samples, settings=settings, events=events)
        a, adot = sol.y
        H = adot / a
        rho, p = self.rho(a), self.p(a)
        acc = self.accel(a)
        with np.errstate(all="ignore"):
            Hdot = acc - H**2
            eps = np.where(np.abs(H) > 1e-8 * np.nanmax(np.abs(H)), -Hdot / H**2, np.nan)
            H2a = self.H2(a) * a**2
            # violacao do vinculo normalizada pela escala GLOBAL de adot^2 (perto do ricochete
            # adot -> 0 e uma normalizacao local divergiria sem significado fisico)
            viol = np.abs(adot**2 - H2a) / max(float(np.max(adot**2)), 1e-300)
        ricci = 6 * (acc + H**2 + self.k / a**2)
        return BackgroundSolution(t=sol.t, a=a, adot=adot, H=H, rho=rho, p=p, rho_eff=self.rho_eff(a), accel=acc,
                                  epsilon=eps, ricci_scalar=ricci, constraint_violation=viol, k=self.k,
                                  model=self.name, parameters=self.parameters(), solver_metadata=sol.as_metadata())


# ----------------------------------------------------------------------------
# Fabricas
# ----------------------------------------------------------------------------
def gr_dust(rho_m0=1.0, k=0.0) -> FriedmannModel:
    return FriedmannModel((Fluid(rho_m0, 0.0, "dust"),), k, GRCorrection())


def einstein_cartan_dust(rho_m0=1.0, sigma=1e-3, rho_r0=0.0, k=0.0) -> FriedmannModel:
    fl = [Fluid(rho_m0, 0.0, "dust")]
    if rho_r0:
        fl.append(Fluid(rho_r0, 1 / 3, "radiation"))
    return FriedmannModel(tuple(fl), k, EinsteinCartanCorrection(sigma))


def lqc_fluid(rho0=1.0, w=0.0, rho_c=None) -> FriedmannModel:
    from ..quantum.lqc import RHO_C_PLANCK
    return FriedmannModel((Fluid(rho0, w, "fluid"),), 0.0, LQCCorrection(RHO_C_PLANCK if rho_c is None else rho_c))


# ----------------------------------------------------------------------------
# Solucao exata (validacao)
# ----------------------------------------------------------------------------
def exact_torsion_dust(t, rho_m0=1.0, sigma=1e-3, a0=1.0):
    """RESULTADO MATEMATICO: poeira + torcao, k = 0.  Em w = a^3:
        w'^2 = 24 pi (rho_m0 w - sigma)  =>  w'' = 12 pi rho_m0,
        a(t)^3 = sigma/rho_m0 + 6 pi rho_m0 (t - t_b)^2,   t_b fixado por a(0) = a0 em contracao."""
    t = np.asarray(t, float)
    w0 = a0**3
    tb = np.sqrt((w0 - sigma / rho_m0) / (6 * np.pi * rho_m0))
    w = sigma / rho_m0 + 6 * np.pi * rho_m0 * (t - tb) ** 2
    return w ** (1 / 3), tb


def classical_dust_collapse(rho_m0=1.0, a0=1.0, n=4000):
    """RESULTADO MATEMATICO: poeira em GR, k = 0: a = a0 (1 - t/t_s)^{2/3}, t_s = 2/(3 H0)."""
    H0 = np.sqrt(8 * np.pi * rho_m0 / 3) * a0 ** (-1.5)
    ts = 2 / (3 * H0)
    t = np.linspace(0, ts * 0.9999, n)
    return dict(t=t, a=a0 * (1 - t / ts) ** (2 / 3), t_singularity=ts)
