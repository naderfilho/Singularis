"""
Colapso gravitacional DINAMICO em simetria esferica (poeira), com dinamica de camada
substituivel: relatividade geral classica ou modelos regularizados com ricochete.

Formulacao (Lemaitre–Tolman–Bondi, coordenadas comoveis (tau, chi), G = c = 1):

    ds^2 = -dtau^2 + R'^2/(1 + 2E(chi)) dchi^2 + R(tau,chi)^2 dOmega^2,      ' = d/dchi

    classico (RG):  Rdot^2 = 2 m(chi)/R + 2 E(chi),   Rddot = -m/R^2                [Lemaitre 1933; Tolman 1934; Bondi 1947]
    densidade:      rho = m' / (4 pi R^2 R')
    massa de Misner-Sharp (definicao geometrica, vale em qualquer modelo):
                    m_MS = (R/2) (1 - grad R . grad R) = (R/2) (1 + Rdot^2 - (1 + 2E))
    aprisionamento: grad R . grad R = 1 + 2E - Rdot^2 <= 0   (superficie marginalmente presa quando = 0)

Dinamica de camada MODIFICADA (ricochete) — MODELO/HIPOTESE, marcado explicitamente:

    Rdot^2 = 2 m/R + 2 E - 3 m^2 / (2 pi rho_* R^4)                                  (BounceShellDynamics)

que e a equacao de Friedmann modificada  H^2 = (8 pi/3) rho_bar (1 - rho_bar/rho_*)  escrita para cada
camada com a densidade MEDIA interior rho_bar = 3m/(4 pi R^3).  Para o caso homogeneo (R = a chi) ela
reproduz exatamente:
  * LQC efetiva com rho_* = rho_c  (Ashtekar–Pawlowski–Singh 2006), e
  * Einstein–Cartan com poeira de fermions, rho_* = rho_b = m_f^2/alpha  (Poplawski 2010),
que para POEIRA sao a mesma familia de um parametro (resultado matematico verificado em testes).
Para camadas nao homogeneas, a extensao "camada a camada com a densidade media" segue a ideia de
Kelly, Santacruz & Wilson-Ewing (2020) para o caso marginalmente ligado (E = 0); com E != 0 e uma
EXTRAPOLACAO deste codigo, e o resultado deve ser lido como tal.

Referencias
-----------
Oppenheimer & Snyder (1939) Phys. Rev. 56, 455.   Misner & Sharp (1964) Phys. Rev. 136, B571.
Kelly, Santacruz & Wilson-Ewing (2020) Phys. Rev. D 102, 106024.   Bojowald, Harada & Tibrewala (2008) PRD 78, 064057.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from ..core.solvers import DEFAULT_ODE, ODESettings, integrate


# ----------------------------------------------------------------------------
# Perfis de materia (m(chi), E(chi))
# ----------------------------------------------------------------------------
@dataclass
class DustProfile:
    """Perfil de poeira: massa de Misner-Sharp conservada m(chi) e energia E(chi) por camada."""

    chi: np.ndarray
    m: np.ndarray
    E: np.ndarray
    R0: np.ndarray          # raio areal inicial de cada camada
    Rdot0: np.ndarray       # velocidade inicial
    name: str = "perfil"
    m_prime: Optional[np.ndarray] = None   # dm/dchi analitico (evita erro de diferencas finitas em rho)

    @property
    def M(self) -> float:
        return float(self.m[-1])

    @classmethod
    def homogeneous(cls, M=1.0, R0=8.0, n_shells=60, bound=True, name="Oppenheimer-Snyder"):
        """Estrela homogenea em repouso (OS). bound=True: E = -m/R0 (fechada, como OS); False: marginalmente ligada
        (E = 0) mas entao a estrela nao pode comecar em repouso — comeca com Rdot0 = -sqrt(2m/R0)."""
        chi = np.linspace(0.02, 1.0, n_shells)
        R = R0 * chi
        m = M * chi**3
        if bound:
            E = -m / R                       # Rdot0 = 0
            Rdot0 = np.zeros_like(R)
        else:
            E = np.zeros_like(R)
            Rdot0 = -np.sqrt(2 * m / R)
        return cls(chi, m, E, R, Rdot0, name=name, m_prime=3 * M * chi**2)

    @classmethod
    def inhomogeneous(cls, M=1.0, R0=8.0, n_shells=80, core_fraction=0.5, contrast=3.0, name="nucleo denso"):
        """Densidade inicial rho(chi) = rho_c [1 + contrast * exp(-(chi/core)^2)], em repouso (E = -m/R0)."""
        chi = np.linspace(0.02, 1.0, n_shells)
        R = R0 * chi
        w = 1 + contrast * np.exp(-((chi / core_fraction) ** 2))
        # m(chi) = 4 pi int rho R^2 dR  com rho propto w  -> normaliza para m(1) = M
        integrand = w * R**2
        m = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(R))])
        norm = M / m[-1]
        m = np.maximum(m * norm, 1e-9 * M)
        mp = integrand * R0 * norm          # dm/dchi = 4 pi rho R^2 dR/dchi  (analitico, mesma normalizacao)
        E = -m / R
        return cls(chi, m, E, R, np.zeros_like(R), name=name, m_prime=mp)


# ----------------------------------------------------------------------------
# Dinamica de camada
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class ShellDynamics:
    """Rdot^2 = F(R; m, E).  Rddot = (1/2) dF/dR.  Classico: F = 2m/R + 2E."""

    name: str = "RG classica"

    def F(self, R, m, E):
        return 2 * m / R + 2 * E

    def accel(self, R, m, E):
        return -m / R**2

    def is_regular(self) -> bool:
        return False


@dataclass(frozen=True)
class BounceShellDynamics(ShellDynamics):
    """F = 2m/R + 2E - 3 m^2 / (2 pi rho_* R^4).  rho_* = rho_c (LQC) ou rho_b (Einstein-Cartan, poeira)."""

    rho_star: float = field(default=1.0, kw_only=True)
    name: str = "ricochete efetivo (rho_*)"

    def F(self, R, m, E):
        return 2 * m / R + 2 * E - 3 * m**2 / (2 * np.pi * self.rho_star * R**4)

    def accel(self, R, m, E):
        return -m / R**2 + 3 * m**2 / (np.pi * self.rho_star * R**5)

    def is_regular(self) -> bool:
        return True

    def mean_density(self, R, m):
        return 3 * m / (4 * np.pi * R**3)


# ----------------------------------------------------------------------------
# Resultado
# ----------------------------------------------------------------------------
@dataclass
class CollapseResult:
    tau: np.ndarray                 # (n_t,)
    chi: np.ndarray                 # (n_chi,)
    R: np.ndarray                   # (n_t, n_chi)
    Rdot: np.ndarray
    rho: np.ndarray                 # densidade local de poeira m'/(4 pi R^2 R')
    rho_mean: np.ndarray            # densidade media interior 3m/(4 pi R^3)
    m_MS: np.ndarray                # massa de Misner-Sharp geometrica (R/2)(1 - grad R.grad R)
    trapped: np.ndarray             # bool: grad R . grad R <= 0
    R_prime: np.ndarray             # dR/dchi (shell crossing quando <= 0)
    ricci_scalar: np.ndarray        # R_Ricci = 8 pi rho_local (poeira, RG) ou 8 pi (rho_eff - 3 p_eff) no modelo efetivo
    kretschmann_homog: Optional[np.ndarray]  # so para perfil homogeneo: K de FRW (exato)
    apparent_horizon_R: np.ndarray  # raio areal do horizonte aparente externo (nan se nao ha regiao presa)
    surface_R: np.ndarray
    model: str
    profile: str
    parameters: dict
    solver_metadata: dict
    shell_crossing_tau: Optional[float] = None
    bounce_tau: np.ndarray = field(default_factory=lambda: np.array([]))   # por camada (nan se nao ricocheteia)
    R_min: np.ndarray = field(default_factory=lambda: np.array([]))         # por camada
    event_horizon: dict = field(default_factory=dict)

    @property
    def M(self) -> float:
        return float(self.m_MS[0, -1])

    def summary(self) -> dict:
        i_s = -1
        s = dict(model=self.model, profile=self.profile, M=self.M,
                 tau_end=float(self.tau[-1]),
                 rho_max=float(np.nanmax(self.rho)), rho_mean_max=float(np.nanmax(self.rho_mean)),
                 ricci_max=float(np.nanmax(np.abs(self.ricci_scalar))),
                 kretschmann_max=None if self.kretschmann_homog is None else float(np.nanmax(self.kretschmann_homog)),
                 R_surface_min=float(np.nanmin(self.surface_R)),
                 trapped_ever=bool(self.trapped.any()),
                 tau_first_trapped=float(self.tau[np.argmax(self.trapped.any(axis=1))]) if self.trapped.any() else None,
                 tau_last_trapped=float(self.tau[self.trapped.any(axis=1)][-1]) if self.trapped.any() else None,
                 shell_crossing_tau=self.shell_crossing_tau,
                 bounce_tau_surface=float(self.bounce_tau[i_s]) if self.bounce_tau.size else None,
                 bounce_tau_center=float(self.bounce_tau[0]) if self.bounce_tau.size else None,
                 mass_deficit_max=float(np.nanmax(1 - self.m_MS[:, -1] / self.M)),
                 event_horizon=self.event_horizon)
        return s


# ----------------------------------------------------------------------------
# Simulador
# ----------------------------------------------------------------------------
@dataclass
class SphericalCollapse:
    profile: DustProfile
    dynamics: ShellDynamics = field(default_factory=ShellDynamics)
    settings: ODESettings = DEFAULT_ODE

    def __post_init__(self):
        # Condicoes iniciais CONSISTENTES com o vinculo Rdot^2 = F(R): a energia E(chi) e ajustada para que
        # Rdot0^2 = F(R0) exatamente (na dinamica de ricochete, E = -m/R0 classico violaria o vinculo por
        # 3 m^2/(2 pi rho_* R0^4) e a trajetoria de 2a ordem nao seria a solucao do modelo).
        p = self.profile
        F0 = self.dynamics.F(p.R0, p.m, p.E)
        delta = (p.Rdot0**2 - F0) / 2.0
        if np.max(np.abs(delta)) > 0:
            self.profile = DustProfile(p.chi, p.m, p.E + delta, p.R0, p.Rdot0, name=p.name, m_prime=p.m_prime)
        self.constraint_adjustment = float(np.max(np.abs(delta)))

    def rhs(self, tau, y):
        n = self.profile.chi.size
        R, Rd = y[:n], y[n:]
        return np.concatenate([Rd, self.dynamics.accel(R, self.profile.m, self.profile.E)])

    def run(self, tau_max: Optional[float] = None, n_samples=1200, stop_at_R=1e-3, stop_after_reexpansion=True) -> CollapseResult:
        p = self.profile
        n = p.chi.size
        if tau_max is None:
            # tempo classico de queda da superficie a partir do repouso: (pi/2) sqrt(R0^3/2M)
            tau_max = 3.0 * (np.pi / 2) * np.sqrt(p.R0[-1] ** 3 / (2 * p.M))
        y0 = np.concatenate([p.R0, p.Rdot0])
        events = []

        def singular(tau, y):   # alguma camada chegou a R ~ 0 (RG classica)
            return float(np.min(y[:n]) - stop_at_R * p.R0[-1])
        singular.terminal = True
        singular.direction = -1
        events.append(singular)
        if stop_after_reexpansion and self.dynamics.is_regular():
            def reexpanded(tau, y):   # superficie voltou ao raio inicial depois do ricochete
                # continua em tau: negativa no inicio (R_s = R0), positiva durante o colapso, volta a
                # negativa (direcao -1) quando a superficie re-expande alem de R0 (1 - 1e-3)
                return float(p.R0[-1] * (1 - 1e-3) - y[n - 1])
            reexpanded.terminal = True
            reexpanded.direction = -1
            events.append(reexpanded)
        sol = integrate(self.rhs, (0.0, tau_max), y0, n_samples=n_samples, settings=self.settings, events=events)
        tau = sol.t
        R = sol.y[:n].T
        Rd = sol.y[n:].T
        Rp = np.gradient(R, p.chi, axis=1)
        mp = np.gradient(p.m, p.chi) if p.m_prime is None else p.m_prime
        with np.errstate(all="ignore"):
            rho = mp / (4 * np.pi * R**2 * Rp)
            rho_mean = 3 * p.m / (4 * np.pi * R**3)
            grad2 = 1 + 2 * p.E - Rd**2
            m_MS = (R / 2) * (1 - grad2)
        trapped = grad2 <= 0
        # curvatura
        if isinstance(self.dynamics, BounceShellDynamics):
            rs = self.dynamics.rho_star
            rho_eff = rho * (1 - rho / rs)
            p_eff = -rho**2 / rs   # p_eff da familia efetiva para poeira (a''/a = -(4pi/3)(rho_eff + 3 p_eff))
            ricci = 8 * np.pi * (rho_eff - 3 * p_eff)
        else:
            ricci = 8 * np.pi * rho
        kret = None
        if p.name.startswith("Oppenheimer") or np.allclose(p.m / p.chi**3, p.m[-1], rtol=1e-6):
            # homogeneo: a = R/chi,  K_FRW = 12 [ (addot/a)^2 + (adot^2/a^2 + k/a^2)^2 ],  k/a^2 = -2E/R^2 (E = -k chi^2/2)
            a = R / p.chi
            H2k = (Rd / R) ** 2 - 2 * p.E / R**2
            acc = self.dynamics.accel(R, p.m, p.E) / R
            kret = 12 * (acc**2 + H2k**2)
        # horizonte aparente externo: camada mais externa presa em cada tau
        ah = np.full(tau.size, np.nan)
        for i in range(tau.size):
            idx = np.where(trapped[i])[0]
            if idx.size:
                ah[i] = R[i, idx.max()]
        # ricochete por camada
        b_tau = np.full(n, np.nan)
        R_min = R.min(axis=0)
        for j in range(n):
            i = int(np.argmin(R[:, j]))
            if 0 < i < tau.size - 1 and Rd[i - 1, j] < 0 < Rd[i + 1, j]:
                b_tau[j] = tau[i]
        # shell crossing
        sc = None
        bad = np.where((Rp <= 0).any(axis=1))[0]
        if bad.size:
            sc = float(tau[bad[0]])
        res = CollapseResult(tau=tau, chi=p.chi, R=R, Rdot=Rd, rho=rho, rho_mean=rho_mean, m_MS=m_MS, trapped=trapped,
                             R_prime=Rp, ricci_scalar=ricci, kretschmann_homog=kret, apparent_horizon_R=ah,
                             surface_R=R[:, -1], model=self.dynamics.name, profile=p.name,
                             parameters=dict(M=p.M, R0=float(p.R0[-1]), n_shells=n,
                                             rho_star=getattr(self.dynamics, "rho_star", None)),
                             solver_metadata=sol.as_metadata(), shell_crossing_tau=sc, bounce_tau=b_tau, R_min=R_min)
        res.event_horizon = self.event_horizon(res)
        return res

    # --- horizonte de eventos por lancamento de raios nulos --------------------------------------
    def event_horizon(self, res: CollapseResult, n_rays=48) -> dict:
        """Raios nulos radiais de saida  dchi/dtau = sqrt(1+2E)/R'  lancados do centro em tau_0.

        Classificacao de cada raio (so com a informacao do INTERIOR simulado):
          escaped                  chega a superficie com R_s > 2M (exterior nao preso): escapa.
          entered_trapped_exterior chega a superficie com R_s <= 2M.  Em RG classica (exterior Schwarzschild)
                                   esta condenado; num modelo com ricochete o exterior NAO e Schwarzschild
                                   perto do ricochete (ver collapse/junction.py) e o destino fica indeterminado.
          unresolved               nao chega a superficie dentro do tempo simulado.
          shell_crossing           atravessa uma regiao com R' <= 0 (metrica degenerada).
        O horizonte de eventos classico e a fronteira escaped -> entered_trapped_exterior (bissecao).
        """
        p = self.profile
        tau, chi = res.tau, res.chi
        M = p.M
        Rp_grid = res.R_prime
        E = p.E

        def ray(t0):
            c = chi[0]
            i0 = int(np.searchsorted(tau, t0))
            for i in range(i0, tau.size - 1):
                dt = tau[i + 1] - tau[i]
                Rp = np.interp(c, chi, Rp_grid[i])
                if Rp <= 0:
                    return "shell_crossing", None
                c = c + np.sqrt(1 + 2 * np.interp(c, chi, E)) / Rp * dt
                if c >= chi[-1]:
                    Rs = float(np.interp(tau[i + 1], tau, res.surface_R))
                    return ("escaped" if Rs > 2 * M else "entered_trapped_exterior"), float(tau[i + 1])
            return "unresolved", None

        t0s = np.linspace(tau[0], tau[-1] * 0.98, n_rays)
        status = [ray(t0)[0] for t0 in t0s]
        info = dict(rays_launch_tau=t0s.tolist(), status=status)
        esc = np.array([s == "escaped" for s in status])
        trap = np.array([s == "entered_trapped_exterior" for s in status])
        if trap.any() and esc.any() and np.where(esc)[0].min() < np.where(trap)[0].min():
            j = int(np.where(trap)[0].min()) - 1
            lo, hi = t0s[j], t0s[j + 1]
            for _ in range(25):
                mid = 0.5 * (lo + hi)
                if ray(mid)[0] == "escaped":
                    lo = mid
                else:
                    hi = mid
            info.update(classical_event_horizon_exists=True, tau_birth_center=float(lo))
            if self.dynamics.is_regular():
                info["note"] = ("raios lancados apos tau_birth_center chegam a superficie na fase presa; com ricochete o "
                                "exterior nao e Schwarzschild ai, logo a existencia de um horizonte de EVENTOS depende da "
                                "geometria exterior alem do simulado (so o horizonte APARENTE e transiente e esta resolvido)")
            else:
                info["note"] = "horizonte de eventos nasce no centro em tau_birth_center (exterior Schwarzschild)"
        elif esc.all():
            info.update(classical_event_horizon_exists=False, note="todos os raios resolvidos escapam")
        else:
            info.update(classical_event_horizon_exists=None, note="sem fronteira escaped->trapped resolvida no tempo simulado")
        # raios que voltam a escapar depois do ricochete: lancados apos o fim da fase presa
        if self.dynamics.is_regular() and res.trapped.any():
            t_last = float(tau[res.trapped.any(axis=1)][-1])
            late = [s for t0, s in zip(t0s, status) if t0 > t_last]
            info["rays_after_trapped_phase"] = dict(n=len(late), escaped=int(sum(s == "escaped" for s in late)),
                                                    unresolved=int(sum(s == "unresolved" for s in late)))
        return info


# ----------------------------------------------------------------------------
# Comparacao classico vs ricochete (interface unica)
# ----------------------------------------------------------------------------
def compare_classical_vs_bounce(M=1.0, R0=8.0, rho_star=None, Rb_over_R0=0.05, n_shells=60, profile="homogeneous",
                                settings: ODESettings = DEFAULT_ODE):
    """Roda o mesmo perfil inicial com RG classica e com a dinamica de ricochete.

    rho_star pode ser dado diretamente ou via o raio de ricochete desejado da superficie (Rb_over_R0):
    para a superficie homogenea, rho_bar(R_b) = rho_*  ->  rho_* = 3M / (4 pi R_b^3) (com E ~ 0);
    com E != 0 o ricochete ocorre onde 2m/R + 2E = 3m^2/(2 pi rho_* R^4).
    """
    if profile == "homogeneous":
        prof = DustProfile.homogeneous(M, R0, n_shells)
    else:
        prof = DustProfile.inhomogeneous(M, R0, n_shells)
    if rho_star is None:
        Rb = Rb_over_R0 * R0
        rho_star = 3 * M / (4 * np.pi * Rb**3)
    classical = SphericalCollapse(prof, ShellDynamics(), settings).run()
    bounce = SphericalCollapse(prof, BounceShellDynamics(rho_star=rho_star), settings).run()
    return classical, bounce
