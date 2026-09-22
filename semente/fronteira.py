"""
Experimentos de fronteira: testes numericos sobre perguntas que a fisica ainda
nao resolveu, feitos DENTRO dos modelos que o projeto implementa.  Nenhum deles
resolve a pergunta; cada um extrai uma consequencia quantitativa que, ate onde
sabemos, nao estava calculada.

A. JUNCAO: o interior de Oppenheimer-Snyder com torcao ricocheteia em R_b.  O
   exterior pode continuar sendo Schwarzschild?  Resposta (calculada): NAO.
   Dentro do horizonte nenhuma superficie pode ter dR/dtau = 0, entao a juncao
   de Israel falha num intervalo finito em torno do ricochete, para qualquer
   camada de energia na superficie.  O unico exterior estatico esfericamente
   simetrico da familia black-bounce que admite a juncao e o de Simpson-Visser
   com  l = R_b  EXATAMENTE (a garganta tem de coincidir com o raio do ricochete).
   Consequencia: o parametro l deixa de ser livre:  l^3 = 3M / (4 pi rho_b).

B. ENTROPIA: se a entropia de Bekenstein-Hawking do pai (4 pi M^2) limita a
   entropia total do filho (hipotese holografica, em aberto), entao a massa do
   pai do NOSSO universo tem um piso, e uma arvore de universos ao estilo Smolin
   tem profundidade finita:  k_max = ln(M_0/m_Pl) / ln(sqrt(N)).

C. SELECAO COM ORCAMENTO: dinamica populacional de Smolin com o orcamento de
   entropia herdado.  Linhagens morrem em poucas geracoes.

D. DADOS: a previsao de Smolin (1992) para a massa maxima de estrelas de
   neutrons contra medidas de pulsares.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from .bounce import L_PLANCK, M_NEUTRON, M_PLANCK, M_SUN, RHO_PLANCK


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


# ----------------------------------------------------------------------------
# B. Limite entropico e profundidade da arvore
# ----------------------------------------------------------------------------
K_B_ENTROPIES = {
    # estimativas de Egan & Lineweaver (2010), ApJ 710, 1825, em unidades de k_B
    "CMB (fotons)": 2.03e89,
    "neutrinos cosmicos": 5.2e89,
    "buracos negros supermassivos": 3.1e104,
    "buracos negros estelares": 5.9e97,
}


def parent_mass_lower_bound(S_child_kB):
    """S_BH = 4 pi M^2 / m_Pl^2  >=  S_child   =>   M >= m_Pl sqrt(S_child / 4 pi).  Retorna kg e M_sol."""
    M = M_PLANCK * np.sqrt(S_child_kB / (4 * np.pi))
    return M, M / M_SUN


def genealogy_depth(M0_kg, N_per_generation, m_min_kg=M_PLANCK):
    """Sob o limite holografico, um filho so cabe se  sum_i m_i^2 <= M_pai^2.
    Com N buracos negros iguais por geracao: m_{k+1} = m_k / sqrt(N).
    A linhagem acaba quando m_k < m_min:  k_max = ln(M0/m_min) / ln(sqrt(N))."""
    return np.log(M0_kg / m_min_kg) / np.log(np.sqrt(N_per_generation))


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
        from .selection import fecundity_landscape
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


# ----------------------------------------------------------------------------
# D. Dados: massa maxima de estrelas de neutrons
# ----------------------------------------------------------------------------
NEUTRON_STAR_MASSES = {
    # (massa em M_sol, incerteza 1 sigma, referencia)
    "PSR J0740+6620": (2.08, 0.07, "Fonseca et al. 2021, ApJL 915, L12"),
    "PSR J0348+0432": (2.01, 0.04, "Antoniadis et al. 2013, Science 340, 448"),
    "PSR J1614-2230": (1.908, 0.016, "Arzoumanian et al. 2018, ApJS 235, 37"),
    "PSR J0952-0607": (2.35, 0.17, "Romani et al. 2022, ApJL 934, L17"),
}


def smolin_neutron_star_test(prediction=1.6, revised=2.0):
    """Smolin (1992) previu massa maxima de estrela de neutrons perto de 1.6 M_sol (mais tarde
    argumentou ~2 M_sol).  Quantos sigma cada pulsar medido esta acima de cada previsao?"""
    rows = {}
    for name, (m, dm, ref) in NEUTRON_STAR_MASSES.items():
        rows[name] = dict(massa=m, erro=dm, sigma_acima_1p6=(m - prediction) / dm, sigma_acima_2p0=(m - revised) / dm, ref=ref)
    return rows


def run_all(verbose=True):
    """Executa os quatro experimentos e devolve um dicionario de resultados."""
    out = {}
    # A
    tc = TorsionCollapse(M=1.0, R0=8.0, Rb_over_R0=0.05)
    sol = tc.solve()
    js = tc.junction_schwarzschild(sol)
    jl = tc.junction_black_bounce(sol, tc.l_predicted())
    ls, frac = tc.scan_l(sol)
    out["A_juncao"] = dict(
        R_b=tc.R_b, R_min_numerico=sol["R_min"],
        fracao_tempo_juncao_impossivel_schwarzschild=js["frac_time_invalid"],
        intervalo_falha_schwarzschild=js["tau_fail"],
        fracao_tempo_juncao_impossivel_black_bounce_l_eq_Rb=jl["frac_time_invalid"],
        massa_camada_max_sobre_M_l_eq_Rb=float(np.nanmax(np.abs(jl["m_shell"])) / tc.M),
        l_minimo_que_funciona_sobre_Rb=float(ls[frac <= 1e-3].min() / tc.R_b) if (frac <= 1e-3).any() else None,
        l_previsto_10Msol_m=l_from_black_hole_mass(10 * M_SUN),
        l_previsto_SgrA_m=l_from_black_hole_mass(4e6 * M_SUN),
        l_previsto_M87_m=l_from_black_hole_mass(6.5e9 * M_SUN),
    )
    # B
    S_total = sum(K_B_ENTROPIES.values())
    Mkg, Msol = parent_mass_lower_bound(S_total)
    Mkg_cmb, Msol_cmb = parent_mass_lower_bound(K_B_ENTROPIES["CMB (fotons)"])
    out["B_entropia"] = dict(
        S_filho_total_kB=S_total, M_pai_min_kg=Mkg, M_pai_min_Msol=Msol,
        M_pai_min_Msol_so_CMB=Msol_cmb,
        profundidade_arvore={f"M0=1e{int(np.log10(M0))} Msol, N=1e{int(np.log10(N))}": genealogy_depth(M0 * M_SUN, N)
                              for M0 in (10, 1e6, 1e10) for N in (1e2, 1e10, 1e20)},
    )
    # C
    hist = BudgetedSelection().run(40)
    out["C_selecao_orcamento"] = dict(geracoes_ate_extincao=len(hist), historico=hist)
    # D
    out["D_estrelas_neutrons"] = smolin_neutron_star_test()
    if verbose:
        import json
        print(json.dumps({k: v for k, v in out.items() if k != "C_selecao_orcamento"}, indent=2, default=float))
        print("C:", [(h["gen"], h["alive"], round(h["N_fit_mean"], 2)) for h in hist])
    return out
