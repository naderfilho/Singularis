"""
Nascimento: o que o modelo SEMENTE preve para um universo que nasce do
ricochete, confrontado com o que medimos no NOSSO universo.

Nada aqui confirma a hipotese.  O que fazemos e o teste honesto: cada
previsao do modelo vira uma linha com "previsto / observado / veredito".
Se alguma linha der "falsificado", a versao minima do modelo morre.

1) ESPECTRO PRIMORDIAL.  O interior de Oppenheimer-Snyder e uma contracao
   dominada por poeira.  Resolvemos numericamente a equacao de Mukhanov-Sasaki
   para um campo de teste sem massa atraves do ricochete de torcao,

        v'' + (k^2 - a''/a) v = 0      (' = d/d eta, tempo conforme),

   com vacuo de Bunch-Davies no passado remoto da contracao, e medimos o
   indice espectral n_s do espectro que sai do outro lado.  Resultado
   conhecido para contracao de poeira (dualidade de Wands 1999): n_s ~ 1,
   como o observado (0.9665 +- 0.0038).  Reproduzimos isso no nosso fundo
   especifico com ricochete NAO singular.

2) CURVATURA.  O filho e um universo FECHADO (k = +1).  Previsao de sinal:
   Omega_k < 0.  Dado o limite observado |Omega_k| < ~0.002, a esfera-3
   precisa ser enorme e isso fixa um piso para a massa do pai.

3) RAZAO TENSOR/ESCALAR.  O ricochete de materia minimo da r = 16 epsilon = 24
   (epsilon = 3/2 para poeira).  Observado: r < 0.036.  Tensao seria, conhecida.

4) IDADE vs RECOLAPSO.  Universo fechado de poeira recolapsa em t = pi a_m/2.
   Tem de ser >> 13.8 Gyr.

Dados usados (Planck 2018 VI, Tabela 2, TT,TE,EE+lowE+lensing+BAO; BICEP/Keck 2021):
   H0 = 67.66 +- 0.42 km/s/Mpc,  Omega_m = 0.3111 +- 0.0056,  Omega_k = 0.0007 +- 0.0019,
   n_s = 0.9665 +- 0.0038,  r < 0.036 (95%).   Planck sozinho (sem BAO): Omega_k = -0.011 +- 0.0065.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from .bounce import C_SI, G_SI, M_PLANCK, M_SUN

# --- dados observacionais -------------------------------------------------------------
PLANCK = dict(
    H0=67.66, H0_err=0.42, Omega_m=0.3111, Omega_m_err=0.0056,
    Omega_k=0.0007, Omega_k_err=0.0019,             # Planck+BAO
    Omega_k_planck_only=-0.011, Omega_k_planck_only_err=0.0065,
    n_s=0.9665, n_s_err=0.0038, r_max_95=0.036, age_Gyr=13.787,
)
MPC = 3.0857e22
GYR = 3.15576e16


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


# ----------------------------------------------------------------------------
# 2-4. Curvatura, massa do pai, idade
# ----------------------------------------------------------------------------
def curvature_constraints(Omega_k=None, Omega_k_err=None, n_sigma=2.0):
    """Do limite observado em Omega_k (fechado: Omega_k < 0) ao piso da massa do pai.

    Universo fechado de poeira:  a_m = (8 pi G / 3 c^2) rho_m0 a0^3,  a0 = c / (H0 sqrt(-Omega_k)).
    O universo observavel (raio comovel D) ocupa chi_obs = D / a0 da esfera-3; a estrela-pai
    cobria chi0 >= chi_obs, e  M_pai = (a_m / 2) sin^3(chi0) c^2/G  >=  (a_m/2) sin^3(chi_obs) c^2/G.
    """
    Ok = PLANCK["Omega_k"] if Omega_k is None else Omega_k
    Oke = PLANCK["Omega_k_err"] if Omega_k_err is None else Omega_k_err
    H0 = PLANCK["H0"] * 1e3 / MPC
    rho_c = 3 * H0**2 / (8 * np.pi * G_SI)
    rho_m0 = PLANCK["Omega_m"] * rho_c
    D_obs = 4.4e26  # m, raio comovel do universo observavel (~46.5 Gly)
    # o valor mais negativo permitido a n_sigma:
    Ok_min = Ok - n_sigma * Oke
    if Ok_min >= 0:
        return dict(compatible_closed=False)
    a0 = C_SI / (H0 * np.sqrt(-Ok_min))  # raio de curvatura MINIMO permitido
    a_m = (8 * np.pi * G_SI / (3 * C_SI**2)) * rho_m0 * a0**3  # em metros (raio maximo do universo fechado)
    chi_obs = D_obs / a0
    M_min_kg = (a_m / 2) * np.sin(chi_obs) ** 3 * C_SI**2 / G_SI
    t_recollapse = np.pi * a_m / 2 / C_SI  # s (poeira pura, sem Lambda)
    # massa da esfera-3 inteira (chi0 = pi/2, se a estrela cobrisse meio universo) para referencia
    M_half = (a_m / 2) * C_SI**2 / G_SI
    return dict(compatible_closed=True, Omega_k_used=Ok_min, a0_min_m=a0, a_m_m=a_m, chi_obs=chi_obs,
                M_parent_min_kg=M_min_kg, M_parent_min_Msol=M_min_kg / M_SUN,
                M_parent_halfsphere_Msol=M_half / M_SUN,
                t_recollapse_Gyr=t_recollapse / GYR, age_Gyr=PLANCK["age_Gyr"])


def matter_bounce_tensor_ratio():
    """Ricochete de materia minimo: r = 16 epsilon, epsilon = 3/2 (poeira)  =>  r = 24."""
    return 16 * 1.5


# ----------------------------------------------------------------------------
# Veredito
# ----------------------------------------------------------------------------
def verdict(ks=None, verbose=True):
    ks = np.logspace(-0.7, 1.5, 36) if ks is None else ks
    bs = BounceSpectrum(a_b=0.05)
    P = bs.spectrum(ks)
    ns = bs.spectral_index(ks, P)
    cc = curvature_constraints()
    cc_planck = curvature_constraints(PLANCK["Omega_k_planck_only"], PLANCK["Omega_k_planck_only_err"], n_sigma=0.0)
    r = matter_bounce_tensor_ratio()
    rows = [
        dict(item="indice espectral n_s (contracao de poeira + ricochete)",
             previsto=f"{ns:.3f} (numerico, plato k << k_b)", observado=f"{PLANCK['n_s']} +- {PLANCK['n_s_err']}",
             veredito="COMPATIVEL em 1a ordem: quase invariante de escala sem inflacao; "
                      "a inclinacao vermelha de ~3% exige correcao (ex.: w ligeiramente negativo)"),
        dict(item="sinal da curvatura espacial", previsto="Omega_k < 0 (fechado)",
             observado=f"Planck+BAO: {PLANCK['Omega_k']} +- {PLANCK['Omega_k_err']};  Planck so: "
                       f"{PLANCK['Omega_k_planck_only']} +- {PLANCK['Omega_k_planck_only_err']}",
             veredito="COMPATIVEL (Planck sozinho ate prefere fechado a ~1.7 sigma); FALSIFICAVEL se Omega_k > 0 for medido"),
        dict(item="massa minima do buraco negro pai (da curvatura + conservacao de massa)",
             previsto=f">= {cc['M_parent_min_Msol']:.1e} M_sol  (Omega_k >= {cc['Omega_k_used']:.4f}, 2 sigma)",
             observado=f"maior buraco negro conhecido ~ 1e10-1e11 M_sol; massa do universo observavel ~ 1e23 M_sol",
             veredito="CONSISTENTE mas exige um pai com a massa de um universo inteiro - nao e evidencia, e conservacao"),
        dict(item="razao tensor/escalar r", previsto=f"r = 16 epsilon = {r:g} (ricochete de materia minimo)",
             observado=f"r < {PLANCK['r_max_95']} (95%, BICEP/Keck 2021)",
             veredito="FALSIFICADO na versao minima; sobrevive so com fisica extra que suprima tensores "
                      "(ex.: ricochete com Lambda/curvatura, campos adicionais)"),
        dict(item="idade vs recolapso do universo fechado", previsto=f"recolapso em {cc['t_recollapse_Gyr']:.1e} Gyr (poeira pura)",
             observado=f"idade {PLANCK['age_Gyr']} Gyr; expansao ACELERADA (Lambda)",
             veredito="CONSISTENTE na idade; a aceleracao observada NAO sai do modelo (precisa de Lambda herdada ou emergente)"),
    ]
    out = dict(n_s_numerico=ns, k_bounce=bs.k_bounce, ks=ks, P=P, curvatura=cc, curvatura_planck_only=cc_planck, r_min=r, rows=rows)
    if verbose:
        for row in rows:
            print(f"- {row['item']}\n    previsto : {row['previsto']}\n    observado: {row['observado']}\n    veredito : {row['veredito']}")
    return out
