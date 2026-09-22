"""
Perturbacoes cosmologicas lineares sobre QUALQUER BackgroundSolution (GR, Einstein-Cartan, LQC).

Equacao de Mukhanov-Sasaki para um campo de teste sem massa / modos tensoriais (mesma equacao):

    v'' + (k^2 - a''/a) v = 0,     ' = d/d eta,     dt = a d eta,
    a''/a = a^2 (H^2 + addot/a)  (calculado do fundo numerico).

Espectros (unidades G = 1, M_Pl^2 = 1/8 pi):
    tensor:   P_t(k)   = 64 pi (k^3 / 2 pi^2) |v_k / a|^2                      (2 polarizacoes, h = 2 sqrt(8 pi) v/a)
    escalar:  P_zeta(k) = (4 pi / epsilon_exit) (k^3 / 2 pi^2) |v_k / a|^2     com z = a sqrt(2 eps) M_Pl
    razao:    r(k) = P_t / P_zeta = 16 epsilon(t_exit(k))                      (Wands 1999; Brandenberger & Peter 2017)

HIPOTESE explicitamente marcada: a perturbacao de curvatura zeta e avaliada com epsilon no instante de saida
do horizonte NA CONTRACAO e assume-se que a amplitude escalar atravessa o ricochete da mesma forma que o
campo de teste (zeta e mal definida em H = 0 porque z ~ a sqrt(2 eps) diverge).  Em contracoes de poeira
esta e a hipotese padrao do "matter bounce"; e ela que da r = 16 * 3/2 = 24.

Nao-gaussianidade: NAO calculada aqui (exige a acao cubica); o valor da literatura para o matter bounce
f_NL^local = -35/8 (Cai, Xue, Brandenberger & Zhang 2009) e reportado em observations/bridge.py como valor
DE LITERATURA, nao deste codigo.

Extracao do modo constante (super-horizonte) depois do ricochete: phi = v/a e amostrado enquanto
k eta_pos < 0.2 e |d ln phi / d ln eta| < tol; o valor no ultimo instante que satisfaz isso e o modo
constante, e a variacao residual e reportada como incerteza numerica.  Para o fundo de poeira exato,
cosmology/perturbations.py faz a extracao analitica (usada como referencia nos testes).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from .friedmann import BackgroundSolution


@dataclass
class SpectrumResult:
    k: np.ndarray
    P_test: np.ndarray            # k^3 |phi_const|^2 / 2 pi^2  (campo de teste; forma = tensor)
    P_tensor: np.ndarray
    P_scalar: np.ndarray
    epsilon_exit: np.ndarray
    r: np.ndarray
    n_s: float
    alpha_s: float
    n_t: float
    k_bounce: float
    extraction_uncertainty: np.ndarray   # variacao relativa de phi no plato (numerica)
    plateau_range: tuple
    background: str
    assumptions: list = field(default_factory=lambda: [
        "campo de teste sem massa (= modos tensoriais)",
        "zeta avaliada com epsilon na saida do horizonte na contracao e transportada pelo ricochete como o campo de teste",
        "vacuo de Bunch-Davies no passado remoto da contracao",
    ])


def bounce_background(model, a0=2.0e4, n_per_side=120000, t_end_factor=1.0) -> BackgroundSolution:
    """Fundo com ricochete amostrado numa grade LOGARITMICA em torno de t_b (resolve o ricochete e alcanca
    |k eta_0| >> 1 para os modos de interesse).  Duas passagens: uma grossa para localizar t_b, outra fina."""
    from ..core.solvers import integrate
    # estimativa do tempo de queda (poeira/radiacao): t_b ~ 2/(3(1+w)) / H0
    H0 = float(np.sqrt(model.H2(a0)))
    w_eff = float(model.p(a0) / model.rho(a0))
    t_b_est = 2.0 / (3.0 * (1.0 + w_eff) * H0)
    y0 = [a0, -H0 * a0]

    def turn(t, y):     # adot cruza zero subindo: instante exato do ricochete (evento, nao amostragem)
        return y[1]
    turn.direction = 1
    turn.terminal = False
    first = integrate(model.rhs, (0.0, 1.6 * t_b_est), y0, n_samples=2000, events=[turn])
    ev = first.events.get("turn", [])
    if not ev:
        raise RuntimeError("nenhum ricochete encontrado no intervalo estimado")
    t_b = float(ev[0])
    t_end = t_b + t_end_factor * t_b
    left = t_b - np.logspace(np.log10(t_b), -8, n_per_side)
    right = t_b + np.logspace(-8, np.log10(t_end - t_b), n_per_side)
    t_eval = np.unique(np.concatenate([[0.0], left, [t_b], right]))
    sol = integrate(model.rhs, (0.0, t_end), y0, t_eval=t_eval)
    a, adot = sol.y
    H = adot / a
    rho, p = model.rho(a), model.p(a)
    acc = model.accel(a)
    with np.errstate(all="ignore"):
        Hdot = acc - H**2
        eps = np.where(np.abs(H) > 1e-8 * np.nanmax(np.abs(H)), -Hdot / H**2, np.nan)
        viol = np.abs(adot**2 - model.H2(a) * a**2) / max(float(np.max(adot**2)), 1e-300)
    ricci = 6 * (acc + H**2 + model.k / a**2)
    return BackgroundSolution(t=sol.t, a=a, adot=adot, H=H, rho=rho, p=p, rho_eff=model.rho_eff(a), accel=acc, epsilon=eps,
                              ricci_scalar=ricci, constraint_violation=viol, k=model.k, model=model.name,
                              parameters=model.parameters(), solver_metadata=sol.as_metadata())


class ModeSolver:
    """Resolve modos sobre um BackgroundSolution amostrado em t (deve cobrir contracao e expansao;
    use `bounce_background` para uma grade que resolve o ricochete e comeca bem dentro do horizonte)."""

    def __init__(self, sol: BackgroundSolution):
        self.sol = sol
        t, a = sol.t, sol.a
        eta = cumulative_trapezoid(1 / a, t, initial=0.0)
        ib = sol.i_bounce if sol.i_bounce is not None else int(np.argmin(a))
        eta -= eta[ib]
        self.eta, self.a, self.t = eta, a, t
        self.app = a**2 * (sol.H**2 + sol.accel)
        self.k_bounce = float(np.sqrt(np.nanmax(self.app)))
        self.ib = ib
        self.eps = sol.epsilon
        self.H = sol.H

    def _app(self, e):
        return np.interp(e, self.eta, self.app)

    def _a(self, e):
        return np.interp(e, self.eta, self.a)

    def _early_offset(self):
        """Ponto singular conformal da lei de potencia da CONTRACAO (a = a1 (eta_s - eta)^p), para o vacuo inicial."""
        m = self.eta < 0
        e, a = self.eta[m], self.a[m]
        sl = slice(0, int(0.5 * e.size))
        eps_early = float(np.nanmedian(self.eps[: int(0.2 * self.eps.size)]))
        p = 1.0 / (eps_early - 1.0)
        c1, c0 = np.polyfit(e[sl], a[sl] ** (1 / p), 1)
        return float(-c0 / c1)

    def mode(self, k, eta_end=None):
        eta0 = self.eta[0]
        x0 = k * (eta0 - self._early_offset())
        # Bunch-Davies com a correcao de primeira ordem em 1/x (exata para poeira; adequada para |x0| >> 1),
        # com eta medido a partir do ponto singular da lei de potencia da contracao
        v0 = np.exp(-1j * x0) * (1 - 1j / x0) / np.sqrt(2 * k)
        dv0 = k * np.exp(-1j * x0) * (-1j * (1 - 1j / x0) + 1j / x0**2) / np.sqrt(2 * k)
        eta_end = self.eta[-1] if eta_end is None else eta_end
        y0 = [v0.real, v0.imag, dv0.real, dv0.imag]

        def rhs(e, y):
            w2 = k * k - self._app(e)
            return [y[2], y[3], -w2 * y[0], -w2 * y[1]]

        sol = solve_ivp(rhs, (eta0, eta_end), y0, method="DOP853", rtol=1e-9, atol=1e-12,
                        max_step=min(0.2 / k, 0.02), dense_output=True)
        return sol

    def late_power(self):
        """Expoente p de a ~ eta^p na expansao tardia (poeira: 2, radiacao: 1), ajustado no ultimo terco."""
        m = self.eta > 0
        e, a = self.eta[m], self.a[m]
        sl = slice(int(0.66 * e.size), e.size)
        p = np.polyfit(np.log(e[sl]), np.log(a[sl]), 1)[0]
        return float(p)

    def _late_powerlaw(self):
        """a = a1 (eta - eta_s)^p no regime tardio: ajusta p, a1 e o deslocamento eta_s."""
        m = self.eta > 0
        e, a = self.eta[m], self.a[m]
        sl = slice(int(0.5 * e.size), e.size)
        # p vem do epsilon tardio do fundo (a ~ t^{1/eps}  =>  a ~ eta^p, p = 1/(eps - 1)); poeira: 2, radiacao: 1
        eps_late = float(np.nanmedian(self.eps[int(0.8 * self.eps.size):]))
        p = 1.0 / (eps_late - 1.0)
        # com p fixo, a^(1/p) e linear em eta: ajuste linear da a1 e o deslocamento eta_s exatamente
        y = a[sl] ** (1 / p)
        c1, c0 = np.polyfit(e[sl], y, 1)
        eta_s = -c0 / c1
        a1 = float(c1**p)
        return float(p), a1, float(eta_s)

    def constant_mode(self, k):
        """phi_const = modo constante de phi = v/a fora do horizonte depois do ricochete.

        No regime tardio a = a1 (eta - eta_s)^p e as solucoes exatas sao  v = sqrt(x) [c1 J_nu(x) + c2 Y_nu(x)],
        x = k (eta - eta_s), nu = p - 1/2 (poeira: nu = 3/2, radiacao: nu = 1/2).  (c1, c2) sao obtidos de
        (v, v') no fim da grade; o modo constante e o ramo J:
            phi_const = c1 k^p 2^(-nu) / (Gamma(nu + 1) a1).
        A incerteza reportada e a diferenca relativa entre as extracoes em dois instantes tardios distintos.
        """
        from scipy.special import gamma, jv, jvp, yv, yvp
        s = self.mode(k)
        p, a1, eta_s = self._late_powerlaw()
        nu = p - 0.5
        vals = []
        for frac in (1.0, 0.7):
            e = self.eta[0] + frac * (self.eta[-1] - self.eta[0])
            y = s.sol(e)
            v = y[0] + 1j * y[1]
            dv = y[2] + 1j * y[3]
            x = k * (e - eta_s)
            sx = np.sqrt(x)
            J, Y = jv(nu, x), yv(nu, x)
            dJ = k * (J / (2 * sx) + sx * jvp(nu, x))
            dY = k * (Y / (2 * sx) + sx * yvp(nu, x))
            A = np.array([[sx * J, sx * Y], [dJ, dY]], dtype=complex)
            c1, c2 = np.linalg.solve(A, np.array([v, dv]))
            vals.append(abs(c1) * k**p * 2 ** (-nu) / (gamma(nu + 1) * a1))
        unc = abs(vals[0] - vals[1]) / max(vals[0], 1e-300)
        return float(vals[0]), float(unc)

    def epsilon_at_exit(self, k):
        """epsilon no instante (na contracao) em que k = a|H| (saida do horizonte)."""
        m = np.arange(self.ib)
        aH = self.a[m] * np.abs(self.H[m])
        # saida do horizonte = PRIMEIRO cruzamento de a|H| = k para cima (a|H| cresce na contracao e
        # volta a zero so no ricochete; o cruzamento descendente perto do ricochete nao e a saida)
        up = np.where((aH[:-1] < k) & (aH[1:] >= k))[0]
        if up.size == 0:
            return float("nan")
        return float(self.eps[up[0]])

    def spectrum(self, ks, plateau_fraction=0.1) -> SpectrumResult:
        """Espectros para os k dados.  O plato (ajuste de n_s, alpha_s) usa k <= plateau_fraction * k_bounce;
        a extracao tem um sistematico de ~3% em P perto de 0.1 k_bounce (validado contra a referencia exata)."""
        ks = np.asarray(ks, float)
        phis, uncs, eps = [], [], []
        for k in ks:
            p, u = self.constant_mode(k)
            phis.append(p)
            uncs.append(u)
            eps.append(self.epsilon_at_exit(k))
        phis, uncs, eps = map(np.array, (phis, uncs, eps))
        P_test = ks**3 * phis**2 / (2 * np.pi**2)
        P_t = 64 * np.pi * P_test
        with np.errstate(all="ignore"):
            P_s = (4 * np.pi / eps) * P_test
            r = P_t / P_s
        m = ks <= plateau_fraction * self.k_bounce
        lk, lP = np.log(ks[m]), np.log(P_test[m])
        c2 = np.polyfit(lk, lP, 2) if m.sum() >= 4 else [0, np.polyfit(lk, lP, 1)[0], 0]
        slope = np.polyfit(lk, lP, 1)[0]
        n_s = 1 + slope
        alpha = 2 * c2[0]
        return SpectrumResult(k=ks, P_test=P_test, P_tensor=P_t, P_scalar=P_s, epsilon_exit=eps, r=r, n_s=float(n_s),
                              alpha_s=float(alpha), n_t=float(slope), k_bounce=self.k_bounce, extraction_uncertainty=uncs,
                              plateau_range=(float(ks[m].min()), float(ks[m].max())), background=self.sol.model)
