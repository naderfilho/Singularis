"""
Ringdown no dominio do tempo: evolucao 1+1 da equacao mestre
    -d^2 Psi/dt^2 + d^2 Psi/dr_*^2 - V(r_*) Psi = 0
por diferencas finitas (leapfrog de 2a ordem) com condicoes de saida (Sommerfeld) nas bordas.

Classificacao: RESULTADO NUMERICO sobre a metrica escolhida (campo de teste de spin 0 ou 1; para vacuo,
spin 2 axial).  Serve para:
  * medir a frequencia e o amortecimento do ringdown diretamente da forma de onda (validacao do WKB);
  * obter ECOS quando o potencial tem duas barreiras (buraco de minhoca de Simpson-Visser, l > 2M):
    o atraso entre ecos e ~ 2 x distancia em r_* entre os picos (Cardoso, Franzin & Pani 2016);
  * gerar espectrogramas e comparar com o molde classico de Schwarzschild.

Unidades geometricas (M = 1); conversao para Hz e ms via core.units.Scale.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import find_peaks, stft

from ..core.units import Scale
from ..stability.perturbations import StaticMetric, master_potential, tortoise_grid


@dataclass
class Waveform:
    t: np.ndarray
    psi: np.ndarray
    metric: str
    spin: int
    ell: int
    r_obs_star: float
    potential_peaks_rs: list
    echo_delay_predicted: float | None
    fitted_omega: complex | None = None
    echoes: list = field(default_factory=list)
    t_direct: float = 0.0            # chegada do pulso direto ao observador
    t_ringdown_start: float = 0.0    # chegada do pulso refletido na barreira (inicio do ringdown)

    def to_si(self, M_solar: float):
        s = Scale.from_solar_masses(M_solar)
        return self.t * s.time_unit, s

    def spectrogram(self, nperseg=256):
        f, tt, Z = stft(self.psi, fs=1.0 / (self.t[1] - self.t[0]), nperseg=nperseg)
        return f, tt, np.abs(Z)


def evolve(metric: StaticMetric, spin=0, ell=2, r_range=None, n_grid=6000, t_max=400.0, courant=0.5,
           pulse_center=None, pulse_width=2.0, r_obs=None) -> Waveform:
    """Integra a equacao mestre com um pulso gaussiano inicial (Psi = 0, dPsi/dt = gaussiana) e registra
    Psi no observador r_obs (em r_*)."""
    if r_range is None:
        if metric.horizons:
            rh = max(metric.horizons)
            r_range = (rh * (1 + 1e-5), rh + 120 * metric.M)
        else:
            r_range = (-120 * metric.M, 120 * metric.M)
    r, rs = tortoise_grid(metric, r_range[0], r_range[1], n_grid)
    V = master_potential(metric, r, spin, ell)
    # grade uniforme em r_*
    x = np.linspace(rs[0], rs[-1], n_grid)
    Vx = np.interp(x, rs, V)
    dx = x[1] - x[0]
    dt = courant * dx
    nt = int(t_max / dt)
    peaks, _ = find_peaks(Vx, height=0.2 * Vx.max())
    peak_rs = [float(x[p]) for p in peaks]
    echo = (2 * abs(peak_rs[-1] - peak_rs[0])) if len(peak_rs) >= 2 else None
    x_peak = x[int(np.argmax(Vx))]
    pc = x_peak + 15.0 if pulse_center is None else pulse_center
    xo = x_peak + 40.0 if r_obs is None else r_obs
    io = int(np.argmin(np.abs(x - xo)))
    psi_prev = np.zeros(n_grid)
    dpsi0 = np.exp(-((x - pc) ** 2) / (2 * pulse_width**2))
    psi = psi_prev + dt * dpsi0  # primeiro passo (Psi(0) = 0)
    out = np.zeros(nt)
    lap = np.zeros(n_grid)
    for i in range(nt):
        lap[1:-1] = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) / dx**2
        psi_next = 2 * psi - psi_prev + dt**2 * (lap - Vx * psi)
        # Sommerfeld: Psi_t = -+ Psi_x nas bordas (onda de saida)
        psi_next[0] = psi[1] + (dt - dx) / (dt + dx) * (psi_next[1] - psi[0])
        psi_next[-1] = psi[-2] + (dt - dx) / (dt + dx) * (psi_next[-2] - psi[-1])
        psi_prev, psi = psi, psi_next
        out[i] = psi[io]
    t = dt * np.arange(1, nt + 1)
    wf = Waveform(t, out, metric.name, spin, ell, float(x[io]), peak_rs, echo)
    wf.t_direct = float(abs(x[io] - pc))
    wf.t_ringdown_start = float(abs(x[io] - pc) + 2 * abs(pc - x_peak))
    wf.fitted_omega = fit_ringdown(wf)
    wf.echoes = detect_echoes(wf)
    return wf


def fit_ringdown(wf: Waveform, t_start_offset=25.0, t_window=40.0) -> complex | None:
    """Ajusta A e^{-gamma t} cos(w t + phi) ao ringdown por minimos quadrados nao lineares (curve_fit),
    com chutes iniciais dos cruzamentos de zero e da envoltoria.

    Janela: comeca 25 M depois da chegada do pulso refletido (a resposta imediata e os sobretons
    contaminam os primeiros ~20 M; estudo de convergencia em tests/test_gw.py: com 25 M o modo
    fundamental de Regge-Wheeler e recuperado a 0.1-0.6% em frequencia e 1-2% em amortecimento) e
    dura 40 M (antes de a cauda de lei de potencia dominar).  Devolve omega = w - i gamma."""
    from scipy.optimize import curve_fit
    t, y = wf.t, wf.psi
    # janela: depois da chegada do pulso refletido (o pulso DIRETO nao e ringdown)
    t0 = wf.t_ringdown_start + t_start_offset
    m = (t > t0) & (t < t0 + t_window)
    if m.sum() < 50:
        return None
    tt, yy = t[m] - t[m][0], y[m]
    zc = np.where(np.diff(np.sign(yy)) != 0)[0]
    if zc.size < 4:
        return None
    w0 = 2 * np.pi / (2 * np.mean(np.diff(tt[zc])))
    pk, _ = find_peaks(np.abs(yy))
    g0 = -np.polyfit(tt[pk], np.log(np.abs(yy[pk])), 1)[0] if pk.size >= 3 else 0.1
    A0 = float(np.max(np.abs(yy)))

    def model(x, A, g, w, ph):
        return A * np.exp(-g * x) * np.cos(w * x + ph)

    try:
        popt, _ = curve_fit(model, tt, yy, p0=[A0, abs(g0), w0, 0.0], maxfev=20000)
    except Exception:
        return complex(w0, -abs(g0))
    A, g, w, ph = popt
    return complex(abs(w), -abs(g))


def detect_echoes(wf: Waveform) -> list:
    """Ecos = trens de pulsos secundarios na envoltoria |Psi|, SO definidos para potenciais com duas barreiras
    (cavidade).  Para uma barreira unica devolve [] por definicao (a reflexao do pulso inicial na barreira
    nao e eco).  A janela de busca comeca meio atraso previsto depois do pico principal e a separacao minima
    entre ecos e metade do atraso previsto  (delay = 2 x distancia em r_* entre os picos)."""
    if wf.echo_delay_predicted is None or len(wf.potential_peaks_rs) < 2:
        return []
    env = np.abs(wf.psi)
    dt = wf.t[1] - wf.t[0]
    i_pk = int(np.argmax(env))
    delay = wf.echo_delay_predicted
    thresh = 0.01 * env[i_pk]
    pk, _ = find_peaks(env, height=thresh, distance=max(1, int(0.5 * delay / dt)))
    out = []
    for p in pk:
        if wf.t[p] > wf.t_ringdown_start + 0.5 * delay:
            out.append((float(wf.t[p]), float(env[p] / env[i_pk])))
    grouped = []
    for tp, a in out:
        if not grouped or tp - grouped[-1][0] > 0.5 * delay:
            grouped.append((tp, a))
    return grouped


def classical_template(M_solar: float, chi: float = 0.0) -> dict:
    """Frequencia e tempo de amortecimento do modo l = m = 2 para um buraco negro classico de massa
    M_solar (Schwarzschild via Leaver; Kerr via ajuste BCW), em Hz e ms."""
    from .qnm import LEAVER_SCHWARZSCHILD, kerr_bcw
    om = LEAVER_SCHWARZSCHILD[(2, 2, 0)] if chi == 0 else kerr_bcw(chi).omega
    s = Scale.from_solar_masses(M_solar)
    f = s.frequency_to_hz(om.real)
    tau_ms = (1.0 / (-om.imag)) * s.time_unit * 1e3
    return dict(M_solar=M_solar, chi=chi, omega_M=om, f_Hz=f, tau_ms=tau_ms, quality=om.real / (2 * -om.imag))
