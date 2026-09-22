"""Testes de ondas gravitacionais: WKB vs Leaver, Kerr BCW, ringdown no dominio do tempo, ecos, templates."""
import numpy as np

from semente.gravitational_waves.qnm import LEAVER_SCHWARZSCHILD, kerr_bcw, wkb3
from semente.gravitational_waves.ringdown import classical_template, evolve
from semente.observations.bridge import compare
from semente.stability.perturbations import schwarzschild_metric, simpson_visser_metric


def test_wkb_matches_leaver_within_half_percent():
    m = schwarzschild_metric(1.0)
    for (spin, ell, n), ref in LEAVER_SCHWARZSCHILD.items():
        q = wkb3(m, spin=spin, ell=ell, n=n)
        assert abs(q.omega - ref) / abs(ref) < 5e-3


def test_black_bounce_qnm_reduces_to_schwarzschild_as_l_to_zero():
    q0 = wkb3(schwarzschild_metric(1.0), spin=0, ell=2)
    ql = wkb3(simpson_visser_metric(1.0, 1e-3), spin=0, ell=2)
    assert abs(ql.omega - q0.omega) < 1e-4
    # a garganta fixada pela juncao (l ~ 1e-14 M) nao muda o ringdown em nada mensuravel
    qj = wkb3(simpson_visser_metric(1.0, 1e-7), spin=0, ell=2)
    assert abs(qj.omega - q0.omega) / abs(q0.omega) < 1e-6


def test_kerr_fit_limits():
    q = kerr_bcw(0.0)
    assert abs(q.omega.real - 0.3737) < 0.01     # ajuste BCW reproduz Schwarzschild a ~1.5%
    assert kerr_bcw(0.9).omega.real > kerr_bcw(0.0).omega.real


def test_time_domain_ringdown_frequency_matches_wkb():
    m = schwarzschild_metric(1.0)
    wf = evolve(m, spin=0, ell=2, n_grid=4000, t_max=300.0)
    q = wkb3(m, spin=0, ell=2)
    assert wf.fitted_omega is not None
    assert abs(wf.fitted_omega.real - q.omega.real) / q.omega.real < 0.03
    assert abs(wf.fitted_omega.imag - q.omega.imag) / abs(q.omega.imag) < 0.15
    assert wf.echoes == []   # barreira unica: sem ecos por definicao


def test_wormhole_has_double_barrier_and_echoes():
    wf = evolve(simpson_visser_metric(1.0, 2.5), spin=0, ell=2, n_grid=4000, t_max=400.0)
    assert len(wf.potential_peaks_rs) == 2
    assert wf.echo_delay_predicted is not None and len(wf.echoes) >= 2
    gaps = np.diff([e[0] for e in wf.echoes[:3]])
    assert np.all(np.abs(gaps / wf.echo_delay_predicted - 1) < 0.5)   # espacamento ~ atraso previsto


def test_gw150914_bridge_shows_rotation_is_needed():
    t0 = classical_template(62.0, 0.0)
    t1 = classical_template(62.0, 0.67)
    c0 = compare("Schwarzschild (sem rotacao)", "GW150914_f_ringdown_Hz", t0["f_Hz"])
    c1 = compare("Kerr BCW (literatura)", "GW150914_f_ringdown_Hz", t1["f_Hz"], literature=True)
    assert c0.status == "excluded"                 # 195 Hz vs 251 +- 8
    assert c1.status.startswith("literature:") and abs(c1.residual_sigma) < 3
