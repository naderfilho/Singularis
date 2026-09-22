"""Testes da Fase 8: termodinamica de horizontes e curvas de Page (exploratorio)."""
import numpy as np

from semente.collapse.dynamic import compare_classical_vs_bounce
from semente.geometry.charged import reissner_nordstrom
from semente.geometry.kerr import Kerr
from semente.geometry.static import StaticSphericalMetric
from semente.information.page_curve import evaporation_time_years, page_curves, timescales
from semente.thermodynamics.horizon import (apparent_horizon_entropy, baby_universe_budget, event_horizon_entropy,
                                            hawking_temperature_kelvin, kerr_thermo, static_thermo)


def _sv(l):
    return StaticSphericalMetric("1 - 2*M/sqrt(r**2 + l**2)", "sqrt(r**2 + l**2)", ("M", "l"), dict(M=1.0, l=l),
                                 name="SV")


def test_schwarzschild_thermodynamics():
    s = StaticSphericalMetric("1 - 2*M/r", "r", ("M",), dict(M=1.0), r_domain=(0, np.inf))
    th = static_thermo(s)
    assert abs(th.temperature - 1 / (8 * np.pi)) < 1e-9
    assert abs(th.entropy - 4 * np.pi) < 1e-9
    assert abs(th.first_law_residual) < 1e-3          # dM = T dS
    assert abs(hawking_temperature_kelvin(1.0) - 6.17e-8) / 6.17e-8 < 0.01


def test_reissner_nordstrom_first_law_at_fixed_charge():
    th = static_thermo(reissner_nordstrom(1.0, 0.5))
    rm, rp = 1 - np.sqrt(0.75), 1 + np.sqrt(0.75)
    assert abs(th.temperature - (rp - rm) / (4 * np.pi * rp**2)) < 1e-9
    assert abs(th.first_law_residual) < 1e-3          # dM = T dS a Q fixo


def test_simpson_visser_entropy_equals_schwarzschild_but_first_law_fails():
    for l in (0.5, 1.0, 1.5):
        th = static_thermo(_sv(l))
        assert abs(th.entropy - 4 * np.pi) < 1e-8
        assert abs(th.temperature - np.sqrt(4 - l**2) / (16 * np.pi)) < 1e-9
        assert abs(th.first_law_residual - (1 - np.sqrt(1 - l**2 / 4))) < 1e-3


def test_kerr_first_law():
    d = kerr_thermo(Kerr(1.0, 0.6))
    assert abs(d["entropy"] - np.pi * (1.8**2 + 0.36)) < 1e-9
    assert abs(d["first_law_residual_fixed_J"]) < 1e-3


def test_event_horizon_area_theorem_and_transient_apparent_horizon():
    c, b = compare_classical_vs_bounce(M=1.0, R0=8.0, Rb_over_R0=0.05, n_shells=30)
    se = event_horizon_entropy(c)
    assert se["defined"] and se["monotonic"]                       # teorema da area para o horizonte de EVENTOS
    assert abs(se["S_final_over_4piM2"] - 1.0) < 0.05              # o gerador chega a superficie em R = 2M
    sa = apparent_horizon_entropy(c)
    assert not sa["monotonic"]                                     # AH de OS se propaga para dentro (classico)
    sb = apparent_horizon_entropy(b)
    assert not sb["monotonic"] and sb["final_S"] == 0.0            # ricochete: regiao presa transiente


def test_page_curve_and_timescales():
    pc = page_curves(M0_planck=1e3, beta=1.5)
    assert abs(pc.t_page_over_tev - (1 - (0.6) ** 1.5)) < 1e-12
    assert np.all(pc.S_page <= pc.S_BH + 1e-9) and np.all(pc.S_page <= pc.S_hawking + 1e-9)
    assert pc.S_page[-1] == 0.0 and pc.S_baby[-1] > 0   # unitario volta a zero; universo-filho congela
    ts = timescales(10.0)
    assert 1e-4 < ts["t_bounce_s"] < 2e-4
    assert 1e69 < ts["t_evap_years"] < 1e71
    assert ts["t_page_over_t_bounce"] > 1e70
    b = baby_universe_budget(10.0, 1.67e-27)
    assert b["ratio_matter_over_horizon"] < 1e-18
