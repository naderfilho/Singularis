"""Testes das geometrias com carga e rotacao (Fase 7): limites analiticos e regularidade."""
import numpy as np

from semente.geometry.charged import (charged_black_bounce, charged_bounce_phase, reissner_nordstrom,
                                      rn_horizons_analytic, rn_kretschmann_analytic)
from semente.geometry.kerr import Kerr
from semente.geometry.spherical import BlackBounce
from semente.geometry.static import StaticSphericalMetric


def test_reissner_nordstrom_horizons_and_kretschmann():
    rn = reissner_nordstrom(1.0, 0.5)
    assert np.allclose(rn.horizons(), rn_horizons_analytic(1.0, 0.5), atol=1e-9)
    r = np.array([2.5, 4.0, 9.0])
    assert np.allclose(rn.kretschmann(r), rn_kretschmann_analytic(r, 1.0, 0.5), rtol=1e-10)
    assert rn.classify()["kind"].startswith("buraco negro com horizonte de Cauchy")
    assert reissner_nordstrom(1.0, 1.0).classify()["n_horizons"] in (0, 1)   # extremal (raiz dupla)
    assert reissner_nordstrom(1.0, 1.3).classify()["kind"] == "singularidade nua"


def test_charged_black_bounce_limits_and_phases():
    # Q -> 0 recupera Simpson-Visser; l -> 0 recupera Reissner-Nordstrom
    cb = charged_black_bounce(1.0, 1e-9, 0.5)
    sv = BlackBounce(1.0, 0.5)
    r = np.array([-3.0, 0.0, 2.5, 6.0])
    assert np.allclose(cb.kretschmann(r), sv.kretschmann(r), rtol=1e-6)
    cb0 = charged_black_bounce(1.0, 0.5, 1e-9)
    rr = np.array([2.5, 4.0])
    assert np.allclose(cb0.kretschmann(rr), rn_kretschmann_analytic(rr, 1.0, 0.5), rtol=1e-6)
    # fases em (Q, l)
    assert charged_bounce_phase(1.0, 0.5, 0.05).startswith("regular BH com horizonte interno")
    assert charged_bounce_phase(1.0, 0.5, 0.5).startswith("black-bounce")
    assert charged_bounce_phase(1.0, 0.5, 2.5).startswith("buraco de minhoca")
    assert charged_black_bounce(1.0, 0.5, 0.05).classify()["n_horizons"] == 4
    assert charged_black_bounce(1.0, 0.5, 0.5).classify()["n_horizons"] == 2
    assert np.isfinite(charged_black_bounce(1.0, 0.5, 0.5).kretschmann(0.0))


def test_kerr_analytic_structure():
    k = Kerr(1.0, 0.6)
    assert np.allclose(k.kerr_horizons, (0.2, 1.8))
    assert abs(k.ergosurface(np.pi / 2) - 2.0) < 1e-12 and abs(k.ergosurface(0.0) - 1.8) < 1e-12
    assert abs(Kerr(1.0, 0.0).isco() - 6.0) < 1e-9
    assert abs(Kerr(1.0, 1 - 1e-12).isco() - 1.0) < 1e-3 and abs(Kerr(1.0, 1 - 1e-12).isco(False) - 9.0) < 1e-3
    assert abs(k.Omega_H - 0.6 / (1.8**2 + 0.36)) < 1e-12
    assert abs(k.horizon_area - 4 * np.pi * (1.8**2 + 0.36)) < 1e-9
    assert abs(k.frame_dragging(1e6)) < 1e-15   # arrasto -> 0 no infinito
    assert Kerr(1.0, 1.2).classify()["kind"].startswith("super-extremal")


def test_rotating_black_bounce_regular_and_reduces_to_kerr():
    """Kretschmann simbolico (custa ~50 s na primeira chamada; cacheado no processo)."""
    k = Kerr(1.0, 0.6)
    r = np.array([2.3, 3.0, 5.0, 9.0])
    assert np.allclose(k.kretschmann(r), k.kretschmann_kerr_analytic(r), rtol=1e-8)
    kb = Kerr(1.0, 0.6, 0.7)
    assert np.allclose(kb.horizons, (np.sqrt(1.8**2 - 0.49),))
    K0 = kb.kretschmann(np.array([0.0, 0.3, 1.0]))
    assert np.all(np.isfinite(K0)) and np.all(np.abs(K0) < 1e4)    # regular na garganta (Mazza et al. 2021)
    kb_small = Kerr(1.0, 0.6, 1e-4)
    assert np.allclose(kb_small.kretschmann(r), k.kretschmann_kerr_analytic(r), rtol=1e-4)   # l -> 0


def test_static_family_reproduces_schwarzschild_photon_sphere():
    s = StaticSphericalMetric("1 - 2*M/r", "r", ("M",), dict(M=1.0), r_domain=(0, np.inf))
    ps = s.photon_spheres()
    assert any(abs(p - 3.0) < 2e-3 for p in ps)
    assert abs(s.surface_gravity(2.0) - 0.25) < 1e-12
