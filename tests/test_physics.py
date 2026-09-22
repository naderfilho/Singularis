"""Testes de consistencia fisica (pytest -q)."""
import numpy as np
import pytest

from semente.bounce import LQCCosmology, SeedUniverse, TorsionCosmology, torsion_dust_exact, M_SUN
from semente.collapse import OppenheimerSnyder
from semente.geodesics import GeodesicSolver, deflection_angle
from semente.geometry import BlackBounce, Schwarzschild
from semente.interior import kantowski_sachs_black_bounce, proper_time_to_singularity, radial_infall_black_bounce
from semente.raytracer import newtonian_thin_disk_flux, page_thorne_flux


def test_kruskal_roundtrip():
    s = Schwarzschild(1.0)
    t = np.array([0.0, 1.5, -2.0, 3.0])
    r = np.array([3.0, 5.0, 1.0, 0.5])
    T, X = s.kruskal_from_tr(t, r)
    assert np.allclose(s.r_from_kruskal(T, X), r, atol=1e-9)
    assert np.allclose(s.t_from_kruskal(T, X), t, atol=1e-9)


def test_kruskal_regions_and_singularity():
    s = Schwarzschild(1.0)
    assert list(s.region(np.array([0, 2, -2, 0]), np.array([1, 0, 0, -1]))) == [1, 2, 3, 4]
    # singularidade: T^2 - X^2 = 1  ->  r = 0
    assert abs(s.r_from_kruskal(np.sqrt(2.0), 1.0)) < 1e-4


def test_kretschmann_black_bounce_reduces_to_schwarzschild():
    r = np.array([3.0, 5.0, 10.0])
    assert np.allclose(BlackBounce(1.0, 0.0).kretschmann(r), 48.0 / r**6)
    assert np.isfinite(BlackBounce(1.0, 0.5).kretschmann(np.array([0.0]))).all()


def test_circular_orbit_is_stable():
    g = GeodesicSolver("schwarzschild", 1.0)
    r0 = 10.0
    dphi = np.sqrt(1.0 / r0**3) / np.sqrt(1 - 3 / r0)
    v0 = g.initial_velocity([0, r0, np.pi / 2, 0], dr=0.0, dphi=dphi)
    sol = g.integrate([0, r0, np.pi / 2, 0], v0, 300.0)
    assert np.allclose(sol.y[1], r0, atol=1e-6)


def test_light_deflection_matches_weak_field_series():
    b = 25.0
    a = deflection_angle(b)
    series = 4 / b + 15 * np.pi / 4 / b**2 + 128 / 3 / b**3
    assert abs(a - series) < 2e-3


def test_proper_time_horizon_to_singularity_is_pi_M():
    assert abs(proper_time_to_singularity(2.0) - 2 * np.pi) < 1e-6


def test_black_bounce_interior_bounces_at_l():
    ks = kantowski_sachs_black_bounce(1.0, 0.7)
    assert abs(ks["a_perp"].min() - 0.7) < 1e-4
    assert np.isfinite(ks["kretschmann"]).all()


def test_infall_crosses_throat_into_other_universe():
    inf = radial_infall_black_bounce(1.0, 0.5, r0=8.0, tau_max=45.0)
    assert inf["r"].min() < -2.0


def test_oppenheimer_snyder_matching():
    os_ = OppenheimerSnyder(M=1.0, R0=10.0)
    assert abs(np.sin(os_.chi0) ** 2 - 2 * 1.0 / 10.0) < 1e-12
    assert abs(os_.R_surface(0.0) - 10.0) < 1e-12
    assert abs(os_.R_surface(os_.eta_h) - 2.0) < 1e-9
    sk = os_.surface_kruskal()
    # a superficie termina na singularidade T^2 - X^2 -> 1
    assert abs((sk["T"][-1] ** 2 - sk["X"][-1] ** 2) - 1.0) < 1e-3


def test_torsion_bounce_matches_exact_solution():
    tc = TorsionCosmology(rho_m0=1.0, sigma=0.05)
    s = tc.solve(a0=1.0, t_span=(0, 1.2))
    exact, _ = torsion_dust_exact(s["t"], 1.0, 0.05)
    assert np.max(np.abs(s["a"] - exact)) < 1e-8
    assert abs(s["a"].min() - (0.05) ** (1 / 3)) < 1e-6


def test_lqc_bounce_density():
    lq = LQCCosmology(rho0=0.01, w=0.0, rho_c=0.41)
    s = lq.solve(a0=1.0, t_span=(0, 3.0))
    assert abs(lq.rho(s["a"].min()) - 0.41) / 0.41 < 1e-3


def test_seed_universe_scales():
    su = SeedUniverse(10 * M_SUN)
    tab = su.table()
    assert tab["R_areal no ricochete / r_s"] < 1e-3  # o ricochete acontece muito dentro do horizonte
    assert tab["fator de expansao a_max/a_bounce"] > 1e3


def test_page_thorne_newtonian_limit():
    R = np.array([5000.0])
    assert abs(page_thorne_flux(R)[0] / newtonian_thin_disk_flux(R)[0] - 1) < 0.03


@pytest.mark.skipif(pytest.importorskip("torch", reason="torch nao instalado") is None, reason="torch")
def test_pinn_learns_bounce():
    from semente.pinn import BouncePINN

    p = BouncePINN(rho_m0=1.0, sigma=0.05, t_max=1.2)
    p.train(epochs=1200, verbose=False, lbfgs_steps=200)
    t = np.linspace(0, 1.2, 200)
    exact, _ = torsion_dust_exact(t, 1.0, 0.05)
    assert np.max(np.abs(p.predict(t) - exact)) < 5e-3
