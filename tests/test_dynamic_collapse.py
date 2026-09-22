"""Testes do colapso dinamico: limites analiticos (OS), conservacao, ricochete, horizontes."""
import numpy as np

from semente.collapse.dynamic import (BounceShellDynamics, DustProfile, ShellDynamics, SphericalCollapse,
                                      compare_classical_vs_bounce)
from semente.collapse.oppenheimer_snyder import OppenheimerSnyder
from semente.cosmology.friedmann import einstein_cartan_dust, lqc_fluid


def test_homogeneous_classical_reproduces_oppenheimer_snyder():
    prof = DustProfile.homogeneous(1.0, 8.0, 40)
    res = SphericalCollapse(prof).run()
    os_ = OppenheimerSnyder(1.0, 8.0)
    # superficie: R(tau) exato de OS
    eta = os_.eta_from_tau(res.tau)
    R_exact = os_.R_surface(eta)
    m = res.tau < 0.98 * os_.tau_s
    assert np.max(np.abs(res.surface_R[m] - R_exact[m]) / R_exact[m]) < 2e-3
    # horizonte de eventos nasce no centro no instante previsto (tolerancia da grade de raios)
    assert abs(res.event_horizon["tau_birth_center"] - os_.event_horizon_interior()["tau_birth"]) < 0.1
    # massa de Misner-Sharp conservada camada a camada em RG (poeira)
    assert np.max(np.abs(res.m_MS - prof.m)) < 1e-6


def test_homogeneous_profile_is_homogeneous_in_density():
    prof = DustProfile.homogeneous(1.0, 8.0, 40)
    res = SphericalCollapse(prof).run()
    i = res.tau.size // 2
    rho = res.rho[i, 3:-3]
    assert np.std(rho) / np.mean(rho) < 1e-3


def test_bounce_matches_target_radius_and_is_transiently_trapped():
    c, b = compare_classical_vs_bounce(M=1.0, R0=8.0, Rb_over_R0=0.05, n_shells=40)
    s = b.summary()
    assert abs(s["R_surface_min"] - 0.4) < 0.02            # ricochete perto de R_b (E != 0 desloca ligeiramente)
    assert s["trapped_ever"] and s["tau_last_trapped"] < b.tau[-1]   # regiao presa TRANSIENTE
    assert np.isfinite(s["bounce_tau_surface"]) and abs(s["bounce_tau_surface"] - s["bounce_tau_center"]) < 1e-6
    assert s["mass_deficit_max"] > 0.9                     # massa de Misner-Sharp na superficie cai ~R_b/R0 no ricochete
    assert c.summary()["trapped_ever"] and c.summary()["tau_last_trapped"] == c.tau[-1]  # classico: preso ate o fim


def test_regular_dynamics_recovers_classical_when_rho_star_to_infinity():
    prof = DustProfile.homogeneous(1.0, 8.0, 30)
    c = SphericalCollapse(prof, ShellDynamics()).run(tau_max=20.0)
    b = SphericalCollapse(prof, BounceShellDynamics(rho_star=1e12)).run(tau_max=20.0)
    assert np.max(np.abs(c.surface_R - b.surface_R)) < 1e-8


def test_shell_dynamics_equals_modified_friedmann_for_homogeneous_interior():
    """Camada homogenea R = a chi com rho_bar = 3m/(4 pi R^3): F/R^2 deve ser H^2 = (8pi/3) rho (1 - rho/rho_*) - k/a^2."""
    rho_star = 5.0
    dyn = BounceShellDynamics(rho_star=rho_star)
    a = np.array([0.3, 0.7, 1.5])
    chi = 1.0
    rho0 = 1.0
    m = (4 * np.pi / 3) * rho0 * chi**3
    E = 0.0
    R = a * chi
    H2_shell = dyn.F(R, m, E) / R**2
    H2_lqc = lqc_fluid(rho0, 0.0, rho_star).H2(a)
    H2_ec = einstein_cartan_dust(rho0, sigma=rho0**2 / rho_star).H2(a)
    assert np.allclose(H2_shell, H2_lqc, rtol=1e-12)
    assert np.allclose(H2_shell, H2_ec, rtol=1e-12)


def test_inhomogeneous_core_collapses_first_and_shell_crossing_detected():
    prof = DustProfile.inhomogeneous(1.0, 8.0, 40)
    res = SphericalCollapse(prof).run()
    assert res.shell_crossing_tau is not None
    assert res.R[-1, 0] < res.R[-1, -1] * 0.05  # o centro chega perto de R=0 antes da superficie


def test_initial_conditions_are_constraint_consistent():
    prof = DustProfile.homogeneous(1.0, 8.0, 20)
    sc = SphericalCollapse(prof, BounceShellDynamics(rho_star=3.0))
    p = sc.profile
    assert np.max(np.abs(sc.dynamics.F(p.R0, p.m, p.E) - p.Rdot0**2)) < 1e-14
