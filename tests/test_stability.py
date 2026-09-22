"""Testes de condicoes de energia e estabilidade."""
import numpy as np

from semente.cosmology.friedmann import einstein_cartan_dust, gr_dust
from semente.stability.bounce import anisotropy_robustness, homogeneous_mode_is_gauge
from semente.stability.energy_conditions import (METRIC_LIBRARY, effective_fluid_functions, friedmann_report,
                                                 static_metric_report)
from semente.stability.perturbations import analyze, master_potential, schwarzschild_metric, simpson_visser_metric


def test_schwarzschild_is_vacuum():
    rep = static_metric_report("schwarzschild", np.linspace(2.1, 40, 200), M=1.0)
    assert np.max(np.abs(rep.rho)) < 1e-14 and np.max(np.abs(rep.p_r)) < 1e-14 and np.max(np.abs(rep.p_t)) < 1e-14
    assert all(v["satisfied_everywhere"] for v in rep.summary().values())


def test_reissner_nordstrom_effective_fluid_is_electromagnetic():
    """rho = -p_r = p_t = Q^2 / (8 pi r^4): resultado classico."""
    r = np.array([3.0, 5.0, 10.0])
    rep = static_metric_report("reissner_nordstrom", r, M=1.0, Q=0.5)
    expected = 0.25 / (8 * np.pi * r**4)
    assert np.allclose(rep.rho, expected, rtol=1e-10)
    assert np.allclose(rep.p_r, -expected, rtol=1e-10)
    assert np.allclose(rep.p_t, expected, rtol=1e-10)
    assert all(v["satisfied_everywhere"] for v in rep.summary().values())


def test_simpson_visser_reduces_to_vacuum_and_violates_nec_outside():
    fn, exprs = effective_fluid_functions(*METRIC_LIBRARY["simpson_visser"])
    r = np.linspace(-8, 8, 400)
    rho0, pr0, pt0 = fn(r, 1.0, 0.0)
    assert np.max(np.abs(np.broadcast_to(rho0, r.shape))) < 1e-14
    rep = static_metric_report("simpson_visser", r, M=1.0, l=0.5)
    outside = np.abs(r) > np.sqrt(4 - 0.25) * 1.01
    assert not rep.NEC[outside].any()            # NEC violada em todo o exterior (Simpson & Visser 2019)
    # forma fechada de rho: -l^2 (R - 4M) / (8 pi R^5)
    R = np.sqrt(r**2 + 0.25)
    assert np.allclose(rep.rho, -0.25 * (R - 4) / (8 * np.pi * R**5), rtol=1e-9)


def test_friedmann_bounce_violates_nec_only_near_bounce():
    sol = einstein_cartan_dust(1.0, 0.05).solve(1.0, (0, 1.2))
    rep = friedmann_report(sol)
    iv = rep.violation_intervals("NEC")
    assert len(iv) == 1 and iv[0][0] < sol.t_bounce < iv[0][1]
    gr = gr_dust(1.0).solve(1.0, (0, 0.2))
    assert friedmann_report(gr).NEC.all()


def test_master_potentials_and_schwarzschild_stability():
    m = schwarzschild_metric(1.0)
    r = np.array([3.0, 6.0, 10.0])
    V2 = master_potential(m, r, spin=2, ell=2)
    assert np.allclose(V2, (1 - 2 / r) * (6 / r**2 - 6 / r**3))
    for spin in (0, 1, 2):
        assert analyze(m, spin=spin, ell=2).classification == "stable"


def test_simpson_visser_test_fields_stable_and_converged():
    for l in (0.5, 2.5):
        res = analyze(simpson_visser_metric(1.0, l), spin=0, ell=2)
        assert res.classification == "stable" and res.convergence["converged"]
        assert res.potential_min > 0


def test_homogeneous_perturbation_is_gauge_and_anisotropy_criterion():
    m = einstein_cartan_dust(1.0, 0.05)
    sol = m.solve(1.0, (0, 1.2))
    assert homogeneous_mode_is_gauge(m, sol) < 1e-2
    rep = anisotropy_robustness(sol, rho_star=20.0, sigma0_sq=1e-6)
    assert rep.Sigma_bounce < 1e-3 and rep.verdict.startswith("ricochete isotropico robusto")
    bad = anisotropy_robustness(sol, rho_star=20.0, sigma0_sq=1.0)
    assert bad.Sigma_bounce > 1
