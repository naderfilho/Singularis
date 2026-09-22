"""Testes da Fase 4: Einstein-Cartan e LQC como estrategias unificadas."""
import numpy as np

from semente.core.units import M_NEUTRON, M_PLANCK, M_SUN
from semente.cosmology.seed_universe import SeedUniverse
from semente.quantum.comparison import compare_models, dust_equivalence_residual, radiation_difference
from semente.quantum.einstein_cartan import ALPHA_SPIN_HALF, bounce_density_planck, sigma_from_matter
from semente.quantum.lqc import RHO_C_PLANCK


def test_einstein_cartan_bounce_density_formula():
    m = M_NEUTRON / M_PLANCK
    assert abs(bounce_density_planck(M_NEUTRON) - 4 * m**2 / np.pi) < 1e-60
    assert abs(sigma_from_matter(1.0, m) * m**2 - ALPHA_SPIN_HALF) < 1e-12  # sigma = alpha (rho/m)^2


def test_lqc_critical_density_value():
    assert abs(RHO_C_PLANCK - 0.41) < 0.01


def test_dust_equivalence_ec_lqc_is_exact():
    h, acc = dust_equivalence_residual(rho0=1.0, rho_star=7.0)
    assert h < 1e-12 and acc < 1e-12


def test_radiation_breaks_equivalence():
    d = radiation_difference(rho0=1.0, rho_star=7.0)
    assert abs(d["a_bounce_EC"] - np.sqrt(1 / 7)) < 1e-9
    assert abs(d["a_bounce_LQC"] - (1 / 7) ** 0.25) < 1e-9


def test_compare_models_equivalent_initial_conditions():
    res = {c.name: c for c in compare_models(rho0=1.0, w=0.0, rho_star=20.0, t_span=(0, 1.5))}
    assert not res["GR"].bounce
    assert res["Einstein-Cartan"].bounce and res["LQC efetiva"].bounce
    assert abs(res["Einstein-Cartan"].a_min - res["LQC efetiva"].a_min) < 1e-9
    assert abs(res["Einstein-Cartan"].rho_max - 20.0) < 1e-3
    assert all(c.constraint_violation_max < 1e-8 for c in res.values())
    rad = {c.name: c for c in compare_models(rho0=1.0, w=1 / 3, rho_star=20.0, t_span=(0, 1.5))}
    assert rad["Einstein-Cartan"].rho_max > 10 * rad["LQC efetiva"].rho_max  # radiacao: EC deixa rho passar de rho_*


def test_seed_universe_mechanisms():
    ec = SeedUniverse(10 * M_SUN, mechanism="einstein_cartan")
    lq = SeedUniverse(10 * M_SUN, mechanism="lqc")
    assert lq.rho_bounce_pl > 1e30 * ec.rho_bounce_pl            # LQC ricocheteia em densidade de Planck
    assert lq.R_bounce_pl < ec.R_bounce_pl                        # e portanto num raio muito menor
    assert abs(ec.tau_horizon_to_bounce_s - lq.tau_horizon_to_bounce_s) / ec.tau_horizon_to_bounce_s < 1e-6  # tempo ~ pi M
