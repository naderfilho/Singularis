"""Testes do fundo de Friedmann unificado: limites classicos, solucao exata, conservacao."""
import numpy as np
import pytest

from semente.cosmology.friedmann import (FriedmannModel, Fluid, GRCorrection, classical_dust_collapse,
                                        einstein_cartan_dust, exact_torsion_dust, gr_dust, lqc_fluid)
from semente.quantum.lqc import LQCCorrection


def test_gr_dust_matches_analytic_collapse():
    m = gr_dust(1.0)
    ref = classical_dust_collapse(1.0)
    s = m.solve(a0=1.0, t_span=(0, 0.95 * ref["t_singularity"]), contracting=True, n_samples=500)
    a_ref = (1 - s.t / ref["t_singularity"]) ** (2 / 3)
    assert np.max(np.abs(s.a - a_ref)) < 1e-7


def test_einstein_cartan_matches_exact_and_conserves_constraint():
    m = einstein_cartan_dust(1.0, sigma=0.05)
    s = m.solve(a0=1.0, t_span=(0, 1.2))
    exact, tb = exact_torsion_dust(s.t, 1.0, 0.05)
    assert np.max(np.abs(s.a - exact)) < 1e-8
    assert abs(s.t_bounce - tb) < 2e-3
    assert np.nanmax(s.constraint_violation) < 1e-7      # vinculo de Friedmann preservado pela integracao
    assert abs(s.a_min - 0.05 ** (1 / 3)) < 1e-6


def test_sigma_to_zero_recovers_gr():
    t = np.linspace(0, 0.2, 200)
    a_gr = gr_dust(1.0).solve(a0=1.0, t_span=(0, 0.2), n_samples=200).a
    a_ec = einstein_cartan_dust(1.0, sigma=1e-12).solve(a0=1.0, t_span=(0, 0.2), n_samples=200).a
    assert np.max(np.abs(a_gr - a_ec)) < 1e-9


def test_lqc_bounce_at_rho_c_and_classical_limit():
    m = lqc_fluid(rho0=0.2, w=0.0, rho_c=0.41)
    s = m.solve(a0=1.0, t_span=(0, 1.2))
    assert abs(m.rho(s.a_min) - 0.41) / 0.41 < 1e-3
    # rho_c -> infinito: H^2 volta a (8pi/3) rho
    big = FriedmannModel((Fluid(0.2, 0.0),), 0.0, LQCCorrection(rho_c=1e12))
    a = np.array([0.3, 1.0, 3.0])
    assert np.allclose(big.H2(a), (8 * np.pi / 3) * 0.2 * a**-3, rtol=1e-9)


def test_lqc_rejects_curvature():
    with pytest.raises(ValueError):
        FriedmannModel((Fluid(1.0, 0.0),), 1.0, LQCCorrection())


def test_closed_universe_turning_points():
    # poeira fechada em GR: raiz unica de H^2 = 0 em a_m = 8 pi rho0 / 3
    m = FriedmannModel((Fluid(1.0, 0.0),), 1.0, GRCorrection())
    roots = m.turning_points(a_lo=1e-3, a_hi=1e3)
    assert len(roots) == 1 and abs(roots[0] - 8 * np.pi / 3) < 1e-6


def test_epsilon_and_efolds_are_reported():
    s = einstein_cartan_dust(1.0, sigma=0.05).solve(a0=1.0, t_span=(0, 1.2))
    assert np.isfinite(s.efolds())
    phases = s.accelerated_phases()
    assert len(phases) >= 1  # perto do ricochete a'' > 0 (super-aceleracao breve)
    assert all(N < 3 for _, _, N in phases)  # e-folds pequenos: nao ha inflacao
