"""Testes do modulo nascimento."""
import numpy as np

from semente.nascimento import BounceSpectrum, curvature_constraints, matter_bounce_tensor_ratio


def test_mode_solver_reproduces_exact_dust_solution_before_bounce():
    """Longe do ricochete (eta << 0) o modo deve seguir f+ = e^{-ix}(1 - i/x)/sqrt(2k) exatamente."""
    bs = BounceSpectrum(a_b=0.05)
    k = 1.0
    v, dv, eta = bs.mode(k, eta_eval=-20.0)
    x = k * eta
    f = np.exp(-1j * x) * (1 - 1j / x) / np.sqrt(2 * k)
    assert abs(v - f) < 2e-3 * abs(f)


def test_spectrum_is_scale_invariant_on_plateau():
    bs = BounceSpectrum(a_b=0.05)
    ks = np.logspace(-0.7, 0.3, 6)
    P = bs.spectrum(ks)
    ns = bs.spectral_index(ks, P)
    assert abs(ns - 1.0) < 0.02
    assert (P.max() / P.min()) < 1.05


def test_curvature_bound_parent_mass_is_universe_scale():
    cc = curvature_constraints()
    assert cc["compatible_closed"]
    assert 1e22 < cc["M_parent_min_Msol"] < 1e25
    assert cc["t_recollapse_Gyr"] > 100 * cc["age_Gyr"]


def test_curvature_bound_when_open_is_incompatible():
    assert curvature_constraints(0.01, 0.001, n_sigma=2.0)["compatible_closed"] is False


def test_matter_bounce_tensor_ratio_is_excluded():
    assert matter_bounce_tensor_ratio() > 0.036
