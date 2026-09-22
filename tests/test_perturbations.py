"""Testes da Fase 5: modos sobre fundo generico, r = 16 eps, expansao emergente, ponte observacional."""
import numpy as np

from semente.cosmology.expansion import analyze_expansion, scan_parameters
from semente.cosmology.friedmann import einstein_cartan_dust, lqc_fluid
from semente.cosmology.modes import ModeSolver, bounce_background
from semente.cosmology.perturbations import BounceSpectrum
from semente.observations.bridge import compare


def _dust_background(a_b=0.05, a0=2.0e4):
    return bounce_background(einstein_cartan_dust(1.0, a_b**3), a0=a0, n_per_side=60000)


def test_generic_solver_recovers_dust_power_law_and_reference_spectrum():
    sol = _dust_background()
    ms = ModeSolver(sol)
    p, a1, eta_s = ms._late_powerlaw()
    assert abs(p - 2.0) < 1e-3 and abs(a1 - 2 * np.pi / 3) < 1e-2
    ks = np.logspace(-0.7, 0.0, 6)
    res = ms.spectrum(ks)
    ref = BounceSpectrum(a_b=0.05).spectrum(ks)
    assert np.max(np.abs(res.P_test / ref - 1)) < 0.05          # concordancia com a referencia exata dentro de 5%
    assert abs(res.n_s - 1.0) < 0.06                            # quase invariante de escala
    assert np.all(np.abs(res.epsilon_exit - 1.5) < 1e-6)        # saida na contracao de poeira: eps = 3/2
    assert np.allclose(res.r, 24.0)                             # r = 16 eps


def test_radiation_background_is_strongly_blue():
    sol = bounce_background(lqc_fluid(1.0, 1 / 3, 20.0), a0=300.0, n_per_side=60000)
    ms = ModeSolver(sol)
    p, _, _ = ms._late_powerlaw()
    assert abs(p - 1.0) < 1e-3
    res = ms.spectrum(np.logspace(-0.5, 0.0, 5))
    assert res.n_s > 2.5   # contracao de radiacao: n_s - 1 = 2 (dualidade de Wands), muito azul


def test_no_emergent_inflation_and_scan_is_reproducible():
    sol = _dust_background()
    rep = analyze_expansion(sol)
    assert not rep.inflation_like and 0 < rep.N_accelerated < 1
    rows = scan_parameters(rho_stars=(2, 200), ws=(0.0,))
    assert all(not r["inflation_like"] and r["N_accelerated"] < 1 for r in rows)
    assert rows[0]["duration"] > rows[1]["duration"]   # duracao decresce com rho_*


def test_bridge_status_rules():
    assert compare("m", "n_s", 0.9665).status == "consistent"
    assert compare("m", "n_s", 0.9665 + 3 * 0.0038).status == "tension"
    assert compare("m", "n_s", 1.05).status == "excluded"
    assert compare("m", "r", 24.0).status == "excluded"
    assert compare("m", "r", 0.01).status == "consistent"
    assert compare("m", "n_s", None).status == "not_predicted"
    assert compare("m", "f_NL_local", -4.375, literature=True).status.startswith("literature:")
