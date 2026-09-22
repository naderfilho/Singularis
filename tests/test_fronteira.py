"""Testes dos experimentos de fronteira."""
import numpy as np

from semente.bounce import M_SUN, M_PLANCK
from semente.fronteira import (TorsionCollapse, genealogy_depth, l_from_black_hole_mass,
                               parent_mass_lower_bound, smolin_neutron_star_test)


def _collapse():
    tc = TorsionCollapse(M=1.0, R0=8.0, Rb_over_R0=0.05)
    return tc, tc.solve()


def test_bounce_radius_matches_design():
    tc, sol = _collapse()
    assert abs(sol["R_min"] - tc.R_b) / tc.R_b < 1e-3


def test_schwarzschild_junction_fails_only_inside_horizon():
    tc, sol = _collapse()
    js = tc.junction_schwarzschild(sol)
    inv = ~js["valid"]
    assert inv.any()
    assert sol["R"][inv].max() < 2 * tc.M  # a falha so ocorre dentro do horizonte
    # e a superficie no ricochete esta dentro da faixa de falha
    ib = int(np.argmin(sol["R"]))
    assert inv[ib]


def test_shell_vanishes_when_torsion_is_off_and_outside():
    """Longe do ricochete (R >> R_b) a torcao e desprezivel e a camada deve ser ~0 (limite OS)."""
    tc, sol = _collapse()
    js = tc.junction_schwarzschild(sol)
    early = sol["tau"] < 0.3 * sol["tau_b"]
    assert np.nanmax(np.abs(js["m_shell"][early])) < 1e-3 * tc.M


def test_only_l_equal_Rb_allows_black_bounce_junction():
    tc, sol = _collapse()
    assert tc.junction_black_bounce(sol, tc.R_b)["frac_time_invalid"] == 0.0
    assert tc.junction_black_bounce(sol, 0.9 * tc.R_b)["frac_time_invalid"] > 0.0
    assert tc.junction_black_bounce(sol, 1.1 * tc.R_b)["frac_time_invalid"] == 1.0


def test_shell_mass_at_bounce_is_Rb_cos_chi0():
    tc, sol = _collapse()
    jl = tc.junction_black_bounce(sol, tc.R_b)
    ib = int(np.argmin(sol["R"]))
    assert abs(jl["m_shell"][ib] - tc.R_b * tc.cos) < 3e-2 * tc.R_b * tc.cos  # limite analitico na garganta: m = R_b cos(chi0)


def test_l_prediction_scales_as_M_to_one_third():
    l1 = l_from_black_hole_mass(10 * M_SUN)
    l2 = l_from_black_hole_mass(10000 * M_SUN)
    assert abs(l2 / l1 - 10.0) < 1e-9
    assert 1e-10 < l1 < 1e-9  # sub-nanometrico para 10 M_sol


def test_entropy_bound_and_depth():
    Mkg, Msol = parent_mass_lower_bound(3.1e104)
    assert 1e13 < Msol < 1e14
    # profundidade: com N = 1 nao ha decaimento (divisao por zero -> inf); com N grande e pequena
    assert genealogy_depth(1e10 * M_SUN, 1e20) < 5
    assert genealogy_depth(1e10 * M_SUN, 1e2) > 40
    assert abs(genealogy_depth(M_PLANCK * 100, 100) - 2.0) < 1e-9  # ln(100)/ln(10) = 2


def test_neutron_star_data_exclude_original_prediction():
    rows = smolin_neutron_star_test()
    assert all(r["sigma_acima_1p6"] > 4 for r in rows.values())
