"""Fachada de compatibilidade dos experimentos de fronteira (docs/02_fronteira.md).

Conteudo real: collapse/junction.py (A), thermodynamics/genealogy.py (B),
cosmology/selection.py (C), observations/pulsars.py (D).
"""
from __future__ import annotations

from .collapse.junction import TorsionCollapse, l_from_black_hole_mass  # noqa: F401
from .core.units import M_SUN  # noqa: F401
from .cosmology.selection import BudgetedSelection  # noqa: F401
from .observations.pulsars import NEUTRON_STAR_MASSES, smolin_neutron_star_test  # noqa: F401
from .thermodynamics.genealogy import K_B_ENTROPIES, genealogy_depth, parent_mass_lower_bound  # noqa: F401


def run_all(verbose=True):
    """Executa os quatro experimentos e devolve um dicionario de resultados."""
    import json

    import numpy as np

    out = {}
    tc = TorsionCollapse(M=1.0, R0=8.0, Rb_over_R0=0.05)
    sol = tc.solve()
    js = tc.junction_schwarzschild(sol)
    jl = tc.junction_black_bounce(sol, tc.l_predicted())
    ls, frac = tc.scan_l(sol)
    out["A_juncao"] = dict(
        R_b=tc.R_b, R_min_numerico=sol["R_min"],
        fracao_tempo_juncao_impossivel_schwarzschild=js["frac_time_invalid"],
        intervalo_falha_schwarzschild=js["tau_fail"],
        fracao_tempo_juncao_impossivel_black_bounce_l_eq_Rb=jl["frac_time_invalid"],
        massa_camada_max_sobre_M_l_eq_Rb=float(np.nanmax(np.abs(jl["m_shell"])) / tc.M),
        l_minimo_que_funciona_sobre_Rb=float(ls[frac <= 1e-3].min() / tc.R_b) if (frac <= 1e-3).any() else None,
        l_previsto_10Msol_m=l_from_black_hole_mass(10 * M_SUN),
        l_previsto_SgrA_m=l_from_black_hole_mass(4e6 * M_SUN),
        l_previsto_M87_m=l_from_black_hole_mass(6.5e9 * M_SUN),
    )
    S_total = sum(K_B_ENTROPIES.values())
    Mkg, Msol = parent_mass_lower_bound(S_total)
    _, Msol_cmb = parent_mass_lower_bound(K_B_ENTROPIES["CMB (fotons)"])
    out["B_entropia"] = dict(
        S_filho_total_kB=S_total, M_pai_min_kg=Mkg, M_pai_min_Msol=Msol, M_pai_min_Msol_so_CMB=Msol_cmb,
        profundidade_arvore={f"M0=1e{int(np.log10(M0))} Msol, N=1e{int(np.log10(N))}": genealogy_depth(M0 * M_SUN, N)
                             for M0 in (10, 1e6, 1e10) for N in (1e2, 1e10, 1e20)},
    )
    hist = BudgetedSelection().run(40)
    out["C_selecao_orcamento"] = dict(geracoes_ate_extincao=len(hist), historico=hist)
    out["D_estrelas_neutrons"] = smolin_neutron_star_test()
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "C_selecao_orcamento"}, indent=2, default=float))
    return out
