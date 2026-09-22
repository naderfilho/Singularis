"""Veredito legado (substituido por observations/bridge.py na fase 6)."""
from __future__ import annotations

import numpy as np

from ..cosmology.curvature import curvature_constraints, matter_bounce_tensor_ratio
from ..cosmology.perturbations import BounceSpectrum
from .constraints import PLANCK


# ----------------------------------------------------------------------------
# Veredito
# ----------------------------------------------------------------------------
def verdict(ks=None, verbose=True):
    ks = np.logspace(-0.7, 1.5, 36) if ks is None else ks
    bs = BounceSpectrum(a_b=0.05)
    P = bs.spectrum(ks)
    ns = bs.spectral_index(ks, P)
    cc = curvature_constraints()
    cc_planck = curvature_constraints(PLANCK["Omega_k_planck_only"], PLANCK["Omega_k_planck_only_err"], n_sigma=0.0)
    r = matter_bounce_tensor_ratio()
    rows = [
        dict(item="indice espectral n_s (contracao de poeira + ricochete)",
             previsto=f"{ns:.3f} (numerico, plato k << k_b)", observado=f"{PLANCK['n_s']} +- {PLANCK['n_s_err']}",
             veredito="COMPATIVEL em 1a ordem: quase invariante de escala sem inflacao; "
                      "a inclinacao vermelha de ~3% exige correcao (ex.: w ligeiramente negativo)"),
        dict(item="sinal da curvatura espacial", previsto="Omega_k < 0 (fechado)",
             observado=f"Planck+BAO: {PLANCK['Omega_k']} +- {PLANCK['Omega_k_err']};  Planck so: "
                       f"{PLANCK['Omega_k_planck_only']} +- {PLANCK['Omega_k_planck_only_err']}",
             veredito="COMPATIVEL (Planck sozinho ate prefere fechado a ~1.7 sigma); FALSIFICAVEL se Omega_k > 0 for medido"),
        dict(item="massa minima do buraco negro pai (da curvatura + conservacao de massa)",
             previsto=f">= {cc['M_parent_min_Msol']:.1e} M_sol  (Omega_k >= {cc['Omega_k_used']:.4f}, 2 sigma)",
             observado=f"maior buraco negro conhecido ~ 1e10-1e11 M_sol; massa do universo observavel ~ 1e23 M_sol",
             veredito="CONSISTENTE mas exige um pai com a massa de um universo inteiro - nao e evidencia, e conservacao"),
        dict(item="razao tensor/escalar r", previsto=f"r = 16 epsilon = {r:g} (ricochete de materia minimo)",
             observado=f"r < {PLANCK['r_max_95']} (95%, BICEP/Keck 2021)",
             veredito="FALSIFICADO na versao minima; sobrevive so com fisica extra que suprima tensores "
                      "(ex.: ricochete com Lambda/curvatura, campos adicionais)"),
        dict(item="idade vs recolapso do universo fechado", previsto=f"recolapso em {cc['t_recollapse_Gyr']:.1e} Gyr (poeira pura)",
             observado=f"idade {PLANCK['age_Gyr']} Gyr; expansao ACELERADA (Lambda)",
             veredito="CONSISTENTE na idade; a aceleracao observada NAO sai do modelo (precisa de Lambda herdada ou emergente)"),
    ]
    out = dict(n_s_numerico=ns, k_bounce=bs.k_bounce, ks=ks, P=P, curvatura=cc, curvatura_planck_only=cc_planck, r_min=r, rows=rows)
    if verbose:
        for row in rows:
            print(f"- {row['item']}\n    previsto : {row['previsto']}\n    observado: {row['observado']}\n    veredito : {row['veredito']}")
    return out
