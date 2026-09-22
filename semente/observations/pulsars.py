"""
Dados de massas de estrelas de neutrons e o teste da previsao de Smolin (1992).
RESULTADO COM CONEXAO OBSERVACIONAL: dados publicados, com fonte e incerteza 1 sigma.
"""
from __future__ import annotations


# ----------------------------------------------------------------------------
# D. Dados: massa maxima de estrelas de neutrons
# ----------------------------------------------------------------------------
NEUTRON_STAR_MASSES = {
    # (massa em M_sol, incerteza 1 sigma, referencia)
    "PSR J0740+6620": (2.08, 0.07, "Fonseca et al. 2021, ApJL 915, L12"),
    "PSR J0348+0432": (2.01, 0.04, "Antoniadis et al. 2013, Science 340, 448"),
    "PSR J1614-2230": (1.908, 0.016, "Arzoumanian et al. 2018, ApJS 235, 37"),
    "PSR J0952-0607": (2.35, 0.17, "Romani et al. 2022, ApJL 934, L17"),
}


def smolin_neutron_star_test(prediction=1.6, revised=2.0):
    """Smolin (1992) previu massa maxima de estrela de neutrons perto de 1.6 M_sol (mais tarde
    argumentou ~2 M_sol).  Quantos sigma cada pulsar medido esta acima de cada previsao?"""
    rows = {}
    for name, (m, dm, ref) in NEUTRON_STAR_MASSES.items():
        rows[name] = dict(massa=m, erro=dm, sigma_acima_1p6=(m - prediction) / dm, sigma_acima_2p0=(m - revised) / dm, ref=ref)
    return rows


