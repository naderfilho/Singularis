"""
Restricoes observacionais com fonte e incerteza.  Toda comparacao modelo-vs-dado do pacote
le daqui; nenhum modulo fisico deve carregar valores observacionais soltos.

Classificacao: RESULTADOS COM CONEXAO OBSERVACIONAL (dados publicados).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Constraint:
    name: str
    value: Optional[float]
    sigma: Optional[float]                 # incerteza 1 sigma (None para limites)
    upper_limit: Optional[float] = None    # limite superior (95%) quando aplicavel
    unit: str = ""
    source: str = ""
    note: str = ""

    def residual(self, prediction: float) -> Optional[float]:
        """(previsto - observado) / sigma, ou (previsto - limite)/limite para limites superiores."""
        if self.value is not None and self.sigma:
            return (prediction - self.value) / self.sigma
        if self.upper_limit is not None:
            return (prediction - self.upper_limit) / self.upper_limit
        return None


# Planck 2018 VI (A&A 641, A6, 2020), Tabela 2, TT,TE,EE+lowE+lensing+BAO salvo indicacao
PLANCK = dict(
    H0=67.66, H0_err=0.42, Omega_m=0.3111, Omega_m_err=0.0056,
    Omega_k=0.0007, Omega_k_err=0.0019,
    Omega_k_planck_only=-0.011, Omega_k_planck_only_err=0.0065,
    n_s=0.9665, n_s_err=0.0038, r_max_95=0.036, age_Gyr=13.787,
    alpha_s=-0.0045, alpha_s_err=0.0067,   # running dn_s/dln k (Planck 2018 X, TT,TE,EE+lowE+lensing)
)

CONSTRAINTS = {
    "n_s": Constraint("indice espectral escalar", 0.9665, 0.0038,
                      source="Planck 2018 VI, Tab. 2 (TT,TE,EE+lowE+lensing+BAO)"),
    "alpha_s": Constraint("running dn_s/dln k", -0.0045, 0.0067, source="Planck 2018 X (TT,TE,EE+lowE+lensing)"),
    "r": Constraint("razao tensor/escalar", None, None, upper_limit=0.036,
                    source="BICEP/Keck 2021, PRL 127, 151301 (95%, k = 0.05 Mpc^-1)"),
    "Omega_k": Constraint("curvatura espacial", 0.0007, 0.0019, source="Planck 2018 VI + BAO"),
    "Omega_k_planck": Constraint("curvatura espacial (Planck so)", -0.011, 0.0065,
                                 source="Planck 2018 VI, TT,TE,EE+lowE+lensing"),
    "f_NL_local": Constraint("nao-gaussianidade local", -0.9, 5.1, source="Planck 2018 IX (A&A 641, A9)"),
    "age_Gyr": Constraint("idade do universo", 13.787, 0.020, unit="Gyr", source="Planck 2018 VI"),
    "GW150914_f_ringdown_Hz": Constraint("frequencia de ringdown GW150914", 251.0, 8.0, unit="Hz",
                                         source="LIGO/Virgo, PRL 116, 221101 (2016)"),
    "GW150914_tau_ms": Constraint("tempo de amortecimento GW150914", 4.0, 0.3, unit="ms",
                                  source="LIGO/Virgo, PRL 116, 221101 (2016)"),
    "GW150914_Mf_Msun": Constraint("massa final GW150914", 62.0, 4.0, unit="M_sun",
                                   source="LIGO/Virgo, PRL 116, 061102 (2016)"),
    "GW150914_chi_f": Constraint("spin final GW150914", 0.67, 0.07, source="LIGO/Virgo, PRL 116, 061102 (2016)"),
}
