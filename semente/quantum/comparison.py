"""
Comparacao GR vs Einstein–Cartan vs LQC efetiva para CONDICOES INICIAIS EQUIVALENTES.

Todas as tres sao estrategias do mesmo FriedmannModel; aqui se fixam o conteudo (fluidos), o fator de
escala inicial e a fase (contracao), e se comparam:
    densidade critica / de ricochete, escala minima, curvatura maxima, tempo do ricochete,
    contribuicao do spin (termo -alpha n^2), robustez a anisotropia (criterio de stability/bounce.py).

RESULTADO MATEMATICO (verificado em testes): para POEIRA, EC com sigma = rho_0^2/rho_b e LQC com
rho_c = rho_b sao identicas (H^2 e a''/a).  Para RADIACAO diferem: a correcao de EC escala com
n^2 ~ a^-6 (spin dos fermions), a de LQC com rho^2 ~ a^-8.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..cosmology.friedmann import BackgroundSolution, Fluid, FriedmannModel, GRCorrection
from ..stability.bounce import anisotropy_robustness
from .einstein_cartan import EinsteinCartanCorrection
from .lqc import LQCCorrection


@dataclass
class ModelComparison:
    name: str
    rho_star: Optional[float]
    a_min: float
    rho_max: float
    ricci_max: float
    t_bounce: Optional[float]
    bounce: bool
    correction_fraction_at_bounce: Optional[float]   # |rho_eff - rho| / rho no ricochete
    anisotropy_Sigma_b: Optional[float]
    efolds_after_bounce: float
    constraint_violation_max: float
    solution: BackgroundSolution = field(repr=False)


def compare_models(rho0=1.0, w=0.0, rho_star=20.0, a0=1.0, t_span=(0.0, 2.0), sigma0_sq=1e-6):
    """Mesmo fluido (rho0, w), mesma escala inicial, mesma fase; tres dinamicas.

    Para EC a correcao e sigma a^-6 com sigma = rho0^2/rho_star (equivale a n0 = rho0/m com m^2/alpha = rho_star).
    """
    fl = (Fluid(rho0, w, "fluid"),)
    models = {
        "GR": FriedmannModel(fl, 0.0, GRCorrection()),
        "Einstein-Cartan": FriedmannModel(fl, 0.0, EinsteinCartanCorrection(sigma=rho0**2 / rho_star)),
        "LQC efetiva": FriedmannModel(fl, 0.0, LQCCorrection(rho_c=rho_star)),
    }
    out = []
    for name, m in models.items():
        try:
            sol = m.solve(a0=a0, t_span=t_span, contracting=True, n_samples=3000)
        except ValueError as e:
            raise RuntimeError(f"{name}: {e}")
        ib = sol.i_bounce
        bounce = ib is not None
        corr = None
        Sig = None
        if bounce:
            ab = sol.a[ib]
            corr = float(abs(m.rho_eff(ab) - m.rho(ab)) / m.rho(ab))
            Sig = anisotropy_robustness(sol, rho_star, sigma0_sq).Sigma_bounce
        out.append(ModelComparison(name=name, rho_star=None if name == "GR" else rho_star, a_min=sol.a_min, rho_max=sol.rho_max,
                                   ricci_max=sol.ricci_max, t_bounce=sol.t_bounce, bounce=bounce,
                                   correction_fraction_at_bounce=corr, anisotropy_Sigma_b=Sig,
                                   efolds_after_bounce=sol.efolds() if bounce else float("nan"),
                                   constraint_violation_max=float(np.nanmax(sol.constraint_violation)), solution=sol))
    return out


def dust_equivalence_residual(rho0=1.0, rho_star=7.0, a=None):
    """max |H^2_EC - H^2_LQC| / H^2 e idem para a''/a, para poeira (deve ser ~1e-15)."""
    a = np.logspace(-0.6, 0.5, 200) if a is None else a
    fl = (Fluid(rho0, 0.0),)
    ec = FriedmannModel(fl, 0.0, EinsteinCartanCorrection(sigma=rho0**2 / rho_star))
    lq = FriedmannModel(fl, 0.0, LQCCorrection(rho_c=rho_star))
    ok = ec.rho(a) < rho_star
    a = a[ok]
    h = np.max(np.abs(ec.H2(a) - lq.H2(a)) / np.maximum(np.abs(lq.H2(a)), 1e-300))
    acc = np.max(np.abs(ec.accel(a) - lq.accel(a)) / np.maximum(np.abs(lq.accel(a)), 1e-300))
    return float(h), float(acc)


def radiation_difference(rho0=1.0, rho_star=7.0):
    """Para radiacao (w = 1/3) EC e LQC diferem: devolve as escalas de ricochete de cada uma."""
    fl = (Fluid(rho0, 1 / 3),)
    ec = FriedmannModel(fl, 0.0, EinsteinCartanCorrection(sigma=rho0**2 / rho_star))
    lq = FriedmannModel(fl, 0.0, LQCCorrection(rho_c=rho_star))
    return dict(a_bounce_EC=ec.bounce_scale(), a_bounce_LQC=lq.bounce_scale(),
                nota="EC: rho a^-4 = sigma a^-6 -> a_b = sqrt(sigma/rho0); LQC: rho a^-4 = rho_c -> a_b = (rho0/rho_c)^(1/4)")
