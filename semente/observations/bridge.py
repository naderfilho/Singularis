"""
Ponte observacional: cada comparacao modelo-vs-dado carrega
    MODEL, PREDICTION, OBSERVATIONAL CONSTRAINT, RESIDUAL, UNCERTAINTY, SOURCE, STATUS.

Regras de status (documentadas; nunca "evidencia positiva" por simples sobreposicao):
    consistent      |residual| <= 2 sigma  (ou previsao abaixo de um limite superior)
    tension         2 < |residual| <= 4 sigma
    excluded        |residual| > 4 sigma  (ou previsao acima de um limite superior por > 2x a incerteza do modelo)
    not_predicted   o modelo nao produz o observavel
    literature      valor tomado da literatura, nao calculado por este codigo (marcado explicitamente)
Um "consistent" significa apenas "nao excluido pelos dados atuais".
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional

import numpy as np

from .constraints import CONSTRAINTS, Constraint


@dataclass
class Comparison:
    model: str
    observable: str
    prediction: Optional[float]
    prediction_uncertainty: Optional[float]
    constraint: Constraint
    residual_sigma: Optional[float]
    status: str
    assumptions: list = field(default_factory=list)
    note: str = ""

    def as_row(self) -> dict:
        d = asdict(self)
        d["constraint"] = dict(value=self.constraint.value, sigma=self.constraint.sigma, upper_limit=self.constraint.upper_limit,
                               unit=self.constraint.unit, source=self.constraint.source)
        return d


def compare(model: str, key: str, prediction: Optional[float], prediction_uncertainty: Optional[float] = None,
            assumptions=None, note="", literature=False) -> Comparison:
    c = CONSTRAINTS[key]
    if prediction is None or (isinstance(prediction, float) and np.isnan(prediction)):
        return Comparison(model, c.name, None, None, c, None, "not_predicted", assumptions or [], note)
    unc_m = prediction_uncertainty or 0.0
    if c.value is not None and c.sigma:
        sig = float(np.hypot(c.sigma, unc_m))
        res = (prediction - c.value) / sig
        status = "consistent" if abs(res) <= 2 else ("tension" if abs(res) <= 4 else "excluded")
    elif c.upper_limit is not None:
        res = (prediction - c.upper_limit) / c.upper_limit
        status = "consistent" if prediction <= c.upper_limit else ("tension" if prediction <= c.upper_limit + 2 * unc_m else "excluded")
    else:
        res, status = None, "not_predicted"
    if literature:
        status = "literature:" + status
    return Comparison(model, c.name, float(prediction), prediction_uncertainty, c, None if res is None else float(res),
                      status, assumptions or [], note)


def table(rows) -> str:
    """Tabela de texto (markdown) das comparacoes."""
    hdr = "| modelo | observavel | previsto | observado | residual (sigma) | status | fonte |\n|---|---|---|---|---|---|---|\n"
    out = hdr
    for r in rows:
        c = r.constraint
        obs = f"{c.value} +- {c.sigma}" if c.value is not None else f"< {c.upper_limit} (95%)"
        pred = "-" if r.prediction is None else (f"{r.prediction:.4g}" + (f" +- {r.prediction_uncertainty:.2g}" if r.prediction_uncertainty else ""))
        res = "-" if r.residual_sigma is None else f"{r.residual_sigma:+.2f}"
        out += f"| {r.model} | {r.observable} | {pred} | {obs} {c.unit} | {res} | {r.status} | {c.source} |\n"
    return out
