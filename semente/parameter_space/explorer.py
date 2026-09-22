"""
Explorador de espaco de parametros: varreduras 1D/2D, Monte Carlo e refinamento adaptativo, com
classificacao de solucoes e persistencia reprodutivel (semente + proveniencia).

Classes de solucao (rotulos padronizados):
    singular, classical_black_hole, regular_black_hole, black_bounce, wormhole, naked_singularity,
    unstable_bounce, stable_bounce, expanding_baby_universe, numerically_unresolved

Classificadores fornecidos (cada um documenta o criterio):
  * classify_static(M, l, Q):  familia black-bounce (carregado) via geometry.static.classify + condicoes
        - Q > M e l = 0 -> naked_singularity; l = 0 -> classical_black_hole (Q<=M); sem horizonte e l > 0 -> wormhole;
        - f(garganta) < 0 com um par de horizontes -> black_bounce; com dois pares -> regular_black_hole.
  * classify_bounce(rho_star, w, sigma0_sq): fundo LQC efetivo (k = 0) com fluido w e anisotropia inicial
        - GR (rho_star = inf) -> singular; solver falha -> numerically_unresolved;
        - Sigma_b >= 1 -> unstable_bounce (criterio BKL de stability/bounce.py);
        - senao, se ha expansao pos-ricochete (H > 0 no fim) -> expanding_baby_universe; caso contrario stable_bounce.
Os rotulos sao CLASSIFICACOES DO MODELO, nao afirmacoes fisicas sobre o universo.
"""
from __future__ import annotations

import itertools
import json
import os
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

from ..core.provenance import RunRecord, save_json

CLASSES = ("singular", "classical_black_hole", "regular_black_hole", "black_bounce", "wormhole", "naked_singularity",
           "unstable_bounce", "stable_bounce", "expanding_baby_universe", "numerically_unresolved")
CLASS_INDEX = {c: i for i, c in enumerate(CLASSES)}


# ----------------------------------------------------------------------------
# Classificadores
# ----------------------------------------------------------------------------
def classify_static(M=1.0, l=0.0, Q=0.0) -> str:
    from ..geometry.charged import charged_black_bounce, reissner_nordstrom
    try:
        if l <= 0:
            info = reissner_nordstrom(M, Q).classify() if Q > 0 else None
            if Q > M:
                return "naked_singularity"
            return "classical_black_hole"
        info = charged_black_bounce(M, Q, l).classify()
    except Exception:
        return "numerically_unresolved"
    n = info["n_horizons"]
    if n == 0:
        return "wormhole"
    if info.get("throat_f", 1.0) < 0:
        return "black_bounce" if n == 2 else "regular_black_hole"
    return "regular_black_hole"


def classify_bounce(rho_star=20.0, w=0.0, sigma0_sq=1e-6, a0=1.0, t_span=(0.0, 3.0)) -> str:
    from ..cosmology.friedmann import Fluid, FriedmannModel, GRCorrection
    from ..quantum.lqc import LQCCorrection
    from ..stability.bounce import anisotropy_robustness
    if not np.isfinite(rho_star):
        return "singular"
    try:
        m = FriedmannModel((Fluid(1.0, w),), 0.0, LQCCorrection(rho_c=rho_star))
        if m.H2(a0) <= 0:
            return "numerically_unresolved"
        sol = m.solve(a0=a0, t_span=t_span, contracting=True, n_samples=1500)
    except Exception:
        return "numerically_unresolved"
    if sol.i_bounce is None:
        return "singular" if sol.a[-1] < 0.05 * a0 else "numerically_unresolved"
    if np.nanmax(sol.constraint_violation) > 1e-5:
        return "numerically_unresolved"
    Sig = anisotropy_robustness(sol, rho_star, sigma0_sq).Sigma_bounce
    if Sig >= 1:
        return "unstable_bounce"
    return "expanding_baby_universe" if sol.H[-1] > 0 else "stable_bounce"


# ----------------------------------------------------------------------------
# Varreduras
# ----------------------------------------------------------------------------
@dataclass
class ScanResult:
    kind: str
    axes: dict                       # nome -> valores (1D/2D) ou limites (MC)
    points: np.ndarray               # (n, d)
    labels: list                     # rotulo por ponto
    seed: int | None
    record: RunRecord
    extra: dict = field(default_factory=dict)

    def as_dict(self):
        return dict(kind=self.kind, axes={k: (np.asarray(v).tolist()) for k, v in self.axes.items()},
                    points=self.points.tolist(), labels=self.labels, seed=self.seed, extra=self.extra,
                    classes=list(CLASSES))

    def save(self, path: str) -> str:
        return save_json(path, self.as_dict(), self.record)

    def grid_labels(self):
        """Para varreduras 2D: matriz (n_y, n_x) de indices de classe."""
        if self.kind != "scan_2d":
            raise ValueError("so para scan_2d")
        (nx_name, xs), (ny_name, ys) = list(self.axes.items())
        G = np.array([CLASS_INDEX[l] for l in self.labels]).reshape(len(ys), len(xs))
        return G


def _record(name, params, seed=None):
    return RunRecord(model=f"parameter_space/{name}", parameters=params, seed=seed)


def scan_1d(classifier: Callable, name: str, values: Sequence[float], fixed: dict | None = None) -> ScanResult:
    fixed = fixed or {}
    labels = [classifier(**{**fixed, name: float(v)}) for v in values]
    pts = np.asarray(values, float).reshape(-1, 1)
    return ScanResult("scan_1d", {name: list(values)}, pts, labels, None, _record("scan_1d", {**fixed, "axis": name}))


def scan_2d(classifier: Callable, x: tuple, y: tuple, fixed: dict | None = None) -> ScanResult:
    fixed = fixed or {}
    xn, xs = x
    yn, ys = y
    pts, labels = [], []
    for yv in ys:
        for xv in xs:
            pts.append((xv, yv))
            labels.append(classifier(**{**fixed, xn: float(xv), yn: float(yv)}))
    return ScanResult("scan_2d", {xn: list(xs), yn: list(ys)}, np.asarray(pts, float), labels, None,
                      _record("scan_2d", {**fixed, "x": xn, "y": yn}))


def monte_carlo(classifier: Callable, bounds: dict, n: int = 200, seed: int = 0, log_scale: Sequence[str] = ()) -> ScanResult:
    """Amostragem uniforme (ou log-uniforme) dentro de `bounds` com semente fixa: reproduzivel."""
    rng = np.random.default_rng(seed)
    names = list(bounds)
    pts = np.zeros((n, len(names)))
    for j, nm in enumerate(names):
        lo, hi = bounds[nm]
        if nm in log_scale:
            pts[:, j] = 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)
        else:
            pts[:, j] = rng.uniform(lo, hi, n)
    labels = [classifier(**dict(zip(names, map(float, p)))) for p in pts]
    return ScanResult("monte_carlo", {k: list(v) for k, v in bounds.items()}, pts, labels, seed,
                      _record("monte_carlo", {"bounds": {k: list(v) for k, v in bounds.items()}, "n": n}, seed))


def adaptive_refine(classifier: Callable, x: tuple, y: tuple, fixed: dict | None = None, levels: int = 2) -> ScanResult:
    """Grade 2D grossa + subdivisao das celulas cuja vizinhanca tem rotulos diferentes (fronteiras de fase)."""
    fixed = fixed or {}
    xn, xs = x
    yn, ys = y
    cache = {}

    def lab(xv, yv):
        key = (round(float(xv), 12), round(float(yv), 12))
        if key not in cache:
            cache[key] = classifier(**{**fixed, xn: float(xv), yn: float(yv)})
        return cache[key]

    xs, ys = list(map(float, xs)), list(map(float, ys))
    for xv in xs:
        for yv in ys:
            lab(xv, yv)
    n_added = 0
    for _ in range(levels):
        new_pts = []
        for i in range(len(xs) - 1):
            for j in range(len(ys) - 1):
                corners = {lab(xs[i], ys[j]), lab(xs[i + 1], ys[j]), lab(xs[i], ys[j + 1]), lab(xs[i + 1], ys[j + 1])}
                if len(corners) > 1:
                    new_pts.append((0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1])))
        for xv, yv in new_pts:
            lab(xv, yv)
            n_added += 1
        # refina a grade base so nas regioes de fronteira: aqui, por simplicidade e reprodutibilidade,
        # subdividimos globalmente a grade (os pontos ja avaliados sao reaproveitados pelo cache)
        xs = sorted(set(xs) | {0.5 * (a + b) for a, b in zip(xs[:-1], xs[1:])})
        ys = sorted(set(ys) | {0.5 * (a + b) for a, b in zip(ys[:-1], ys[1:])})
    pts = np.array(list(cache.keys()), float)
    labels = [cache[k] for k in cache]
    return ScanResult("adaptive", {xn: xs, yn: ys}, pts, labels, None,
                      _record("adaptive_refine", {**fixed, "x": xn, "y": yn, "levels": levels}), extra=dict(n_evaluations=len(cache)))
