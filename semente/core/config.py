"""
Arquivos de configuracao para estudos reprodutiveis.

Um config e um JSON com:
    {"study": "<nome>", "version": "1", "seed": 0, "parameters": {...}, "solver": {...}, "outputs": "output/<dir>"}
`run_study` despacha para a funcao registrada em STUDIES, grava os resultados com RunRecord (proveniencia)
e devolve o caminho do JSON de saida.  Cada estudo e uma funcao pura (parametros -> dict serializavel).
"""
from __future__ import annotations

import json
import os
from typing import Callable

import numpy as np

from .provenance import RunRecord, save_json
from .solvers import ODESettings


def _study_collapse(p: dict, seed: int) -> dict:
    from ..collapse.dynamic import compare_classical_vs_bounce
    c, b = compare_classical_vs_bounce(M=p.get("M", 1.0), R0=p.get("R0", 8.0), Rb_over_R0=p.get("Rb_over_R0", 0.05),
                                       n_shells=p.get("n_shells", 40), profile=p.get("profile", "homogeneous"))
    return dict(classico=c.summary(), ricochete=b.summary())


def _study_friedmann(p: dict, seed: int) -> dict:
    from ..quantum.comparison import compare_models
    comps = compare_models(rho0=p.get("rho0", 1.0), w=p.get("w", 0.0), rho_star=p.get("rho_star", 20.0),
                           t_span=tuple(p.get("t_span", (0.0, 1.5))), sigma0_sq=p.get("sigma0_sq", 1e-6))
    return {c.name: dict(a_min=c.a_min, rho_max=c.rho_max, ricci_max=c.ricci_max, t_bounce=c.t_bounce, bounce=c.bounce,
                         Sigma_b=c.anisotropy_Sigma_b, efolds=c.efolds_after_bounce, constraint=c.constraint_violation_max)
            for c in comps}


def _study_spectrum(p: dict, seed: int) -> dict:
    from ..cosmology.friedmann import einstein_cartan_dust
    from ..cosmology.modes import ModeSolver, bounce_background
    sol = bounce_background(einstein_cartan_dust(1.0, p.get("a_b", 0.05) ** 3), a0=p.get("a0", 2e4), n_per_side=p.get("n_per_side", 60000))
    ks = np.logspace(*p.get("log10_k", (-0.7, 0.0)), p.get("n_k", 6))
    r = ModeSolver(sol).spectrum(ks)
    return dict(k=ks.tolist(), P_test=r.P_test.tolist(), n_s=r.n_s, alpha_s=r.alpha_s, r=r.r.tolist(), k_bounce=r.k_bounce)


def _study_qnm(p: dict, seed: int) -> dict:
    from ..gravitational_waves.qnm import wkb3
    from ..stability.perturbations import simpson_visser_metric
    out = {}
    for l in p.get("l_values", [0.0, 0.5, 1.0]):
        q = wkb3(simpson_visser_metric(p.get("M", 1.0), max(l, 1e-9)), spin=p.get("spin", 0), ell=p.get("ell", 2))
        out[str(l)] = dict(omega_re=q.frequency, omega_im=-q.damping_rate)
    return out


def _study_parameter_scan(p: dict, seed: int) -> dict:
    from ..parameter_space.explorer import classify_bounce, classify_static, monte_carlo
    clf = {"static": classify_static, "bounce": classify_bounce}[p.get("classifier", "static")]
    res = monte_carlo(clf, {k: tuple(v) for k, v in p["bounds"].items()}, n=p.get("n", 100), seed=seed,
                      log_scale=tuple(p.get("log_scale", ())))
    return res.as_dict()


STUDIES: dict[str, Callable] = {
    "collapse": _study_collapse, "friedmann": _study_friedmann, "spectrum": _study_spectrum,
    "qnm": _study_qnm, "parameter_scan": _study_parameter_scan,
}


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        cfg = json.load(fh)
    for key in ("study", "parameters"):
        if key not in cfg:
            raise ValueError(f"config sem campo obrigatorio '{key}'")
    if cfg["study"] not in STUDIES:
        raise ValueError(f"estudo desconhecido: {cfg['study']} (disponiveis: {list(STUDIES)})")
    return cfg


def run_study(cfg: dict, out_dir: str | None = None) -> str:
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)
    result = STUDIES[cfg["study"]](cfg["parameters"], seed)
    rec = RunRecord(model=f"study/{cfg['study']}", model_version=str(cfg.get("version", "1")), parameters=cfg["parameters"],
                    solver=ODESettings(**cfg.get("solver", {})).as_dict(), seed=seed, notes=cfg.get("notes", ""))
    out_dir = out_dir or cfg.get("outputs", "output/studies")
    path = os.path.join(out_dir, f"{cfg['study']}.json")
    save_json(path, result, rec)
    return path
