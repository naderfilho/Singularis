"""
Registro de validacao cientifica automatizada.  `python -m semente.validation` executa todas as
verificacoes e grava output/validation_report.json (com proveniencia).

Cada verificacao declara: nome, categoria (dimensional / limiting_case / conservation / convergence /
stability / analytical_benchmark), criterio numerico e resultado.  A suite pytest cobre mais casos; este
registro e o resumo auditavel "de uma pagina".
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass

import numpy as np


@dataclass
class Check:
    name: str
    category: str
    criterion: str
    value: float
    tolerance: float
    passed: bool
    doc: str = ""


def _run_all_checks() -> list[Check]:
    from .collapse.dynamic import BounceShellDynamics, DustProfile, ShellDynamics, SphericalCollapse
    from .collapse.oppenheimer_snyder import OppenheimerSnyder
    from .core.units import C_SI, G_SI, L_PLANCK, M_PLANCK, Scale
    from .cosmology.friedmann import einstein_cartan_dust, exact_torsion_dust, gr_dust
    from .geometry.spherical import BlackBounce
    from .gravitational_waves.qnm import LEAVER_SCHWARZSCHILD, wkb3
    from .quantum.comparison import dust_equivalence_residual
    from .stability.energy_conditions import static_metric_report
    from .stability.perturbations import schwarzschild_metric
    checks = []

    def add(name, cat, crit, value, tol, doc=""):
        checks.append(Check(name, cat, crit, float(value), float(tol), bool(value <= tol), doc))

    # dimensional
    add("l_Pl = G m_Pl / c^2", "dimensional", "|G m_Pl/c^2 / l_Pl - 1|", abs(G_SI * M_PLANCK / C_SI**2 / L_PLANCK - 1), 1e-12, "core/units")
    add("r_s(M_sol) = 2953 m", "dimensional", "|r_s/2953.25 - 1|", abs(Scale.from_solar_masses(1.0).schwarzschild_radius_m() / 2953.25 - 1), 1e-3)
    # limiting cases
    r = np.array([3.0, 5.0, 10.0])
    add("black-bounce l->0 recupera Kretschmann de Schwarzschild", "limiting_case", "max rel err",
        np.max(np.abs(BlackBounce(1.0, 0.0).kretschmann(r) / (48 / r**6) - 1)), 1e-12, "docs/02")
    a = np.linspace(0.3, 1.0, 50)
    add("sigma->0 recupera GR (H^2)", "limiting_case", "max rel err", np.max(np.abs(einstein_cartan_dust(1.0, 1e-14).H2(a) / gr_dust(1.0).H2(a) - 1)), 1e-10, "docs/05")
    h, acc = dust_equivalence_residual()
    add("EC == LQC para poeira", "limiting_case", "max rel err (H^2, a''/a)", max(h, acc), 1e-12, "docs/04, 05")
    rep = static_metric_report("schwarzschild", np.linspace(2.1, 30, 100), M=1.0)
    add("Schwarzschild e vacuo (G_mu_nu = 0)", "limiting_case", "max |rho|,|p|", max(np.max(np.abs(rep.rho)), np.max(np.abs(rep.p_r))), 1e-14, "docs/07")
    # conservation
    sol = einstein_cartan_dust(1.0, 0.05).solve(1.0, (0, 1.2))
    add("vinculo de Friedmann preservado pelo integrador", "conservation", "max violacao normalizada", np.nanmax(sol.constraint_violation), 1e-8, "docs/05")
    prof = DustProfile.homogeneous(1.0, 8.0, 30)
    res = SphericalCollapse(prof).run()
    add("massa de Misner-Sharp conservada camada a camada (RG)", "conservation", "max |m_MS - m|", np.max(np.abs(res.m_MS - prof.m)), 1e-6, "docs/03")
    # analytical benchmarks
    exact, _ = exact_torsion_dust(sol.t, 1.0, 0.05)
    add("ricochete de torcao vs solucao exata", "analytical_benchmark", "max |a - a_exato|", np.max(np.abs(sol.a - exact)), 1e-8, "docs/05")
    os_ = OppenheimerSnyder(1.0, 8.0)
    eta = os_.eta_from_tau(res.tau)
    m = res.tau < 0.98 * os_.tau_s
    add("colapso dinamico vs Oppenheimer-Snyder", "analytical_benchmark", "max rel err R_s(tau)", np.max(np.abs(res.surface_R[m] / os_.R_surface(eta)[m] - 1)), 2e-3, "docs/03")
    q = wkb3(schwarzschild_metric(1.0), spin=2, ell=2)
    add("QNM WKB vs Leaver (Regge-Wheeler l=2)", "analytical_benchmark", "|omega - omega_Leaver|/|omega|", abs(q.omega - LEAVER_SCHWARZSCHILD[(2, 2, 0)]) / abs(LEAVER_SCHWARZSCHILD[(2, 2, 0)]), 5e-3, "docs/09")
    # convergence
    b1 = SphericalCollapse(DustProfile.homogeneous(1.0, 8.0, 20), BounceShellDynamics(rho_star=3.7)).run()
    b2 = SphericalCollapse(DustProfile.homogeneous(1.0, 8.0, 60), BounceShellDynamics(rho_star=3.7)).run()
    add("convergencia do raio minimo com n_shells (20 vs 60)", "convergence", "|R_min(20)/R_min(60) - 1|", abs(b1.summary()["R_surface_min"] / b2.summary()["R_surface_min"] - 1), 1e-3, "docs/03")
    q1 = wkb3(schwarzschild_metric(1.0), spin=0, ell=2, n_grid=3000)
    q2 = wkb3(schwarzschild_metric(1.0), spin=0, ell=2, n_grid=12000)
    add("convergencia do QNM com a grade (3000 vs 12000)", "convergence", "|omega1 - omega2|/|omega|", abs(q1.omega - q2.omega) / abs(q2.omega), 5e-4, "docs/09")
    # stability
    from .stability.perturbations import analyze, simpson_visser_metric
    st = analyze(simpson_visser_metric(1.0, 0.5), spin=0, ell=2)
    add("black-bounce l=0.5M: sem estado ligado (estavel, campo de teste)", "stability", "-(menor autovalor)", -st.lowest_eigenvalue, 0.0, "docs/07")
    add("classico: dinamica de ricochete rho_*->inf recupera RG", "limiting_case", "max |R_s - R_s(GR)|",
        np.max(np.abs(SphericalCollapse(prof, BounceShellDynamics(rho_star=1e12)).run(tau_max=20.0).surface_R
                      - SphericalCollapse(prof, ShellDynamics()).run(tau_max=20.0).surface_R)), 1e-8, "docs/04")
    return checks


def run(out_path="output/validation_report.json") -> dict:
    from .core.provenance import RunRecord, save_json
    t0 = time.time()
    checks = _run_all_checks()
    report = dict(n_checks=len(checks), n_passed=sum(c.passed for c in checks), all_passed=all(c.passed for c in checks),
                  elapsed_s=time.time() - t0, checks=[asdict(c) for c in checks])
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    save_json(out_path, report, RunRecord(model="validation/registry"))
    return report


if __name__ == "__main__":
    rep = run()
    for c in rep["checks"]:
        print(f"[{'ok' if c['passed'] else 'FALHOU'}] {c['category']:22s} {c['name']}: {c['criterion']} = {c['value']:.2e} (tol {c['tolerance']:.0e})")
    print(f"{rep['n_passed']}/{rep['n_checks']} verificacoes passaram em {rep['elapsed_s']:.0f}s")
