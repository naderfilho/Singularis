<div align="center">

# SEMENTE · blackhole-genesis

**A computational laboratory for gravitation and cosmology.**<br>
Guiding question: *under what physical conditions can gravitational collapse transition into a nonsingular cosmological bounce, and what observable signatures could distinguish such a scenario from classical black-hole formation?*

[![Português](https://img.shields.io/badge/idioma-Português-009c3b?style=for-the-badge&logo=googletranslate&logoColor=white)](README.md)
[![English](https://img.shields.io/badge/language-English-1f6feb?style=for-the-badge&logo=googletranslate&logoColor=white)](README.en.md)

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![NumPy](https://img.shields.io/badge/NumPy-013243?logo=numpy&logoColor=white)](requirements.txt)
[![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?logo=scipy&logoColor=white)](requirements.txt)
[![SymPy](https://img.shields.io/badge/SymPy-3B5526?logo=sympy&logoColor=white)](semente/stability/energy_conditions.py)
[![PyTorch](https://img.shields.io/badge/PyTorch-EE4C2C?logo=pytorch&logoColor=white)](semente/ml/pinn.py)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c?logo=plotly&logoColor=white)](semente/figures/)
[![WebGL2 / GLSL](https://img.shields.io/badge/WebGL2-GLSL-990000?logo=webgl&logoColor=white)](web/index.html)
[![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?logo=javascript&logoColor=black)](web/index.html)

[![Tests](https://img.shields.io/badge/tests-89%20passing-2ea043?logo=pytest&logoColor=white)](tests/)
[![Validation](https://img.shields.io/badge/validation%20registry-15%2F15-2ea043)](semente/validation.py)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/naderfilho/blackhole-genesis)](https://github.com/naderfilho/blackhole-genesis/commits/main)
[![Stars](https://img.shields.io/github/stars/naderfilho/blackhole-genesis?style=social)](https://github.com/naderfilho/blackhole-genesis/stargazers)

<img src="docs/img/bounce_sweep.gif" width="720" alt="Sweep of the ℓ parameter: from Schwarzschild to a wormhole">

*Computed, not drawn: the Simpson–Visser parameter ℓ runs from 0 (Schwarzschild) to 2.6M (traversable wormhole). In between, the sky of the other universe appears inside the "shadow" (valid for the eternal geometry).*

</div>

---

## Contents

- [In 30 seconds](#in-30-seconds)
- [The five categories of claim](#the-five-categories-of-claim)
- [What the laboratory does (by study)](#what-the-laboratory-does-by-study)
- [Main results and negative results](#main-results-and-negative-results)
- [The model against the data](#the-model-against-the-data)
- [Get started](#get-started)
- [Reproducibility and validation](#reproducibility-and-validation)
- [Layout](#layout)
- [Documentation](#documentation-portuguese)

---

## In 30 seconds

| | |
|---|---|
| **What it is** | A modular platform (Python + WebGL) to investigate gravitational collapse, black holes, regular geometries, black-bounces, white-hole transitions, Einstein–Cartan, effective LQC, cosmological bounces, baby universes, perturbations, gravitational waves, stability, thermodynamics, information and PINNs — sharing metrics, units, solvers, results, validations and figures. |
| **What it has concluded so far** | Dust contraction inside a black hole yields a nearly scale-invariant spectrum **without inflation**, but the minimal version predicts $r = 24$ (excluded), produces **no** emergent inflation ($N<1$), is **not** robust to anisotropies and leaves **no** detectable ringdown signature. The black-bounce throat compatible with the junction is $\ell = R_b$. For dust, Einstein–Cartan and effective LQC are the same dynamics. |
| **What it is not** | Proof of anything. Every claim carries one of the five categories below; negative results are kept. |
| **Try it** | `python scripts/run_all.py` produces 38 figures; `python -m semente.validation` runs 15 checks; `web/index.html` is the real-time ray-tracer. |

---

## The five categories of claim

| category | examples in this repository |
|---|---|
| **mathematical/numerical result** | Kruskal, Oppenheimer–Snyder, $\tau=\pi M$, Regge–Wheeler, Leaver QNMs |
| **theoretical model from the literature** | Einstein–Cartan (Popławski), effective LQC (Ashtekar–Pawlowski–Singh), Simpson–Visser, charged and rotating black-bounces |
| **speculative hypothesis** | holographic bound on the child, cosmological natural selection, baby-universe Page curve |
| **result produced by this code** | $\ell = R_b$ from the Israel junction; EC ≡ LQC for dust; $n_s = 1.00$ through the bounce; no emergent inflation; first-law failure for SV with $S=A/4$ |
| **possible observational connection** | $n_s$, $r$, $\alpha_s$, $\Omega_k$, GW150914, pulsar masses |

Mandatory language: "the model predicts", "under these assumptions", "the simulation is consistent with", "this remains speculative". Never "this proves that black holes create universes". See [docs/01_foundations.md](docs/01_foundations.md) (Portuguese).

---

## What the laboratory does (by study)

| # | study | module | doc | figures |
|---|---|---|---|---|
| 1 | Dynamic gravitational collapse (LTB; density, mass, apparent and event horizons, curvature, shell crossing), classical **or** regularized interior | `collapse/dynamic.py` | 03 | 22, 23 |
| 2 | Dynamic black-bounce: collapse → horizon → high curvature → bounce → re-expansion; Israel junction with the exterior | `collapse/dynamic.py`, `collapse/junction.py` | 04 | 16, 22, 23 |
| 3 | Einstein–Cartan (torsion/spin): critical density, minimum scale, curvature, GR vs EC | `quantum/einstein_cartan.py`, `quantum/comparison.py` | 05 | 06, 27, 28 |
| 4 | Effective LQC ($H^2 = \tfrac{8\pi}{3}\rho(1-\rho/\rho_c)$), progenitor → child initial conditions | `quantum/lqc.py`, `cosmology/seed_universe.py` | 06 | 27, 28 |
| 5 | Energy conditions NEC/WEC/SEC/DEC (symbolic Einstein tensor) for any solution | `stability/energy_conditions.py` | 07 | 24, 33 |
| 6 | Dynamical stability: master potentials, bound-state criterion, bounce robustness to shear | `stability/perturbations.py`, `stability/bounce.py` | 07 | 25, 26 |
| 7 | Cosmological perturbations: scalar/tensor, $n_s$, $r$, running, bridge to Planck/BICEP | `cosmology/modes.py`, `cosmology/perturbations.py` | 08 | 20, 29 |
| 8 | Emergent expansion: $H$, $\epsilon$, e-folds, entry/exit — the code decides whether inflation occurs | `cosmology/expansion.py` | 08 | 30 |
| 9 | Rotating geometries: Kerr (horizons, ergosphere, frame dragging, ISCO) and the rotating black-bounce | `geometry/kerr.py` | 02 | 34 |
| 10 | Charged geometries: Reissner–Nordström and the charged black-bounce (phases in $(Q,\ell)$) | `geometry/charged.py`, `geometry/static.py` | 02 | 33 |
| 11 | Gravitational waves: QNMs (WKB), time-domain ringdown, echoes, spectrograms, GW150914 | `gravitational_waves/` | 09 | 31, 32 |
| 12 | Thermodynamics: $S=A/4$, $T$, first law (Schwarzschild, RN, Kerr, SV), apparent vs event horizon in collapse | `thermodynamics/horizon.py` | 10 | 35 |
| 13 | Information (exploratory): evaporation, Page curves, baby-universe scenario, timescales | `information/page_curve.py` | 11 | 36 |
| 14 | PINNs: analytical vs solver vs network; error, conservation, constraint, seed stability, hybrid mode | `ml/pinn.py`, `ml/benchmarks.py` | 12 | 11, 38 |
| 15 | Parameter space: 1D/2D scans, seeded Monte Carlo, adaptive refinement, classification | `parameter_space/explorer.py` | 13 | 37 |
| 16 | Observational bridge: MODEL / PREDICTION / CONSTRAINT / RESIDUAL / UNCERTAINTY / SOURCE | `observations/` | 14 | tables |
| 17 | Automated validation: dimensional, limits, conservation, convergence, stability, benchmarks | `semente/validation.py` | 01 | JSON report |
| 18 | Reproducibility: configs, seeds, provenance (commit, versions, solver, timestamp) on every figure and JSON | `core/provenance.py`, `core/config.py`, `configs/` | 01 | footers |

The original studies (Kruskal, Schwarzschild interior, Oppenheimer–Snyder, SEMENTE model, cosmological selection, ray-tracing) remain in docs/01_matematica, 02_fronteira and 03_nascimento (Portuguese).

<details>
<summary><b>Dynamic collapse: classical vs bounce</b></summary>
<br>

![collapse](docs/img/22_colapso_dinamico.png)

</details>

---

## Main results and negative results

| result | category | where |
|---|---|---|
| The interior of a collapsing star is a contracting Friedmann universe; Kruskal's white hole is erased by the collapse | mathematical result | 01, fig. 01, 05 |
| Only $\ell = R_b$ (throat = bounce radius) admits the Israel junction; the wall weighs $R_b c^2/G$ | result of this code | 02_fronteira, fig. 16 |
| For dust, Einstein–Cartan and effective LQC are identical ($\rho_c \leftrightarrow \rho_b$); they differ for radiation | result of this code | 05, 06, fig. 27 |
| Regularized collapse: **transient** trapped region (τ ∈ [23.7, 26.6] M), quasi-local mass drops to 5% at the bounce, event horizon undecided | result of this code | 03, 04, fig. 22 |
| Test-field spectrum through the bounce: $n_s = 1.00$ without inflation (Wands mechanism) | result of this code | 08, fig. 29 |
| **Negative:** $r = 24$ vs $r < 0.036$: minimal scenario excluded | code + observation | 08, 14 |
| **Negative:** no emergent inflation ($N<1$ e-fold, independent of $\rho_*$) | result of this code | 08, fig. 30 |
| **Negative:** isotropic bounce not robust to shear for stellar contraction ratios | result of this code | 07, fig. 26 |
| **Negative:** with $\ell = R_b$ the QNM shift is $<10^{-6}$; echoes only for $\ell > 2M$ | result of this code | 09, fig. 31, 32 |
| SV: $S = A/4$ is independent of $\ell$ and the first law fails by $1-\sqrt{1-\ell^2/4M^2}$ | result of this code | 10, fig. 35 |
| Test fields on black-bounces are linearly stable (no bound state) | numerical result | 07, fig. 25 |
| Under the holographic bound the tree of universes has < 5 generations and selection becomes drift | speculative hypothesis → computed consequence | 02_fronteira, fig. 17, 18 |

---

## The model against the data

| model | prediction | observed | residual | status |
|---|---|---|---|---|
| dust + torsion | $n_s = 1.01 \pm 0.02$ | 0.9665 ± 0.0038 | +2.4σ | tension |
| dust + torsion | $r = 24$ | < 0.036 | ×670 | **excluded** |
| radiation + LQC | $n_s = 2.9$ | 0.9665 | +20σ | **excluded** |
| closed interior | $\Omega_k < 0$ | 0.0007 ± 0.0019 / −0.011 ± 0.0065 | ≤ 1.7σ | consistent, falsifiable |
| Schwarzschild, non-rotating | $f_{\rm ringdown}$ = 179 Hz | 251 ± 8 Hz (GW150914) | −5σ | excluded (rotation) |
| Kerr χ = 0.67 (literature) | 250 Hz, 4.0 ms | 251 ± 8 Hz, 4.0 ± 0.3 ms | 0.1σ | literature:consistent |
| Smolin (1992) | $M_{\max,NS}$ ≈ 1.6 M☉ | 2.08 ± 0.07 M☉ | > 4σ | **excluded** |

"consistent" only means "not excluded". Full table in docs/14.

---

## Get started

```bash
git clone https://github.com/naderfilho/blackhole-genesis.git
cd blackhole-genesis
pip install -r requirements.txt
python -m pytest -q tests                # 89 tests (limits, conservation, convergence, benchmarks)
python -m semente.validation             # registry of 15 checks -> output/validation_report.json
python scripts/run_all.py                # 38 figures + JSONs with provenance in ./output
python scripts/reproduce.py "configs/*.json"   # re-run the studies from configuration files
```

Real-time viewer: [`web/index.html`](web/index.html) (WebGL2).

---

## Reproducibility and validation

- Every figure carries a footer with model, parameters, package version, commit and date; every result JSON has a `.meta.json` (`RunRecord`: parameters, solver, seed, numpy/scipy/sympy/torch versions, platform, timestamp).
- `configs/*.json` describe studies; `scripts/reproduce.py` re-runs them.
- `semente/validation.py` is the auditable registry: units, limits ($\ell\to0$, $\sigma\to0$, $\rho_c\to\infty$, $\rho_*\to\infty$), conservation (Friedmann constraint, Misner–Sharp mass), convergence (shells, grid), stability, analytical benchmarks (exact solution, Oppenheimer–Snyder, Leaver).
- One integrator (`core/solvers.py`) with recorded tolerances.

---

## Layout

```
semente/
  core/            units, provenance, single solver, study configs
  geometry/        generic static family, Schwarzschild/Kruskal/Penrose, SV, RN, charged, Kerr, geodesics, interior
  collapse/        Oppenheimer-Snyder, dynamic collapse (LTB + bounce), Israel junction
  quantum/         Einstein-Cartan, effective LQC, comparison
  cosmology/       unified Friedmann, modes, expansion, curvature, seed universe, selection
  stability/       energy conditions, master potentials, bounce robustness
  gravitational_waves/  QNM (WKB), time-domain ringdown, echoes
  thermodynamics/  horizons, genealogy
  information/     Page curves (EXPLORATORY)
  ml/              PINNs, benchmarks
  raytracing/      images (CPU)          web/index.html: real-time GPU
  observations/    sourced constraints, observational bridge, pulsars
  parameter_space/ scans, Monte Carlo, adaptive, classification
  figures/         all visualisation
  validation.py    check registry
configs/           reproducible studies      scripts/   run_all, reproduce
tests/             89 tests                  docs/      00-15 + the three original documents
```

The old modules (`semente.geometry`, `semente.bounce`, `semente.fronteira`, `semente.nascimento`, ...) remain importable as facades.

---

## Documentation (Portuguese)

00 audit · 01 foundations · 02 geometries · 03 collapse · 04 black-bounce · 05 Einstein–Cartan · 06 LQC · 07 stability · 08 perturbations · 09 gravitational waves · 10 thermodynamics · 11 information · 12 PINN · 13 parameter space · 14 observational constraints · 15 limitations — all under [`docs/`](docs/). Each document has: objective, mathematical formulation, assumptions, numerical method, parameters, validation, results, limitations, references.

<div align="center">

MIT License · Author: **Nader Filho**

</div>
