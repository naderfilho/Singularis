<div align="center">

# 🌌 SEMENTE · blackhole-genesis

**Every black hole carries the seed of a white hole. Every white hole is a Big Bang.**<br>
This repository turns that sentence into solved equations, trained neural networks and computed images, without hiding where physics ends and speculation begins.

[![Português](https://img.shields.io/badge/idioma-Português-009c3b?style=for-the-badge&logo=googletranslate&logoColor=white)](README.md)
[![English](https://img.shields.io/badge/language-English-1f6feb?style=for-the-badge&logo=googletranslate&logoColor=white)](README.en.md)

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![NumPy](https://img.shields.io/badge/NumPy-013243?logo=numpy&logoColor=white)](requirements.txt)
[![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?logo=scipy&logoColor=white)](requirements.txt)
[![SymPy](https://img.shields.io/badge/SymPy-3B5526?logo=sympy&logoColor=white)](semente/geodesics.py)
[![PyTorch](https://img.shields.io/badge/PyTorch-EE4C2C?logo=pytorch&logoColor=white)](semente/pinn.py)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c?logo=plotly&logoColor=white)](semente/figures.py)
[![WebGL2 / GLSL](https://img.shields.io/badge/WebGL2-GLSL-990000?logo=webgl&logoColor=white)](web/index.html)
[![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?logo=javascript&logoColor=black)](web/index.html)

[![Tests](https://img.shields.io/badge/tests-27%20passing-2ea043?logo=pytest&logoColor=white)](tests/)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Top language](https://img.shields.io/github/languages/top/naderfilho/blackhole-genesis?color=3776AB)](https://github.com/naderfilho/blackhole-genesis)
[![Last commit](https://img.shields.io/github/last-commit/naderfilho/blackhole-genesis)](https://github.com/naderfilho/blackhole-genesis/commits/main)
[![Stars](https://img.shields.io/github/stars/naderfilho/blackhole-genesis?style=social)](https://github.com/naderfilho/blackhole-genesis/stargazers)

<img src="docs/img/bounce_sweep.gif" width="720" alt="Sweep of the ℓ parameter: from Schwarzschild to a wormhole">

*Computed, not drawn: the Simpson–Visser parameter ℓ goes from 0 (Schwarzschild, black shadow) to 2.6M (traversable wormhole). In between, the sky of the other universe appears inside the "shadow", seen through the white hole.*

</div>

---

## 📑 Contents

- [In 30 seconds](#-in-30-seconds)
- [The idea, as a chain of facts](#-the-idea-as-a-chain-of-facts)
- [What is new: open questions](#-what-is-new-experiments-on-open-questions)
- [The model against our universe](#-the-model-against-our-universe)
- [Get started](#-get-started)
- [Gallery](#%EF%B8%8F-gallery)
- [Layout](#%EF%B8%8F-layout)
- [Scientific honesty](#%EF%B8%8F-scientific-honesty-read-before-citing)
- [References](#-references)

---

## ⚡ In 30 seconds

| | |
|---|---|
| **What it is** | A Python package plus a WebGL viewer that solve, numerically and with neural networks, the chain *collapse → black hole → bounce → white hole → Big Bang of another universe*. |
| **What is new** | Three results computed here: the throat of a regular black hole is **fixed** by its interior (ℓ = R_b); under the holographic bound a tree of universes has **fewer than 5 generations**; the primordial spectrum produced by the bounce comes out with **n_s ≈ 1 without inflation**. |
| **What it is not** | A proof. Nothing here has been observed. Every claim is labelled *theorem*, *model* or *hypothesis*. |
| **Try it** | `python scripts/run_all.py` produces 21 figures. Open `web/index.html` and drag ℓ. |

---

## 🧭 The idea, as a chain of facts

| # | claim | status | where in the code |
|---|---|---|---|
| 1 | The maximal extension of Schwarzschild necessarily contains a **white** hole and a second universe. | **theorem** (Kruskal 1960) | `geometry.py`, fig. 01, 02 |
| 2 | Inside the horizon $r$ is time: the interior of a black hole **is** a cosmology (Kantowski–Sachs) lasting exactly $\pi M$. | **theorem** | `interior.py`, fig. 04 |
| 3 | The interior of a collapsing star **is** a closed Friedmann universe in contraction (Oppenheimer–Snyder). | **theorem** (1939) | `collapse.py`, fig. 05 |
| 4 | In a real collapse the white hole and the other universe **do not exist**: the star takes their place. | **theorem** | fig. 01 (hatched region) |
| 5 | If the singularity is replaced by a bounce, the contraction of item 3 becomes an expansion: a Big Bang. Einstein–Cartan torsion and loop quantum cosmology do this. | **published model** (Popławski 2010; Ashtekar et al. 2006) | `bounce.py`, fig. 06 |
| 6 | There is an exact metric (Simpson–Visser 2019) in which the black hole crosses a regular throat and exits as a white hole into **another universe**; the interior contracts to $R=\ell$ and re-expands. | **exact solution** with effective exotic matter | `geometry.py`, `raytracer.py`, `web/`, fig. 03, 04, 10, 12–15 |
| 7 | For a parent black hole of mass $M$ one can compute the density, size and time of the child's "Big Bang". | **SEMENTE model** (this project's synthesis) | `bounce.SeedUniverse`, fig. 07 |
| 8 | If each child inherits the parent's constants with mutations, the population of universes evolves to maximise black hole production (Smolin). | **falsifiable hypothesis** | `selection.py`, fig. 08 |
| 9 | A neural network that only sees the equations (PINN) reconstructs the bounce with error $10^{-5}$. | **result of this project** | `pinn.py`, fig. 11 |

<details>
<summary><b>📐 Kruskal diagram with the star erasing the white hole</b></summary>
<br>

![kruskal](docs/img/01_kruskal.png)

The hatched region is the star's interior: there the metric is Friedmann, not Schwarzschild. In a real collapse everything to the left of the surface (regions III and IV included) is replaced by the star. The white hole exists only in the eternal vacuum solution, or if the singularity is replaced by a bounce.

</details>

<details>
<summary><b>🕳️ The universe inside the black hole bouncing</b></summary>
<br>

![interior](docs/img/04_interior_universo.png)

</details>

Full derivations in [`docs/01_matematica.md`](docs/01_matematica.md) (Portuguese).

---

## 🔬 What is new: experiments on open questions

Four numerical tests run inside the models above, with results that, as far as we know, had not been computed in this form (`semente/fronteira.py`, figures 16–19, [`docs/02_fronteira.md`](docs/02_fronteira.md)).

| # | open question | what the test found |
|---|---|---|
| A | What happens **outside** the star when the interior bounces? | The Israel junction with Schwarzschild is impossible in a finite interval around the bounce. Within the black-bounce family **only ℓ = R_b** works: $\ell^3 = 3M/4\pi\rho_b$. For 10 M☉, ℓ ≈ 5×10⁻¹⁰ m; the wall between the universes weighs $R_b c^2/G \approx 7\times10^{17}$ kg. |
| B | Does the parent's horizon bound the child's entropy? | If so, the parent of our universe had ≥ 5×10¹³ M☉ and a tree of universes with our fecundity has **fewer than 5 generations**. |
| C | Does cosmological natural selection work under that bound? | No: effective fecundity drops to ~1 child per universe and the dynamics becomes neutral drift. **CNS and holography are nearly incompatible.** |
| D | Does Smolin's prediction (neutron star M_max ≈ 1.6 M☉) survive the data? | Excluded at > 4σ by four pulsars; the revised version (2 M☉) is at the edge. |

<details>
<summary><b>📊 The junction figure (why only ℓ = R_b works)</b></summary>
<br>

![juncao](docs/img/16_fronteira_juncao.png)

</details>

Each result depends on an explicit assumption (effective treatment of torsion; holographic bound). That is where physics is undecided, which is why these are tests, not theorems.

---

## 🌍 The model against our universe

Nothing confirms that we were born from a black hole. What can be done is the consistency test: the model predicts properties of the child universe; we measure ours (Planck 2018, BICEP/Keck 2021); we compare (`semente/nascimento.py`, figures 20–21, [`docs/03_nascimento.md`](docs/03_nascimento.md)).

| model prediction | predicted | observed | verdict |
|---|---|---|---|
| Primordial spectrum (Mukhanov–Sasaki through the bounce, Bunch–Davies vacuum in the dust contraction) | $n_s = 1.004$, < 1% variation over 1.2 decades | $n_s = 0.9665 \pm 0.0038$ | 🟢 **compatible**: near scale invariance **without inflation** |
| Curvature: the child is closed | $\Omega_k < 0$ | Planck+BAO $0.0007\pm0.0019$; Planck alone $-0.011\pm0.0065$ | 🟢 **compatible**, falsifiable |
| Parent mass (curvature + conservation) | $\ge 4.7\times10^{23}\,M_\odot$ | mass of the observable universe ~ $10^{23}\,M_\odot$ | 🟡 consistent (not evidence) |
| Tensor-to-scalar ratio of the minimal matter bounce | $r = 24$ | $r < 0.036$ | 🔴 **falsified** in the minimal version |
| Accelerated expansion | not predicted | $\Lambda$ dominates | ⚪ not predicted |

<details>
<summary><b>📈 The primordial spectrum coming out of the bounce</b></summary>
<br>

![espectro](docs/img/20_nascimento_espectro.png)

</details>

**Score:** four compatible, one falsified, one not predicted. The model is **alive and constrained**. The natural next test (open in this repository) is to redo the spectrum with $k=+1$ and with tensor modes, to see whether $r$ drops.

---

## 🚀 Get started

```bash
git clone https://github.com/naderfilho/blackhole-genesis.git
cd blackhole-genesis
pip install -r requirements.txt
python -m pytest -q tests          # 27 physical-consistency tests
python scripts/run_all.py          # 21 figures + renders in ./output (~4 min on CPU)
```

**Real-time viewer:** open [`web/index.html`](web/index.html) in a WebGL2 browser and drag **ℓ** from 0 to 3. Every pixel integrates the exact geodesic on the GPU.

<details>
<summary><b>🐍 Use as a library</b></summary>
<br>

```python
from semente.geometry import BlackBounce
from semente.interior import radial_infall_black_bounce
from semente.bounce import SeedUniverse, M_SUN
from semente.raytracer import render, Scene, Camera, save_png
from semente.nascimento import BounceSpectrum
import numpy as np

bb = BlackBounce(M=1.0, l=0.8)
print(bb.kind, bb.horizons)                      # horizons at ±sqrt(4M² − ℓ²)
q = radial_infall_black_bounce(l=0.8, r0=8.0)    # r(τ) crosses r = 0 and emerges on the other side
print(SeedUniverse(10 * M_SUN).table())          # the child universe of a 10 solar-mass black hole
save_png(render(Scene(l=0.8), Camera(width=960, height=540)), "my_black_hole.png")

bs = BounceSpectrum(a_b=0.05)                    # primordial spectrum through the bounce
ks = np.logspace(-0.7, 0.3, 8)
print(bs.spectral_index(ks, bs.spectrum(ks)))    # ≈ 1.00
```

</details>

<details>
<summary><b>🧪 Run only part of it</b></summary>
<br>

```bash
python scripts/run_all.py --sem-render      # scientific figures only (~1 min)
python scripts/run_all.py --rapido          # renders at 640x360
python -m semente.figures_fronteira         # open-question experiments only
python -m semente.figures_nascimento        # verdict vs Planck/BICEP only
```

</details>

---

## 🖼️ Gallery

<details>
<summary><b>Full list of the 21 figures</b></summary>
<br>

| file | content |
|---|---|
| `01_kruskal.png` | The four Kruskal regions, the Oppenheimer–Snyder star (which erases the white hole) and an infalling observer. |
| `02_penrose.png` | Penrose diagram of eternal Schwarzschild and the black-bounce "ladder" of universes. |
| `03_flamm_wormhole.png` | Flamm's paraboloid (Einstein–Rosen bridge) and the smooth throat of a wormhole. |
| `04_interior_universo.png` | Scale factors of the universe inside the black hole; with bounce; finite Kretschmann; infall through the throat. |
| `05_colapso_oppenheimer_snyder.png` | Collapse: surface, event horizon born at the centre, FRW $a(\tau)$ and its time-reversed version. |
| `06_ricochete.png` | Pure GR vs torsion vs LQC: $a(t)$, $H(t)$, effective density. Exact solution overlaid. |
| `07_tabela_universo_filho.png/.json` | Child-universe properties for parents from 3 to 6.5×10⁹ solar masses. |
| `08_selecao_cosmologica.png` | Cosmological natural selection: mean fecundity and population migration on the landscape. |
| `09_deflexao.png` | Light deflection angle vs impact parameter; divergence at the photon sphere. |
| `10_geodesicas.png` | Photon and particle orbits (general integrator with symbolic Christoffel symbols). |
| `11_pinn.png` | PINN vs integrator: bounce and photon orbit. |
| `12–15_render_*.png` | Ray-tracer: Schwarzschild, black-bounce, wormhole, lensing seen from above. |
| `16_fronteira_juncao.png` | Israel junction of the bounce: where Schwarzschild fails, why only ℓ = R_b works, predicted ℓ vs mass. |
| `17_fronteira_entropia.png` | Parent-mass floor from the entropy bound; maximum depth of the tree of universes. |
| `18_fronteira_selecao_orcamento.png` | Smolin's selection with a holographic budget: effective fecundity collapses to ~1. |
| `19_fronteira_estrelas_neutrons.png` | Pulsar masses vs Smolin's prediction. |
| `20_nascimento_espectro.png` | Mukhanov–Sasaki through the bounce: potential, spectrum with n_s ≈ 1, freezing modes. |
| `21_nascimento_veredito.png` | Parent mass vs observed curvature and the verdict table. |

</details>

<div align="center">
<img src="docs/img/13_render_black_bounce.png" width="720" alt="Black-bounce render">

*Regular black hole (ℓ = 1M) with a Page–Thorne accretion disk. Exact geodesics, gravitational + Doppler redshift, blackbody colour. The only aesthetic freedom is the disk temperature.*
</div>

---

## 🗂️ Layout

```
semente/
  geometry.py          Schwarzschild, Kruskal (Lambert W), Penrose, Flamm, black-bounce, Kretschmann
  geodesics.py         Christoffel via sympy -> DOP853 integrator; orbital form u(phi); deflection
  interior.py          Kantowski-Sachs inside the BH; interior bounce; infall through the throat
  collapse.py          Oppenheimer-Snyder: FRW/Schwarzschild junction, event horizon, Kruskal surface
  bounce.py            Friedmann with torsion (+ exact solution), LQC, SEMENTE model, Pathria
  selection.py         Smolin's population dynamics
  pinn.py              PINNs (torch): bounce in the volume variable; photon orbit
  raytracer.py         numpy ray-tracer: Page-Thorne (closed form), redshift, lensing, two skies
  fronteira.py         experiments on open questions  (+ figures_fronteira.py)
  nascimento.py        child-universe predictions vs Planck/BICEP  (+ figures_nascimento.py)
  figures.py           all figures
web/index.html         real-time WebGL2 ray-tracer (GLSL)
tests/                 27 tests
docs/                  derivations (01), open questions (02), birth (03)
```

---

## ⚖️ Scientific honesty (read before citing)

- Items 1–4 of the chain are **theorems** of general relativity. Item 4 is what most popular accounts omit: Kruskal's white hole is erased by the collapse.
- Items 5–6 are **models**: they depend on physics beyond classical relativity. They are exact solutions *of those models' equations*, not observations.
- The argument "the Schwarzschild radius of the observable universe equals its size" is an **identity** of the Friedmann equation for a flat universe ($r_s = c/H_0$ exactly), not evidence.
- The minimal SEMENTE model produces a **small, short-lived** child universe and predicts $r = 24$, excluded by data. Turning it into a universe like ours needs extra physics. The verdict table says so with numbers.
- Nothing here has been observed. White holes have never been detected.

---

## 📚 References

<details>
<summary><b>Show references</b></summary>
<br>

- Oppenheimer, J. R. & Snyder, H. (1939). *On Continued Gravitational Contraction.* Phys. Rev. 56, 455.
- Kruskal, M. (1960). Phys. Rev. 119, 1743. Szekeres, G. (1960). Publ. Math. Debrecen 7, 285.
- Novikov, I. (1964). *Delayed explosion of a part of the Fridman universe and quasars.* Astron. Zh. 41, 1075.
- Pathria, R. K. (1972). *The Universe as a Black Hole.* Nature 240, 298.
- Page, D. N. & Thorne, K. S. (1974). *Disk-Accretion onto a Black Hole.* ApJ 191, 499.
- Frolov, V., Markov, M. & Mukhanov, V. (1990). *Black holes as possible sources of closed and semiclosed worlds.* Phys. Rev. D 41, 383.
- Smolin, L. (1992). *Did the universe evolve?* Class. Quantum Grav. 9, 173.
- Wands, D. (1999). *Duality invariance of cosmological perturbation spectra.* Phys. Rev. D 60, 023507.
- Ashtekar, A., Pawlowski, T. & Singh, P. (2006). *Quantum nature of the Big Bang.* Phys. Rev. Lett. 96, 141301.
- Popławski, N. (2010). *Cosmology with torsion: An alternative to cosmic inflation.* Phys. Lett. B 694, 181.
- Popławski, N. (2012). *Nonsingular, big-bounce cosmology from spinor-torsion coupling.* Phys. Rev. D 85, 107502.
- Haggard, H. & Rovelli, C. (2015). *Black hole fireworks.* Phys. Rev. D 92, 104020.
- Ashtekar, A., Olmedo, J. & Singh, P. (2018). *Quantum Transfiguration of Kruskal Black Holes.* Phys. Rev. Lett. 121, 241301.
- Simpson, A. & Visser, M. (2019). *Black-bounce to traversable wormhole.* JCAP 02, 042.
- Raissi, M., Perdikaris, P. & Karniadakis, G. (2019). *Physics-informed neural networks.* J. Comput. Phys. 378, 686.
- Planck Collaboration (2020). *Planck 2018 results. VI.* A&A 641, A6. BICEP/Keck (2021). Phys. Rev. Lett. 127, 151301.

</details>

<div align="center">

MIT License · Author: **Nader Filho**

</div>
