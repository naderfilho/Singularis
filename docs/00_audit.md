# 00 — Auditoria do repositório (Fase 1)

Data: 2026-09-22. Estado auditado: commit `9ed8fce` (27 testes passando, 21 figuras geradas por `scripts/run_all.py --sem-render`).

## 1. O que já existe (inventário por implementação matemática)

| tema | módulo original | o que implementa | status |
|---|---|---|---|
| Schwarzschild | `geometry.py` | f(r), tartaruga, Kretschmann 48M²/r⁶, Flamm | ok |
| Kruskal–Szekeres | `geometry.py` | (t,r)→(T,X), inversa via Lambert W, regiões I–IV, curvas de r const | ok, testado (round-trip) |
| Penrose | `geometry.py`, `figures.py` | compactificação arctan de U,V; diagrama de Schwarzschild eterno; escada do black-bounce (esquemática) | ok |
| Flamm / Einstein–Rosen | `geometry.py`, `figures.py` | z(r) = √(8M(r−2M)); embedding do buraco de minhoca SV | ok |
| Simpson–Visser (black-bounce) | `geometry.py` | f, R=√(r²+ℓ²), horizontes, esfera de fótons, Kretschmann (base ortonormal), aceleração nula/tipo-tempo | ok, testado (ℓ→0) |
| Geodésicas | `geodesics.py` | Christoffel simbólico (sympy) → DOP853; forma orbital u(φ); deflexão | ok, testado (órbita circular, série de Einstein) |
| Oppenheimer–Snyder | `collapse.py` | junção FRW fechado/Schwarzschild; horizonte de eventos interior; superfície em Kruskal (via v de EF) | ok, testado |
| Interior de Schwarzschild | `interior.py` | Kantowski–Sachs, τ = πM, versão black-bounce, queda pela garganta | ok, testado |
| Ricochete (Friedmann modificada) | `bounce.py` | Einstein–Cartan (termo σa⁻⁶), LQC efetiva, solução exata a³ = σ/ρ₀ + 6πρ₀(t−t_b)², referência GR | ok, testado contra a exata (1e-11) |
| Universo-filho (SEMENTE) | `bounce.py` | ρ_b, R_b, a_b, tempo, fator de expansão vs massa do pai | ok; **tempo horizonte→ricochete usa aproximação ad hoc** |
| Seleção cosmológica | `selection.py`, `fronteira.py` | replicador + mutação; versão com orçamento holográfico | brinquedo declarado |
| PINN | `pinn.py` | ricochete na variável w=a³; órbita de fóton; Adam+L-BFGS+currículo | ok, testado (5e-3) |
| Ray-tracing | `raytracer.py`, `web/index.html` | geodésicas nulas SV, Page–Thorne (forma fechada + integral), redshift, dois céus; GPU | ok, testado (limite newtoniano) |
| Junção de Israel | `fronteira.py` | interior FRW+torção vs exterior Schwarzschild/SV; ℓ = R_b | ok, testado |
| Limite entrópico / genealogia | `fronteira.py` | M_pai ≥ m_Pl√(S/4π); k_max | ok, testado |
| Perturbações | `nascimento.py` | Mukhanov–Sasaki campo de teste através do ricochete; n_s; extração exata do modo constante | ok, testado |
| Curvatura / massa do pai | `nascimento.py` | Ω_k → a₀ → a_m → M_pai; idade vs recolapso | ok, testado |
| Dados observacionais | `nascimento.py`, `fronteira.py` | Planck 2018, BICEP/Keck 2021, 4 pulsares | hard-coded, sem camada própria |

## 2. Duplicações encontradas

1. **Equação de Friedmann com torção** aparece quatro vezes: `bounce.TorsionCosmology` (k=0, poeira+radiação, σ), `fronteira.TorsionCollapse` (k=+1, poeira, s ligado a R_b), `nascimento.BounceSpectrum` (fundo exato k=0) e `pinn.BouncePINN` (resíduo próprio). → unificar em `cosmology/friedmann.py` com estratégias de densidade efetiva (GR, Einstein–Cartan, LQC) e curvatura k.
2. **Constantes físicas**: SI em `bounce.py`; `MPC`, `GYR` em `nascimento.py`; importações cruzadas. → `core/units.py`.
3. **Linha de mundo em Kruskal via v de Eddington–Finkelstein** construída três vezes (`collapse.surface_kruskal`, `figures.fig_kruskal` duas vezes). → `geometry/kruskal.py: worldline_to_kruskal`.
4. **Órbita de fóton por RK4** em `geodesics.deflection_angle` e em `figures.fig_geodesic_gallery`. → `geodesics.photon_orbit`.
5. **Integração de τ em Kantowski–Sachs** com grades ad hoc em duas funções. → parametrização exata única.
6. Oito chamadas a `solve_ivp` com tolerâncias diferentes e não registradas. → `core/solvers.py`.

## 3. Inconsistências e fragilidades

- `SeedUniverse.tau_horizon_to_bounce_s` usa fórmula aproximada `πM(1 − (R_b/r_s)^{3/2})`; a forma exata existe (η de OS). **Corrigir.**
- Retornos em `dict` sem tipo em todos os módulos físicos; nomes de chaves inconsistentes (`a_bounce` vs `R_min`). → dataclasses de resultado.
- `figures.py` contém física inline (queda livre, horizonte) — visualização acoplada ao modelo. → mover para a biblioteca.
- Dados observacionais hard-coded em dois módulos, sem fonte estruturada nem incerteza padronizada. → `observations/constraints.py`.
- Números mágicos: `BudgetedSelection.m_ref = 1e39`, paisagem de fecundidade com larguras fixas; `raytracer` com `T_peak_K`, `R_out`, `R_far` (estes são escolhas declaradas de cena; manter como parâmetros de `Scene`).
- Interpretação do "céu do outro universo" para 0<ℓ<2M vale só na geometria eterna; está documentado, mas deve constar da camada de limitações.
- Sem registro de proveniência: nenhuma figura sabe com que commit, versão de pacote, semente ou tolerâncias foi gerada.

## 4. Ausências (relativas ao objetivo da plataforma)

Colapso dinâmico genérico (além de OS homogêneo); black-bounce **dinâmico** (evolução, não métrica assumida); condições de energia automatizadas; análise de estabilidade; ondas gravitacionais/ringdown; termodinâmica de horizonte como módulo; informação/Page; geometrias com rotação/carga; explorador de espaço de parâmetros; ponte observacional estruturada; infraestrutura de validação (limites, convergência); reprodutibilidade.

## 5. Arquitetura de evolução (decidida nesta fase)

Pacote continua `semente` (nome público do repositório: blackhole-genesis). Subpacotes:

```
semente/
  core/            units, provenance (metadados de execução), solvers (configuração única de ODE), results
  geometry/        metrics (base estática esfericamente simétrica + Schwarzschild, Simpson–Visser,
                   Reissner–Nordström, charged black-bounce), kruskal/penrose, geodesics, kerr (fase 7)
  collapse/        oppenheimer_snyder, junction (Israel), dynamic (fase 2: colapso dinâmico + bounce dinâmico)
  quantum/         einstein_cartan, lqc (estratégias de densidade efetiva, com referências)
  cosmology/       friedmann (fundo unificado), seed_universe, perturbations, expansion, selection
  stability/       energy_conditions, perturbations (potenciais mestres, estados ligados), criteria
  gravitational_waves/  qnm (WKB), ringdown (formas de onda, espectrogramas, ecos)
  thermodynamics/  horizon (área, S, T), genealogy (limite entrópico)
  information/     page_curve (EXPLORATÓRIO)
  ml/              pinn, benchmarks
  raytracing/      raytracer
  observations/    constraints (dados com fonte e incerteza), bridge (MODEL/PREDICTION/CONSTRAINT/RESIDUAL/UNCERTAINTY/SOURCE)
  parameter_space/ explorer (varreduras 1D/2D, Monte Carlo, adaptativo, classificação, persistência)
  figures/         toda a visualização (separada da física)
```

Compatibilidade: os módulos antigos (`semente.geometry`, `semente.bounce`, `semente.fronteira`, `semente.nascimento`, ...) continuam importáveis como fachadas que reexportam a nova localização. Testes antigos permanecem válidos.

## 6. Plano por fases (o que cada fase entrega e testa)

| fase | entrega | validação |
|---|---|---|
| 1 | auditoria, reorganização, `core` (unidades, proveniência, solvers), correção do tempo de queda | suíte antiga + testes de unidades |
| 2 | colapso dinâmico (OS/LTB) com interior substituível; black-bounce dinâmico | limites: LTB homogêneo = OS; bounce→GR quando correção→0 |
| 3 | condições de energia; estabilidade (potenciais mestres, estados ligados; perturbação homogênea do bounce) | Schwarzschild satisfaz todas; SV viola NEC perto da garganta (resultado de SV 2019) |
| 4 | Einstein–Cartan e LQC como estratégias unificadas; comparação GR vs EC vs LQC | solução exata; ρ_c → ∞ recupera GR |
| 5 | perturbações (escalar/tensorial, n_s, running, r), expansão emergente (ε, N) | n_s ≈ 1 na poeira; N e-folds reportado sem ajuste |
| 6 | QNMs (WKB), ringdown, ecos; ponte observacional | Schwarzschild ℓ=2: ω conhecido; residuais com fonte |
| 7 | Reissner–Nordström, charged black-bounce, Kerr (horizontes, ergosfera) | limites Q→0, a→0 |
| 8 | termodinâmica de horizonte (S, T, 2ª lei durante colapso/bounce), informação (Page) | S_BH de Schwarzschild; Page time conhecido |
| 9 | benchmarks PINN; explorador de parâmetros com classificação e persistência | reprodutibilidade por semente |
| 10 | documentação 01–15, limitações, integração final | suíte completa |
