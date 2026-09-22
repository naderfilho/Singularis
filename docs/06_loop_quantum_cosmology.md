# 06 — Loop Quantum Cosmology (dinâmica efetiva)

Módulos: `semente/quantum/lqc.py`, `semente/cosmology/friedmann.py`, `semente/quantum/comparison.py`,
`semente/cosmology/seed_universe.py` (`mechanism="lqc"`). Figuras 06, 27, 28.

## Distinção obrigatória

| | o que é | implementado? |
|---|---|---|
| **LQC efetiva** | equações de Friedmann modificadas por correções de holonomia, válidas para estados semiclássicos (APS 2006; Taveras 2008; Singh 2009) | **sim** (k = 0) |
| **LQC completa** | dinâmica quântica do operador de Hamilton no espaço de Hilbert de cosmologia de laço | não |
| **Gravidade quântica de laço (LQG)** | teoria completa | não |

Tudo que este repositório chama de "LQC" é o **modelo efetivo**.

## Objective

Implementar as equações efetivas, verificar o ricochete em $\rho_c$, e ligar a massa do buraco negro
progenitor às condições iniciais do universo-filho.

## Mathematical formulation

$$H^2 = \frac{8\pi}{3}\rho\Big(1-\frac{\rho}{\rho_c}\Big),\qquad
\dot H = -4\pi(\rho+p)\Big(1-\frac{2\rho}{\rho_c}\Big),\qquad
\rho_c = \frac{\sqrt3}{32\pi^2\gamma^3}\rho_{Pl}\approx0.41\,\rho_{Pl}\ (\gamma = 0.2375).$$

Implementação como estratégia `LQCCorrection` com $\ddot a/a = \dot H + H^2$; recusa $k\ne0$ (as equações
efetivas com curvatura têm outra forma, APSV 2007).

## Assumptions

Estado semiclássico (correções efetivas válidas); k = 0; fluido barotrópico; sem inomogeneidades quânticas.

## Numerical method

Como em 05. Nota: para $\rho_0>\rho_c$ o ponto inicial está além do ricochete e o modelo levanta erro.

## Parameters

`rho0`, `w`, `rho_c` (padrão 0.41 ρ_Pl).

## Validation

| teste | resultado |
|---|---|
| ricochete exatamente em $\rho=\rho_c$ | 10⁻³ relativo |
| $\rho_c\to\infty$ recupera $H^2 = 8\pi\rho/3$ | 10⁻⁹ |
| identidade com Einstein–Cartan para poeira | 10⁻¹⁴ |
| k ≠ 0 rejeitado | passa |

## Results (fig. 27, 28)

Ver tabela em docs/05: para poeira, LQC ≡ EC; para radiação, LQC mantém $\rho\le\rho_c$ (a_min = 0.473
contra 0.224 em EC, com ρ_* = 20ρ₀). Progenitor → filho com $\rho_c = 0.41\rho_{Pl}$: $R_b = 1.3\times10^{-22}$ m
para 10 M☉ (contra $5\times10^{-10}$ m em EC com nêutrons); o tempo horizonte → ricochete é o mesmo
($\approx\pi M$). Consequência: **a escala do ricochete é a única diferença observável entre os dois
mecanismos no cenário de poeira**, e ela é inacessível (sub-Planckiana ou sub-nanométrica).

## Limitations

Sem a transição BN → BB de Ashtekar–Olmedo–Singh (2018) em Kantowski–Sachs quântico; a relação
progenitor → filho usa o interior FRW homogêneo (OS) e não o interior de Schwarzschild anisotrópico;
sem k = +1.

## References

Ashtekar, Pawlowski & Singh (2006) PRL 96, 141301; PRD 74, 084003. Ashtekar, Pawlowski, Singh & Vandersloot
(2007) PRD 75, 024035. Taveras (2008) PRD 78, 064072. Singh (2009) CQG 26, 125005.
Ashtekar, Olmedo & Singh (2018) PRL 121, 241301.
