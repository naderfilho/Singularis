# 04 — Black-bounce dinâmico

Módulos: `semente/collapse/dynamic.py` (evolução), `semente/collapse/junction.py` (exterior),
`semente/geometry/spherical.py` (métrica estática de Simpson–Visser, usada como candidata a exterior e
para ray-tracing). Figuras 22, 23 (dinâmica), 16 (junção), 04 (interior de SV), 12–15 (imagens).

## Objective

Não assumir que uma métrica black-bounce existe: simular a sequência
colapso → formação de horizonte → regime de alta curvatura → ricochete → re-expansão, e extrair tempo do
ricochete, raio mínimo, densidade e curvatura máximas, evolução da massa quasi-local e dos horizontes, e a
estrutura causal.

## Mathematical formulation

Dinâmica de camada regularizada (`BounceShellDynamics`):

$$\dot R^2 = \frac{2m}{R} + 2E - \frac{3m^2}{2\pi\rho_*R^4},\qquad
\ddot R = -\frac{m}{R^2} + \frac{3m^2}{\pi\rho_*R^5}.$$

Origem: a equação de Friedmann modificada $H^2 = \tfrac{8\pi}{3}\bar\rho(1-\bar\rho/\rho_*)$ escrita para cada
camada com a densidade **média** interior $\bar\rho = 3m/(4\pi R^3)$ e $R = a\chi$. Para o interior
homogêneo ela é **idêntica** a

- LQC efetiva (k = 0): $\rho_* = \rho_c$ (Ashtekar–Pawlowski–Singh 2006);
- Einstein–Cartan com poeira de férmions: $\rho_* = \rho_b = m_f^2/\alpha$ (Popławski 2010).

**Resultado matemático deste projeto (verificado em `test_shell_dynamics_equals_modified_friedmann_for_homogeneous_interior`):**
para poeira, LQC efetiva e Einstein–Cartan são a mesma família de um parâmetro, inclusive na equação de
aceleração ($\ddot a/a = -\tfrac{4\pi}{3}\rho(1-4\rho/\rho_*)$ em ambas). Elas diferem para outras equações
de estado (a correção de EC escala com $n^2\propto a^{-6}$; a de LQC com $\rho^2$).

## Assumptions

1. Poeira; camadas independentes (sem pressão, sem fluxo de energia entre camadas).
2. A correção age camada a camada com a densidade média interior: exato no caso homogêneo; para
   $E = 0$ segue Kelly–Santacruz–Wilson-Ewing (2020); para $E \ne 0$ é uma **extrapolação** deste código.
3. Condições iniciais ajustadas ao vínculo (documentado em 03).
4. O exterior não é simulado; é analisado por junção de Israel (`collapse/junction.py`).

## Numerical method

Como em 03. Parâmetro de ricochete escolhido por `Rb_over_R0` ($\rho_* = 3M/4\pi R_b^3$, caso $E=0$) ou
por `rho_star` diretamente.

## Parameters

`M`, `R0`, `rho_star` ou `Rb_over_R0`, `n_shells`, perfil.

## Validation

- $\rho_*\to\infty$ → RG (10⁻⁸).
- Interior homogêneo → Friedmann modificada (10⁻¹²).
- Raio mínimo da superfície = $R_b$ dentro de 2% (deslocamento por $E\ne0$).
- Ricochete simultâneo em todas as camadas no perfil homogêneo.

## Results (M = 1, R₀ = 8M, R_b/R₀ = 0.05)

| grandeza | valor |
|---|---|
| tempo do ricochete (superfície e centro) | 25.16 M |
| raio mínimo da superfície | 0.407 M (R_b nominal 0.4 M) |
| densidade máxima | 3.53 M⁻² (= ρ_*) |
| curvatura máxima (Kretschmann) | 2.0×10⁴ M⁻⁴ |
| horizonte aparente | existe de τ = 23.7 M a 26.6 M (transiente); raio máximo 1.98 M |
| massa de Misner–Sharp na superfície | 1 → 0.05 M no ricochete → 1 depois |
| raios nulos lançados depois da fase presa | escapam (15 de 22 resolvidos; 7 não alcançam a superfície no tempo simulado) |
| estrutura causal | fig. 23: região presa é uma faixa fina em torno do ricochete nas camadas externas |

**Consistência com a junção de Israel (docs/02, seção A):** o déficit de massa quasi-local na superfície
(95% no ricochete) é exatamente o que a junção exige que fique numa camada de fronteira ou numa geometria
exterior não-vácuo; o único exterior estático da família Simpson–Visser compatível tem $\ell = R_b$.

## Limitations

- Não há dinâmica do exterior: o "black-bounce dinâmico" aqui é o interior + a condição de junção. Uma
  simulação completa exigiria o exterior com matéria efetiva (violando a condição de energia nula perto
  da garganta, ver 05/stability) evoluído junto.
- A geometria de Simpson–Visser é usada como exterior *candidato*, não derivada da evolução.
- Para $0<\ell<2M$ a interpretação "luz do universo anterior" nas imagens vale só na geometria eterna.

## References

Simpson & Visser (2019) JCAP 02, 042; Kelly, Santacruz & Wilson-Ewing (2020) PRD 102, 106024;
Ashtekar, Olmedo & Singh (2018) PRL 121, 241301; Haggard & Rovelli (2015) PRD 92, 104020;
Popławski (2010) PLB 694, 181; Ashtekar, Pawlowski & Singh (2006) PRL 96, 141301.
