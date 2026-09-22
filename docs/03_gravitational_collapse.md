# 03 — Colapso gravitacional dinâmico

Módulo: `semente/collapse/dynamic.py`. Figuras: `22_colapso_dinamico.png`, `23_estrutura_causal_colapso.png`
(resultados numéricos em `output/22_colapso_dinamico.json` com proveniência em `.meta.json`).

## Objective

Evoluir dinamicamente o colapso de poeira em simetria esférica, extrair densidade, raio areal, massa
interna, curvatura, formação de horizonte aparente e de horizonte de eventos, e permitir substituir a
dinâmica clássica por uma dinâmica regularizada (ricochete) mantendo **as mesmas condições iniciais**.

## Mathematical formulation

Coordenadas comóveis de Lemaître–Tolman–Bondi (G = c = 1):

$$ds^2 = -d\tau^2 + \frac{R'^2}{1+2E(\chi)}d\chi^2 + R^2(\tau,\chi)\,d\Omega^2 .$$

**RG clássica (resultado matemático, LTB 1933–1947):**
$\dot R^2 = 2m(\chi)/R + 2E(\chi)$, $\ddot R = -m/R^2$, $\rho = m'/(4\pi R^2R')$.

**Diagnósticos geométricos (válidos em qualquer modelo):**
massa de Misner–Sharp $m_{MS} = \tfrac{R}{2}(1-\nabla R\cdot\nabla R) = \tfrac{R}{2}(1+\dot R^2-(1+2E))$;
superfície aprisionada onde $\nabla R\cdot\nabla R = 1+2E-\dot R^2 \le 0$; horizonte aparente externo = camada
presa mais externa; raios nulos radiais $d\chi/d\tau = \pm\sqrt{1+2E}/R'$.

**Dinâmica de camada regularizada (modelo, ver 04_black_bounce.md):**
$\dot R^2 = 2m/R + 2E - 3m^2/(2\pi\rho_* R^4)$.

## Assumptions

Poeira sem pressão; simetria esférica; sem rotação; sem radiação; condições iniciais em repouso
($E = -m/R_0$) ou marginalmente ligadas; para a dinâmica regularizada, as condições iniciais são ajustadas
para satisfazer exatamente o vínculo $\dot R^2 = F(R)$ (correção $\le 3m^2/(4\pi\rho_*R_0^4)$, registrada em
`SphericalCollapse.constraint_adjustment`).

## Numerical method

Sistema de $2N$ EDOs (uma por camada) integrado com DOP853 (`core.solvers.integrate`, rtol $10^{-10}$,
atol $10^{-12}$), forma de 2ª ordem $\ddot R = \tfrac12\,dF/dR$ (regular nos pontos de retorno).
Eventos: alguma camada atinge $R < 10^{-3}R_0$ (singularidade clássica); superfície re-expande até
$R_0(1-10^{-3})$ (fim de um ciclo de ricochete). Derivadas em $\chi$: $m'$ analítico para os perfis
embutidos, $R'$ por diferenças centrais. Horizonte de eventos: lançamento de raios nulos do centro com
bissecção no instante de partida.

## Parameters

`M`, `R0`, `n_shells`, perfil (`homogeneous` = Oppenheimer–Snyder; `inhomogeneous` = núcleo denso com
`core_fraction`, `contrast`), `rho_star` (ou `Rb_over_R0`), tolerâncias do solver.

## Validation (tests/test_dynamic_collapse.py)

| teste | critério | resultado |
|---|---|---|
| superfície homogênea vs OS exato | erro relativo em $R(\tau)$ < 2×10⁻³ | passa |
| nascimento do horizonte de eventos vs OS ($\eta_h-\chi_0$) | < 0.1 M | passa (20.61 vs 20.57) |
| conservação de $m_{MS}$ camada a camada (RG) | < 10⁻⁶ | passa |
| homogeneidade de $\rho$ no perfil OS | dispersão < 10⁻³ | passa |
| $\rho_*\to\infty$ recupera RG | < 10⁻⁸ | passa |
| dinâmica de camada = Friedmann modificada (LQC e EC) | 10⁻¹² | passa |
| perfil inomogêneo: shell crossing detectado e centro colapsa primeiro | — | passa |

## Results (M = 1, R₀ = 8M, 40 camadas, R_b/R₀ = 0.05)

| grandeza | RG clássica | ricochete efetivo (ρ_* = 3.73 M⁻²) |
|---|---|---|
| tempo de colapso da superfície até R = 2M | 23.68 M | 23.72 M |
| singularidade / ricochete | τ = 25.13 M (OS: $\pi\sqrt{R_0^3/8M}$ = 25.13) | ricochete simultâneo em τ = 25.16 M, R_min = 0.407 M |
| densidade máxima | diverge (cortada em 10 M⁻²) | 3.53 M⁻² = ρ_* (saturação) |
| Kretschmann máximo | diverge | 2.0×10⁴ M⁻⁴ |
| região presa | permanente a partir de τ = 23.7 M | **transiente**: τ ∈ [23.7, 26.6] M |
| massa de Misner–Sharp na superfície | conservada (= M) | cai a 5% de M no ricochete (déficit = energia da correção) |
| horizonte de eventos | nasce no centro em τ = 20.6 M | fronteira escaped/trapped em 20.6 M, mas destino dos raios indeterminado sem a geometria exterior (ver 04) |

**Leitura:** no modelo regularizado a "estrela" atravessa uma fase presa de ~3 M e volta a expandir; o
buraco negro clássico é substituído por um objeto com horizonte aparente temporário. Isso é uma
consequência do modelo, não uma observação.

## Limitations

- Poeira: sem pressão, sem radiação, sem rotação. Estrelas reais têm equação de estado.
- A dinâmica regularizada para camadas com $E\ne0$ é uma extrapolação (ver 04).
- O horizonte de eventos é um conceito global; só o interior é simulado. O exterior durante o ricochete
  não é Schwarzschild (`collapse/junction.py`), então a existência de um horizonte de eventos no
  cenário com ricochete **não está decidida** por este módulo.
- Perfis inomogêneos desenvolvem shell crossing (singularidade fraca de coordenadas), após o qual $\rho$
  deixa de ser confiável; o instante é registrado e a simulação continua com essa ressalva.

## References

Lemaître (1933) Ann. Soc. Sci. Bruxelles A53, 51; Tolman (1934) PNAS 20, 169; Bondi (1947) MNRAS 107, 410;
Oppenheimer & Snyder (1939) Phys. Rev. 56, 455; Misner & Sharp (1964) Phys. Rev. 136, B571;
Kelly, Santacruz & Wilson-Ewing (2020) Phys. Rev. D 102, 106024.
