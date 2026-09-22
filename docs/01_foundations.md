# 01 — Fundamentos: como ler este laboratório

## A pergunta

> Sob que condições físicas um colapso gravitacional pode transitar para um ricochete cosmológico não
> singular, e que assinaturas observáveis distinguiriam esse cenário da formação clássica de um buraco negro?

Tudo no repositório é um passo para responder a isso: geometrias (02), colapso dinâmico (03), black-bounce
dinâmico (04), os mecanismos de ricochete (05, 06), estabilidade e condições de energia (07), perturbações e
expansão (08), ondas gravitacionais (09), termodinâmica (10), informação (11), ML (12), varreduras (13),
restrições observacionais (14) e limitações (15).

## As cinco categorias de afirmação (obrigatórias em todo módulo, figura e doc)

| categoria | o que significa | exemplo |
|---|---|---|
| **resultado matemático/numérico** | teorema ou cálculo exato/numérico dentro de uma teoria dada | Kruskal, Oppenheimer–Snyder, $\tau = \pi M$, QNMs de Schwarzschild |
| **modelo teórico da literatura** | equações publicadas, implementadas e citadas | Einstein–Cartan (Popławski), LQC efetiva (APS), Simpson–Visser, Kerr black-bounce |
| **hipótese especulativa** | conjectura explicitamente marcada | limite holográfico no filho, seleção cosmológica, curva "universo-filho" |
| **resultado produzido por este código** | consequência calculada aqui, com testes | $\ell = R_b$ pela junção; EC ≡ LQC para poeira; $n_s = 1.00$ pelo ricochete; ausência de inflação emergente |
| **conexão observacional possível** | comparação com dado real, com residual e fonte | $n_s$, $r$, $\Omega_k$, GW150914, pulsares |

Linguagem: "o modelo prevê", "sob estas hipóteses", "a simulação é consistente com", "permanece
especulativo". Nunca "isto prova que buracos negros criam universos". Resultados negativos são preservados.

## Unidades

- Módulos geométricos: $G = c = 1$, com $M$ como escala. `core.units.Scale` converte para SI.
- Módulos cosmológicos/quânticos: unidades de Planck ($G = c = \hbar = 1$).
- Observáveis: SI (Hz, ms, kg, m) com fonte.

## Arquitetura

```
semente/
  core/            unidades, proveniencia (RunRecord), solver unico (ODESettings), config
  geometry/        metricas estaticas genericas, Schwarzschild/Kruskal/Penrose, SV, RN, carregado, Kerr, geodesicas
  collapse/        Oppenheimer-Snyder, colapso dinamico (LTB + ricochete), juncao de Israel
  quantum/         Einstein-Cartan, LQC efetiva (estrategias de densidade efetiva), comparacao
  cosmology/       Friedmann unificado, modos (Mukhanov-Sasaki), expansao, curvatura, universo-filho, selecao
  stability/       condicoes de energia (simbolico), potenciais mestres, robustez do ricochete
  gravitational_waves/  QNM (WKB), ringdown no tempo, ecos
  thermodynamics/  horizontes, genealogia (hipotese)
  information/     curvas de Page (EXPLORATORIO)
  ml/              PINNs e benchmarks
  raytracing/      imagens (CPU); web/index.html (GPU)
  observations/    restricoes com fonte, ponte MODEL/PREDICTION/CONSTRAINT/RESIDUAL/UNCERTAINTY/SOURCE
  parameter_space/ varreduras, Monte Carlo, adaptativo, classificacao
  figures/         toda a visualizacao (separada da fisica)
  validation.py    registro de verificacoes (limites, conservacao, convergencia, benchmarks)
```

## Reprodutibilidade

Cada figura leva um rodapé com modelo, parâmetros, versão e commit; cada JSON de resultado tem um
`.meta.json` com `RunRecord` (parâmetros, solver, semente, versões, plataforma, timestamp).
`scripts/reproduce.py configs/<nome>.json` reexecuta um estudo a partir de um arquivo de configuração.
`python -m semente.validation` roda o registro de validações e grava um relatório.
