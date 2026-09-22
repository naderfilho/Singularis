# 13 — Explorador de espaço de parâmetros

Módulo: `semente/parameter_space/explorer.py`. Figura 37. Resultados: `output/37_scan_estatico.json`,
`output/37_monte_carlo_ricochete.json` (com `.meta.json` de proveniência). Testes: `tests/test_phase9.py`.

## Objective

Varrer automaticamente parâmetros ($M$, $\ell$, $Q$, $\rho_*$, $w$, anisotropia inicial) e classificar as
soluções de forma reprodutível, com varreduras 1D/2D, Monte Carlo com semente e refinamento adaptativo
nas fronteiras de fase.

## Classes (rótulos do modelo, não afirmações físicas)

`singular`, `classical_black_hole`, `regular_black_hole`, `black_bounce`, `wormhole`, `naked_singularity`,
`unstable_bounce`, `stable_bounce`, `expanding_baby_universe`, `numerically_unresolved`.

## Classifiers

- `classify_static(M, ℓ, Q)`: via `geometry.static.classify` (número de horizontes, sinal de $f$ na garganta,
  regularidade do centro). Critérios em docs/02.
- `classify_bounce(ρ_*, w, σ₀²)`: fundo LQC efetivo (k = 0): GR → `singular`; falha do solver ou violação do
  vínculo > 10⁻⁵ → `numerically_unresolved`; $\Sigma_b \ge 1$ → `unstable_bounce` (critério BKL, docs/07);
  senão `expanding_baby_universe` se $H>0$ no fim, ou `stable_bounce`.

## Methods

`scan_1d`, `scan_2d` (grade), `monte_carlo` (uniforme/log-uniforme, `seed`), `adaptive_refine` (subdivide
células cujos cantos têm rótulos diferentes; cache de avaliações). `ScanResult.save` grava pontos, rótulos,
semente e `RunRecord` (commit, versões, timestamp).

## Results (fig. 37)

- Família estática $(Q,\ell)$: as cinco fases de docs/02 aparecem na varredura 40×40; o refinamento
  adaptativo reproduz as fronteiras com ~150 avaliações em vez de 1600.
- Ricochete LQC (poeira), Monte Carlo 400 pontos, semente 11: a fronteira `expanding_baby_universe` /
  `unstable_bounce` é a curva $\sigma_0^2 = \rho_* a_b^6$ do critério BKL; nenhuma amostra é
  `numerically_unresolved` (vínculo preservado a 10⁻⁹).

## Reproducibility

Mesma semente → mesmos pontos e rótulos (teste); `RunRecord` no sidecar responde "como este mapa foi
produzido" (parâmetros, limites, semente, commit).

## Limitations

Classificação por heurísticas explícitas (não por invariantes globais como diagramas de Penrose calculados);
o custo cresce com classificadores dinâmicos (o de ricochete integra uma EDO por ponto).
