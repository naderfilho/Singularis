# 02 — Geometrias

Módulos: `semente/geometry/spherical.py` (Schwarzschild, Kruskal, Penrose, Flamm, Simpson–Visser),
`geometry/static.py` (família estática genérica), `geometry/charged.py` (Reissner–Nordström, black-bounce
carregado), `geometry/kerr.py` (Kerr, black-bounce rotativo), `geometry/geodesics.py`, `geometry/interior.py`.
Figuras 01–04, 10, 33, 34. Testes: `test_physics.py`, `test_geometries.py`.

## Classificação de cada geometria

| geometria | status | fonte |
|---|---|---|
| Schwarzschild, Kruskal–Szekeres, Penrose, Flamm | resultado matemático | Kruskal 1960; Szekeres 1960 |
| Simpson–Visser (black-bounce) | solução exata com matéria efetiva (viola NEC) | Simpson & Visser 2019 |
| Reissner–Nordström | resultado matemático (eletrovácuo) | — |
| black-bounce carregado | modelo da literatura | Franzin, Liberati, Mazza, Simpson & Visser 2021 |
| Kerr | resultado matemático (vácuo) | Kerr 1963; BPT 1972 |
| black-bounce rotativo ($r\to\sqrt{r^2+\ell^2}$ em Kerr) | modelo da literatura | Mazza, Franzin & Liberati 2021 |

## Família estática genérica (`StaticSphericalMetric`)

$ds^2 = -f\,dt^2 + dr^2/f + R(r)^2d\Omega^2$ com $f$, $R$ como expressões sympy. Fornece horizontes
(raízes de $f$), esferas de fótons (extremos de $f/R^2$), gravidade superficial $|f'(r_h)|/2$, Kretschmann
em base ortonormal com derivadas simbólicas, e classificação causal (singularidade nua / BH / BH com
horizonte de Cauchy / black-bounce / minhoca) pelo sinal de $f$ na garganta e pelo número de horizontes.

## Carga

RN: $r_\pm = M\pm\sqrt{M^2-Q^2}$; $K = 8(6M^2r^2 - 12MQ^2r + 7Q^4)/r^8$; fluido efetivo
$\rho = -p_r = p_t = Q^2/8\pi r^4$ (todas as condições de energia satisfeitas).

Black-bounce carregado: $f = 1 - 2M/R + Q^2/R^2$, $R = \sqrt{r^2+\ell^2}$. Horizontes em $R = R_\pm$,
i.e. $r_h = \pm\sqrt{R_\pm^2 - \ell^2}$. Fases (fig. 33): $\ell < R_-$: dois pares de horizontes (garganta
tipo-tempo entre horizontes de Cauchy); $R_- < \ell < R_+$: black-bounce (garganta espacial, ricochete);
$\ell > R_+$: buraco de minhoca; $Q > M$: sem horizontes (regular para $\ell>0$). Kretschmann finito na
garganta para todo $\ell>0$ (verificado). NEC violada em faixas (fig. 33c).

## Rotação

Kerr: horizontes $r_\pm$, ergosuperfície $r_E(\theta)$, arrasto $\omega = 2Mar/[(r^2+a^2)^2 - a^2\Delta\sin^2\theta]$,
$\Omega_H$, $\kappa$, área, ISCO de BPT (testado: 6M, 1M, 9M), Kretschmann fechado (testado vs simbólico, 10⁻⁸).

Black-bounce rotativo: mesma métrica com $r\to\sqrt{r^2+\ell^2}$ em $\Sigma$ e $\Delta$. Horizontes só se
$r_\pm > \ell$ ($r_h = \pm\sqrt{r_\pm^2-\ell^2}$); ergosuperfície $\sqrt{r^2+\ell^2} = M+\sqrt{M^2-a^2\cos^2\theta}$.
O Kretschmann é obtido **simbolicamente** (Riemann completo em Boyer–Lindquist, ~50 s na primeira chamada)
e é finito na garganta $r=0$ (testado para $a=0.6M$, $\ell=0.7M$), reduzindo-se a Kerr para $\ell\to0$.

## O que NÃO está implementado (limitações declaradas)

- Estabilidade de Kerr/black-bounce rotativo (equação de Teukolsky): não implementada; o ringdown de Kerr
  usa o ajuste publicado de Berti–Cardoso–Will (2006).
- Perturbações gravitacionais acopladas à matéria efetiva em black-bounces: não implementadas.
- Condições de energia do black-bounce rotativo: não calculadas (o tensor de Einstein da métrica rotativa
  não está lambdificado; Mazza et al. reportam violação da NEC).
- O ISCO do black-bounce rotativo é o de Kerr em $\sqrt{r^2+\ell^2}$: aproximação documentada, não derivada.

## References

Simpson & Visser (2019) JCAP 02, 042; Franzin, Liberati, Mazza, Simpson & Visser (2021) JCAP 07, 036;
Mazza, Franzin & Liberati (2021) JCAP 04, 082; Bardeen, Press & Teukolsky (1972) ApJ 178, 347;
Kerr (1963) PRL 11, 237.
