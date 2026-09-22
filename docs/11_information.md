# 11 — Informação e curva de Page (EXPLORATÓRIO)

Módulo: `semente/information/page_curve.py`. Figura 36. Testes: `tests/test_thermo_info.py`.

**Status: exploratório / teórico.** Nada aqui resolve o paradoxo da informação. O módulo produz curvas
$S_{\rm ent}(t)$ sob hipóteses explícitas para comparar cenários e escalas de tempo.

## Objective

Colocar lado a lado, com as mesmas unidades e o mesmo fundo de evaporação, os cenários
formação → evaporação → tempo de Page → (ricochete) → universo-filho.

## Mathematical formulation

Evaporação semiclássica: $dM/dt = -\alpha/M^2$, $\alpha = 1/15360\pi$; $t_{ev} = 5120\pi M^3$;
$M(t) = M_0(1-t/t_{ev})^{1/3}$; $S_{BH} = 4\pi M^2$.

Curvas (hipóteses marcadas):
- Hawking: $S_{\rm rad} = \beta[S_{BH}(0) - S_{BH}(t)]$, $\beta = 1.5$ (Page 1983/2013), cresce sempre.
- Page: $S_{\rm ent} = \min(S_{\rm rad}, S_{BH})$; $t_{\rm Page}/t_{ev} = 1 - (\beta/(1+\beta))^{3/2} = 0.535$.
- Universo-filho: $S_{\rm ent}$ segue Hawking até $t_*$ e congela (parceiros de emaranhamento no filho,
  inacessíveis ao pai). $t_* = 0.3\,t_{ev}$ é **ilustrativo**; no cenário SEMENTE o ricochete interior ocorre
  em $\tau\sim\pi M$, i.e. $t_*\approx0$ em unidades de $t_{ev}$.

## Results

| grandeza (10 M☉) | valor |
|---|---|
| ricochete interior | $1.5\times10^{-4}$ s |
| evaporação | $2\times10^{70}$ anos |
| $t_{\rm Page}/t_{\rm ricochete}$ | $\sim10^{78}$ |

Consequência (resultado deste código, dentro das hipóteses): se o interior forma um universo-filho no
ricochete, isso acontece com o buraco negro praticamente com toda a massa. A curva de Page do horizonte
remanescente é então uma questão *independente* do destino do interior, e a entropia de emaranhamento da
radiação com o filho fica congelada do ponto de vista do pai (curva "universo-filho" na fig. 36). Isso é
uma reformulação do cenário de Frolov–Markov–Mukhanov (1990), não uma resolução.

## Limitations

Sem contagem de espécies na evaporação; sem correções de greybody; sem ilhas/réplicas (a fórmula de Page
mínima é um proxy da unitariedade); sem dinâmica do remanescente. Tudo permanece especulativo.

## References

Hawking (1975) CMP 43, 199; Page (1976) PRD 13, 198; Page (1993) PRL 71, 3743; Page (2013) JCAP 09, 028;
Almheiri et al. (2021) RMP 93, 035002; Frolov, Markov & Mukhanov (1990) PRD 41, 383.
