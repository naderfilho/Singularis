# 09 — Ondas gravitacionais: ringdown, QNMs e ecos

Módulos: `semente/gravitational_waves/qnm.py`, `semente/gravitational_waves/ringdown.py`,
`semente/observations/bridge.py`. Figuras 31, 32. Tabela: `output/32_gw_ponte.md`. Testes: `tests/test_gw.py`.

## Objective

Quantificar as assinaturas de ondas gravitacionais que distinguiriam um black-bounce / buraco de minhoca de
um buraco negro clássico: frequências características, tempos de amortecimento, ecos, e comparar com um
evento real (GW150914) sem afirmar detectabilidade sem análise.

## Mathematical formulation

- **QNMs (WKB 3ª ordem, Iyer & Will 1987)** para os potenciais mestres de `stability/perturbations.py`;
  derivadas até 6ª ordem em $r_*$ por ajuste polinomial no pico. Validado contra Leaver (1985):
  erro 0.2% para spin 0, 1, 2 ($\ell=2$) e $\ell=3$.
- **Kerr** ($\ell=m=2$, $n=0$): ajuste publicado de Berti–Cardoso–Will (2006), marcado como literatura.
- **Domínio do tempo**: $-\partial_t^2\Psi + \partial_{r_*}^2\Psi - V\Psi = 0$, leapfrog de 2ª ordem, Sommerfeld
  nas bordas, pulso gaussiano inicial. O ringdown é ajustado por $Ae^{-\gamma t}\cos(\omega t+\phi)$ numa
  janela que começa 25 M após a chegada do pulso refletido (estudo de convergência: 0.1–0.6% em $\omega$,
  1–2% em $\gamma$ para Regge–Wheeler).
- **Ecos**: só definidos para potenciais com duas barreiras; atraso previsto $\Delta t = 2\,\Delta r_*$ entre
  os picos (Cardoso, Franzin & Pani 2016), verificado na forma de onda.
- Conversão para Hz/ms por `core.units.Scale`; massas de GW150914 no referencial do detector, $M(1+z)$.

## Assumptions

Campos de teste (spin 0/1) nas métricas com matéria efetiva; spin 2 axial só em vácuo; perturbações
lineares; sem rotação (até a fase 7); templates de Kerr vêm de ajuste da literatura.

## Validation

| teste | resultado |
|---|---|
| WKB vs Leaver (4 modos) | < 0.5% |
| $\ell\to0$ recupera Schwarzschild; $\ell = 10^{-7}M$: desvio $<10^{-6}$ | passa |
| ajuste BCW em $\chi=0$ vs Schwarzschild | 1.5% |
| ringdown no domínio do tempo vs WKB (freq. 3%, amort. 15%) | passa |
| barreira única: sem ecos; buraco de minhoca $\ell=2.5M$: duas barreiras e trem de ecos com espaçamento ≈ atraso previsto | passa |
| GW150914: Schwarzschild sem rotação excluído; Kerr (literatura) consistente | passa |

## Results

| grandeza | valor |
|---|---|
| desvio de frequência do QNM fundamental (spin 0/1, $\ell_{\rm mult}=2$) para black-bounce $\ell/M\in[0, 1.95]$ | $\lesssim 0.5\%$ |
| desvio do amortecimento | até −35% perto de $\ell\to2M$ |
| $\ell = R_b$ fixado pela junção ($\sim10^{-14}M$ para 10 M☉) | desvio $<10^{-6}$: **indistinguível** do clássico |
| buraco de minhoca $\ell=2.5M$: atraso de eco | 28.7 M (previsto) vs ≈ 20–37 M (medido, trem irregular por interferência) |
| GW150914 ($M_{\rm det}=67.6$ M☉): Schwarzschild sem rotação | $f = 179$ Hz vs $251\pm8$: **excluded** (falta rotação) |
| GW150914: Kerr $\chi=0.67$ (BCW, literatura) | $f = 250$ Hz, $\tau = 4.0$ ms: literature:consistent |

**Leitura.** (i) A rotação domina o ringdown; qualquer teste do black-bounce com dados reais exige a fase 7.
(ii) Para a garganta fixada pela física do interior, o modelo **não prevê** nenhuma assinatura de ringdown
detectável — um resultado negativo. (iii) Ecos só aparecem no regime de buraco de minhoca ($\ell>2M$), que
a junção de Israel não seleciona; se observados, apontariam para outra física. Não afirmamos detectabilidade:
não há análise de razão sinal/ruído aqui.

## Limitations

Sem perturbações gravitacionais polares/axiais acopladas à matéria efetiva; sem rotação; sem análise de
SNR/detectabilidade; o trem de ecos do buraco de minhoca de SV é irregular (o potencial tem picos assimétricos
e interferência), e a identificação automática é heurística.

## References

Iyer & Will (1987) PRD 35, 3621; Leaver (1985) Proc. R. Soc. A 402, 285; Berti, Cardoso & Will (2006) PRD 73,
064030; Cardoso, Franzin & Pani (2016) PRL 116, 171101; Churilova & Stuchlík (2020) CQG 37, 075014;
LIGO/Virgo (2016) PRL 116, 061102; 221101.
