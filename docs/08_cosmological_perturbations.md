# 08 — Perturbações cosmológicas e expansão emergente

Módulos: `semente/cosmology/modes.py` (solver genérico), `semente/cosmology/perturbations.py` (referência
exata para poeira), `semente/cosmology/expansion.py`, `semente/observations/bridge.py`.
Figuras 20, 29, 30. Tabela de comparação: `output/29_perturbacoes_ponte.md`. Testes: `tests/test_perturbations.py`.

## Objective

Depois do ricochete, determinar (i) o espectro primordial de um campo de teste (= modos tensoriais),
(ii) a razão tensor/escalar sob a hipótese padrão do matter bounce, (iii) o índice espectral e o running,
(iv) se emerge uma fase de expansão acelerada — e comparar com Planck/BICEP com residuais e fontes.

## Mathematical formulation

Equação de Mukhanov–Sasaki $v'' + (k^2 - a''/a)v = 0$ em tempo conforme, sobre qualquer
`BackgroundSolution` (GR/EC/LQC), com $a''/a = a^2(H^2 + \ddot a/a)$ do fundo numérico.
Vácuo de Bunch–Davies no passado remoto da contração (com $\eta$ medido a partir do ponto singular da lei de
potência da contração, obtido por ajuste).

Extração do modo constante após o ricochete: no regime tardio $a = a_1(\eta-\eta_s)^p$, $p = 1/(\epsilon-1)$
(poeira 2, radiação 1), as soluções exatas são $v = \sqrt{x}[c_1J_\nu(x) + c_2Y_\nu(x)]$, $\nu = p-\tfrac12$;
$(c_1,c_2)$ de $(v,v')$ no fim da grade e $\phi_{\rm const} = c_1 k^p 2^{-\nu}/(\Gamma(\nu+1)a_1)$
(o ramo $J$ é o modo constante). Validado contra a extração analítica de poeira (`BounceSpectrum`):
**concordância de 0.05% a 1.3%** em $P(k)$ para $k\le0.06\,k_b$; $a_1 = 2\pi/3$ recuperado com 10⁻⁶.

Espectros ($G=1$, $M_{Pl}^2 = 1/8\pi$): $P_t = 64\pi\,(k^3/2\pi^2)|v/a|^2$;
$P_\zeta = (4\pi/\epsilon_{\rm exit})(k^3/2\pi^2)|v/a|^2$ com $z = a\sqrt{2\epsilon}M_{Pl}$;
$r = P_t/P_\zeta = 16\,\epsilon(t_{\rm exit})$.

Running: $\alpha_s = d n_s/d\ln k$ do ajuste quadrático no platô. Não-gaussianidade: **não calculada**;
valor de literatura $f_{NL}^{\rm local} = -35/8$ (Cai et al. 2009) marcado como tal.

Expansão: $\epsilon = -\dot H/H^2$, $N = \int H\,dt$, fases com $\ddot a>0$ após $t_b$; inflação sse $N\ge60$.

## Assumptions (marcadas)

1. Campo de teste sem massa = tensores; $\zeta$ estimada com $\epsilon$ na saída do horizonte na contração
   e assumida transportada pelo ricochete como o campo de teste ($\zeta$ é mal definida em $H=0$).
2. Bunch–Davies na contração.
3. Fundo isotrópico (ver docs/07 para o problema das anisotropias).

## Numerical method

`bounce_background`: localiza $t_b$ por evento ($\dot a = 0$) e amostra o fundo numa grade logarítmica
em torno de $t_b$ (resolve o ricochete e alcança $|k\eta_0|\gg1$: $\eta_0 = -97.6$ para $a_0 = 2\times10^4$).
Modos: DOP853, rtol 10⁻⁹, passo ≤ 0.2/k.

## Validation

| teste | resultado |
|---|---|
| lei de potência tardia $p=2$, $a_1 = 2\pi/3$ | 10⁻⁶ |
| $P(k)$ vs referência exata (poeira) | < 1.3% |
| $\epsilon_{\rm exit} = 3/2$ e $r = 24$ para poeira | exato |
| radiação: $n_s - 1 \approx 2$ (dualidade de Wands) | $n_s = 2.9$ |
| ausência de inflação emergente; duração ∝ $\rho_*^{-1/2}$ | passa |
| regras de status da ponte (consistent/tension/excluded/not_predicted/literature) | passa |

## Results

| modelo | observável | previsto | observado | residual | status |
|---|---|---|---|---|---|
| poeira + torção (SEMENTE) | $n_s$ | 1.01 ± 0.02 (sist.) | 0.9665 ± 0.0038 | +2.4σ | tension |
| poeira + torção | $\alpha_s$ | ≈ 0.03 ± 0.05 | −0.0045 ± 0.0067 | +0.6σ | consistent |
| poeira + torção | $r$ | 24 | < 0.036 | ×670 | **excluded** |
| matter bounce (literatura) | $f_{NL}^{\rm local}$ | −4.4 | −0.9 ± 5.1 | −0.7σ | literature:consistent |
| radiação + LQC | $n_s$ | 2.9 | 0.9665 | +19σ | **excluded** |
| qualquer ricochete (k=0) | inflação emergente | $N_{\rm acel}\approx0.3$–$0.5$ | — | — | not predicted |

Leitura: a contração de poeira produz quase-invariância de escala **sem inflação** (mecanismo de Wands),
mas (i) a inclinação vermelha observada não sai da poeira pura, (ii) a versão mínima prevê $r=24$, excluído
por quase três ordens de grandeza, (iii) uma contração dominada por radiação está excluída, e (iv) não há
fase inflacionária emergente em nenhum dos modelos ($N<1$, independente de $\rho_*$). O cenário mínimo
"colapso → ricochete → nosso universo" **não sobrevive** aos tensores sem física adicional na contração ou
no ricochete. Este é um resultado negativo e é preservado como tal.

## Limitations

$\zeta$ através do ricochete depende do modelo (assunção 1); sem não-gaussianidade calculada; sem
perturbações de matéria/entropia; sem curvatura $k=+1$ nos modos (a contração de OS é fechada; o efeito de
$k$ nos modos de grande escala não está incluído); modos escalares em EC com torção não derivados
(o termo de spin pode alterar a equação de $\zeta$).

## References

Mukhanov, Feldman & Brandenberger (1992) Phys. Rep. 215, 203; Wands (1999) PRD 60, 023507;
Brandenberger & Peter (2017) Found. Phys. 47, 797; Cai, Xue, Brandenberger & Zhang (2009) JCAP 05, 011;
Planck 2018 VI, X; BICEP/Keck (2021) PRL 127, 151301.
