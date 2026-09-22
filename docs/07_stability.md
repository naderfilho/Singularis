# 07 — Condições de energia e estabilidade

Módulos: `semente/stability/energy_conditions.py`, `semente/stability/perturbations.py`,
`semente/stability/bounce.py`. Figuras 24–26. Testes: `tests/test_stability.py`.

## Objective

Automatizar (i) a verificação das condições de energia de qualquer solução e (ii) a classificação de
estabilidade linear de métricas estáticas e do ricochete homogêneo, com critérios documentados.

## Mathematical formulation

### Condições de energia

Para $ds^2 = -f\,dt^2 + dr^2/f + R^2 d\Omega^2$ o tensor de Einstein é derivado simbolicamente (sympy)
para $f(r)$, $R(r)$ genéricos; o fluido efetivo em base ortonormal é
$\rho = -G^t{}_t/8\pi$, $p_r = G^r{}_r/8\pi$, $p_t = G^\theta{}_\theta/8\pi$.
NEC: $\rho+p_r\ge0$ e $\rho+p_t\ge0$; WEC: NEC e $\rho\ge0$; SEC: NEC e $\rho+p_r+2p_t\ge0$;
DEC: $\rho\ge|p_r|$, $\rho\ge|p_t|$.

Para fundos FRW: $\rho_{\rm eff} = 3(H^2+k/a^2)/8\pi$, $p_{\rm eff} = -(2\ddot a/a + H^2 + k/a^2)/8\pi$.

Forma fechada obtida pelo código para Simpson–Visser (coincide com SV 2019):
$\rho = -\dfrac{\ell^2(R-4M)}{8\pi R^5}$, $p_r = -\dfrac{\ell^2}{8\pi R^4}$, com $R=\sqrt{r^2+\ell^2}$.
Logo $\rho + p_r = -\ell^2(2R-4M)/(8\pi R^5) < 0$ para $R > 2M$: **NEC violada em todo o exterior**.

### Estabilidade linear de métricas estáticas

Equação mestre $-\partial_t^2\Psi + \partial_{r_*}^2\Psi - V\Psi = 0$ com
$V_0 = f\,\ell(\ell+1)/R^2 + (f/R)\,(fR')'$ (escalar), $V_1 = f\,\ell(\ell+1)/R^2$ (EM),
$V_2 = f[\ell(\ell+1)/r^2 - 6M/r^3]$ (axial gravitacional, **só vácuo**, Regge–Wheeler).

**Critério:** modos $e^{-i\omega t}$ com $\omega^2<0$ crescem exponencialmente; existem sse o operador
$H = -d^2/dr_*^2 + V$ tem autovalor negativo (estado ligado). Se $V\ge0$ em todo lugar, $H\ge0$.
Numericamente: menor autovalor por diferenças finitas em grade uniforme em $r_*$ (Dirichlet nas bordas,
que só pode *subestimar* instabilidades de longo alcance; por isso o domínio é $\pm60M$), com convergência
em três resoluções.

Classes: `stable`; `weakly_unstable` (tempo de crescimento $>100M$); `strongly_unstable` ($\le100M$);
`numerically_unresolved` (sem convergência).

### Ricochete homogêneo

- Perturbação homogênea $\delta a$: é a translação temporal $\delta a = \dot a\,\delta t$ (gauge). O código
  verifica que ela resolve a equação linearizada (resíduo $10^{-4}$, ruído de diferenças finitas) e
  **não** a usa como critério.
- Anisotropia (cisalhamento de Bianchi I, $\sigma^2\propto a^{-6}$): $\Sigma_b = \sigma_0^2 a_b^{-6}/\rho_*$.
  Robusto se $\Sigma_b\ll1$; anisotropia inicial máxima $\sigma_0^2 < \rho_* a_b^6/a_0^6$. Critério
  conservador: o cisalhamento **não** é regularizado no modelo isotrópico (em LQC de Bianchi I ele é
  limitado, Gupt & Singh 2012).

## Assumptions

Campos de teste desacoplados da matéria efetiva (para SV, o setor gravitacional axial depende da resposta
da matéria e não está implementado). Perturbações lineares. Condições de contorno de Dirichlet em $\pm60M$.

## Validation

| teste | resultado |
|---|---|
| Schwarzschild: $\rho=p_r=p_t=0$ (10⁻¹⁴) e todas as condições satisfeitas | passa |
| Reissner–Nordström: $\rho=-p_r=p_t=Q^2/8\pi r^4$ (10⁻¹⁰) | passa |
| SV com $\ell\to0$ → vácuo; $\rho$ fechada de SV 2019 (10⁻⁹) | passa |
| $V_2$ de Schwarzschild = Regge–Wheeler | passa |
| Schwarzschild spin 0, 1, 2: `stable` | passa |
| SV spin 0/1, $\ell/M\in\{0.5, 2.5\}$: `stable`, convergido, $V_{\min}>0$ | passa |
| ricochete FRW: NEC violada só num intervalo em torno de $t_b$; GR puro satisfaz | passa |

## Results

- **Simpson–Visser (fig. 24):** para todo $\ell>0$, NEC/WEC/SEC/DEC violadas em todo o exterior ($|r|>r_h$);
  dentro do horizonte $\rho+p_r>0$. Isso é o conteúdo de matéria exótica exigido pela regularização — não é
  um defeito numérico, é o preço da geometria (consistente com o teorema de Penrose).
- **Campos de teste em SV (fig. 25):** potenciais escalares e EM positivos em todo o domínio para
  $\ell/M\in[0.1, 3.5]$; menor autovalor $\approx +0.008$ (limitado pela caixa de $\pm60M$, i.e. sem
  estado ligado). Classificação `stable` para buraco negro regular e para buraco de minhoca. Consistente
  com a literatura de campos de teste em black-bounces (Churilova & Stuchlík 2020).
- **Ricochete de torção (fig. 25c):** $\rho_{\rm eff}+p_{\rm eff} < 0$ para $t\in[0.173, 0.276]$ em torno
  de $t_b = 0.224$ (unidades de Planck, $\sigma=0.05$): o ricochete exige violação da NEC pelo fluido
  efetivo total, como todo ricochete com $k=0$.
- **Anisotropia (fig. 26):** para $a_0/a_b = 2.7$ e $\rho_*=20\rho_0$, o ricochete isotrópico é robusto
  se $\sigma_0^2 < 5\times10^{-2}\rho_0$. Em colapsos reais (contração por fatores $10^{14}$ para uma estrela
  de 10 M☉) o ganho $a^{-6}$ torna esse limite **extremamente** restritivo: $\sigma_0^2/\rho_0 < 10^{-84}$.
  Resultado negativo importante: **sem um mecanismo que limite o cisalhamento, o ricochete isotrópico do
  modelo SEMENTE não é robusto** para as razões de contração de um colapso estelar.

## Limitations

Sem perturbações gravitacionais polares/axiais acopladas à matéria efetiva; sem análise não linear;
caixa finita em $r_*$; cisalhamento tratado como fluido rígido não regularizado.

## References

Hawking & Ellis (1973); Visser (1995) *Lorentzian Wormholes*; Regge & Wheeler (1957) PR 108, 1063;
Chandrasekhar (1983); Bronnikov, Konoplya & Zhidenko (2012) PRD 86, 024028; Churilova & Stuchlík (2020)
CQG 37, 075014; Simpson & Visser (2019) JCAP 02, 042; Gupt & Singh (2012) PRD 85, 044011;
Erickson et al. (2004) PRD 69, 063514.
