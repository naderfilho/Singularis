# Experimentos de fronteira — o diferencial do projeto

Quatro testes numéricos sobre perguntas em aberto, feitos *dentro* dos modelos
que o projeto implementa. Nenhum resolve a pergunta. Cada um extrai uma
consequência quantitativa que, até onde sabemos, não estava calculada nesta
forma. Código: `semente/fronteira.py`. Figuras 16–19. Resultados brutos em
`output/fronteira_resultados.json`.

Os números abaixo foram produzidos por `python -m semente.figures_fronteira`
e são reproduzíveis.

---

## A. O buraco negro pai não pode continuar sendo Schwarzschild durante o ricochete

**Pergunta aberta:** quando a singularidade é substituída por um ricochete, o que
acontece com a geometria *fora* da estrela? Os modelos de ricochete (Popławski,
LQC) tratam só o interior homogêneo; o exterior fica implícito.

**Teste.** Pegamos o interior exato de Oppenheimer–Snyder com torção
(FRW fechado, $H^2 + 1/a^2 = \tfrac{8\pi}{3}(\rho - s\,a^{-6})$), integramos o
colapso, o ricochete em $R_b$ e a re-expansão, e impusemos as condições de
junção de Israel com dois exteriores candidatos.

Para uma superfície FRW, a curvatura extrínseca interior é exatamente
$K^\theta{}_\theta = \cos\chi_0 / R$ (segue da equação de Friedmann). Para o
exterior Schwarzschild, $K^\theta{}_\theta = \sqrt{1+\dot R^2 - 2M/R}\,/R$. A
densidade de energia da camada fina necessária é

$$\sigma(\tau) = \frac{\cos\chi_0 - \sqrt{1+\dot R^2 - 2M/R}}{4\pi R}.$$

No colapso clássico $\sigma\equiv0$ (teste automatizado). Com torção:

| resultado | valor (M = 1, R₀ = 8M, R_b = 0.4M) |
|---|---|
| intervalo em que $1+\dot R^2-2M/R<0$ (junção **impossível** para qualquer σ) | τ ∈ [24.92, 25.34] M, em torno do ricochete em τ = 25.13 M |
| raios onde falha | R ∈ [0.40, 0.63] M, tudo dentro do horizonte |

Dentro do horizonte não existe superfície com $\dot R = 0$ em Schwarzschild.
Logo **um ricochete dentro de um buraco negro obriga o exterior a deixar de ser
vácuo de Schwarzschild** perto do ricochete.

**Segundo candidato:** o black-bounce de Simpson–Visser, com
$K^\theta{}_\theta = (r/R^2)\sqrt{f + \dot r^2}$, $R = \sqrt{r^2+\ell^2}$. A
varredura em ℓ (fig. 16c) dá:

| ℓ / R_b | fração do tempo com junção impossível |
|---|---|
| 0.5 | 0.72 % |
| 0.9 | 0.32 % |
| 0.99 | 0.10 % |
| **1.00** | **0** |
| > 1 | 100 % (ℓ maior que o raio mínimo do interior) |

**Só $\ell = R_b$ funciona.** O parâmetro "livre" de Simpson–Visser fica
determinado pelo interior: a garganta *tem* de coincidir com o raio areal do
ricochete. Com $\rho_b = 4m^2/\pi$ (torção de férmions de massa $m$):

$$\ell^3 = \frac{3M}{4\pi\rho_b}\qquad\Rightarrow\qquad
\ell(10\,M_\odot) = 5.0\times10^{-10}\ \mathrm{m},\quad
\ell(\text{Sgr A*}) = 3.7\times10^{-8}\ \mathrm{m},\quad
\ell(\text{M87*}) = 4.3\times10^{-7}\ \mathrm{m}.$$

Uma garganta sub-nanométrica para um buraco negro estelar, crescendo como
$M^{1/3}$. E a camada de fronteira no instante do ricochete tem massa
$m_{\rm shell} = R_b\cos\chi_0$ (limite analítico na garganta, onde
$K^\theta{}_\theta$ exterior se anula): para 10 M☉, cerca de $7\times10^{17}$ kg,
ou $10^{-13}$ da massa do pai. Esse é o "preço" energético de colar um
universo-filho a um buraco negro regular.

*Ressalvas.* Tratamento efetivo (a torção entra como fluido efetivo em GR; as
condições de junção próprias de Einstein–Cartan podem adicionar termos de spin
na superfície). Não é uma prova de que o exterior *é* Simpson–Visser; é a prova
de que, dentro dessa família de exteriores estáticos, só um membro é compatível.

---

## B. Se o horizonte limita a entropia do filho, o pai do nosso universo pesa ≥ 5×10¹³ M☉

**Pergunta aberta:** a entropia de Bekenstein–Hawking $S = 4\pi M^2$ (em
unidades de Planck) é o máximo de informação que cabe dentro de um buraco negro?
(Hipótese holográfica; problema da informação.)

**Teste.** Se sim, tudo que nasce dentro do buraco negro — o universo-filho
inteiro — carrega no máximo $S_{\rm pai}$. Usando as entropias do universo
observável de Egan & Lineweaver (2010):

| componente do filho | $S$ (k_B) | $M_{\rm pai} \ge m_{Pl}\sqrt{S/4\pi}$ |
|---|---|---|
| fótons da CMB | 2.0×10⁸⁹ | 1.4×10⁶ M☉ |
| neutrinos cósmicos | 5.2×10⁸⁹ | 2.2×10⁶ M☉ |
| buracos negros estelares | 5.9×10⁹⁷ | 4.7×10¹⁰ M☉ |
| buracos negros supermassivos | 3.1×10¹⁰⁴ | 5.4×10¹³ M☉ |

Se o nosso universo nasceu num buraco negro e o limite vale, **o pai tinha pelo
menos ~5×10¹³ massas solares** — mil vezes o maior buraco negro conhecido.
(Nota: a entropia total do universo observável hoje é dominada pelos buracos
negros supermassivos; a CMB sozinha exige só um pai de ~10⁶ M☉. O piso vem do
que o filho *produziu depois*, o que já é uma tensão interessante: um filho
não pode gerar mais entropia de horizonte do que herdou.)

**Consequência para a genealogia.** Se cada filho respeita
$\sum_i m_i^2 \le M_{\rm pai}^2$, com $N$ buracos negros iguais por geração a
massa cai como $m_{k+1} = m_k/\sqrt N$ e a linhagem acaba na massa de Planck em

$$k_{\max} = \frac{\ln(M_0/m_{Pl})}{\ln\sqrt N}.$$

| $M_0$ | N = 10² | N = 10¹⁰ | N = 10²⁰ (≈ nosso universo) |
|---|---|---|---|
| 10 M☉ | 39 | 7.8 | **3.9** |
| 10⁶ M☉ | 44 | 8.8 | 4.4 |
| 10¹⁰ M☉ | 48 | 9.6 | 4.8 |

Uma árvore de universos com fecundidade parecida com a nossa tem **menos de
cinco gerações** de profundidade sob o limite holográfico.

---

## C. Seleção natural cosmológica com orçamento herdado perde a alavanca

**Pergunta aberta:** a seleção natural cosmológica de Smolin funciona? Ela
precisa de muitas gerações com fecundidade diferencial grande.

**Teste.** Repetimos a dinâmica populacional de `selection.py`, mas cada
universo herda capacidade $C = m_{\rm pai}^2$ (o quadrado da massa do buraco
negro que o gerou) e só os filhos que cabem ($N m^2 \le C$) se reproduzem.

Resultado (fig. 18): a fecundidade *efetiva* cai de ~10³ para **~1 filho por
universo já na primeira geração** e fica lá por 300 gerações. Sem diferença de
fecundidade não há seleção — só deriva neutra (a massa típica dos buracos
negros cai lentamente porque universos com buracos mais leves cabem mais
filhos). O oposto do que a hipótese de Smolin precisa.

Ou seja: **CNS e limite holográfico são quase incompatíveis.** Quem defende a
CNS precisa que a informação dentro do buraco negro *não* seja limitada pela
área — o que é uma posição sobre o problema da informação, não só sobre cosmologia.

---

## D. Dados: a previsão original de Smolin está excluída

Smolin (1992) previu que a massa máxima de estrelas de nêutrons ficaria perto
de 1.6 M☉ (o mínimo que a física permite, porque isso maximiza a produção de
buracos negros); depois argumentou por ~2 M☉.

| pulsar | massa (M☉) | σ acima de 1.6 | σ acima de 2.0 | ref. |
|---|---|---|---|---|
| PSR J0740+6620 | 2.08 ± 0.07 | 6.9 | 1.1 | Fonseca+ 2021 |
| PSR J0348+0432 | 2.01 ± 0.04 | 10.3 | 0.25 | Antoniadis+ 2013 |
| PSR J1614−2230 | 1.908 ± 0.016 | 19.3 | −5.8 | Arzoumanian+ 2018 |
| PSR J0952−0607 | 2.35 ± 0.17 | 4.4 | 2.1 | Romani+ 2022 |

A previsão original está excluída por > 4σ em todos; a revisada está no limite
(J0952−0607 a 2σ). Isso é conhecido na comunidade; o valor aqui é ter o teste
reproduzível ao lado das simulações.

---

## O que é novo, em uma frase cada

1. **ℓ não é livre**: a garganta de um buraco negro regular formado por colapso
   com ricochete tem de ter o raio do ricochete, $\ell^3 = 3M/4\pi\rho_b$, e a
   parede entre os universos pesa $R_b c^2/G$.
2. **Piso de massa para o pai do nosso universo** (5×10¹³ M☉) e **profundidade
   máxima da árvore** (< 5 gerações) sob o limite holográfico.
3. **CNS + holografia ≈ incompatíveis**: com orçamento herdado a seleção vira deriva.

Cada uma depende de uma hipótese explícita (tratamento efetivo; limite
holográfico). É exatamente onde a física está indecisa — e é por isso que são
testes, não teoremas.
