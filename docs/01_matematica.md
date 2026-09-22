# A matemática por trás do SEMENTE

Este documento acompanha o código. Cada seção aponta para o módulo que
implementa a equação e para a figura que a mostra. Unidades geométricas
$G = c = 1$; quando $\hbar$ entra, unidades de Planck.

---

## 1. O buraco branco vem "de brinde": Kruskal–Szekeres

A métrica de Schwarzschild

$$ds^2 = -\Big(1-\frac{2M}{r}\Big)dt^2 + \frac{dr^2}{1-2M/r} + r^2 d\Omega^2$$

é singular em $r=2M$ apenas por culpa das coordenadas. Kruskal (1960) e
Szekeres (1960) trocaram $(t,r)$ por $(T,X)$ com

$$X^2 - T^2 = \Big(\frac{r}{2M}-1\Big)e^{r/2M}, \qquad \frac{T}{X}=\tanh\frac{t}{4M}\ \ (\text{fora}),$$

e a métrica vira $ds^2 = \frac{32M^3}{r}e^{-r/2M}(-dT^2+dX^2)+r^2d\Omega^2$,
perfeitamente regular em $r=2M$. O preço: o plano $(T,X)$ tem **quatro**
regiões, não duas:

| região | onde | o que é |
|---|---|---|
| I | $X>\lvert T\rvert$ | nosso exterior |
| II | $T>\lvert X\rvert$ | buraco **negro** (tudo cai para $r=0$ no futuro) |
| III | $T<-\lvert X\rvert$ | buraco **branco** (tudo sai de $r=0$ no passado) |
| IV | $X<-\lvert T\rvert$ | um segundo exterior, causalmente desconectado |

Inverter $(T,X)\to r$ exige a função W de Lambert: $r/2M - 1 = W\big((X^2-T^2)/e\big)$.
A singularidade $r=0$ é a hipérbole $T^2-X^2=1$ (o ponto de ramificação $W(-1/e)=-1$).

*Código:* `semente/geometry.py` (`Schwarzschild.kruskal_from_tr`, `r_from_kruskal`).
*Figura:* `output/01_kruskal.png`, `output/02_penrose.png`.

**Honestidade obrigatória.** Num colapso real (seção 3) a estrela ocupa tudo
"à esquerda" da sua superfície no diagrama, e as regiões III e IV
simplesmente não existem. O buraco branco só existe na solução de vácuo
*eterna*, ou — e esse é o ponto do projeto — se a singularidade for
substituída por um ricochete (seções 4 e 5).

---

## 2. O interior de um buraco negro **é** um universo (Kantowski–Sachs)

Para $r<2M$, $g_{tt}>0$ e $g_{rr}<0$: $r$ virou tempo e $t$ virou espaço.
Escrevendo $T\equiv r$:

$$ds^2 = -\frac{dT^2}{2M/T-1} + \Big(\frac{2M}{T}-1\Big)dt^2 + T^2 d\Omega^2 .$$

É uma cosmologia homogênea e anisotrópica (Kantowski–Sachs) com dois fatores de escala

$$a_\parallel(T)=\sqrt{2M/T-1},\qquad a_\perp(T)=T,$$

e tempo próprio $d\tau = dT/\sqrt{2M/T-1}$. Parametrizando $T = M(1+\cos\eta)$,
$\tau = M(\eta+\sin\eta)$: o tempo do horizonte até a singularidade é
**exatamente $\pi M$** (para um buraco de 10 massas solares, $1.5\times10^{-4}$ s).
Quando $T\to0$: $a_\parallel\to\infty$ (estica), $a_\perp\to0$ (esmaga) e o
escalar de Kretschmann $K = 48M^2/r^6\to\infty$: singularidade física.

*Código:* `semente/interior.py`. *Figura:* `output/04_interior_universo.png`.

---

## 3. Colapso de Oppenheimer–Snyder: o interior da estrela é FRW

Uma estrela de poeira homogênea que colapsa a partir do repouso com raio $R_0$
tem (Oppenheimer & Snyder 1939):

* **interior**: métrica de Friedmann fechada
  $ds^2 = -d\tau^2 + a(\eta)^2[d\chi^2+\sin^2\chi\, d\Omega^2]$, $0\le\chi\le\chi_0$,
  $a(\eta)=\tfrac{a_m}{2}(1+\cos\eta)$, $\tau(\eta)=\tfrac{a_m}{2}(\eta+\sin\eta)$;
* **exterior**: Schwarzschild com massa $M$;
* **junção** (continuidade da métrica induzida e da curvatura extrínseca):
  $\sin^2\chi_0 = 2M/R_0$, $a_m = R_0/\sin\chi_0$.

O horizonte de eventos nasce **no centro** da estrela em $\eta = \eta_h-\chi_0$
(antes de a superfície cruzar $2M$) e cresce até $r=2M$: `OppenheimerSnyder.event_horizon_interior`.

Conclusão exata, não metáfora: *dentro de um buraco negro formado por colapso há um
universo de Friedmann fechado em contração.* Inverta o tempo e é um Big Bang.

*Código:* `semente/collapse.py`. *Figura:* `output/05_colapso_oppenheimer_snyder.png`.

---

## 4. O ricochete

### 4a. Torção de Einstein–Cartan (Popławski 2010, 2012)

Na teoria de Einstein–Cartan–Sciama–Kibble o spin dos férmions gera torção,
e a torção entra nas equações de Friedmann como uma densidade **negativa**
proporcional ao quadrado da densidade numérica $n$:

$$H^2 + \frac{k}{a^2} = \frac{8\pi}{3}\big(\rho - \alpha n^2\big),\qquad
\frac{\ddot a}{a} = -\frac{4\pi}{3}\big(\rho + 3p - 4\alpha n^2\big),$$

com $\alpha = \kappa/32 = \pi/4$ para spin $\tfrac12$ (usando $s^2 = n^2/8$, $\kappa=8\pi$).
Como $n\propto a^{-3}$, o termo $\alpha n^2\propto a^{-6}$ domina em $a$ pequeno e
**inverte a contração**. Para poeira ($\rho\propto a^{-3}$, $k=0$), na variável de
volume $w=a^3$ a equação vira

$$\dot w^2 = 24\pi(\rho_0 w - \sigma)\quad\Rightarrow\quad \ddot w = 12\pi\rho_0,$$

uma **parábola**: $a(t)^3 = \sigma/\rho_0 + 6\pi\rho_0 (t-t_b)^2$. Essa solução exata
(`bounce.torsion_dust_exact`) é usada para testar o integrador (erro $10^{-11}$) e a PINN.
Ricochete em $\rho_b = 4m^2/\pi$ (unidades de Planck); para nêutrons,
$\rho_b \approx 3.9\times10^{58}\ \mathrm{kg/m^3}$.

### 4b. Cosmologia quântica de laço (Ashtekar–Pawlowski–Singh 2006)

$$H^2 = \frac{8\pi}{3}\rho\Big(1-\frac{\rho}{\rho_c}\Big),\qquad
\dot H = -4\pi(\rho+p)\Big(1-\frac{2\rho}{\rho_c}\Big),\qquad \rho_c\simeq0.41\,\rho_{Pl}.$$

Ricochete quando $\rho=\rho_c$. Mesma fenomenologia, mecanismo diferente.

*Código:* `semente/bounce.py`. *Figura:* `output/06_ricochete.png`.

### 4c. A geometria completa: black-bounce de Simpson–Visser (2019)

$$ds^2 = -f\,dt^2 + \frac{dr^2}{f} + (r^2+\ell^2)\,d\Omega^2,\qquad
f = 1-\frac{2M}{\sqrt{r^2+\ell^2}},\qquad r\in(-\infty,\infty).$$

* $\ell=0$: Schwarzschild.
* $0<\ell<2M$: buraco negro **regular**. A singularidade vira uma hipersuperfície
  espacial regular em $r=0$ (raio areal mínimo $R=\ell$). Tudo que cai atravessa
  $r=0$ e sai de um **buraco branco** no universo $r<0$. O diagrama de Penrose é
  uma escada infinita de universos.
* $\ell>2M$: buraco de minhoca atravessável.

Dentro do horizonte a análise da seção 2 se repete com $R=\sqrt{r^2+\ell^2}$:
$a_\perp$ contrai até $\ell$ e **re-expande** — o universo interno ricocheteia.
O escalar de Kretschmann é finito em todo lugar (calculado em base ortonormal em
`BlackBounce.kretschmann`, testado contra $48M^2/r^6$ em $\ell=0$).

Geodésicas nulas: $(dr/d\lambda)^2 = E^2 - L^2 f/(r^2+\ell^2)$, logo

$$\frac{d^2r}{d\lambda^2} = \frac{L^2\, r\,(R-3M)}{R^5},\qquad \frac{d\phi}{d\lambda}=\frac{L}{R^2},$$

regular nos horizontes e na garganta. É isso que o ray-tracer integra (CPU e GPU).

*Código:* `semente/geometry.py`, `semente/interior.py`, `semente/raytracer.py`, `web/index.html`.

---

## 5. Modelo SEMENTE (síntese autoral)

Junte 3 + 4a: o interior FRW exato do colapso, feito de $N=M/m$ férmions, com
torção. Então, para um buraco negro pai de massa $M$:

* raio areal do ricochete: $R_b^3 = 3M/(4\pi\rho_b)$ — sempre **muito** dentro do
  horizonte ($R_b/r_s\sim10^{-14}$ para $10\,M_\odot$);
* tempo próprio horizonte → ricochete $\approx \pi M$ (segundos para buracos estelares,
  ~1 dia para o de M87);
* o universo-filho nasce com raio próprio $a_b\chi_0$ (fração de nanômetro para
  massas estelares) e, **se for só poeira**, é um universo fechado que re-expande
  até $a_m$ e recolapsa — pequeno e efêmero.

Esse último ponto é o que as equações realmente dizem, e por isso qualquer versão
séria da hipótese precisa de um mecanismo extra de expansão (produção de partículas
pela torção em Popławski; inflação; constante cosmológica herdada). O projeto
deixa isso explícito em vez de esconder.

*Código:* `bounce.SeedUniverse`. *Figura/tabela:* `output/07_tabela_universo_filho.png/.json`.

---

## 6. Seleção natural cosmológica (Smolin 1992)

Se os filhos herdam os parâmetros do pai com pequenas mutações, a dinâmica do
replicador com fecundidade $N_{bh}(\theta)$ leva a população ao máximo de
$N_{bh}$. A paisagem usada é um brinquedo; a dinâmica é exata.
Previsão falsificável de Smolin: mudar qualquer constante deve **reduzir** o
número de buracos negros. *Código:* `semente/selection.py`. *Figura:* `output/08_selecao_cosmologica.png`.

---

## 7. A rede neural (PINN)

A rede $N_\theta$ nunca vê dados: minimiza o resíduo da equação diferencial em
pontos de colocação. Com "hard constraints" $w(t) = w_0 + \dot w_0 t + t^2 N_\theta(t)$
as condições iniciais valem exatamente. Perda:

$$\mathcal L = \big\langle (w\ddot w - \tfrac23\dot w^2 + 4\pi\rho_0 w - 16\pi\sigma)^2\big\rangle
+ \lambda\big\langle(\dot w^2 - 24\pi(\rho_0 w-\sigma))^2\big\rangle .$$

Curriculum temporal (o intervalo cresce até $t_{max}$) + Adam + L-BFGS.
Erro final $\sim10^{-5}$ contra a solução exata. Segunda PINN: órbita de fóton
$u''+u=3Mu^2$, erro $\sim10^{-4}$.

Por que a variável de volume: em $a$ a equação é rígida ($\ddot a\sim a^{-5}$) e a
rede "trapaceia" levando $a\to0$; em $w=a^3$ tudo é polinomial e bem condicionado.
(Esse foi um resultado empírico deste projeto.)

*Código:* `semente/pinn.py`. *Figura:* `output/11_pinn.png`.

---

## 8. O ray-tracer

* Cada raio vive num plano pela origem (simetria esférica): integra-se $r(\lambda),\phi(\lambda)$
  no plano e reconstrói-se a posição 3D.
* Parâmetro de impacto a partir da direção no referencial estático da câmera:
  $b = R\sin\alpha/\sqrt{f}$.
* Disco fino de Novikov–Thorne/Page–Thorne entre $R_{ISCO}=6M$ e $R_{out}$. Forma fechada
  derivada neste projeto (e testada contra a integral e o limite newtoniano), com $x=\sqrt{R/M}$:
  $$F(R)=\frac{3\dot M}{8\pi M^2}\,\frac{(x-\sqrt6)-\frac{\sqrt3}{2}\Big[\ln\frac{x-\sqrt3}{x+\sqrt3}-\ln\frac{\sqrt6-\sqrt3}{\sqrt6+\sqrt3}\Big]}{x^5(x^2-3)} .$$
  Curiosidade: em termos do raio areal $R$, toda a física de órbitas circulares do
  black-bounce ($\Omega$, $E$, $L$, ISCO) é idêntica à de Schwarzschild.
* Desvio para o vermelho/azul: $g = \sqrt{1-3M/R}\,/\,(1-\Omega\lambda_z)$, brilho $\propto g^4 F$,
  cor de corpo negro com $T_{obs}=gT_{em}$.
* Para $\ell>0$, os raios que atravessam $r=0$ mostram o céu do **outro** universo;
  para $0<\ell<2M$ isso é, na geometria eterna, luz vinda do universo anterior pelo buraco branco.

*Código:* `semente/raytracer.py` (numpy, ~1 s em 320×180) e `web/index.html` (WebGL2, tempo real).
