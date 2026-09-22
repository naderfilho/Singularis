# SEMENTE — Cosmologia de Buracos Negros

*English name for the repository:* **`black-hole-genesis`** (SEMENTE = "seed": the white-hole seed inside every black hole).

> *Todo buraco negro guarda a semente de um buraco branco. Todo buraco branco é um Big Bang.*
> Este projeto pega essa frase e a transforma em equações resolvidas, redes neurais treinadas
> e imagens calculadas — sem esconder onde a física termina e a especulação começa.

![render](docs/img/13_render_black_bounce.png)

*Imagem calculada, não desenhada: buraco negro regular (black-bounce, ℓ = 1M) com disco de
acreção de Page–Thorne. Dentro da "sombra" aparece o céu do outro universo, visto através do
buraco branco. Geodésicas exatas, redshift gravitacional + Doppler, cor de corpo negro.*

---

## A ideia, em uma cadeia de fatos

| # | afirmação | status | onde no código |
|---|---|---|---|
| 1 | A extensão maximal de Schwarzschild contém, obrigatoriamente, um buraco **branco** e um segundo universo. | **teorema** (Kruskal 1960) | `geometry.py`, fig. 01, 02 |
| 2 | Dentro do horizonte, $r$ é tempo: o interior de um buraco negro **é** uma cosmologia (Kantowski–Sachs) que dura exatamente $\pi M$. | **teorema** | `interior.py`, fig. 04 |
| 3 | O interior de uma estrela em colapso **é** um universo de Friedmann fechado em contração (Oppenheimer–Snyder). | **teorema** (1939) | `collapse.py`, fig. 05 |
| 4 | Num colapso real, o buraco branco e o outro universo **não existem**: a estrela ocupa o lugar deles. | **teorema** | fig. 01 (região hachurada) |
| 5 | Se a singularidade for substituída por um ricochete, a contração do item 3 vira uma expansão — um Big Bang. Torção de Einstein–Cartan e cosmologia quântica de laço fazem isso. | **modelo publicado** (Popławski 2010; Ashtekar et al. 2006) | `bounce.py`, fig. 06 |
| 6 | Existe uma métrica exata (Simpson–Visser 2019) em que o buraco negro atravessa uma garganta regular e sai como buraco branco em **outro universo**; o interior contrai até $R=\ell$ e re-expande. | **solução exata** com matéria efetiva exótica | `geometry.py`, `raytracer.py`, `web/`, fig. 03, 04, 10, 12–15 |
| 7 | Para um buraco negro pai de massa $M$, dá para calcular densidade, tamanho e tempo do "Big Bang" do filho. | **modelo SEMENTE** (síntese deste projeto) | `bounce.SeedUniverse`, fig. 07 |
| 8 | Se cada filho herda as constantes do pai com mutações, a população de universos evolui para maximizar a produção de buracos negros (Smolin). | **hipótese falsificável** | `selection.py`, fig. 08 |
| 9 | Uma rede neural que só vê as equações (PINN) reconstrói o ricochete com erro $10^{-5}$. | **resultado deste projeto** | `pinn.py`, fig. 11 |

Tudo isso está derivado em [`docs/01_matematica.md`](docs/01_matematica.md).

![kruskal](docs/img/01_kruskal.png)

---

## O diferencial: experimentos sobre perguntas em aberto

Quatro testes numéricos feitos dentro dos modelos acima, com resultados que, até onde sabemos,
não estavam calculados nesta forma (`semente/fronteira.py`, figuras 16–19,
[`docs/02_fronteira.md`](docs/02_fronteira.md)):

| # | pergunta aberta | o que o teste encontrou |
|---|---|---|
| A | O que acontece **fora** da estrela quando o interior ricocheteia? | A junção de Israel com Schwarzschild é impossível num intervalo finito em torno do ricochete (nenhuma camada de energia salva). Na família black-bounce, **só ℓ = R_b** funciona: o parâmetro de Simpson–Visser fica fixado pelo interior, $\ell^3 = 3M/4\pi\rho_b$. Para 10 M☉, ℓ ≈ 5×10⁻¹⁰ m; a parede entre os universos pesa $R_b c^2/G \approx 7\times10^{17}$ kg. |
| B | O horizonte do pai limita a entropia do filho? | Se sim, o pai do nosso universo tinha ≥ 5×10¹³ M☉, e uma árvore de universos com a nossa fecundidade tem **menos de 5 gerações** ($k_{max} = \ln(M_0/m_{Pl})/\ln\sqrt N$). |
| C | A seleção natural cosmológica funciona com esse limite? | Não: a fecundidade efetiva cai para ~1 filho por universo na primeira geração e a dinâmica vira deriva neutra. **CNS e holografia são quase incompatíveis.** |
| D | A previsão de Smolin (M_max de estrelas de nêutrons ≈ 1.6 M☉) sobrevive aos dados? | Excluída por > 4σ em quatro pulsares; a versão revisada (2 M☉) está no limite. |

Cada resultado depende de uma hipótese explícita (tratamento efetivo da torção; limite
holográfico). É onde a física está indecisa, e por isso são testes, não teoremas.

![juncao](docs/img/16_fronteira_juncao.png)

---

## O modelo contra o nosso universo

Nada confirma que nascemos de um buraco negro. O que dá para fazer é o teste de consistência:
o modelo prevê propriedades do universo-filho; medimos o nosso (Planck 2018, BICEP/Keck 2021);
comparamos (`semente/nascimento.py`, figuras 20–21, [`docs/03_nascimento.md`](docs/03_nascimento.md)):

| previsão do modelo | previsto | observado | veredito |
|---|---|---|---|
| Espectro primordial: resolvemos Mukhanov–Sasaki através do ricochete de torção com vácuo de Bunch–Davies na contração de poeira | $n_s = 1.003$ (platô, < 1% de variação em 1.2 décadas) | $n_s = 0.9665 \pm 0.0038$ | **compatível**: quase-invariância de escala **sem inflação**; a inclinação de 3% precisa de correção |
| Curvatura: o filho é fechado | $\Omega_k < 0$ | Planck+BAO $0.0007\pm0.0019$; Planck só $-0.011\pm0.0065$ | **compatível**, falsificável |
| Massa do pai a partir da curvatura + conservação | $\ge 4.7\times10^{23}\,M_\odot$ | massa do universo observável ~ $10^{23}\,M_\odot$ | consistente (não é evidência) |
| Razão tensor/escalar do ricochete de matéria mínimo | $r = 24$ | $r < 0.036$ | **falsificado** na versão mínima |
| Expansão acelerada | não prevista | $\Lambda$ domina | não previsto |

![espectro](docs/img/20_nascimento_espectro.png)

Placar: quatro compatíveis, um falsificado, um não previsto. O modelo está **vivo e restrito**.
O próximo teste natural (aberto no repositório) é refazer o espectro com $k=+1$ e tensores para ver se $r$ cai.

---

## Instalação e uso

```bash
pip install -r requirements.txt
python -m pytest -q tests          # 14 testes de consistencia fisica
python scripts/run_all.py          # gera tudo em ./output (~3 min em CPU)
```

Abra `web/index.html` num navegador com WebGL2 para o ray-tracer em tempo real:
arraste **ℓ** de 0 (Schwarzschild) até 3 e veja a singularidade virar uma janela para outro universo.

Uso como biblioteca:

```python
from semente.geometry import BlackBounce
from semente.interior import radial_infall_black_bounce
from semente.bounce import SeedUniverse, M_SUN
from semente.raytracer import render, Scene, Camera, save_png

bb = BlackBounce(M=1.0, l=0.8)
print(bb.kind, bb.horizons)                     # horizontes em +-sqrt(4M^2 - l^2)
q = radial_infall_black_bounce(l=0.8, r0=8.0)   # r(tau) atravessa r=0 e emerge no outro lado
print(SeedUniverse(10 * M_SUN).table())         # o universo-filho de um BN de 10 massas solares
save_png(render(Scene(l=0.8), Camera(width=960, height=540)), "meu_buraco.png")
```

---

## O que cada figura mostra

| arquivo | conteúdo |
|---|---|
| `01_kruskal.png` | As quatro regiões de Kruskal, a estrela de Oppenheimer–Snyder (que apaga o buraco branco) e um observador em queda. |
| `02_penrose.png` | Penrose de Schwarzschild eterno e a "escada" de universos do black-bounce. |
| `03_flamm_wormhole.png` | Paraboloide de Flamm (ponte de Einstein–Rosen) e a garganta lisa de um buraco de minhoca. |
| `04_interior_universo.png` | Os fatores de escala do universo dentro do buraco negro; com ricochete; Kretschmann finito; queda através da garganta. |
| `05_colapso_oppenheimer_snyder.png` | Colapso: superfície, horizonte de eventos nascendo no centro, $a(\tau)$ FRW e sua versão com tempo invertido (Big Bang). |
| `06_ricochete.png` | GR pura vs torção vs LQC: $a(t)$, $H(t)$, densidade efetiva. Solução exata sobreposta. |
| `07_tabela_universo_filho.png/.json` | Propriedades do universo-filho para pais de 3 a 6.5×10⁹ massas solares. |
| `08_selecao_cosmologica.png` | Seleção natural cosmológica: fecundidade média e migração da população na paisagem. |
| `09_deflexao.png` | Ângulo de deflexão da luz vs parâmetro de impacto; divergência na esfera de fótons. |
| `10_geodesicas.png` | Órbitas de fótons e partículas (integrador geral com Christoffel simbólico). |
| `11_pinn.png` | PINN vs integrador: ricochete e órbita de fóton. |
| `12–15_render_*.png` | Ray-tracer: Schwarzschild, black-bounce, buraco de minhoca, lente vista de cima. |
| `16_fronteira_juncao.png` | Junção de Israel do ricochete: onde Schwarzschild falha, por que só ℓ = R_b funciona, ℓ previsto vs massa. |
| `17_fronteira_entropia.png` | Piso de massa do pai pelo limite entrópico; profundidade máxima da árvore de universos. |
| `18_fronteira_selecao_orcamento.png` | Seleção de Smolin com orçamento holográfico: a fecundidade efetiva colapsa para ~1. |
| `19_fronteira_estrelas_neutrons.png` | Massas de pulsares vs previsão de Smolin. |
| `20_nascimento_espectro.png` | Mukhanov–Sasaki através do ricochete: potencial, espectro com n_s ≈ 1, modos congelando. |
| `21_nascimento_veredito.png` | Massa do pai vs curvatura observada e a tabela de vereditos. |

---

## Estrutura

```
semente/
  geometry.py    Schwarzschild, Kruskal (Lambert W), Penrose, Flamm, black-bounce, Kretschmann
  geodesics.py   Christoffel via sympy -> integrador DOP853; forma orbital u(phi); deflexao
  interior.py    Kantowski-Sachs dentro do BN; ricochete do interior; queda pela garganta
  collapse.py    Oppenheimer-Snyder: juncao FRW/Schwarzschild, horizonte de eventos, Kruskal da superficie
  bounce.py      Friedmann com torcao (+ solucao exata), LQC, modelo SEMENTE, Pathria
  selection.py   dinamica populacional de Smolin
  pinn.py        PINNs (torch): ricochete na variavel de volume; orbita de foton
  raytracer.py   ray-tracer numpy: Page-Thorne, redshift, lente, dois ceus
  figures.py     todas as figuras
  fronteira.py   experimentos sobre perguntas em aberto (juncao de Israel, limite entropico,
                 selecao com orcamento, dados de pulsares) + figures_fronteira.py
  nascimento.py  previsoes para o universo-filho vs Planck/BICEP (Mukhanov-Sasaki atraves do
                 ricochete, curvatura, massa do pai, r) + figures_nascimento.py
web/index.html   ray-tracer WebGL2 em tempo real (GLSL)
tests/           27 testes (round-trips, limites classicos, solucoes exatas, PINN, fronteira, nascimento)
docs/            derivacoes
```

---

## Honestidade científica (leia antes de citar)

* Itens 1–4 da tabela são **teoremas** da relatividade geral. O item 4 é o que a maioria das
  divulgações omite: o buraco branco de Kruskal é apagado pelo colapso.
* Itens 5–6 são **modelos**: dependem de física além da relatividade clássica (torção, gravidade
  quântica de laço, matéria que viola a condição de energia nula perto da garganta). São
  soluções exatas *das equações desses modelos*, não observações.
* O argumento "o raio de Schwarzschild do universo observável é igual ao seu tamanho" é uma
  **identidade** da equação de Friedmann para universo plano ($r_s = c/H_0$ exatamente), não uma
  evidência. `bounce.observable_universe_as_black_hole` calcula e explica.
* O modelo SEMENTE, levado à risca com só poeira, gera um universo-filho **pequeno e efêmero**
  (fechado, recolapsa). Para que vire um universo como o nosso é preciso um mecanismo de expansão
  extra (produção de partículas por torção, inflação). A tabela 07 diz isso com números.
* Nada aqui foi observado. Buracos brancos nunca foram detectados. O que o projeto oferece é a
  cadeia completa, calculada e testada, entre a relatividade que conhecemos e a hipótese.

---

## Referências

* Oppenheimer, J. R. & Snyder, H. (1939). *On Continued Gravitational Contraction.* Phys. Rev. 56, 455.
* Kruskal, M. (1960). Phys. Rev. 119, 1743. Szekeres, G. (1960). Publ. Math. Debrecen 7, 285.
* Novikov, I. (1964). *Delayed explosion of a part of the Fridman universe and quasars.* Astron. Zh. 41, 1075.
* Pathria, R. K. (1972). *The Universe as a Black Hole.* Nature 240, 298.
* Frolov, V., Markov, M. & Mukhanov, V. (1990). *Black holes as possible sources of closed and semiclosed worlds.* Phys. Rev. D 41, 383.
* Smolin, L. (1992). *Did the universe evolve?* Class. Quantum Grav. 9, 173.
* Ashtekar, A., Pawlowski, T. & Singh, P. (2006). *Quantum nature of the Big Bang.* Phys. Rev. Lett. 96, 141301.
* Popławski, N. (2010). *Cosmology with torsion: An alternative to cosmic inflation.* Phys. Lett. B 694, 181.
* Popławski, N. (2012). *Nonsingular, big-bounce cosmology from spinor-torsion coupling.* Phys. Rev. D 85, 107502.
* Haggard, H. & Rovelli, C. (2015). *Black hole fireworks: quantum-gravity effects outside the horizon spark black to white hole tunneling.* Phys. Rev. D 92, 104020.
* Ashtekar, A., Olmedo, J. & Singh, P. (2018). *Quantum Transfiguration of Kruskal Black Holes.* Phys. Rev. Lett. 121, 241301.
* Simpson, A. & Visser, M. (2019). *Black-bounce to traversable wormhole.* JCAP 02, 042.
* Page, D. N. & Thorne, K. S. (1974). *Disk-Accretion onto a Black Hole.* ApJ 191, 499.
* Raissi, M., Perdikaris, P. & Karniadakis, G. (2019). *Physics-informed neural networks.* J. Comput. Phys. 378, 686.

Licença MIT. Autor: Nader Filho. Construído com auxílio de Claude (Anthropic).
