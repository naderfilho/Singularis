# Nascimento — o modelo contra o nosso universo

Pergunta do usuário: "algo que confirme o surgimento do mundo nascendo disso".
Resposta honesta: **nada confirma**. O que existe é o teste de consistência:
o modelo faz previsões para um universo que nasce do ricochete; medimos o
nosso; comparamos linha por linha. Código: `semente/nascimento.py`.
Figuras 20–21. Resultados em `output/nascimento_veredito.json`.

Dados: Planck 2018 (VI, Tabela 2, TT,TE,EE+lowE+lensing+BAO), BICEP/Keck 2021.

---

## 1. O espectro primordial sai do ricochete — sem inflação

O interior de Oppenheimer–Snyder é uma **contração dominada por poeira**. Há
um resultado clássico (Wands 1999; "matter bounce", Brandenberger e colaboradores)
de que uma contração de poeira gera um espectro invariante de escala para
perturbações que saem do horizonte. Testamos isso no *nosso* fundo específico,
com o ricochete de torção não singular $a^3 = a_b^3 + 6\pi\rho_0 t^2$:

* Equação de Mukhanov–Sasaki para um campo de teste sem massa, $v'' + (k^2 - a''/a)v = 0$,
  em tempo conforme, com $a''/a$ finito na garganta (fig. 20a).
* Vácuo de Bunch–Davies no passado remoto da contração ($k|\eta_i| \ge 22$).
* Extração exata da amplitude que sobrevive fora do horizonte depois do ricochete:
  decompomos $v$ nas soluções exatas do regime de poeira,
  $f_\pm = e^{\mp ix}(1 \mp i/x)$, e usamos $f_+ + f_- = -\tfrac23 x^2 + O(x^4)$
  para isolar o modo constante de $\phi = v/a$ (isso remove a contaminação do modo
  decrescente $\propto\eta^{-3}$, que enviesava o resultado em 17% quando avaliávamos
  ingenuamente em $k\eta = 0.1$).

Resultado numérico:

| | |
|---|---|
| $n_s$ no platô ($k \ll k_{\rm ricochete}$) | **1.003** (a_b = 0.05), **1.002** (a_b = 0.02) |
| variação de $P(k)$ ao longo de 1.2 décadas | < 1% |
| Planck 2018 | $n_s = 0.9665 \pm 0.0038$ |

O modelo produz, sem inflação, a quase-invariância de escala observada. A
inclinação vermelha de 3.4% **não** sai da poeira pura; precisa de uma
componente com $w$ ligeiramente negativo ou de correções na fase de contração.
Isso é uma previsão qualitativa correta, não um ajuste.

## 2. O sinal da curvatura

O filho é um universo **fechado** ($k=+1$): $\Omega_k < 0$, sem exceção. Isso é
falsificável.

| medida | $\Omega_k$ | leitura |
|---|---|---|
| Planck+BAO | $0.0007 \pm 0.0019$ | compatível com fechado dentro de 1σ |
| Planck sozinho | $-0.011 \pm 0.0065$ | prefere fechado a ~1.7σ (o debate "Planck evidence for a closed Universe", Di Valentino, Melchiorri & Silk 2019) |

Se DESI/Euclid medirem $\Omega_k > 0$ com significância, a versão
Oppenheimer–Snyder do modelo morre.

## 3. O tamanho do pai, pela curvatura

Universo fechado de poeira: $a_m = \frac{8\pi G}{3c^2}\rho_{m0}\,a_0^3$, com
$a_0 = c/(H_0\sqrt{-\Omega_k})$. Nosso universo observável ocupa
$\chi_{\rm obs} = D/a_0$ da esfera-3, e a estrela-mãe cobria $\chi_0 \ge \chi_{\rm obs}$.
Logo $M_{\rm pai} \ge \tfrac{a_m}{2}\sin^3\chi_{\rm obs}\,c^2/G$.

| $\Omega_k$ usado | $a_0$ mínimo | $\chi_{\rm obs}$ | $M_{\rm pai}$ mínimo | recolapso |
|---|---|---|---|---|
| $-0.0031$ (Planck+BAO, 2σ) | 2.5×10²⁷ m | 0.18 rad | **4.7×10²³ M☉** | 4×10⁴ Gyr |
| $-0.011$ (Planck só) | 1.3×10²⁷ m | 0.34 rad | 4.5×10²³ M☉ | 6×10³ Gyr |

O pai tem de ter, no mínimo, a massa do universo observável (~10²³ M☉).
Isso é conservação de massa, não evidência — mas mostra que o piso entrópico
do experimento B (5×10¹³ M☉) é irrelevante perto deste. E o universo fechado de
poeira não recolapsaria antes de milhares de Gyr: consistente com a idade.

## 4. Onde o modelo mínimo falha

| item | previsto | observado | veredito |
|---|---|---|---|
| razão tensor/escalar | $r = 16\epsilon = 24$ (ricochete de matéria mínimo) | $r < 0.036$ | **falsificado** na versão mínima |
| expansão acelerada | não há (poeira pura) | $\Lambda$ domina desde z ≈ 0.3 | não previsto; $\Lambda$ tem de ser herdada/emergente |

O problema dos tensores é o calcanhar de Aquiles de todo ricochete de matéria
e é conhecido (Cai, Easson & Brandenberger 2012; revisão Brandenberger & Peter 2017).
Escapes propostos na literatura: ricochete com curvatura positiva (o nosso caso!
o interior é $k=+1$), campos extras, ou fase de contração com $w \ne 0$. Fica como
o próximo teste natural deste projeto: refazer o cálculo do item 1 com $k=+1$ e
com os tensores, e ver se $r$ cai.

---

## Placar

| previsão | veredito |
|---|---|
| espectro quase invariante de escala sem inflação | compatível |
| universo fechado ($\Omega_k<0$) | compatível, falsificável |
| pai com massa de um universo | consistente (conservação) |
| idade vs recolapso | consistente |
| $r = 24$ | **falsificado** (versão mínima) |
| aceleração | não previsto |

Um modelo que passa em quatro testes, falha em um e não fala de outro não
está confirmado. Está **vivo e restrito** — que é o máximo que qualquer hipótese
sobre a origem do universo pode dizer de si mesma hoje.
