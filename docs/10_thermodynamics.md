# 10 — Termodinâmica de buracos negros

Módulos: `semente/thermodynamics/horizon.py`, `semente/thermodynamics/genealogy.py`. Figuras 17, 35.
Testes: `tests/test_thermo_info.py`.

## Objective

Calcular área, entropia, temperatura e a primeira lei para as geometrias do repositório; acompanhar a
entropia do horizonte aparente durante o colapso e o ricochete; relacionar a massa do progenitor aos graus
de liberdade do universo-filho. Separar o que é derivado do modelo do que é conjectura.

## Mathematical formulation

$A = 4\pi R(r_h)^2$, $S = A/4$, $T = \kappa/2\pi$ com $\kappa = |f'(r_h)|/2$ (família $-f\,dt^2+dr^2/f+R^2d\Omega^2$).
Primeira lei testada numericamente: $1 - T\,dS/dM$ a parâmetros restantes fixos. Kerr: $S = \pi(r_+^2+a^2)$,
$dM = T\,dS + \Omega_H dJ$. Horizonte aparente dinâmico: $S_{AH}(\tau) = \pi R_{AH}^2$ (Hayward 1994).

## Results

| geometria | $S$ | $T$ | $1 - T\,dS/dM$ | status |
|---|---|---|---|---|
| Schwarzschild | $4\pi M^2$ | $1/8\pi M$ | 0 | resultado matemático |
| Reissner–Nordström (Q fixo) | $\pi r_+^2$ | $(r_+-r_-)/4\pi r_+^2$ | 0 | resultado matemático |
| Kerr (J fixo) | $\pi(r_+^2+a^2)$ | $\kappa/2\pi$ | 0 | resultado matemático |
| **Simpson–Visser** ($\ell$ fixo) | $4\pi M^2$ (**independe de** $\ell$: $R(r_h)=2M$) | $\sqrt{4M^2-\ell^2}/16\pi M^2$ | $1-\sqrt{1-\ell^2/4M^2}$ | resultado deste código |

O resíduo de SV é exato (verificado a 10⁻³): com $S = A/4$ a primeira lei falha. Duas leituras possíveis,
não decididas aqui: (i) há um termo de trabalho $\Phi_\ell\,d\ell$ no estilo de Smarr; (ii) a entropia
de um black-bounce não é $A/4$. Para $\ell = R_b\sim10^{-14}M$ o resíduo é $\sim10^{-29}$: irrelevante na
prática, mas conceitualmente presente.

**Colapso dinâmico (fig. 35b):** o teorema da área (Hawking 1971) é verificado para o horizonte de
**eventos** clássico: $S_{EH}(\tau) = \pi R_{EH}^2$ cresce monotonicamente ao longo do gerador lançado do
centro e chega a $4\pi M^2$ na superfície (teste: 5%). O horizonte **aparente** de Oppenheimer–Snyder nasce
na superfície e se propaga para dentro ($R_{AH} = 2m(\chi_{AH})$ decresce): sua área não obedece ao
teorema, e isso é clássico e conhecido — o código o registra em vez de esconder. No ricochete a região
presa é transiente e $S_{AH}$ volta a zero; o horizonte de eventos clássico deixa de estar definido
(docs/03), consistentemente com a violação da NEC (docs/07). Resultado do modelo, não conjectura.

**Escalas de tempo (fig. 35c):** para 10 M☉, ricochete interior $\sim\pi M = 1.5\times10^{-4}$ s; evaporação
$\sim10^{70}$ anos; tempo de Page $\sim0.54\,t_{ev}$. O destino do interior é decidido $\sim10^{78}$ vezes
antes do tempo de Page.

**Graus de liberdade progenitor → filho:** $S_{BH} = 4\pi M^2$ vs entropia de ordem $N$ da matéria:
para 10 M☉, $N/S_{BH}\sim10^{-20}$. A matéria que forma o filho carrega uma fração ínfima dos graus de
liberdade do horizonte; a genealogia sob o limite holográfico está em docs/02_fronteira (seção B) e
`thermodynamics/genealogy.py` — **hipótese especulativa** marcada.

## Second law and speculation boundary

- Segunda lei generalizada no ricochete: **não** verificada nem verificável aqui (exigiria a entropia da
  matéria efetiva que viola a NEC e a entropia do universo-filho).
- Relação "massa do progenitor ↔ graus de liberdade do filho": só o balanço de contagem acima é derivado;
  qualquer identificação $S_{\rm filho} \le S_{BH}$ é conjectura sobre gravidade quântica.

## Limitations

Entropia $A/4$ assumida (pode não valer fora de RG); sem termodinâmica do black-bounce rotativo; sem
contribuições quânticas (correções logarítmicas).

## References

Bekenstein (1973) PRD 7, 2333; Bardeen, Carter & Hawking (1973) CMP 31, 161; Hawking (1975) CMP 43, 199;
Hayward (1994) PRD 49, 6467; Simpson & Visser (2019) JCAP 02, 042.
