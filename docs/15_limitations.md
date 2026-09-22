# 15 — Limitações e resultados negativos

Este documento é a lista honesta do que o laboratório **não** faz, do que assumiu, e do que descobriu que
**não** funciona. Resultados negativos são resultados.

## Resultados negativos produzidos pelo código

1. **Tensores excluem o cenário mínimo.** A contração de poeira dá $n_s\approx1$ sem inflação, mas
   $r = 16\epsilon = 24$ contra $r < 0.036$: a versão mínima está excluída por quase três ordens de
   grandeza (docs/08, 14).
2. **Não há inflação emergente.** Em todo o espaço de parâmetros varrido (ρ_*, w) a fase acelerada
   pós-ricochete rende $N < 1$ e-fold (docs/08).
3. **Anisotropias.** Sem um mecanismo que limite o cisalhamento, o ricochete isotrópico não é robusto
   para razões de contração estelares ($\sigma_0^2/\rho_0 < 10^{-84}$) (docs/07).
4. **Nenhuma assinatura de ringdown.** Para $\ell = R_b$ fixado pela junção, o desvio de QNM é $<10^{-6}$
   (docs/09).
5. **Universo-filho pequeno.** Só com poeira o filho é fechado, pequeno e efêmero (docs/03_nascimento,
   `SeedUniverse`).
6. **Primeira lei falha para SV com $S = A/4$** (docs/10).
7. **Seleção cosmológica + holografia ≈ incompatíveis** (docs/02_fronteira C) e a previsão original de
   Smolin está excluída pelos pulsares (docs/14).

## Hipóteses que sustentam os resultados positivos

- Tratamento efetivo: torção e LQC entram como fluidos efetivos em RG; as condições de junção próprias de
  Einstein–Cartan não estão incluídas.
- Dinâmica de camada com densidade média (extrapolação para $E\ne0$ em LTB).
- Exterior de Simpson–Visser é *candidato*, não derivado da evolução.
- $\zeta$ transportada pelo ricochete como o campo de teste (matter bounce).
- Campos de teste desacoplados da matéria efetiva nas análises de estabilidade e QNM.
- $S = A/4$; evaporação semiclássica sem espécies.
- Limite holográfico no filho (hipótese especulativa marcada).

## O que não está implementado

- Estabilidade de geometrias rotativas (Teukolsky); perturbações gravitacionais com matéria efetiva.
- Perturbações escalares em Einstein–Cartan com spin (a equação de $\zeta$ pode mudar).
- Curvatura $k=+1$ nos modos e na LQC efetiva.
- Não-gaussianidade (só valor de literatura).
- Análise de detectabilidade (SNR) de ecos/QNMs.
- Dinâmica completa do exterior durante o ricochete (matéria efetiva evoluída).
- Correções quânticas à entropia; ilhas/réplicas na curva de Page.
- Equação de estado realista (a poeira é irrealista no regime de ricochete).

## Interpretações que o repositório recusa

- "O raio de Schwarzschild do universo é igual ao seu tamanho" é uma identidade de Friedmann.
- "Vemos o outro universo dentro da sombra" vale só na geometria eterna.
- "Consistent" na ponte observacional significa "não excluído", nunca "confirmado".

## O que ainda poderia salvar o cenário (próximos testes)

Refazer os modos com $k=+1$ e com matéria relativística na contração; incluir um campo extra ou uma fase
com $w\ne0$ para baixar $r$; incluir cisalhamento regularizado (Bianchi I em LQC); derivar as condições de
junção em Einstein–Cartan; evoluir o exterior. Cada um é um estudo, não uma afirmação.
