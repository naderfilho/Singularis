# 12 — Physics-Informed Neural Networks

Módulos: `semente/ml/pinn.py`, `semente/ml/benchmarks.py`. Figuras 11, 38. Testes: `tests/test_physics.py`,
`tests/test_phase9.py`.

## Objective

Resolver equações diferenciais gravitacionais por redes neurais informadas por física e **medir** o
resultado contra a solução analítica e o solver numérico. O ML é ferramenta de verificação; a validação
matemática é feita pelo solver e pela solução exata, nunca pela rede.

## Mathematical formulation

Perda geral $\mathcal L = \mathcal L_{\rm physics} + \lambda\mathcal L_{\rm boundary} + \mu\mathcal L_{\rm data}$.
Nas PINNs deste pacote $\mathcal L_{\rm boundary} = 0$ por construção (condições iniciais impostas na
parametrização: $w(t) = w_0 + \dot w_0 t + t^2 N_\theta(t)$); $\mu > 0$ liga o modo híbrido com pontos do solver.

Problema 1 (ricochete de torção, variável de volume $w = a^3$):
$w\ddot w - \tfrac23\dot w^2 + 4\pi\rho_0 w - 16\pi\sigma = 0$, vínculo $\dot w^2 = 24\pi(\rho_0 w - \sigma)$,
solução exata $w = \sigma/\rho_0 + 6\pi\rho_0(t-t_b)^2$. Problema 2 (órbita de fóton): $u'' + u = 3Mu^2$.

Descoberta metodológica (docs/03_nascimento): em $a$ a rede "trapaceia" levando $a\to0$ (multiplicar o
resíduo por $a^6$ cria um atrator espúrio); em $w = a^3$ a equação é polinomial e bem condicionada.

## Numerical method

MLP tanh 4×64, Adam com currículo temporal (cosseno) + L-BFGS; pontos de colocação uniformes; autograd para
derivadas. Métricas (`benchmarks.py`): erro absoluto e relativo vs exato, erro de conservação (vínculo de
Friedmann normalizado), violação do resíduo da EDO, estabilidade entre sementes, tempo.

## Results (ρ₀ = 1, σ = 0.05, t ≤ 1.2, 2000 épocas + 200 L-BFGS, 4 sementes)

| métrica | só física | híbrido (μ = 1, 16 pontos) |
|---|---|---|
| max |a_PINN − a_exato| (média ± dp entre sementes) | 1.18×10⁻⁵ ± 6×10⁻⁷ | 1.20×10⁻⁵ ± 6×10⁻⁷ |
| erro de conservação (vínculo) | 4×10⁻⁷ | 4.5×10⁻⁷ |
| violação do resíduo da EDO | 3.5×10⁻⁶ | 3.5×10⁻⁶ |
| tempo (CPU) | 16 s | 17 s |
| solver DOP853 vs exato | 9×10⁻¹¹ em < 0.1 s | — |

Leitura: a PINN reproduz o ricochete a 10⁻⁵ de forma estável entre sementes, e os dados do solver **não
melhoram** a solução (a física já a determina). O solver numérico é 10⁵ vezes mais preciso e 100× mais
rápido: a PINN tem valor como método independente de verificação e como base para problemas inversos, não
como substituto.

## Limitations

Só EDOs (não EDPs); sem quantificação de incerteza bayesiana; sem problemas inversos ainda (estimar σ a
partir de dados é a extensão natural).

## References

Raissi, Perdikaris & Karniadakis (2019) J. Comput. Phys. 378, 686; Krishnapriyan et al. (2021) NeurIPS
(falhas de PINNs e currículo); Wang, Sankaran & Perdikaris (2022) (causalidade em PINNs).
