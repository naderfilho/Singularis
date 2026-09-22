# 05 — Einstein–Cartan (torção e spin)

Módulos: `semente/quantum/einstein_cartan.py` (física), `semente/cosmology/friedmann.py` (dinâmica),
`semente/quantum/comparison.py` (GR vs EC vs LQC), `semente/cosmology/seed_universe.py` (progenitor → filho).
Figuras 06, 27, 28. Testes: `tests/test_quantum.py`, `tests/test_friedmann.py`.

Classificação: **modelo teórico existente na literatura**, implementado como estratégia de densidade
efetiva; os números são resultados numéricos do modelo.

## Objective

Quantificar a contribuição da torção associada ao spin da matéria na dinâmica de Friedmann e comparar com
a relatividade geral e com a LQC para condições iniciais equivalentes.

## Mathematical formulation

Teoria ECSK com fluido de spin de Weyssenhoff (spins não polarizados), unidades de Planck:

$$H^2 + \frac{k}{a^2} = \frac{8\pi}{3}(\rho - \alpha n^2),\qquad
\frac{\ddot a}{a} = -\frac{4\pi}{3}(\rho + 3p - 4\alpha n^2),\qquad \alpha = \frac{\kappa}{32} = \frac{\pi}{4}.$$

$n\propto a^{-3}$ é a densidade numérica dos férmions; $\langle s^2\rangle = n^2/8$. O termo $-\alpha n^2 \propto a^{-6}$
domina para $a$ pequeno e produz o ricochete. Para poeira de férmions de massa $m$:
$\rho_b = m^2/\alpha = 4m^2/\pi$; para nêutrons $\rho_b = 3.9\times10^{58}$ kg/m³ ($7.5\times10^{-39}\rho_{Pl}$).

Solução exata (poeira, k = 0): $a^3 = \sigma/\rho_0 + 6\pi\rho_0(t-t_b)^2$, $\sigma = \alpha n_0^2$.

## Assumptions

Fluido de spin não polarizado; torção tratada como fonte efetiva em RG (condições de junção próprias de EC
não incluídas); equação de estado de poeira até o ricochete, embora $\rho_b$ esteja muito acima da densidade
nuclear — o valor é uma ordem de grandeza, não uma previsão de precisão.

## Numerical method

`FriedmannModel` com `EinsteinCartanCorrection(sigma)`; DOP853 (rtol 10⁻¹⁰); vínculo de Friedmann
monitorado (`constraint_violation` < 10⁻⁹).

## Parameters

`rho0`, `w`, `sigma` (ou `rho_star = rho0²/sigma`), `k`; para o universo-filho: `M_kg`, `m_fermion_kg`,
`R0_over_rs`.

## Validation

| teste | resultado |
|---|---|
| integrador vs solução exata (poeira) | 10⁻⁸ (a), 10⁻¹¹ no módulo base |
| $\sigma\to0$ recupera RG | 10⁻⁹ |
| $\rho_b = 4m^2/\pi$ e $\sigma = \alpha(\rho/m)^2$ | exato |
| **equivalência EC ≡ LQC para poeira** ($H^2$ e $\ddot a/a$) | 10⁻¹⁴ |

## Results (ρ₀ = 1, ρ_* = 20, a₀ = 1, contração, k = 0; fig. 27)

| | GR | Einstein–Cartan | LQC efetiva |
|---|---|---|---|
| ricochete (poeira) | não (singular) | a_min = 0.368, ρ_max = 20 = ρ_*, t_b = 0.2246 | idêntico a EC |
| ricochete (radiação) | não | a_min = 0.224, **ρ_max = 400 = 20ρ_\*** | a_min = 0.473, ρ_max = 20 |
| curvatura máxima (|R|) poeira / radiação | ∞ | 1.5×10³ / 2.0×10⁴ | 1.5×10³ / 2.0×10³ |
| e-folds pós-ricochete (até t = 1.5) | — | 2.14 / 2.50 | 2.14 / 1.77 |
| Σ_b (σ₀² = 10⁻⁶) | — | 2×10⁻⁵ / 4×10⁻⁴ | 2×10⁻⁵ / 4×10⁻⁶ |

**Resultado matemático deste projeto:** para poeira, EC e LQC efetiva são a *mesma* dinâmica de um
parâmetro ($\rho_c \leftrightarrow \rho_b$). A distinção física aparece só para matéria relativística: em
EC a correção escala com $n^2\propto a^{-6}$ e a densidade de energia pode exceder $\rho_b$ (radiação:
$\rho_{\max} = \rho_0^2/\rho_b\cdot$…), enquanto em LQC $\rho\le\rho_c$ sempre.

**Progenitor → filho (fig. 28):** $R_b = (3M/4\pi\rho_*)^{1/3}$; para 10 M☉, $R_b = 5\times10^{-10}$ m (EC,
nêutrons) ou $1.3\times10^{-22}$ m (LQC). O tempo próprio horizonte → ricochete é $\approx\pi M$ em ambos
($6.6\times10^{-5}$ s para 10 M☉), pois $R_b\ll 2M$.

## Limitations

Fluido de spin idealizado; sem polarização; sem acoplamento spin–torção além do termo quadrático médio;
equação de estado irrealista no regime de ricochete; a robustez a anisotropias (docs/07) é um problema
sério para razões de contração estelares.

## References

Hehl, von der Heyde, Kerlick & Nester (1976) RMP 48, 393; Gasperini (1986) PRL 56, 2873;
Popławski (2010) PLB 694, 181; Popławski (2012) PRD 85, 107502; Popławski (2016) ApJ 832, 96.
