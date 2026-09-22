# 14 — Restrições observacionais (ponte consolidada)

Módulos: `semente/observations/constraints.py` (dados com fonte e incerteza),
`semente/observations/bridge.py` (regras de status), tabelas geradas em `output/29_perturbacoes_ponte.md`
e `output/32_gw_ponte.md`, pulsares em `observations/pulsars.py`, curvatura em `cosmology/curvature.py`.

Regras: `consistent` = |residual| ≤ 2σ (**apenas "não excluído"**); `tension` = 2–4σ; `excluded` > 4σ ou acima
de limite superior; `not_predicted`; `literature:` = valor da literatura, não deste código.

| modelo | previsão | valor previsto | restrição observacional | residual | incerteza (modelo) | fonte | status |
|---|---|---|---|---|---|---|---|
| poeira + torção (SEMENTE) | $n_s$ | 1.01 | 0.9665 ± 0.0038 | +2.4σ | ±0.02 (sist.) | Planck 2018 VI | tension |
| poeira + torção | $\alpha_s$ | ≈ +0.03 | −0.0045 ± 0.0067 | +0.6σ | ±0.05 | Planck 2018 X | consistent |
| poeira + torção | $r$ | 24 | < 0.036 (95%) | ×670 | — | BICEP/Keck 2021 | **excluded** |
| matter bounce | $f_{NL}^{\rm local}$ | −4.4 | −0.9 ± 5.1 | −0.7σ | — | Planck 2018 IX | literature:consistent |
| radiação + LQC | $n_s$ | 2.9 | 0.9665 | +20σ | ±0.1 | Planck 2018 VI | **excluded** |
| qualquer ricochete (k=0) | inflação emergente | $N \approx 0.3$–0.5 | — | — | — | — | not_predicted (N ≥ 60 nunca ocorre) |
| interior OS | sinal de $\Omega_k$ | < 0 | 0.0007 ± 0.0019 (Planck+BAO); −0.011 ± 0.0065 (Planck) | ≤ 1.7σ | — | Planck 2018 VI | consistent, falsificável |
| interior OS | $M_{\rm pai}$ | ≥ 4.7×10²³ M☉ | massa do universo observável ~10²³ M☉ | — | — | — | consistent (conservação, não evidência) |
| Schwarzschild sem rotação | $f_{\rm ringdown}$ GW150914 | 179 Hz | 251 ± 8 Hz | −5.2σ | ±7% (massa) | LIGO/Virgo 2016 | **excluded** (falta rotação) |
| Kerr χ = 0.67 (BCW) | $f$, τ | 250 Hz, 4.0 ms | 251 ± 8 Hz, 4.0 ± 0.3 ms | −0.1σ, +0.1σ | ±7% | LIGO/Virgo 2016 | literature:consistent |
| black-bounce ℓ = R_b (junção) | desvio de QNM | < 10⁻⁶ | — | — | — | — | not_predicted (indistinguível) |
| black-bounce ℓ > 2M | ecos | atraso 2Δr* | nenhum eco reportado com significância | — | — | Abedi et al. 2017 (contestado) | not tested here |
| seleção cosmológica (Smolin) | $M_{\max}$ estrelas de nêutrons ≈ 1.6 M☉ | 1.6 | 2.08 ± 0.07 (J0740+6620) e outros | > 4σ | — | Fonseca+ 2021 etc. | **excluded** (versão original) |

## Síntese

Nenhuma linha constitui evidência positiva. O cenário mínimo "colapso → ricochete → nosso universo" está
excluído pelos tensores ($r = 24$) e não produz inflação; as assinaturas de ringdown do black-bounce com a
garganta fixada pelo interior são indistinguíveis do clássico; a única previsão falsificável que sobrevive
é o sinal da curvatura ($\Omega_k < 0$). Ver docs/15 para o que isso implica.
