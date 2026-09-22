"""
EXPLORATORIO / TEORICO.  Modulo de informacao quantica: evaporacao de Hawking, curva de Page e o
cenario de universo-filho.  NAO e uma solucao do paradoxo da informacao; e um conjunto de curvas
S_ent(t) sob hipoteses explicitas, para comparar cenarios.

Evaporacao (RESULTADO MATEMATICO semiclassico, Hawking 1975; Page 1976, so fotons/gravitons e sem
dependencia detalhada de especies):
    dM/dt = -alpha / M^2,  alpha = 1/(15360 pi)  (unidades de Planck)
    t_ev = 5120 pi M^3,    M(t) = M0 (1 - t/t_ev)^{1/3},   S_BH(t) = 4 pi M(t)^2

Curvas de entropia de emaranhamento da radiacao (HIPOTESES explicitamente marcadas):
    Hawking:   S_rad(t) = beta [S_BH(0) - S_BH(t)]   (radiacao termica; beta ~ 1.5 e o excesso de entropia
               da radiacao sobre a perda de entropia do horizonte, Page 1983/2013; cresce sempre)
    Page:      S_ent(t) = min( S_rad_Hawking(t), S_BH(t) )   (evaporacao unitaria: a entropia de
               emaranhamento nao pode exceder a do buraco negro; t_Page onde as curvas se cruzam)
    Universo-filho: S_ent segue Hawking ate t_*, quando o interior deixa de estar acessivel (remanescente
               ou universo-filho desconectado), e depois fica CONSTANTE em S_rad(t_*): os parceiros
               de emaranhamento estao no filho; do ponto de vista do universo-pai a informacao nao volta.

Escalas de tempo (RESULTADO DESTE CODIGO): para um buraco negro estelar o ricochete interior ocorre em
tau ~ pi M ~ 1e-4 s depois da formacao, e a evaporacao leva ~1e67 anos.  Logo, no cenario SEMENTE, o
destino do interior e decidido ~1e78 vezes antes do tempo de Page: o universo-filho nasce enquanto o buraco
negro ainda tem essencialmente toda a massa, e a curva de Page (se existe) diz respeito ao horizonte que fica.

Referencias: Hawking (1975) CMP 43, 199; Page (1976) PRD 13, 198; Page (1993) PRL 71, 3743; Page (2013)
JCAP 09, 028; Almheiri, Hartman, Maldacena, Shaghoulian & Tajdini (2021) RMP 93, 035002 (revisao); Frolov,
Markov & Mukhanov (1990) PRD 41, 383 (universos-filhos dentro de buracos negros).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.units import M_PLANCK, M_SUN, T_PLANCK

ALPHA_EVAP = 1.0 / (15360 * np.pi)


def evaporation_time_planck(M_planck: float) -> float:
    return 5120 * np.pi * M_planck**3


def evaporation_time_years(M_solar: float) -> float:
    Mp = M_solar * M_SUN / M_PLANCK
    return evaporation_time_planck(Mp) * T_PLANCK / 3.15576e7


def mass_of_time(t, M0, t_ev):
    return M0 * np.clip(1 - np.asarray(t, float) / t_ev, 0.0, None) ** (1 / 3)


@dataclass
class PageCurves:
    t_over_tev: np.ndarray
    S_BH: np.ndarray
    S_hawking: np.ndarray
    S_page: np.ndarray
    S_baby: np.ndarray
    t_page_over_tev: float
    t_star_over_tev: float
    S0: float
    beta: float
    assumptions: tuple = ("radiacao termica com beta = 1.5", "unitariedade: S_ent <= S_BH (curva de Page)",
                          "universo-filho: S_ent congela em t_* (parceiros inacessiveis)")
    status: str = "EXPLORATORIO"


def page_curves(M0_planck: float = 1e3, beta: float = 1.5, t_star_over_tev: float | None = None, n=2000) -> PageCurves:
    t_ev = evaporation_time_planck(M0_planck)
    x = np.linspace(0, 1, n)
    M = mass_of_time(x * t_ev, M0_planck, t_ev)
    S_BH = 4 * np.pi * M**2
    S0 = 4 * np.pi * M0_planck**2
    S_h = beta * (S0 - S_BH)
    S_p = np.minimum(S_h, S_BH)
    t_page = 1 - (beta / (1 + beta)) ** 1.5
    if t_star_over_tev is None:
        # por padrao, o instante em que o ricochete interior (tau ~ pi M) e "resolvido" e imediato: t_* ~ 0.
        # Para VISUALIZAR o cenario usamos t_* = 0.3 t_ev como caso ilustrativo (marcado).
        t_star_over_tev = 0.3
    S_b = np.where(x < t_star_over_tev, S_h, S_h[np.searchsorted(x, t_star_over_tev)])
    return PageCurves(x, S_BH, S_h, S_p, S_b, float(t_page), float(t_star_over_tev), float(S0), beta)


def timescales(M_solar: float) -> dict:
    """Ricochete interior (~pi M) vs tempo de Page vs evaporacao, em segundos e anos."""
    from ..core.units import Scale
    s = Scale.from_solar_masses(M_solar)
    t_bounce = np.pi * s.time_unit
    t_ev_s = evaporation_time_years(M_solar) * 3.15576e7
    t_page_s = (1 - (1.5 / 2.5) ** 1.5) * t_ev_s
    return dict(M_solar=M_solar, t_bounce_s=t_bounce, t_page_s=t_page_s, t_evap_s=t_ev_s,
                t_page_over_t_bounce=t_page_s / t_bounce, t_evap_years=evaporation_time_years(M_solar))
