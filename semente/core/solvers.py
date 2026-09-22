"""
Configuracao unica de integradores de EDO.  Toda integracao do pacote passa por
`integrate`, para que tolerancias e metodo sejam registrados na proveniencia.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class ODESettings:
    method: str = "DOP853"
    rtol: float = 1e-10
    atol: float = 1e-12
    max_step: float = np.inf
    dense_output: bool = True

    def as_dict(self) -> dict:
        d = asdict(self)
        d["max_step"] = None if not np.isfinite(self.max_step) else self.max_step
        return d


DEFAULT_ODE = ODESettings()
FAST_ODE = ODESettings(rtol=1e-8, atol=1e-10)


@dataclass
class ODESolution:
    """Solucao amostrada numa grade uniforme, com metadados do solver."""

    t: np.ndarray
    y: np.ndarray                      # shape (n_state, n_t)
    settings: ODESettings
    n_rhs_evals: int
    status: int
    message: str
    events: dict = field(default_factory=dict)

    def as_metadata(self) -> dict:
        return dict(solver=self.settings.as_dict(), n_rhs_evals=self.n_rhs_evals, status=self.status,
                    message=self.message, t_span=[float(self.t[0]), float(self.t[-1])], n_samples=int(self.t.size))


def integrate(rhs: Callable, t_span: Sequence[float], y0: Sequence[float], n_samples: int = 4000,
              settings: ODESettings = DEFAULT_ODE, events: Optional[list] = None, t_eval=None) -> ODESolution:
    """Integra dy/dt = rhs(t, y) e amostra em `n_samples` pontos uniformes (ou em `t_eval`).

    Eventos terminais encurtam o intervalo; a amostragem respeita o fim efetivo.
    """
    kw = dict(method=settings.method, rtol=settings.rtol, atol=settings.atol, dense_output=True)
    if np.isfinite(settings.max_step):
        kw["max_step"] = settings.max_step
    sol = solve_ivp(rhs, tuple(t_span), np.asarray(y0, float), events=events, **kw)
    t_end = sol.t[-1]
    t = np.asarray(t_eval, float) if t_eval is not None else np.linspace(t_span[0], t_end, n_samples)
    t = t[(t >= min(t_span[0], t_end)) & (t <= max(t_span[0], t_end))] if t_eval is not None else t
    y = sol.sol(t)
    ev = {}
    if events:
        for i, e in enumerate(events):
            name = getattr(e, "__name__", f"event{i}")
            ev[name] = [float(x) for x in sol.t_events[i]]
    return ODESolution(t=t, y=y, settings=settings, n_rhs_evals=int(sol.nfev), status=int(sol.status),
                       message=str(sol.message), events=ev)
