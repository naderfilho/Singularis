"""
Benchmarks quantitativos das PINNs:  ANALITICO  vs  SOLVER NUMERICO  vs  PINN.

Perda geral  L = L_physics + lambda L_boundary + mu L_data.  Nas PINNs deste pacote as condicoes
iniciais sao impostas por construcao (L_boundary = 0 exatamente); `data_weight > 0` liga um termo de
dados com pontos do solver numerico (modo hibrido), para medir quanto a "ajuda" de dados muda o erro.

Metricas (todas adimensionais, definidas aqui):
    abs_error        max_t |a_PINN - a_ref|
    rel_error        max_t |a_PINN - a_ref| / a_ref
    conservation     max_t |w'^2 - 24 pi (rho0 w - sigma)| / max_t w'^2      (vinculo de Friedmann em w = a^3)
    constraint       max_t |residuo da EDO de 2a ordem| / escala               (a EDO que a rede minimiza)
    stability        media e desvio-padrao do erro final entre sementes (treino repetido)
    time_s           tempo de treino
Referencia analitica: exact_torsion_dust (RESULTADO MATEMATICO).  Referencia numerica: FriedmannModel.solve.

Classificacao: RESULTADO NUMERICO (ML como ferramenta, validada contra matematica; nunca o contrario).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from ..cosmology.friedmann import einstein_cartan_dust, exact_torsion_dust


@dataclass
class PINNBenchmark:
    problem: str
    seeds: list
    abs_error: list
    rel_error: list
    conservation: list
    constraint: list
    final_loss: list
    time_s: list
    numeric_vs_exact_abs: float
    data_weight: float
    epochs: int
    parameters: dict = field(default_factory=dict)

    def summary(self) -> dict:
        f = lambda x: dict(mean=float(np.mean(x)), std=float(np.std(x)), min=float(np.min(x)), max=float(np.max(x)))  # noqa: E731
        return dict(problem=self.problem, n_seeds=len(self.seeds), abs_error=f(self.abs_error), rel_error=f(self.rel_error),
                    conservation=f(self.conservation), constraint=f(self.constraint), final_loss=f(self.final_loss),
                    time_s=f(self.time_s), numeric_vs_exact_abs=self.numeric_vs_exact_abs, data_weight=self.data_weight,
                    epochs=self.epochs, parameters=self.parameters)


def benchmark_bounce_pinn(rho_m0=1.0, sigma=0.05, t_max=1.2, seeds=(0, 1, 2), epochs=2000, lbfgs_steps=200,
                          data_weight=0.0, n_data=16, verbose=False) -> PINNBenchmark:
    try:
        import torch
        from .pinn import BouncePINN
    except ImportError as e:  # pragma: no cover
        raise ImportError("torch necessario para os benchmarks de PINN") from e
    t = np.linspace(0, t_max, 400)
    exact, _ = exact_torsion_dust(t, rho_m0, sigma)
    num = einstein_cartan_dust(rho_m0, sigma).solve(a0=1.0, t_span=(0, t_max), n_samples=400).a
    nve = float(np.max(np.abs(num - exact)))
    out = dict(abs=[], rel=[], cons=[], cstr=[], loss=[], time=[])
    for seed in seeds:
        t0 = time.time()
        p = BouncePINN(rho_m0=rho_m0, sigma=sigma, t_max=t_max, seed=seed)
        data = None
        if data_weight > 0:
            td = np.linspace(0, t_max, n_data)
            data = (torch.tensor(td, dtype=torch.float32).reshape(-1, 1),
                    torch.tensor(np.interp(td, t, num) ** 3, dtype=torch.float32).reshape(-1, 1))
        hist = p.train(epochs=epochs, verbose=verbose, lbfgs_steps=lbfgs_steps, data=data, data_weight=data_weight)
        dt = time.time() - t0
        a_p = p.predict(t)
        err = np.abs(a_p - exact)
        # vinculo e residuo da EDO avaliados com autograd nos pontos t
        tt = torch.tensor(t, dtype=torch.float32).reshape(-1, 1)
        res, fr, w = p.residual(tt)
        w_np = w.detach().numpy().ravel()
        cons = float(np.max(np.abs(fr.detach().numpy().ravel())) / max(np.max(24 * np.pi * (rho_m0 * w_np)), 1e-300))
        cstr = float(np.max(np.abs(res.detach().numpy().ravel())) / max(np.max(np.abs(4 * np.pi * rho_m0 * w_np)), 1e-300))
        out["abs"].append(float(err.max()))
        out["rel"].append(float((err / exact).max()))
        out["cons"].append(cons)
        out["cstr"].append(cstr)
        out["loss"].append(float(hist[-1]))
        out["time"].append(float(dt))
    return PINNBenchmark("ricochete de torcao (w = a^3)", list(seeds), out["abs"], out["rel"], out["cons"], out["cstr"],
                         out["loss"], out["time"], nve, data_weight, epochs,
                         dict(rho_m0=rho_m0, sigma=sigma, t_max=t_max, lbfgs_steps=lbfgs_steps))
