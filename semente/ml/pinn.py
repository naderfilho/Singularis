"""
Rede neural informada por fisica (PINN) para o ricochete cosmologico.

A rede aprende a(t) diretamente das EQUACOES (nao de dados): a perda e o
residuo da equacao de Friedmann com torcao (Einstein-Cartan) em pontos de
colocacao, mais as condicoes iniciais.  Depois comparamos com o integrador
numerico (semente.bounce.TorsionCosmology) - a rede nunca viu a solucao.

    a'' = -(4 pi/3) [ rho + 3p - 4 sigma a^-6 ] a      (escrita na variavel de volume w = a^3)

Parametrizacao com "hard constraints": w = w0 + w0' t + t^2 N(t), o que impoe as
condicoes iniciais exatamente.

Tambem inclui uma PINN para a orbita de fotons  u'' + u = 3 M u^2  (a mesma
equacao que gera as imagens do ray-tracer), como segunda demonstracao.

Requer torch (pip install torch).  Roda em CPU em ~1 min.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    import torch
    import torch.nn as nn

    TORCH_OK = True
except Exception:  # pragma: no cover
    TORCH_OK = False


def _require_torch():
    if not TORCH_OK:
        raise ImportError("Instale torch para usar a PINN: pip install torch")


class MLP(nn.Module if TORCH_OK else object):
    def __init__(self, width=64, depth=4, inputs=1, outputs=1):
        super().__init__()
        layers = [nn.Linear(inputs, width), nn.Tanh()]
        for _ in range(depth - 1):
            layers += [nn.Linear(width, width), nn.Tanh()]
        layers += [nn.Linear(width, outputs)]
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)
        nn.init.zeros_(self.net[-1].weight)  # saida inicial = 0 -> chute inicial = solucao "livre"

    def forward(self, x):
        return self.net(x)


@dataclass
class BouncePINN:
    """PINN na variavel de VOLUME w = a^3 (a variavel natural da cosmologia quantica de laco).

    Para poeira + torcao (k = 0), a equacao de Friedmann de 2a ordem vira

        w w'' - (2/3) w'^2 + 4 pi rho_m0 w - 16 pi sigma = 0,        (residuo polinomial, sem rigidez)
        w'^2 - 24 pi (rho_m0 w - sigma) = 0                          (vinculo de Friedmann)

    Condicoes iniciais impostas exatamente: w(t) = w0 + w0' t + t^2 N(t/T).
    A rede so ve as equacoes; a solucao numerica (e a exata) servem apenas para conferir.
    """

    rho_m0: float = 1.0
    sigma: float = 0.05
    a0: float = 1.0
    t_max: float = 1.2
    width: int = 64
    depth: int = 4
    seed: int = 0

    def __post_init__(self):
        _require_torch()
        torch.manual_seed(self.seed)
        self.net = MLP(self.width, self.depth)
        self.w0 = self.a0**3
        H2 = (8 * np.pi / 3) * (self.rho_m0 * self.a0**-3 - self.sigma * self.a0**-6)
        if H2 <= 0:
            raise ValueError("a0 abaixo do ricochete")
        self.dw0 = -3 * self.w0 * np.sqrt(H2)  # comeca contraindo
        self.T = self.t_max

    def w(self, t):
        return self.w0 + self.dw0 * t + t**2 * self.net(t / self.T)

    def a(self, t):
        return torch.clamp(self.w(t), min=1e-9) ** (1.0 / 3.0)

    def residual(self, t):
        t = t.requires_grad_(True)
        w = self.w(t)
        dw = torch.autograd.grad(w, t, torch.ones_like(w), create_graph=True)[0]
        ddw = torch.autograd.grad(dw, t, torch.ones_like(dw), create_graph=True)[0]
        res = w * ddw - (2.0 / 3.0) * dw**2 + 4 * np.pi * self.rho_m0 * w - 16 * np.pi * self.sigma
        fried = dw**2 - 24 * np.pi * (self.rho_m0 * w - self.sigma)
        return res, fried, w

    def _loss(self, t, w_fried, data=None, data_weight=0.0):
        """L = L_physics + w_fried * L_vinculo (+ data_weight * L_data, modo hibrido).  L_boundary = 0 por
        construcao (condicoes iniciais impostas exatamente na parametrizacao)."""
        res, fr, _ = self.residual(t)
        loss = torch.mean(res**2) + w_fried * torch.mean(fr**2)
        if data is not None and data_weight > 0:
            td, wd = data
            loss = loss + data_weight * torch.mean((self.w(td) - wd) ** 2)
        return loss

    def train(self, epochs=3000, n_col=512, lr=3e-3, w_fried=0.1, verbose=True, lbfgs_steps=300, curriculum=True,
              data=None, data_weight=0.0):
        """Adam com curriculo temporal (o intervalo [0, t_k] cresce ate t_max) + refinamento L-BFGS.
        `data = (t_pontos, w_pontos)` e `data_weight > 0` ligam o termo de dados (hibrido)."""
        opt = torch.optim.Adam(self.net.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
        hist = []
        for ep in range(epochs):
            frac = min(1.0, 0.3 + ep / (0.6 * epochs)) if curriculum else 1.0
            t = torch.rand(n_col, 1) * self.t_max * frac
            loss = self._loss(t, w_fried, data, data_weight)
            opt.zero_grad()
            loss.backward()
            opt.step()
            sched.step()
            hist.append(loss.item())
            if verbose and ep % 500 == 0:
                print(f"  epoca {ep}: perda {loss.item():.3e}")
        if lbfgs_steps:
            lb = torch.optim.LBFGS(self.net.parameters(), max_iter=lbfgs_steps, line_search_fn="strong_wolfe",
                                   tolerance_grad=1e-12, tolerance_change=1e-14)
            t = torch.linspace(0, self.t_max, 2048).reshape(-1, 1)

            def closure():
                lb.zero_grad()
                loss = self._loss(t.clone(), w_fried, data, data_weight)
                loss.backward()
                hist.append(loss.item())
                return loss

            lb.step(closure)
        return np.array(hist)

    def predict(self, t: np.ndarray) -> np.ndarray:
        with torch.no_grad():
            tt = torch.tensor(t, dtype=torch.float32).reshape(-1, 1)
            return self.a(tt).numpy().ravel()


@dataclass
class PhotonOrbitPINN:
    """u'' + u = 3 M u^2 com u(0) = 0, u'(0) = 1/b: a orbita de um foton que vem do infinito."""

    b: float = 6.0
    M: float = 1.0
    phi_max: float = 3.0
    width: int = 48
    depth: int = 4
    seed: int = 0

    def __post_init__(self):
        _require_torch()
        torch.manual_seed(self.seed)
        self.net = MLP(self.width, self.depth)

    def u(self, phi):
        return phi / self.b + phi**2 * self.net(phi / self.phi_max)

    def residual(self, phi):
        phi = phi.requires_grad_(True)
        u = self.u(phi)
        du = torch.autograd.grad(u, phi, torch.ones_like(u), create_graph=True)[0]
        ddu = torch.autograd.grad(du, phi, torch.ones_like(du), create_graph=True)[0]
        return ddu + u - 3 * self.M * u**2

    def train(self, epochs=2000, n_col=400, lr=3e-3, verbose=False, lbfgs_steps=300):
        opt = torch.optim.Adam(self.net.parameters(), lr=lr)
        hist = []
        for ep in range(epochs):
            phi = torch.rand(n_col, 1) * self.phi_max
            loss = torch.mean(self.residual(phi) ** 2)
            opt.zero_grad()
            loss.backward()
            opt.step()
            hist.append(loss.item())
        if lbfgs_steps:
            lb = torch.optim.LBFGS(self.net.parameters(), max_iter=lbfgs_steps, line_search_fn="strong_wolfe")
            phi = torch.linspace(0, self.phi_max, 1024).reshape(-1, 1)

            def closure():
                lb.zero_grad()
                loss = torch.mean(self.residual(phi.clone()) ** 2)
                loss.backward()
                hist.append(loss.item())
                return loss

            lb.step(closure)
        return np.array(hist)

    def predict(self, phi):
        with torch.no_grad():
            return self.u(torch.tensor(phi, dtype=torch.float32).reshape(-1, 1)).numpy().ravel()


def reference_photon_orbit(b, M=1.0, phi_max=3.0, n=600):
    """Solucao numerica de referencia (RK45) de u'' + u = 3Mu^2, u(0)=0, u'(0)=1/b."""
    from scipy.integrate import solve_ivp

    sol = solve_ivp(lambda p, y: [y[1], -y[0] + 3 * M * y[0] ** 2], (0, phi_max), [0.0, 1 / b],
                    rtol=1e-10, atol=1e-12, dense_output=True)
    phi = np.linspace(0, phi_max, n)
    return phi, sol.sol(phi)[0]
