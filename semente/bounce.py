"""Fachada de compatibilidade: o ricochete vive em cosmology/friedmann.py; o universo-filho em
cosmology/seed_universe.py; as constantes em core/units.py."""
from __future__ import annotations

from .core.units import (C_SI, G_SI, HBAR_SI, L_PLANCK, M_NEUTRON, M_PLANCK, M_SUN, RHO_PLANCK, T_PLANCK)  # noqa: F401
from .cosmology.friedmann import (Fluid, FriedmannModel, GRCorrection, classical_dust_collapse,  # noqa: F401
                                  einstein_cartan_dust, exact_torsion_dust, lqc_fluid)
from .cosmology.seed_universe import SeedUniverse, observable_universe_as_black_hole  # noqa: F401
from .quantum.einstein_cartan import EinsteinCartanCorrection  # noqa: F401
from .quantum.lqc import LQCCorrection  # noqa: F401

torsion_dust_exact = exact_torsion_dust
classical_collapse_reference = classical_dust_collapse


class TorsionCosmology:
    """Adaptador legado sobre FriedmannModel + EinsteinCartanCorrection (mesma API de antes)."""

    def __init__(self, rho_m0=1.0, rho_r0=0.0, sigma=1e-3, k=0.0):
        self.rho_m0, self.rho_r0, self.sigma, self.k = rho_m0, rho_r0, sigma, k
        self.model = einstein_cartan_dust(rho_m0, sigma, rho_r0, k)

    def rho(self, a):
        return self.model.rho(a)

    def p(self, a):
        return self.model.p(a)

    def rho_eff(self, a):
        return self.model.rho_eff(a)

    def H2(self, a):
        return self.model.H2(a)

    def a_bounce(self):
        return self.model.bounce_scale()

    def solve(self, a0=1.0, t_span=(0.0, 6.0), contracting=True, n=4000, **kw):
        s = self.model.solve(a0=a0, t_span=t_span, contracting=contracting, n_samples=n)
        return dict(t=s.t, a=s.a, adot=s.adot, H=s.H, rho=s.rho, rho_eff=s.rho_eff, a_bounce=self.a_bounce(), solution=s)


class LQCCosmology:
    """Adaptador legado sobre FriedmannModel + LQCCorrection (k = 0)."""

    def __init__(self, rho0=1.0, w=0.0, rho_c=0.41):
        self.rho0, self.w, self.rho_c = rho0, w, rho_c
        self.model = lqc_fluid(rho0, w, rho_c)

    def rho(self, a):
        return self.model.rho(a)

    def H2(self, a):
        return self.model.H2(a)

    def a_bounce(self):
        return (self.rho0 / self.rho_c) ** (1.0 / (3 * (1 + self.w)))

    def solve(self, a0=1.0, t_span=(0.0, 6.0), n=4000):
        if self.H2(a0) < 0:
            raise ValueError("rho(a0) > rho_c: a0 ja esta alem do ricochete; use rho0 menor.")
        s = self.model.solve(a0=a0, t_span=t_span, contracting=True, n_samples=n)
        return dict(t=s.t, a=s.a, H=s.H, rho=s.rho, a_bounce=self.a_bounce(), solution=s)
