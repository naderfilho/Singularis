"""Fachada de compatibilidade (docs/03_nascimento.md).

Conteudo real: cosmology/perturbations.py (espectro), cosmology/curvature.py (curvatura e massa
do pai), observations/constraints.py (dados), observations/verdict_legacy.py (veredito).
"""
from __future__ import annotations

from .core.units import GYR, MPC  # noqa: F401
from .cosmology.curvature import curvature_constraints, matter_bounce_tensor_ratio  # noqa: F401
from .cosmology.perturbations import BounceSpectrum  # noqa: F401
from .observations.constraints import PLANCK  # noqa: F401
from .observations.verdict_legacy import verdict  # noqa: F401
