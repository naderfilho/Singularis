"""Testes da Fase 9: explorador de parametros (reprodutibilidade, classificacao) e benchmark de PINN."""
import numpy as np
import pytest

from semente.parameter_space.explorer import (CLASSES, adaptive_refine, classify_bounce, classify_static, monte_carlo,
                                              scan_1d, scan_2d)


def test_static_classifier_known_cases():
    assert classify_static(1.0, 0.0, 0.0) == "classical_black_hole"
    assert classify_static(1.0, 0.0, 1.2) == "naked_singularity"
    assert classify_static(1.0, 0.5, 0.0) == "black_bounce"
    assert classify_static(1.0, 0.05, 0.5) == "regular_black_hole"
    assert classify_static(1.0, 2.5, 0.3) == "wormhole"


def test_bounce_classifier_known_cases():
    assert classify_bounce(np.inf) == "singular"
    assert classify_bounce(20.0, 0.0, 1e-6) == "expanding_baby_universe"
    assert classify_bounce(20.0, 0.0, 10.0) == "unstable_bounce"


def test_scans_shapes_and_labels_are_valid():
    s1 = scan_1d(classify_static, "l", np.linspace(0, 2.5, 6), fixed=dict(M=1.0, Q=0.0))
    assert len(s1.labels) == 6 and all(l in CLASSES for l in s1.labels)
    s2 = scan_2d(classify_static, ("Q", np.linspace(0, 1.3, 5)), ("l", np.linspace(0, 2.6, 4)), fixed=dict(M=1.0))
    assert s2.grid_labels().shape == (4, 5)


def test_monte_carlo_is_reproducible_with_seed():
    kw = dict(bounds=dict(rho_star=(2, 1e3), sigma0_sq=(1e-6, 10)), n=8, log_scale=("rho_star", "sigma0_sq"))
    a = monte_carlo(classify_bounce, seed=5, **kw)
    b = monte_carlo(classify_bounce, seed=5, **kw)
    c = monte_carlo(classify_bounce, seed=6, **kw)
    assert np.allclose(a.points, b.points) and a.labels == b.labels
    assert not np.allclose(a.points, c.points)
    assert a.record.seed == 5


def test_adaptive_refinement_adds_points_near_boundaries(tmp_path):
    ad = adaptive_refine(classify_static, ("Q", np.linspace(0, 1.3, 5)), ("l", np.linspace(0, 2.6, 5)), fixed=dict(M=1.0), levels=1)
    assert ad.extra["n_evaluations"] > 25
    meta = ad.save(str(tmp_path / "scan.json"))
    assert meta.endswith(".meta.json")


@pytest.mark.skipif(pytest.importorskip("torch", reason="torch nao instalado") is None, reason="torch")
def test_pinn_benchmark_metrics():
    from semente.ml.benchmarks import benchmark_bounce_pinn
    b = benchmark_bounce_pinn(seeds=(0,), epochs=600, lbfgs_steps=100)
    s = b.summary()
    assert s["numeric_vs_exact_abs"] < 1e-8
    assert np.isfinite(s["abs_error"]["mean"]) and s["abs_error"]["mean"] < 5e-2
    assert s["conservation"]["mean"] < 1e-1 and s["constraint"]["mean"] < 1e-1
