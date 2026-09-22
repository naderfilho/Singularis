"""Testes do nucleo: unidades, solver unico, proveniencia."""
import json
import os

import numpy as np

from semente.core.provenance import RunRecord, save_json
from semente.core.solvers import ODESettings, integrate
from semente.core.units import (C_SI, G_SI, L_PLANCK, M_PLANCK, M_SUN, RHO_PLANCK, T_PLANCK, Scale)


def test_planck_units_consistent():
    # l_Pl = G m_Pl / c^2 e t_Pl = l_Pl / c
    assert abs(G_SI * M_PLANCK / C_SI**2 / L_PLANCK - 1) < 1e-12
    assert abs(T_PLANCK * C_SI / L_PLANCK - 1) < 1e-12
    assert abs(RHO_PLANCK * L_PLANCK**3 / M_PLANCK - 1) < 1e-12


def test_scale_conversions_for_solar_mass():
    s = Scale.from_solar_masses(1.0)
    assert abs(s.schwarzschild_radius_m() - 2953.25) < 1.0          # 2 G M_sun / c^2 = 2953 m
    assert abs(s.time_unit - 4.925e-6) / 4.925e-6 < 1e-3            # G M_sun / c^3
    assert abs(s.M_planck - M_SUN / M_PLANCK) < 1e-6 * s.M_planck
    # frequencia de uma oscilacao omega = 0.5/M para 10 M_sol ~ 1.6 kHz
    f = Scale.from_solar_masses(10).frequency_to_hz(0.5)
    assert 1500 < f < 1700


def test_integrate_harmonic_oscillator_and_event():
    def rhs(t, y):
        return [y[1], -y[0]]

    def zero_crossing(t, y):
        return y[0]
    zero_crossing.terminal = True
    zero_crossing.direction = -1
    sol = integrate(rhs, (0, 10), [1.0, 0.0], n_samples=200, settings=ODESettings(rtol=1e-10, atol=1e-12), events=[zero_crossing])
    assert abs(sol.t[-1] - np.pi / 2) < 1e-6
    assert np.max(np.abs(sol.y[0] - np.cos(sol.t))) < 1e-7
    assert sol.as_metadata()["solver"]["method"] == "DOP853"


def test_provenance_sidecar(tmp_path):
    rec = RunRecord(model="teste", parameters=dict(M=1.0, l=0.5), solver=ODESettings().as_dict(), seed=7)
    path = os.path.join(tmp_path, "r.json")
    meta = save_json(path, dict(x=np.arange(3)), rec)
    d = json.load(open(meta, encoding="utf-8"))
    assert d["model"] == "teste" and d["parameters"]["l"] == 0.5 and d["seed"] == 7
    assert "numpy" in d["dependencies"] and d["timestamp_utc"]
    assert json.load(open(path, encoding="utf-8"))["x"] == [0, 1, 2]
