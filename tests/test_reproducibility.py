"""Testes da Fase 10: configuracoes reprodutiveis, proveniencia e registro de validacao."""
import json
import os

from semente.core.config import STUDIES, load_config, run_study
from semente.validation import run as run_validation


def test_config_study_roundtrip(tmp_path):
    cfg = dict(study="qnm", version="1", seed=3, parameters=dict(M=1.0, spin=0, ell=2, l_values=[0.0, 0.5]),
               outputs=str(tmp_path))
    p = tmp_path / "qnm.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")
    loaded = load_config(str(p))
    out = run_study(loaded)
    res = json.load(open(out, encoding="utf-8"))
    meta = json.load(open(out[:-5] + ".meta.json", encoding="utf-8"))
    assert abs(res["0.0"]["omega_re"] - 0.4834) < 2e-3
    assert meta["seed"] == 3 and meta["model"] == "study/qnm" and meta["parameters"]["l_values"] == [0.0, 0.5]
    assert "git_commit" in meta and "timestamp_utc" in meta and meta["solver"]["method"] == "DOP853"


def test_unknown_study_rejected(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(dict(study="nope", parameters={})), encoding="utf-8")
    try:
        load_config(str(p))
    except ValueError as e:
        assert "desconhecido" in str(e)
    else:
        raise AssertionError("config invalida aceita")


def test_all_registered_studies_have_configs():
    names = {json.load(open(os.path.join("configs", f), encoding="utf-8"))["study"] for f in os.listdir("configs") if f.endswith(".json")}
    assert names == set(STUDIES)


def test_validation_registry_passes(tmp_path):
    rep = run_validation(str(tmp_path / "validation.json"))
    failed = [c["name"] for c in rep["checks"] if not c["passed"]]
    assert rep["all_passed"], f"verificacoes falhando: {failed}"
    cats = {c["category"] for c in rep["checks"]}
    assert {"dimensional", "limiting_case", "conservation", "convergence", "stability", "analytical_benchmark"} <= cats
