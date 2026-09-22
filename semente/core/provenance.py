"""
Proveniencia: cada resultado/figura salvo responde "como exatamente isto foi produzido?".

`RunRecord` captura: nome do modelo e versao, parametros (snapshot), configuracao do
solver, semente aleatoria, versao do pacote e das dependencias, commit git, timestamp.
`save_json` grava o resultado junto com o registro num arquivo sidecar `<nome>.meta.json`.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import subprocess
from dataclasses import asdict, dataclass, field, is_dataclass
from typing import Any, Optional

import numpy as np

from .. import __version__ as PACKAGE_VERSION


def _git_commit() -> Optional[str]:
    try:
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:
        return None


def _dependency_versions() -> dict:
    vers = {}
    for name in ("numpy", "scipy", "sympy", "matplotlib", "torch"):
        try:
            mod = __import__(name)
            vers[name] = getattr(mod, "__version__", "?")
        except Exception:
            vers[name] = None
    return vers


def _jsonable(obj: Any):
    if is_dataclass(obj):
        return {k: _jsonable(v) for k, v in asdict(obj).items()}
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist() if obj.size <= 20000 else {"__array__": list(obj.shape), "summary": [float(np.nanmin(obj)), float(np.nanmax(obj))]}
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (bool, int, float, str)) or obj is None:
        return obj
    return str(obj)


@dataclass
class RunRecord:
    model: str
    model_version: str = "1"
    parameters: dict = field(default_factory=dict)
    solver: dict = field(default_factory=dict)
    seed: Optional[int] = None
    assumptions: list = field(default_factory=list)
    notes: str = ""
    package_version: str = PACKAGE_VERSION
    git_commit: Optional[str] = field(default_factory=_git_commit)
    dependencies: dict = field(default_factory=_dependency_versions)
    python: str = field(default_factory=lambda: platform.python_version())
    platform: str = field(default_factory=lambda: platform.platform())
    timestamp_utc: str = field(default_factory=lambda: _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"))

    def as_dict(self) -> dict:
        return _jsonable(self)


def save_json(path: str, result: Any, record: RunRecord) -> str:
    """Grava `result` (dict/dataclass/arrays) em `path` e o registro de proveniencia em `<path>.meta.json`."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(_jsonable(result), fh, indent=2, ensure_ascii=False)
    meta = path[:-5] + ".meta.json" if path.endswith(".json") else path + ".meta.json"
    with open(meta, "w", encoding="utf-8") as fh:
        json.dump(record.as_dict(), fh, indent=2, ensure_ascii=False)
    return meta


def stamp_figure(fig, record: RunRecord, extra: str = ""):
    """Rodape discreto com modelo, parametros principais e commit, para a figura ser auditavel."""
    params = ", ".join(f"{k}={v:.3g}" if isinstance(v, (int, float, np.floating)) else f"{k}={v}"
                       for k, v in list(record.parameters.items())[:6])
    commit = (record.git_commit or "no-git")[:7]
    txt = f"{record.model} v{record.model_version} | {params} | semente v{record.package_version} @ {commit} | {record.timestamp_utc[:10]}"
    if extra:
        txt += " | " + extra
    fig.text(0.005, 0.005, txt, fontsize=6.5, color="#8a93a6", ha="left", va="bottom")
