"""Reexecuta um estudo a partir de um arquivo de configuracao e grava resultado + proveniencia.

    python scripts/reproduce.py configs/collapse_bounce.json
    python scripts/reproduce.py configs/*.json
"""
import glob
import sys

sys.path.insert(0, ".")
from semente.core.config import load_config, run_study  # noqa: E402

if __name__ == "__main__":
    paths = [p for arg in sys.argv[1:] for p in glob.glob(arg)] or glob.glob("configs/*.json")
    for p in paths:
        cfg = load_config(p)
        out = run_study(cfg)
        print(f"{p} -> {out} (+ .meta.json)")
