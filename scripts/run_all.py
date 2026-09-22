"""Gera todas as figuras, tabelas e renders do projeto em ./output.

    python scripts/run_all.py              # tudo (renders 1600x900, ~3 min em CPU)
    python scripts/run_all.py --sem-render # so as figuras cientificas (~1 min)
    python scripts/run_all.py --rapido     # renders em 640x360
"""
import sys

sys.path.insert(0, ".")
from semente.figures import make_all  # noqa: E402

if __name__ == "__main__":
    size = (640, 360) if "--rapido" in sys.argv else (1600, 900)
    make_all(with_renders="--sem-render" not in sys.argv, render_size=size)
