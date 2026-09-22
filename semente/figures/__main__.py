import sys

from .core import make_all

make_all(with_renders="--sem-render" not in sys.argv)
