"""Python wrappers around Cell-Tree 2D spatial index."""

# Just bring in the CellTree into this namespace.
from .cell_tree2d import CellTree

try:
    from ._version import __version__
except ImportError:  # pragma: nocover
    # package is not installed
    __version__ = "0.0.0.dev0"

__all__ = ["__version__", "CellTree"]
