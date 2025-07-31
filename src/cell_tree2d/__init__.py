#!/usr/bin/env python

from __future__ import absolute_import

# Just bring in the CellTree into this namespace.
# using as to make the linters happy
from .cell_tree2d import CellTree as CellTree
from .cell_tree2d import sanity_check as sanity_check


__version__ = "1.0.0dev"
