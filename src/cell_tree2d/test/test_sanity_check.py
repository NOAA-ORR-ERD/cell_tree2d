"""
Tests of the sanity checking code.
"""

import numpy as np
import warnings

from cell_tree2d import CellTree, sanity_check

import pytest

def test_sanity_check_good():
    nodes = np.array([[0.0, 0.0],
                      [2.0, 0.0],
                      [2.0, 2.0],
                      [0.0, 2.0],
                      [1.0, 1.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sanity_check(nodes)


def test_sanity_check_duplicate_node():
    nodes = np.array([[0.0, 0.0],
                      [2.0, 0.0],
                      [2.0, 2.0],
                      [0.0, 2.0],
                      [1.0, 1.0],
                      [0.0, 0.0]])
    with pytest.warns(UserWarning, match="duplicate nodes"):
        sanity_check(nodes)

    with pytest.raises(ValueError):
        sanity_check(nodes, error='error')

@pytest.mark.filterwarnings("ignore:There are 11 duplicate nodes")
def test_zero_size_bb():
    """
    If a grid has zero-size bounding boxes, it shouldn't build,
    and should raise an Error instead.

    This is only testing the really pathological case

    almost all zeros, but it's something.
    """
    # A simple quad grid
    nodes = np.array([[0.0, 0.0], #0
                      [0.0, 0.0], #1
                      [0.0, 0.0], #2
                      [0.0, 0.0], #3
                      [0.0, 0.0], #4
                      [0.0, 0.0], #5
                      [0.0, 0.0], #6
                      [0.0, 0.0], #7
                      [0.0, 0.0], #8
                      [0.0, 0.0], #9
                      [4.0, 4.0], #10
                      [6.0, 4.0]  #11
                      ])
#    nodes = np.zeros_like(nodes)
    nodes = np.ones_like(nodes)

    faces = np.array([0, 8, 9, 5, 2,
                      9, 11, 7, 5,
                      4, 7, 6]
                      , dtype=np.int32)

    n_verts_arr = np.array([5, 4, 3], dtype=np.ubyte)

    with pytest.raises(ValueError):
        tree = CellTree(nodes,
                        faces,
                        len_arr=n_verts_arr,
                        num_buckets=2,
                        cells_per_leaf=1,
                        check_inputs='warn')


def test_bad_grid_raise():
    """
    If a grid is really bad, it shouldn't crash

    See:
    https://github.com/NOAA-ORR-ERD/gridded/issues/79

    this is only testing the really pathological case

    almost all zeros, but it's something.
    """
    # A simple quad grid
    nodes = np.array([[0.0, 0.0], #0
                      [0.0, 0.0], #1
                      [0.0, 0.0], #2
                      [0.0, 0.0], #3
                      [0.0, 0.0], #4
                      [0.0, 0.0], #5
                      [0.0, 0.0], #6
                      [0.0, 0.0], #7
                      [0.0, 0.0], #8
                      [0.0, 0.0], #9
                      [4.0, 4.0], #10
                      [6.0, 4.0] #11
                      ])
#    nodes = np.zeros_like(nodes)
    nodes = np.ones_like(nodes)

    faces = np.array([0, 8, 9, 5, 2,
                      9, 11, 7, 5,
                      4, 7, 6]
                      , dtype=np.int32)

    n_verts_arr = np.array([5, 4, 3], dtype=np.ubyte)

    with pytest.raises(ValueError):
        tree = CellTree(nodes,
                        faces,
                        len_arr=n_verts_arr,
                        num_buckets=2,
                        cells_per_leaf=1,
                        check_inputs='raise')




