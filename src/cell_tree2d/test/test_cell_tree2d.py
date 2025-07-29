#!/usr/bin/env python

"""
unit tests for cell_tree_2d

uses pytest

to run:

$ py.test

in this dir.
"""

import pytest
import numpy as np

from cell_tree2d import CellTree

# some very simple test data:
# two triangles:
nodes2 = [[0.0, 0.0],
          [2.0, 0.0],
          [1.0, 2.0],
          [3.0, 2.0],
          ]

faces2 = [[0, 1, 2],
          [1, 3, 2],
          ]
nodes = np.array(nodes2, dtype=np.float64)
faces = np.array(faces2, dtype=np.intc)

# A basic triangle grid with a hole and a "tail"
# see 21_tri_mesh.png to see what it looks like

nodes21 = [(5, 1),
           (10, 1),
           (3, 3),
           (7, 3),
           (9, 4),
           (12, 4),
           (5, 5),
           (3, 7),
           (5, 7),
           (7, 7),
           (9, 7),
           (11, 7),
           (5, 9),
           (8, 9),
           (11, 9),
           (9, 11),
           (11, 11),
           (7, 13),
           (9, 13),
           (7, 15),
           ]

faces21 = [(0, 1, 3),
           (0, 2, 6),
           (0, 3, 6),
           (1, 3, 4),
           (1, 4, 5),
           (2, 6, 7),
           (6, 7, 8),
           (7, 8, 12),
           (6, 8, 9),
           (8, 9, 12),
           (9, 12, 13),
           (4, 5, 11),
           (4, 10, 11),
           (9, 10, 13),
           (10, 11, 14),
           (10, 13, 14),
           (13, 14, 15),
           (14, 15, 16),
           (15, 16, 18),
           (15, 17, 18),
           (17, 18, 19),
           ]


def test_init():
    """
    can  a tree be initialized
    """

    # with defaults
    tree = CellTree(nodes, faces)

    # with everything specified
    tree = CellTree(nodes, faces, num_buckets=2, cells_per_leaf=1)

    # with num_buckets
    tree = CellTree(nodes, faces, num_buckets=4)

    # with cells_per_leaf
    tree = CellTree(nodes, faces, cells_per_leaf=2)

    assert True


def test_lists():
    """
    python lists should get converted to numpy arrays
    """
    tree = CellTree(nodes2, faces2)
    assert True


def test_types():
    """
    It should auto-cast the types to the right types for you
    """
    nodes = np.array(nodes2, dtype=np.float32)
    faces = np.array(faces2, dtype=np.int64)

    tree = CellTree(nodes, faces)
    assert True


def test_shape_error():
    nodes = [(1, 2, 3),
             (3, 4, 5),
             (4, 5, 6)]
    faces = [[0, 1, 2],
             [1, 3, 2],
             ]

    with pytest.raises(ValueError):
        # nodes is wrong shape
        tree = CellTree(nodes, faces)

        tree = CellTree(nodes2, (2, 3, 4, 5))

        tree = CellTree(
            nodes2, ((2, 3, 4, 5), (1, 2, 3, 4, 5), (1, 2, 3, 4, 5)))


def test_bounds_errors():
    with pytest.raises(ValueError):
        tree = CellTree(nodes, faces, cells_per_leaf=-1)

    with pytest.raises(ValueError):
        tree = CellTree(nodes, faces, num_buckets=0)


def test_triangle_lookup():
    tree = CellTree(nodes, faces)
    point = np.array([1., 1.])  # in triangle 1
    result = tree.locate(point)
    assert result == 0
    point[0] = 2.0  # tri 2
    result = tree.locate(point)
    assert result == 1
    point[0] = -1.0  # out of grid
    result = tree.locate(point)
    assert result == -1

def test_poly_lookup():
    # A simple quad grid
    nodes = np.array([[0.0, 0.0], #0
                      [0.0, 2.0], #1
                      [2.0, 0.0], #2
                      [2.0, 2.0], #3
                      [4.0, 0.0], #4
                      [4.0, 2.0], #5
                      [6.0, 0.0], #6
                      [6.0, 2.0], #7
                      [0.0,4.0], #8
                      [2.0,4.0], #9
                      [4.0,4.0], #10
                      [6.0,4.0] #11
                      ])

    #quads
    faces1 = np.array([[0, 2, 3, 1],
                      [4, 6, 7, 5],
                      ], dtype=np.intc)
    #Pentas
    faces2 = np.array([[0, 8, 9, 5, 2],
                      [9, 11, 6, 2, 5],
                      ], dtype=np.intc)
    tree1 = CellTree(nodes, faces1, num_buckets = 2, cells_per_leaf = 1)
    point = np.array([1., 1.])  # in triangle 1
    result = tree1.locate(point)
    assert result == 0
    point[0] = 5.0  # tri 2
    result = tree1.locate(point)
    assert result == 1
    point[0] = -1.0  # out of grid
    result = tree1.locate(point)
    assert result == -1

    tree2 = CellTree(nodes, faces2, num_buckets=2, cells_per_leaf=1)
    point = np.array([1., 2.])  # in triangle 1
    result = tree2.locate(point)
    assert result == 0
    point[0] = 5.0  # tri 2
    result = tree2.locate(point)
    assert result == 1
    point[0] = -1.0  # out of grid
    result = tree2.locate(point)
    assert result == -1

def test_multi_poly_lookup():
    # A simple quad grid
    nodes = np.array([[0.0, 0.0], #0
                      [0.0, 2.0], #1
                      [2.0, 0.0], #2
                      [2.0, 2.0], #3
                      [4.0, 0.0], #4
                      [4.0, 2.0], #5
                      [6.0, 0.0], #6
                      [6.0, 2.0], #7
                      [0.0,4.0], #8
                      [2.0,4.0], #9
                      [4.0,4.0], #10
                      [6.0,4.0] #11
                      ])

    faces = np.array([[0, 8, 9, 5, 2],
                      [9, 11, 7, 5, -1],
                      [4, 7 ,6, -1, -1]
                      ], dtype=np.intc)
    tree = CellTree(nodes, faces, num_buckets = 2, cells_per_leaf = 1)
    point = np.array([1., 1.])  # in POLY 1
    result = tree.locate(point)
    assert result == 0
    point[0] = 5.0  # tri 2
    result = tree.locate(point)
    assert result == 2

    point[1] = 3.0
    result = tree.locate(point)
    assert result == 1

    point[0] = -1.0  # out of grid
    result = tree.locate(point)
    assert result == -1

def test_multi_poly_1D_lookup():
    # A simple quad grid
    nodes = np.array([[0.0, 0.0], #0
                      [0.0, 2.0], #1
                      [2.0, 0.0], #2
                      [2.0, 2.0], #3
                      [4.0, 0.0], #4
                      [4.0, 2.0], #5
                      [6.0, 0.0], #6
                      [6.0, 2.0], #7
                      [0.0,4.0], #8
                      [2.0,4.0], #9
                      [4.0,4.0], #10
                      [6.0,4.0] #11
                      ])

    faces = np.array([0, 8, 9, 5, 2, 9, 11, 7, 5, 4, 7, 6]
                      , dtype=np.intc)

    n_verts_arr = np.array([5,4,3],dtype=np.ubyte)
    tree = CellTree(nodes, faces, len_arr = n_verts_arr, num_buckets = 2, cells_per_leaf = 1)
    point = np.array([1., 1.])  # in POLY 1
    result = tree.locate(point)
    assert result == 0
    point[0] = 5.0  # tri 2
    result = tree.locate(point)
    assert result == 2

    point[1] = 3.0
    result = tree.locate(point)
    assert result == 1

    point[0] = -1.0  # out of grid
    result = tree.locate(point)
    assert result == -1


def test_edge_cases():
    nodes = np.array([[0.0, 0.0],
                      [2.0, 0.0],
                      [2.0, 2.0],
                      [0.0, 2.0],
                      [1.0, 1.0]])

    faces1 = np.array([[0, 1, 3],
                       [1, 2, 3]], dtype=np.intc)
    faces2 = np.array([[0, 1, 3],
                       [1, 2, 4],
                       [2, 4, 3]], dtype=np.intc)
    tree1 = CellTree(nodes, faces1)
    tree2 = CellTree(nodes, faces2)


def test_multipoint():
    tree = CellTree(nodes21, faces21)
    points = [(4.2,  3.0),
              (7.7, 13.5),
              (3.4, 7.000000001),
              (7.0,  5.0),  # out of bounds points
              (8.66, 10.99),
              (7.3, 0.74),
              (2.5, 5.5),
              (9.8, 12.3),
              ]
    correct_indexes = (1, 20, 7, -1, -1, -1, -1, -1)

    ind = tree.locate(points)
    assert np.array_equal(ind, correct_indexes)


if __name__ == "__main__":
    test_triangle_lookup()
    test_multi_poly_lookup()
    test_poly_lookup()
    test_multipoint()
