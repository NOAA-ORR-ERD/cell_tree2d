#!/usr/bin/env python

"""
tests for using celltree with a quad grid

also serves as an example for how to use it with curvilear grids
"""

import numpy as np

from cell_tree2d import CellTree

def curv_grid(n_x=4, n_y=12, center=(20.0, 30.0), min_radius=1.0, max_radius=1.1,angle = np.pi / 2.0):
    """
    example of quad grid of a partial circle to use for tests
    """

    lon = np.zeros((n_y, n_x), dtype=np.float64)
    lat = np.zeros((n_y, n_x), dtype=np.float64)

    for i, theta in enumerate(np.linspace(0, angle, n_x)):
        for j, r in enumerate(np.linspace(min_radius, max_radius, n_y)):
            lon[j, i] = center[1] + r * np.cos(theta)
            lat[j, i] = center[0] + r * np.sin(theta)
    return lon, lat
 
def nodes_from_coords(x, y):
    """
    generates nodes and indeces for the cells from arrays of x and y nodes_from_coords
    """
    nodes = np.ascontiguousarray(np.column_stack((x[:].reshape(-1),
                                                          y[:].reshape(-1)))).astype(np.float64)
    y_size, x_size = x.shape
    faces = np.array([np.array([[x, x + 1, x + x_size + 1, x + x_size]
                                        for x in range(0, x_size - 1, 1)]) + y * x_size for y in range(0, y_size - 1)])
    faces = np.ascontiguousarray(faces.reshape(-1, 4).astype(np.int32))

    return nodes, faces

def test_build_tree_from_coords():

    x, y = curv_grid()
    nodes, faces = nodes_from_coords(x, y)
    #these range from (30.,  20.) to (31.1,  21.1)

    print(nodes.min(axis=0), nodes.max(axis=0))

    tree = CellTree(nodes, faces)

    # try to find some points

    # points outside the domain:
    result = tree.multi_locate([(10.0, 12.0),
                                (40.0, 50.0),
                                ],
                               )
    assert np.array_equal(result, (-1, -1))


    






if __name__ == "__main__":

    test_build_tree_from_coords()


