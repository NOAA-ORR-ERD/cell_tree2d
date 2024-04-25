"""Python wrappers around Cell-Tree 2D spatial index."""

import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

# cython extensions
CY_MODULES = [
    Extension(
        name=f"cell_tree2d.cell_tree2d",
        sources=[
            os.path.join("src", "cell_tree2d", "cell_tree2d.pyx"),
            os.path.join("src", "cell_tree2d_cpp", "cell_tree2d.cpp"),
        ],
        include_dirs=[np.get_include(), os.path.join("src", "cell_tree2d_cpp")],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        language="c++",
    )
]

setup(ext_modules=cythonize(CY_MODULES), include_package_data=False)
