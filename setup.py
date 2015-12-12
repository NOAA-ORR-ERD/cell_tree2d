#!/usr/bin/env python

try:
    from setuptools import setup, Extension
except ImportError:
    print "You need setuptools to build this module"

from Cython.Build import cythonize

import numpy as np #for the include dirs...
import os, sys


include_dirs = [np.get_include(),
                os.path.join('.', 'src')]

if sys.platform.startswith('win'):
    # need the stdint header for Windows (VS2008)
    include_dirs.append(os.path.join('.','src', 'win_headers'))


ext_modules=[ Extension("cell_tree2d.cell_tree2d",
                        ["cell_tree2d/cell_tree2d.pyx", "src/cell_tree2d.cpp"],
                        include_dirs = include_dirs,
                        language = "c++",
                         )]

setup(
    name = "cell_tree2d",
    version='0.1.1',
    description = "python wrappers around Cell-Tree 2D spatial index",
    long_description=open('README.md').read(),
    author = "Jay Hennen",
    author_email = "jay.hennen@noaa.gov",
    url="https://github.com/NOAA-ORR-ERD",
    license = "Public Domain",
    #keywords = "",
    ext_modules = cythonize(ext_modules),
    packages = ["cell_tree2d", "cell_tree2d/test"],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: Public Domain",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 2 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Utilities",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
      )
