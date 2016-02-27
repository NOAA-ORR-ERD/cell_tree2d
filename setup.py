#!/usr/bin/env python

try:
    from setuptools import setup, Extension
    from setuptools.command.test import test as TestCommand
except ImportError:
    print("You need setuptools to build this module.")

from Cython.Build import cythonize

import os
import sys
import numpy as np  # For the include directory.


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.verbose = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


include_dirs = [np.get_include(),
                os.path.join('.', 'src')]

# Need the stdint header for Windows (VS2008).
if sys.platform.startswith('win') and sys.version_info.major <= 2:
    include_dirs.append(os.path.join('.', 'src', 'win_headers'))


ext_modules = [Extension("cell_tree2d.cell_tree2d",
                         ["cell_tree2d/cell_tree2d.pyx",
                          "src/cell_tree2d.cpp"],
                         include_dirs=include_dirs,
                         language="c++",
                         )]

setup(
    name="cell_tree2d",
    version='0.1.1',
    description="Python wrappers around Cell-Tree 2D spatial index",
    long_description=open('README.rst').read(),
    author="Jay Hennen",
    author_email="jay.hennen@noaa.gov",
    url="https://github.com/NOAA-ORR-ERD",
    license="Public Domain",
    # keywords = "",
    ext_modules=cythonize(ext_modules),
    packages=["cell_tree2d", "cell_tree2d/test"],
    cmdclass=dict(test=PyTest),
    install_requires=['numpy'],
    setup_requires=['cython>0.23', 'setuptools'],
    tests_require=['pytest'],    
    classifiers=[
        "Development Status :: 4 - Beta",
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
