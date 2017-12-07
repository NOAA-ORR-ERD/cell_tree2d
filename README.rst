============
cell_tree2d
============

A module that provides the CellTree data structure as described by Garth and Joy in their 2010 paper:

http://graphics.cs.ucdavis.edu/~joy/NSF-IIS-0916289/Papers/GarthTVCG2011.pdf  

This implementation is 2D specific and includes some additions useful to answering one question:

"What is the index of the polygon that contains this point?"

Installation
------------

cell_tree2d should be buildable on any system with a properly set up compiler to work with python with the usual::

    python setup.py build
    python setup.py install


It is also available on PyPi::

    pip install cell_tree2d

But only as source, so the same compiler setup will be required.


conda packages
..............

For people using the Anaconda or miniconda python distributions, there are pre-built conda packges available on the conda-forge channel::

    conda install -c conda-forge cell_tree2d



Algorithm Notes
---------------

There are two major benefits to this algorithm over other types of BVHs. First is that overlaps in volumes
bounded by nodes do not create duplicates, decreasing the memory footprint. Secondly, the tree is balanced
by splitting a node along the plane that minimizes a cost function that 'weighs' each half. The result is
a tree with no duplicates and that becomes increasingly balanced as the number of buckets used when building
the tree increases (though this linearly increases build time).

This is a 2D version of the algorithm that is oriented towards maximum lookup speed and immediate answer
checking. The basic CellTree does not hold enough information to give conclusive point-in-polygon answers; 
the best it can do is provide the (short) list of cell bounding boxes that contain the point. Since a point
can be within the bounds of two different cells, and it is very possible both children of a parent node will 
need to be searched, implementing immediate point-in-polygon checks on each cell as they are encountered is
highly beneficial, as an early success will avoid all further tree traversal.

Usage Notes
-----------

The tree needs certain information to be built:

1. 'verts' - A 2xV numpy array containing x/y coordinates of the V vertices   

2. 'faces' - A PxN numpy array containing N arrays of P indices of vertices that describe one
   'face' or polygon of degree P  

3. 'num_buckets' - The # of buckets desired. Must be >= 2 The default is 4. Values higher than
   8 begin to provide diminishing returns.  

4. 'cells\_per\_leaf' - The # of polygons per leaf node. The default is 2. Using 1 is possible,
   but doubles memory footprint for only slightly. faster lookup. If memory is a concern, this
   value can be increased, but lookup performance will quickly be impacted  

**IMPORTANT:** 'verts' and 'faces' MUST describe a *properly formed* unstructured grid. Assume
that degenerate (0 area) or overlapping polygons **will** cause a build failure. If the construction
of the tree causes a segfault, this is probably the cause.


