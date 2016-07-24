import cython
import numpy as np
cimport numpy as cnp
from libcpp.vector cimport vector

cdef extern from "cell_tree2d.h" :

    cdef cppclass CellTree2D:
        CellTree2D(double*, int, int*, int, int, int, int) except +
        int locate_points(double*, int*, int)
        int size()
        int num_buckets, boxes_per_leaf, poly, v_len, n_polys
        
        cppclass node:
            node()
            int child
            double Lmax
            double Rmin
            int ptr
            int size
            bint dim
        
        vector[node] nodes
        vector[int] bb_indices
        
#     ctypedef struct CellTree2D.node:
#         int child
#         double Lmax
#         double Rmin
#         int ptr
#         int size
#         bint dim
        

cdef class CellTree:

    cdef object verts, faces
    cdef CellTree2D* thisptr

    def __cinit__(self,
                  verts,
                  faces,
                  int num_buckets=4,
                  int cells_per_leaf=2):

        """
        Initilize a CellTree

        :param verts: The vertices of the nodes of the mesh
        :type verts: Nx2 numpy array of float64, or something that can be turned into one

        :param faces: The indices of the nodes of each face
        :type faces: Nx3 numpy array of np.intc (for triangles -- Nx4 for quads)

        :param num_buckets=4: Number of buckets to divide cells when building tree.
                              More buckets results in a better balanced tree --
                              fewer buckets results in a faster tree build time.
        :type num_buckets: integer

        :param cells_per_leaf=2: number of cells in the final leaf -- more cells in the leaf
                                 results  in a smaller tree (less memory), but then it takes
                                 longer to find which cell you are in at the end.
        :type cells_per_leaf: integer
        """


        cdef cnp.ndarray[double, ndim=2, mode="c"] verts_arr
        cdef cnp.ndarray[int, ndim=2, mode="c"] faces_arr
        cdef int num_verts, num_faces, num_poly_vert

        # convert to numpy arrays:
        verts = np.asarray(verts).astype(np.float64)
        verts = np.ascontiguousarray(verts)
        if len(verts.shape) <> 2 or verts.shape[1] <> 2:
            raise ValueError("verts must be convertible to a Nx2 numpy array of float64")
        verts_arr = verts
        self.verts=verts

        faces = np.asarray(faces).astype(np.int32)
        faces = np.ascontiguousarray(faces)
        if  len(faces.shape)<>2 or not (faces.shape[1] == 3 or faces.shape[1] == 4):
            raise ValueError("faces must be convertible to a Nx3 (for triangles)"
                             "or Nx4 (for quads) numpy array of int32")
        faces_arr = faces
        self.faces=faces

        ## a bit more error checking:
        if num_buckets < 2:
            raise ValueError("num_buckets must be an integer greater than 2")
        if cells_per_leaf < 1:
            raise ValueError("cells_per_leaf must be >= 1")

        num_verts = verts.shape[0]
        num_faces = faces.shape[0]
        num_poly_vert = faces.shape[1]
        self.thisptr = new CellTree2D(&verts_arr[0,0],
                                      num_verts,
                                      &faces_arr[0,0],
                                      num_faces,
                                      num_poly_vert,
                                      num_buckets,
                                      cells_per_leaf)

    def __del__(self):
        del self.thisptr

    @property
    def size(self):
        return self.thisptr.size()
    
    @property
    def num_buckets(self):
        return self.thisptr.num_buckets
    
    @property
    def boxes_per_leaf(self):
        return self.thisptr.boxes_per_leaf
    
    @property
    def bb_indices(self):
        l = []
        for i in range(0, self.thisptr.n_polys):
            l.append(self.thisptr.bb_indices[i])
        return l
    
    @property
    def nodes(self):
        l = []
        for i in range(0,self.size):
            l.append((self.thisptr.nodes[i].child,
                      self.thisptr.nodes[i].Lmax,
                      self.thisptr.nodes[i].Rmin,
                      self.thisptr.nodes[i].ptr,
                      self.thisptr.nodes[i].size,
                      self.thisptr.nodes[i].dim))
        return l

    cdef c_locate(self, double* points, int* results, int len):
        self.thisptr.locate_points(points, results, len)

    def locate(self, points_in):
        cdef int size
        cdef cnp.ndarray[int, ndim=1, mode="c"] locations
        
        # convert to memoryview:
        cdef cnp.ndarray[double, ndim=2, mode="c"] points
        points_in = np.ascontiguousarray(points_in, dtype=np.float64)
        if len(points_in.shape) < 2: #single [x,y]
            points_in = np.expand_dims(points_in, axis=0)
        if points_in.shape[1] <> 2:
            raise ValueError("points must be convertible to a Nx2 numpy array of float64")
        points = points_in


        size = points.shape[0]
        locations = np.zeros((size,), dtype=np.intc)

        self.c_locate(&points[0, 0], &locations[0], size)

        return locations

