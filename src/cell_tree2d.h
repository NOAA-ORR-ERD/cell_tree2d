/*
 * .h
 *
 *  Created on: Aug 20, 2015
 *      Author: jay.hennen
 */

#ifndef CELLTREE2D_H_
#define CELLTREE2D_H_

#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <iterator>
#include <vector>
#include <iostream>
#include "grid_data.h"

typedef std::vector<double> BB;

struct poly{
    int len;
    int* verts_ptr;
    bool contains_point(double*);
    void bounding_box(double*);

};

class CellTree2D {
public:
    class node{
    public:
        node(int pt, int sz, int d):
            child(-1), Lmax(-1), Rmin(-1), ptr(pt), size(sz), dim(d){};
        ~node(){};
        int child;          //index of left child...right child is child+1
        double Lmax;
        double Rmin;
        int ptr;            //index into the bounding box index array
        int size;           //number of bbs in this node
        bool dim;           //false = 0 = x, true = 1 = y;
    };
    class bucket;
    class centroid_test; //for partition
    ~CellTree2D();
    //the double* and int* are pointers to the first element of a 1D multidimensional array
    //inside it builds a pointer array to re-impose the multidimensional structure
    CellTree2D( double*,int, int*,int, int, int, int);
    void build(int, int);

    void locate_points(double*, int*, int);
    int size();

    bool point_in_poly (int, double*);
    void build_BB_vector();
    int locate_point_helper(double*, int);
    void get_bounds(bucket&, int);
    void sort_bbs(std::vector<bucket>&, node&, int);

    double* vertices;
    unsigned int v_len;

    int* poly_data
    unsigned int len_poly_data

    int* faces;
    int num_buckets;
    int boxes_per_leaf;
//    std::vector<vertex>& vertices;
//    std::vector<face>& faces;
    int f_len;
    bool mixed_polys
    int* n_verts; // # of vertices per face (affects bounding box calculation)
    int max_verts
    std::vector<std::vector<double> > dataset;
    //bb_indices initially appears as : [0,1,2,3,...,dataset.size()-1]
    //as the tree is built, it will be partitioned and sub-partitioned, until the order of elements
    //represents the in-order traversal of every leaf node
    std::vector<poly> polys;
    std::vector<node> nodes;
};

#endif /* RELEASE_CellTree2D_H_ */
