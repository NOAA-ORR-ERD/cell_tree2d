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

typedef std::vector<double> BB;



class CellTree2D {
public:
    ~CellTree2D();
    //the double* and int* are pointers to the first element of a 1D multidimensional array
    //inside it builds a pointer array to re-impose the multidimensional structure
    CellTree2D( double*,int, int*,int, int, int, int);
    void build(int, int);

    int FindBoxLeaf(double*);
    int size();
protected:
    class node;
    class bucket;
    class centroid_test; //for partition
private:
    bool point_in_poly (int, double*);
    void build_BB_vector();
    int FindBoxLeafHelper(double*, int);
    void get_bounds(bucket&, int);
    void sort_bbs(std::vector<bucket>&, node&, int);
    double** vertices;
    int** faces;
    int num_buckets;
    int boxes_per_leaf;
//    std::vector<vertex>& vertices;
//    std::vector<face>& faces;
    int poly; // # of vertices per face (affects bounding box calculation)
    int v_len;
    int f_len;
    std::vector<std::vector<double> > dataset;
    //bb_indices initially appears as : [0,1,2,3,...,dataset.size()-1]
    //as the tree is built, it will be partitioned and sub-partitioned, until the order of elements
    //represents the in-order traversal of every leaf node
    std::vector<int> bb_indices;
    std::vector<node> nodes;
};

#endif /* RELEASE_CellTree2D_H_ */
