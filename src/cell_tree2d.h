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
#include <stdint.h>

typedef std::vector<double> BB;



class CellTree2D {
public:
    ~CellTree2D();
    //the double* and int32_t* are pointers to the first element of a 1D multidimensional array
    //inside it builds a pointer array to re-impose the multidimensional structure
    CellTree2D( double*,int32_t, int32_t*,int32_t, int32_t, int32_t, int32_t);
    void build(int32_t, int32_t);

    int32_t FindBoxLeaf(double*);
    int32_t size();
protected:
    class node;
    class bucket;
    class centroid_test; //for partition
private:
    bool point_in_poly (int32_t, double*);
    void build_BB_vector();
    int32_t FindBoxLeafHelper(double*, int32_t);
    void get_bounds(bucket&, int32_t);
    void sort_bbs(std::vector<bucket>&, node&, int32_t);
    double** vertices;
    int32_t** faces;
    int32_t num_buckets;
    int32_t boxes_per_leaf;
//    std::vector<vertex>& vertices;
//    std::vector<face>& faces;
    int32_t poly; // # of vertices per face (affects bounding box calculation)
    int32_t v_len;
    int32_t f_len;
    std::vector<std::vector<double> > dataset;
    //bb_indices initially appears as : [0,1,2,3,...,dataset.size()-1]
    //as the tree is built, it will be partitioned and sub-partitioned, until the order of elements
    //represents the in-order traversal of every leaf node
    std::vector<int32_t> bb_indices;
    std::vector<node> nodes;
};

#endif /* RELEASE_CellTree2D_H_ */
