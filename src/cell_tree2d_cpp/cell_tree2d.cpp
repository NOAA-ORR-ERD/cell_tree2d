/*
 * CellTree2D.cpp
 *
 *  Created on: Aug 20, 2015
 *      Author: jay.hennen
 */
#include "cell_tree2d.h"

using namespace std;

class CellTree2D::bucket {
public:
    //range of the bucket
    double max;
    double min;
    //range of the bounding boxes inside the bucket
    double Rmin;
    double Lmax;
    int index;          //pointer into the bounding box index array
    int size;           //number of bbs in this bucket
};

//tests whether the centroid of the bounding box in the selected dimension falls within this bucket
class CellTree2D::centroid_test {
public:
        unsigned int dim;
        vector<vector<double> >& ds;
        bucket b;

        centroid_test( unsigned int _d, vector<vector<double> >& _ds, bucket _b ) : dim(_d), ds(_ds), b(_b) {}

        bool operator()( const int bb ){ //format is xmin,xmax,ymin,ymax
            vector<double>& box = ds.at(bb);
            double point = (box.at(2*dim+1) - box.at(2*dim)) / 2  + box.at(2*dim);
            return point >= b.min && point < b.max;
        }
    };

CellTree2D::CellTree2D() {
    this->vertices = NULL;
    this->faces = NULL;
    this->n_verts = 0;
    this->poly_data = NULL;
    this->n_polys = 0;
    this->num_buckets = 4;
    this->boxes_per_leaf = 2;

}

void CellTree2D::add_vertices(double* verts, unsigned int v_len) {
    this->vertices = new double*[v_len];
    this->vertices[0] = verts;
    for (unsigned int i = 1; i < v_len; i++)
        this->vertices[i] = this->vertices[i-1] + 2; //points always have 2 coordinates
    this->v_len = v_len;
}

void CellTree2D::add_polys(int* poly_data, unsigned int n_polys, unsigned short n_verts) {
    this->faces = new int*[n_polys];
    this->faces[0] = poly_data;
    for (unsigned int i = 1; i < n_polys; i++)
        this->faces[i] = this->faces[i-1] + n_verts; //you can have n vertices in your polygon
    this->n_verts = n_verts;
    this->max_n_verts = 0;
    this->n_polys = n_polys;
    this->n_verts_arr = NULL;
}

void CellTree2D::add_polys(int* poly_data, unsigned short* n_verts_arr, unsigned int n_polys) {
    this->faces = new int*[n_polys];
    this->faces[0] = poly_data;
    for (unsigned int i = 1; i < n_polys; i++) {
        this->faces[i] = this->faces[i-1] + n_verts_arr[i-1];
    }
    this->n_verts = 0;
    this->max_n_verts = 0;
    this->n_polys = n_polys;
    this->n_verts_arr = n_verts_arr;
}

void CellTree2D::add_polys(int* poly_data, unsigned short max_n_verts, unsigned int n_polys) {
    this->faces = new int*[n_polys];
    this->faces[0] = poly_data;
    for (unsigned int i = 1; i < n_polys; i++) {
        this->faces[i] = this->faces[i-1] + max_n_verts;
    }
    this->n_verts = 0;
    this->max_n_verts = max_n_verts;
    this->n_polys = n_polys;
    this->n_verts_arr = new unsigned short[n_polys];
    for (unsigned int i = 0; i < n_polys; i++) {
        int* poly = this->faces[i];
        int* poly_end = poly;
        while (*poly_end != -1 && poly_end != poly + max_n_verts) {
            poly_end++;
        }
        this->n_verts_arr[i] = (unsigned short)(poly_end - poly);
    }
}

void CellTree2D::finalize(int n_buckets, int bb_per_leaf) {
    num_buckets = n_buckets;
    boxes_per_leaf = bb_per_leaf;
    build_BB_vector();
    if (nodes.size() == 0) { //special case for root node
        nodes.push_back(node(0,bb_indices.size(),0));
    }
    build(0,0);
    for (unsigned int i = 0; i < dataset.size(); i++){
        std::vector<double>(dataset[i]).swap(dataset[i]);
    }
    dataset.clear();
    std::vector<std::vector<double> > (dataset).swap(dataset);
    return;
}

CellTree2D::CellTree2D(double* vertices, unsigned int v_len, int* faces, unsigned int n_polys, unsigned short n_verts, int n_buckets, int bb_per_leaf) {
    //in comes the bounding box vector...this is the underlying data store
    //tree is built on top of it
    this->add_vertices(vertices, v_len);
    this->add_polys(faces, n_verts, n_polys);
    this->finalize(n_buckets, bb_per_leaf);
}

CellTree2D::~CellTree2D() {delete[] faces; delete[] vertices;}

//finds the bounding box of each triangle and inserts it into the dataset and adds it's index to bb_indices
void CellTree2D::build_BB_vector() {
    bb_indices.resize(n_polys);
    dataset.resize(n_polys);
    double v[4];
    for (unsigned int i = 0; i < n_polys; i++) {
        int* poly = faces[i];
        unsigned short poly_len = 0;
        if (n_verts == 0) {
            poly_len = this->n_verts_arr[i];
        } else {
            poly_len = n_verts;
        }
        double x_min, x_max, y_min, y_max;
        double* vt = vertices[poly[0]];
        x_min = x_max = vt[0];
        y_min = y_max = vt[1];
        for(unsigned short j = 1; j < poly_len; j++) {
             vt = vertices[poly[j]];
             x_min = min(x_min, vt[0]);
             x_max = max(x_max, vt[0]);
             y_min = min(y_min, vt[1]);
             y_max = max(y_max, vt[1]);
        }
        v[0] = x_min;
        v[1] = x_max;
        v[2] = y_min;
        v[3] = y_max;
        dataset[i] = vector<double>(v, v + sizeof v / sizeof v[0]);
        bb_indices[i]=i;
    }
}

//Finds the range of the bounding boxes contained by a bucket in dimension d
void CellTree2D::get_bounds(bucket& buk, int d) {
    buk.Rmin = std::numeric_limits<double>::max();
    buk.Lmax = -std::numeric_limits<double>::max();
    for (int i = buk.index; i < (buk.index + buk.size); i++) {
        int data_index = bb_indices.at(i);
        if (dataset.at(data_index).at(2*d) < buk.Rmin) {
            buk.Rmin = dataset.at(data_index).at(2*d);
        }
        if (dataset.at(data_index).at(2*d+1) > buk.Lmax) {
            buk.Lmax = dataset.at(data_index).at(2*d+1);
        }
    }
}



//Sorts all the bounding boxes specified by node into the buckets. Does alter order of bb_indices!
//This should partition all the bounding boxes within the current node into the buckets.
//The idea is that when you split the node, the bb_indices are already arranged
//and you can retrieve them with an integer+size rather than a vector
void CellTree2D::sort_bbs(vector<bucket>& buks, node& node, int dim) {
    vector<int>::iterator current = bb_indices.begin() + node.ptr;
    vector<int>::iterator end = bb_indices.begin() + node.ptr + node.size;
    buks.at(0).index = node.ptr;
    for (unsigned int i = 1; current != end; i++) {
        bucket& b = buks.at(i-1);
        current = std::stable_partition(current,end,centroid_test(dim,dataset,b));
        vector<int>::iterator start = bb_indices.begin()+b.index;
        buks.at(i-1).size = current - start;
        if(i < buks.size()){
            buks.at(i).index = buks.at(i-1).index + buks.at(i-1).size;
        }
        start = current;
    }
}

void CellTree2D::build(int root_ind, int dim) {
    int dim_flag = dim;
        if (dim < 0)
            dim+=2;
    node& root = nodes.at(root_ind);
    //is it a leaf? if so, we're done, otherwise split
    if (root.size <= boxes_per_leaf) {
        return;
    }

    //find bounding range of node's entire dataset in dimension 0 (x-axis)
    bucket range;
    range.index = root.ptr;
    range.size = root.size;
    get_bounds(range, dim);

    vector<bucket> buks;
    double bucket_len = (range.Lmax - range.Rmin) / num_buckets;
    //create buckets and specify their ranges
    for (int b_ind=0;b_ind < num_buckets; b_ind++) {
        bucket b = {(b_ind+1)*bucket_len + range.Rmin, b_ind*bucket_len + range.Rmin};
        buks.push_back(b);
    }

    //Now that bucket ranges are setup, sort bounding boxes contained in node into buckets in the dimension specified.
    sort_bbs(buks, root, dim);

    //Determine Lmax & Rmin for each bucket
    for (int b_ind=0;b_ind < num_buckets; b_ind++) {
        get_bounds(buks.at(b_ind),dim);
    }

    //Special case: 2 bounding boxes share the same centroid, but boxes_per_leaf is 1
    //This will break most of the usual bucketing code. Unless the grid has overlapping
    //triangles (which it shouldnt!) this is the only case to deal with
    if (boxes_per_leaf == 1 && root.size == 2) {
        root.Lmax = range.Lmax;
        root.Rmin = range.Rmin;
        node left_child(root.ptr, 1, !dim);
        node right_child(root.ptr+1, 1, !dim);
        root.child = nodes.size();
        nodes.push_back(left_child);
        nodes.push_back(right_child);
        return;
    }

    while (buks[0].size == 0) {
        buks[1].min = buks[0].min;
        buks.erase(buks.begin());
    }
    for (unsigned int b = 1;b < buks.size(); b++) {
        bucket& cur_buk = buks.at(b-1);
        bucket& next_buk = buks.at(b);
        //if a empty bucket is encountered, merge it with the previous one and continue as normal. As long as the
        //ranges of the merged buckets are still proper, calcualting cost for empty buckets can be avoided, and
        //the split will still happen in the right place
        if (next_buk.size == 0) {
            cur_buk.max = next_buk.max;
            buks.erase(buks.begin() + b);
            b--;
        }
    }
    //CHECK....are all the triangles in one bucket? If so, restart and switch dimension
    for(int b_ind = 0; b_ind < buks.size(); b_ind++) {
        if (buks.at(b_ind).size == root.size){
            if (dim_flag >= 0) { //dim_flag will be negative after one switch
                dim_flag = !dim - 2;
                root.dim = !root.dim;
                build(root_ind,dim_flag);
            } else { //Already split once
                //can't split...convert to leaf
                root.Lmax = -1;
                root.Rmin = -1;
            }
            return;
        }
    }

    //plane is the separation line to split on...0 [bucket0] 1 [bucket1] 2 [bucket2] 3 [bucket3] 4
    double plane_cost;
    double plane_min_cost = std::numeric_limits<double>::max();
    int plane = std::numeric_limits<int>::max();

    //if we split here, lmax is from bucket 0, and rmin is from bucket 1
    //after computing those, we can compute the cost to split here, and if this is the minimum, we split here.
    int bbs_in_left = 0;
    int bbs_in_right = 0;
    for (unsigned int b = 1;b < buks.size(); b++) {
        bucket& cur_buk = buks.at(b-1);
        bucket& next_buk = buks.at(b);
        bbs_in_left += cur_buk.size;
        bbs_in_right = root.size-bbs_in_left;
        //compute volumes bounded by left & right.
        double lvol = (cur_buk.Lmax - range.Rmin) / bucket_len;
        double rvol = (range.Lmax - next_buk.Rmin) / bucket_len;
        plane_cost = lvol*bbs_in_left + rvol*(bbs_in_right);
        if (plane_cost < plane_min_cost) {
            plane_min_cost = plane_cost;
            plane = b;
        }
    }
    //we've found the plane, now simply split by creating the children nodes
    //by assigning the appropriate buckets. The underlying information is already
    //rearranged.
    int right_index = buks.at(plane).index;
    int right_size = root.ptr + root.size - right_index;
    int left_index = root.ptr;
    int left_size = root.size - right_size;
    root.Lmax = buks.at(plane-1).Lmax;
    root.Rmin = buks.at(plane).Rmin;
    node left_child(left_index, left_size, !dim);
    node right_child(right_index, right_size, !dim);
    root.child = nodes.size();
    int child_ind = root.child;
    nodes.push_back(left_child);
    nodes.push_back(right_child);
    build(child_ind, left_child.dim);
    build(child_ind+1, right_child.dim);
}

//returns true if test node is within the polygon
//checking the point immediately is important....potentially saves a lot of traversal if found.
bool CellTree2D::point_in_poly (int bb, double* test){
    int* f = faces[bb];
    unsigned short poly_len = 0;
    if (n_verts_arr != NULL) {
        poly_len = n_verts_arr[bb];
    } else {
        poly_len = n_verts;
    }
    unsigned short i, j = 0;
    bool c = 0; /*really need a bool here...*/
    for (i = 0, j = poly_len-1; i < poly_len; j = i++) {
        double* v1 = vertices[f[i]];
        double* v2 = vertices[f[j]];
        if ( ((v1[1]>test[1]) != (v2[1]>test[1])) &&
                (test[0] < (v2[0]-v1[0]) * (test[1]-v1[1]) / (v2[1]-v1[1]) + v1[0]) )
            c = !c;
    }
    return c;
}

int CellTree2D::locate_point_helper(double* point, int node) {
    //where does point lie in this node's dimension?
    //if l && r, go both, otherwise go one or the other, or leave results empty
    CellTree2D::node& current = nodes[node];
    if (current.size <= boxes_per_leaf || current.Lmax == -1) {
        for (int i = current.ptr; i < current.ptr + current.size; i++) {
            if (point_in_poly(bb_indices[i],point)){
                return bb_indices[i];
            }
        }
        return -1;
    }
    int d = current.dim ? 1 : 0;
    bool l = point[d] <= current.Lmax;
    bool r = point[d] >= current.Rmin;
    int ret = -1;
    if (l && r) {
        if( current.Lmax-point[d] < point[d]-current.Rmin ){
            // go left first
            ret = locate_point_helper(point,current.child);
            if(ret == -1)
                return ret = locate_point_helper(point,current.child+1);
        } else {
            // go right first
            ret = locate_point_helper(point,current.child+1);
            if(ret == -1)
                 return ret = locate_point_helper(point,current.child);
        }
    } else if (l) {
        return ret = locate_point_helper(point,current.child);
    } else if (r) {
        return ret = locate_point_helper(point,current.child+1);
    }
    return ret;
}

// Finds the index of the triangle for every point in pts, storing result in res
void CellTree2D::locate_points(double* pts, int* res, int len) {
    for (int i = 0; i < len; i++) {
        res[i] = locate_point_helper(&pts[2*i],0);
    }
}

int CellTree2D::size() {
    return nodes.size();
}
