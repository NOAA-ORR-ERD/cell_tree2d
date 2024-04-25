/*

 * testmain.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: jay.hennen
 *
 *  Description:
 *  Build and run this program with a test file to ensure the CellTree works
 *  and to measure the effects of performance tweaks.
 *
 *  Test file format:
 *  0   nodes:
 *  1   0,0
 *  2   2,0
 *  3   1,2
 *  4
 *  5   faces:
 *  6   0,1,2
 *  7
 */

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iostream>

#include "cell_tree2d.h"
using namespace std;
using namespace std::chrono;


int parse_nodes_to_coords_vector(ifstream& reader, double*** nodes) {
    string buf;
    getline(reader, buf);
    if (string(buf).compare("nodes") != 0) {
        cerr << "not at nodes section of file" << endl;
        return -1;
    }
    streampos start = reader.tellg();
    int len = 0;
    for(;buf[0] != '\0';len++) {
        getline(reader,buf);
    }
    len--;
    reader.seekg(start);
    *nodes = new double*[len];
    (*nodes)[0] = new double[len * 2];
    for (int i = 1; i < len; i++)
        (*nodes)[i] = (*nodes)[i-1] + 2;


    getline(reader, buf);
    for(int i = 0 ; buf[0] != '\0'; i++) {
        string coord;
        istringstream s(buf);
        getline(s,coord,',');
        double x = strtod(coord.c_str(),NULL);
        getline(s,coord);
        double y = strtod(coord.c_str(),NULL);
        (*nodes)[i][0] = x;
        (*nodes)[i][1] = y;
        getline(reader, buf);
    }
    cout << "Processed nodes" << endl;
    return len;
}

int parse_faces_to_vector(ifstream& reader, int*** faces) {


    string buf;
    getline(reader, buf);
    if (string(buf).compare("faces") != 0) {
        cerr << "not at faces section of file" << endl;
        return -1;
    }
    streampos start = reader.tellg();
    int len = 0;
    for(;buf[0] != '\0';len++) {
        getline(reader,buf);
    }
    len--;
    reader.clear(ios::goodbit);
    reader.seekg(start);
    *faces = new int*[len];
    (*faces)[0] = new int[len * 3];
    for (int i = 1; i < len; i++)
        (*faces)[i] = (*faces)[i-1] + 3;

    getline(reader, buf);
    for(int i = 0;buf[0] != '\0';i++) {
        string coord;
        istringstream s(buf);
        getline(s,coord,',');
        int a = atoi(coord.c_str());
        getline(s,coord,',');
        int b = atoi(coord.c_str());
        getline(s,coord);
        int c = atoi(coord.c_str());
        (*faces)[i][0] = a;
        (*faces)[i][1] = b;
        (*faces)[i][2] = c;
        getline(reader, buf);
    }
    cout << "Processed faces" << endl;
    return len;
}

void build_BB_vector(double** nodes, int** faces, int faces_len, vector<vector<double> >& BB_vec) {
    for (int i = 0; i < faces_len; i++) {
        int * f = faces[i];
        double* a = nodes[f[0]];
        double* b = nodes[f[1]];
        double* c = nodes[f[2]];
        double xmin = min(a[0],min(b[0],c[0]));
        double xmax = max(a[0],max(b[0],c[0]));
        double ymin = min(a[1],min(b[1],c[1]));
        double ymax = max(a[1],max(b[1],c[1]));
        BB_vec.push_back({xmin,xmax,ymin,ymax});
    }

}

//guaranteed to return a point within the grid
double* get_contained_point(double** nodes, int nodes_len, int** faces, int faces_len) {
    int r = rand() % faces_len;
    int* f = faces[r];
    double* a = nodes[f[0]];
    double* b = nodes[f[1]];
    double* c = nodes[f[2]];
    double* ret = new double[2];
    ret[0] = ((a[0]+b[0]+c[0]) / 3);
    ret[1] = ((a[1]+b[1]+c[1]) / 3);
    return ret;
}

//determinant of 2 vectors
double det(double* u, double* v) {
    return u[0]*v[1]-u[1]*v[0];
}

//returns true if test node is within the triangle abc
bool point_in_triangle (double* a, double* b, double* c, double* test)
{
    double v[] = {test[0] - a[0], test[1]-a[1]};
    double v0[] = {0,0};
    double v1[] = {b[0] - a[0], b[1]-a[1]};
    double v2[] = {c[0] - a[0], c[1]-a[1]};
    double c1 = (det(v,v2) - det(v0,v2))/det(v1,v2);
    double c2 = -(det(v,v1) - det(v0,v1))/det(v1,v2);

    return c1 > 0 && c2 > 0 && c1 + c2 < 1;
}


// Checks the output of the CellTree2D search against the result of a linear search.
bool result_validation_test(int n, double** nodes, int nodes_len, int** faces, int faces_len, CellTree2D& bbt) {
    cout << "Running result validation test..." << endl;
    for (int i = 0; i < n; i++) {
        double* test_point = get_contained_point(nodes, nodes_len, faces, faces_len);
        vector<int> boxes;
        int bbt_result = bbt.FindBoxLeaf(test_point);
        int* res = faces[bbt_result];
        if (point_in_triangle(nodes[res[0]],nodes[res[0]],nodes[res[0]],test_point)) {
            cerr << "INCONSISTENT RESULTS!" << endl;
            cerr << "Test point: " << test_point[0] << "," << test_point[1] << endl;
            cerr << "Test index: " << i << endl;
            cerr << "CellTree2D Result: " << bbt_result << endl;
            free(test_point);
            return false;
        }
        free(test_point);
    }
    cout << "\tResults validation test successful" << endl;
    cout << "\tSuccessfully found " << n << " points within grid" << endl << endl;
    return true;
}

// tests a number of points off the grid, including several at the limits of the coordinate system
//bool point_off_grid_test(vector<node>& nodes, vector<face>& faces, CellTree2D& bbt) {
//    cout << "Running point_off_grid test..." << endl;
//    vector<node> test_points ={
//            {std::numeric_limits<double>::max(),std::numeric_limits<double>::max()},
//            {std::numeric_limits<double>::max(), 0},
//            {std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()},
//            {0,std::numeric_limits<double>::max()},
//            {0,0},
//            {0,-std::numeric_limits<double>::max()},
//            {-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()},
//            {-std::numeric_limits<double>::max(),0},
//            {-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()}};
//    vector<int> boxes;
//    int result = -1;
//    int brute_result = 0;
//    for (unsigned int i = 0; i < test_points.size(); i++) {
//        brute_result = linear_node_search(nodes,faces, test_points[i]);
//        if (brute_result == -1) {
//            node test_point = test_points.at(i);
//            result = bbt.FindBoxLeaf(test_point.data());
//            if (result != -1) {
//                cerr << "INCONSISTENT RESULTS!" << endl;
//                cerr << "Test point: " << test_point[0] << "," << test_point[1] << endl;
//                cerr << "Test index: " << i << endl;
//                cerr << "CellTree2D Result: " << result << endl;
//                cerr << "Linear Result: " << brute_result << endl;
//                return false;
//            }
//        }
//    }
//    cout << "\tPoints off grid test successful" << endl << endl;
//    return true;
//}

void performance_test(double** nodes, int nodes_len, int** faces, int faces_len, CellTree2D& bbt) {
    cout << "Starting performance test...populating test point array" << endl;
    int n = 1000000;
    double* test_points = new double[2*n];
    cout << "\tStarting clock..." << endl;

    //just pick a point close to a node
    for(int i = 0;i < n;i++) {
        double x = nodes[i%nodes_len][0]+0.0000000001*i;
        double y = nodes[i%nodes_len][1]+0.0000000001*i;
        test_points[i] = x;
        test_points[i+1] = y;
    }
    int res;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    for (int i = 0; i < n ; i++) {
        res = bbt.FindBoxLeaf(&test_points[i*2]);
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto time_span = duration_cast<milliseconds>(end - start);
    std::cout << "\tIt took me " << time_span.count() << " milliseconds to check " << n << " points with CellTree2D";
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    string filename =string(argv[1]);
    ifstream reader(argv[1]);
    if(reader.fail()) {
        return 1;
    }
    int nodes_len;
    double** nodes;
    int faces_len;
    int** faces;
    nodes_len = parse_nodes_to_coords_vector(reader, &nodes);
    faces_len = parse_faces_to_vector(reader, &faces);
    high_resolution_clock::time_point start = high_resolution_clock::now();

    CellTree2D bbt(nodes[0],nodes_len,faces[0],faces_len, 3, 8,2);
    high_resolution_clock::time_point build = high_resolution_clock::now();
    auto time_span = duration_cast<milliseconds>(build - start);
    std::cout << "It took me " << time_span.count() << " milliseconds to build a CellTree2D of "<< bbt.size() << " nodes" << endl;

    //test calls
    result_validation_test(1000,nodes, nodes_len, faces, faces_len, bbt);
//    point_off_grid_test(nodes,faces,bbt);
    performance_test(nodes, nodes_len, faces, faces_len,bbt);

    return 0;
}

