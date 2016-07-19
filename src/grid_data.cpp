#include grid_data.h

int structure_data::add_poly(int* p, int len) {
    this->polys.push_back(new grid_data::poly(len, p));
}


