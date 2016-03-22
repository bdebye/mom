

#ifndef GMSH_H_
#define GMSH_H_


#include <vector>
#include <array>
#include <string>

#include "type.h"

struct mesh_data {

    int triang_num;
    int segment_num;
    int node_num;

    std::vector<Point> node;
    std::vector<blitz::TinyVector<int, 2>> segment;
    std::vector<blitz::TinyVector<int, 3>> triang;

};

extern mesh_data mesh;

void load_gmsh_file(std::string filename);


Point get_point(int n);

TrianIndex get_triangle_index(int n);

Point get_triangle_center(int n);

double triangle_area(int n);

Triangle get_triangle(int n);

double geometry_size();

void test_mesh();

#endif /* GMSH_H_ */
