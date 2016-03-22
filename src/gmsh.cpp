
#include "gmsh.h"
#include "mathlib.h"
#include <fstream>
#include <sstream>

using namespace std;
mesh_data mesh;

void clearMesh() {
    mesh.node_num = 0;
    mesh.triang_num = 0;
    mesh.segment_num = 0;

    mesh.node.clear();
    mesh.segment.clear();
    mesh.triang.clear();
}

void readNode(string line) {
    istringstream iss(line);
    Point node;
    int n;
    iss >> n;
    iss >> node(0);
    iss >> node(1);
    iss >> node(2);
    node = node / 100.0; // centi meter
    mesh.node.push_back(node);
}

const int SEGMENT_TYPE = 1;
const int TRIANG_TYPE = 2;
const int TETRAH_TYPE = 4;

void readElement(string line) {
    istringstream iss(line);
    blitz::TinyVector<int, 2> segment;
    blitz::TinyVector<int, 3> triang;
    int n;
    int code;
    iss >> n;
    iss >> code;
    iss >> n; iss >> n; iss >> n;
    if(code == SEGMENT_TYPE) {
        iss >> segment(0);
        iss >> segment(1);
        mesh.segment.push_back(segment);
        mesh.segment_num++;
    }
    else if(code == TRIANG_TYPE) {
        iss >> triang(0);
        iss >> triang(1);
        iss >> triang(2);
        mesh.triang.push_back(triang);
        mesh.triang_num++;
    }
}

void load_gmsh_file(string filename) {
    fstream fs(filename);
    string line;
    istringstream iss;
    int element_num = 0;
    clearMesh();

    if(!fs.good()) {
        cout << "Broken mesh file or file doesn't exist..." << endl;
    }

    mesh.node.push_back(Point()); // The index of nodes starts from 1
    for(int i = 0; i < 4; i++)
        getline(fs, line);
    getline(fs, line);

    iss.str(line);
    iss >> mesh.node_num;
    mesh.node.reserve(mesh.node_num);
    for(int i = 0; i < mesh.node_num; i++) {
        getline(fs, line);
        readNode(line);
    }

    getline(fs, line);
    getline(fs, line);
    getline(fs, line);

    iss.clear();
    iss.str(line);
    iss >> element_num;

    for (int i = 0; i < element_num; i++) {
        getline(fs, line);
        readElement(line);
    }
}


Point get_point(int n) {
	return mesh.node[n];
}

TrianIndex get_triangle_index(int n) {
	return mesh.triang[n];
}

double triangle_area(int n) {
	Triangle trian;
	TrianIndex index = get_triangle_index(n);
	trian(0) = get_point(index(0));
	trian(1) = get_point(index(1));
	trian(2) = get_point(index(2));
	return trian_area(trian);
}


Triangle get_triangle(int n) {
	Triangle trian;
	TrianIndex index = get_triangle_index(n);
	trian(0) = get_point(index(0));
	trian(1) = get_point(index(1));
	trian(2) = get_point(index(2));
	return trian;
}

Point get_triangle_center(int n) {
	TrianIndex trian = get_triangle_index(n);
	return trian_center(get_triangle(n));
}

double geometry_size() {
	double max_x = (get_point(1))(0);
	double max_y = (get_point(1))(1);
	double max_z = (get_point(1))(2);
	double min_x = (get_point(1))(0);
	double min_y = (get_point(1))(1);
	double min_z = (get_point(1))(2);
	for(int i = 2; i <= mesh.node_num; i++) {
		if((get_point(i))(0) > max_x)
			max_x = (get_point(i))(0);
		if((get_point(i))(1) > max_y)
			max_y = (get_point(i))(1);
		if((get_point(i))(2) > max_z)
			max_z = (get_point(i))(2);
		if((get_point(i))(0) < min_x)
			min_x = (get_point(i))(0);
		if((get_point(i))(1) < min_y)
			min_y = (get_point(i))(1);
		if((get_point(i))(2) < min_z)
			min_z = (get_point(i))(2);
	}
	return sqrt(pow(max_x - min_x, 2) + pow(max_y - min_y, 2) + pow(max_z - min_z, 2));
}

bool triangle_has_point(int trian_index, Point p) {
    const double eps = 1e-7;
	Triangle trian = get_triangle(trian_index);
	for(int i = 0; i < 3; i++) {
        double xd = trian(i)(0) - p(0);
        double yd = trian(i)(1) - p(1);
        double zd = trian(i)(2) - p(2);
		if(norm(Point(xd, yd, zd)) < eps)
			return true;
	}
	return false;
}

void print_triangle(int trian_index) {
	Triangle trian = get_triangle(trian_index);
	cout << "--------- triangle " << trian_index << "-----------" << endl;
	cout << trian(0) << endl;
	cout << trian(1) << endl;
	cout << trian(2) << endl;
}

void test_mesh() {
	for(int i = 0; i < mesh.triang_num; i++) {
		if(triangle_has_point(i, Point(0.01 / 100, 0,  -3.2 / 100)) &&
				triangle_has_point(i, Point(-0.01 / 100, 0, -3.2 / 100))) {
			print_triangle(i);
		}
	}
}

