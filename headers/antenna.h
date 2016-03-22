

#ifndef ANTENNA_H_
#define ANTENNA_H_

#include "type.h"
#include "rwg_edge.h"
#include "gmsh.h"
#include "dataio.h"
#include "mathlib.h"
#include "graphic.h"
#include "feed.h"

#include <map>

extern double omega;

void generate_edge_elements();

void set_radiate_model(const std::vector<Feed> &fd, double omega);

void set_scatter_model(Vector_c E_inc, double omega);

Complex get_input_impedance(int feed_edge);

Complex get_input_power(int feed_edge);

const RWGedge &get_RWGedge(int n);

int get_RWGedge_num();

const Array<complex<double>, 2> &get_impedance_matrix();

double get_omega();

double near_field_radius();

double U_density(const Point &r);

void save_impedance_matrix();

void print_shared_edge();

void directional_radiate_pattern_XZ();

void solve_moment_equation();

void directional_radiate_pattern_XY(int feed_edge);

void print_feed_information();

void special_output();

#endif /* ANTENNA_H_ */
