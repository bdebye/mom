
#include <iostream>
#include <map>
#include <complex>

using namespace std;


#include "antenna.h"

#define HORN


double f = 35.22e9;
double omega = 2 * pi_ * f;

#ifdef HORN

int main()
{
	load_gmsh_file("horn.msh");
	test_mesh();
	vector<Feed> fd;
	fd.push_back(Feed(8843, 8871, Complex(1, 0)));
	//fd.push_back(Feed(5670, 1400, Complex(1, 0)));
	set_radiate_model(fd, omega);
    //print_feed_information();
	//print_shared_edge();
	solve_moment_equation();
	//cout << near_field_radius() << endl;

	directional_radiate_pattern_XZ();
    //special_output();

}

#endif
#ifdef DIPOLE

int main()
{
    load_gmsh_file("untitled.msh");
    test_mesh();
    vector<Feed> fd;
    fd.push_back(Feed(8, 11, Complex(1, 0)));
    set_radiate_model(fd, omega);
    //print_feed_information();
    //print_shared_edge();
    solve_moment_equation();
    //cout << near_field_radius() << endl;
    
    directional_radiate_pattern_XZ();
    special_output();
    
}
#endif
