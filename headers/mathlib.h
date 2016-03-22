

#ifndef MATHLIB_H_
#define MATHLIB_H_

#include <blitz/array.h>
#include "type.h"
using namespace blitz;

Array<double, 1> linspace(double a, double b, int n);

Array<double, 2> transpose(const Array<double, 2> & matrix);

Array<double, 2> transpose(const Array<double, 1> & vect);
Array<Complex, 2> transpose(const Array<Complex, 1> & vect);

double norm(const Vector & vect);

double trian_area(const Triangle & trian);

Point trian_center(const Triangle & trian);
Point trian_center(const Point &p1, const Point &p2, const Point &p3);

double distance(const Point & point1, const Point & point2);

void gauss_elimination(Array<Complex, 1> & I, Array<Complex, 2> & Z, Array<Complex, 1> & V);

Point edge_center(const Point & point1,
				   const Point & point2);



Vector_c exp(Vector_c comp);

Complex operator/(const Complex & c1, const Complex & c2);

inline double square(double v) {
	return v * v;
}

double abs(Complex v);

Complex transv(const Vector_c & v1, const Vector & v2);
Complex transv(const Vector_c & v1, const Vector_c & v2);

Vector_c expand(const Vector &v);

Vector_c conj(const Vector_c &vect_c);

Vector real(const Vector_c &vect_c);


#endif /* MATHLIB_H_ */
