

#ifndef TYPE_H_
#define TYPE_H_


#include <complex>
#include <blitz/array.h>
using namespace blitz;

typedef blitz::TinyVector<double, 3> Point ;
typedef blitz::TinyVector<double, 3> Vector ;
typedef blitz::TinyVector<blitz::complex<double>, 3> Vector_c ;
typedef blitz::TinyVector<int, 3> TrianIndex ;
typedef blitz::TinyVector<int, 2> EdgeIndex ;
typedef blitz::TinyVector<Point, 3> Triangle;

typedef std::complex<double> Complex ;

typedef blitz::TinyVector<int, 3> Index3 ;
typedef blitz::TinyVector<int, 2> Index2 ;

//typedef mxArray *matArray;

extern const double epsilon_ ;
extern const double mu_ ;
extern const double c_ ;	// speed of light
extern const double eta_ ;  // free-space impedance
//extern const double epsilon_R ;
extern const double pi_ ;

extern const Complex j ;

extern const std::string GRAPH_DATA_FILE ;



#endif /* TYPE_H_ */
