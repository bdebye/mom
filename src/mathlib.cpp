
#include "mathlib.h"

//
//用于产生等差数列，类似于MATLAB住的linspace函数，但不同的是第三个参数规定区间划分的段数
//也就是说实际产生的等差数列长度为n+1。
//参数： a 等差数列的起始元素。
//参数： b 等差数列的终止元素。
//参数： n 等差数列的间隔段数。
//返回： 产生的等差数列数组。
//
Array<double, 1> linspace(double a, double b, int n) {
	Array<double, 1> vect(n + 1);
	double h = (b - a) / n;
	double pos = a;
	for(int i = 0; i <= n; i++) {
		vect(i) = pos;
		pos += h;
	}
	return vect;
}

//
//提供矩阵的转置操作，得益于blitz函数库的数组索引的封装，转置过程中不发生矩阵存储形状的改变，
//而只改变数组索引的策略。
//参数： matrix 需要转置的矩阵。
//返回： 转置后的矩阵。
//
Array<double, 2> transpose(const Array<double, 2> & matrix) {
	return matrix.transpose(1, 0);
}

//
//将一维行向量转置为1乘n的矩阵。
//参数： vect 用于转置的行向量。
//返回：转置之后的列向量。
//
Array<double, 2> transpose(const Array<double, 1> & vect) {
	int leng = vect.extent(0);
	Array<double, 2> vvec(leng, 1);
	for(int i = 0; i < leng; i++)
		vvec(i, 0) = vect(i);
	return vvec;
}

//
//将一维行向量转置为1乘n的矩阵，此函数为针对复数向量的重载形式。
//参数： vect 用于转置的行向量。
//返回：转置之后的列向量。
//
Array<Complex, 2> transpose(const Array<Complex, 1> & vect) {
	int leng = vect.extent(0);
	Array<Complex, 2> vvec(leng, 1);
	for(int i = 0; i < leng; i++)
		vvec(i, 0) = vect(i);
	return vvec;
}

//
//求一空间向量的二范数。
//参数：vect 待求的空间三维向量。
//返回：空间向量的二范数。
//
double norm(const Vector & vect) {
	int n = vect.extent(0);
	double result = 0.0;
	for(int i = 0; i < n; i++)
		result += vect(i) * vect(i);
	return sqrt(result);
}

//
//给出三角形三个顶点的坐标数据则算出三角形的面积。
//参数：trian 三角形顶点数据。
//返回：三角形的面积。
//
double trian_area(const Triangle & trian) {

	Vector vect1(trian(0) - trian(1));
	Vector vect2(trian(2) - trian(1));
	return sqrt(norm(cross(vect1, vect2))) / 2;
}

//
//给出三角形三个顶点的坐标数据则算出三角形的中心。
//参数：trian 三角形顶点数据。
//返回：三角形的中心。
//
Point trian_center(const Triangle & trian) {
	return Point((trian(0) + trian(1) + trian(2)) / 3);
}

//
//给出三角形三个顶点的坐标数据则算出三角形的中心。
//参数：p1, p2, p3 三角形的三个顶点。
//返回：三角形的中心。
//
Point trian_center(const Point &p1, const Point &p2, const Point &p3) {
	return  Point((p1 + p2 + p3) / 3);
}

//
//求空间内两个点的距离。
//参数：point1， point2 空间两个点的坐标。
//返回：两个点之间的距离。
//
double distance(const Point & point1, const Point & point2) {
	return norm(Point(point1 - point2));
}

/*
void gauss_elimination_pivot(Array<double, 2> &extended, int n) {
    int max_i = 0;
    double max_r = 0;
    double temp;
    for(int i = 0; i < extended.extent(0); i++) {
        if(abs(extended(i, n)) > max_r) {
            max_r = abs(extended(i, n));
            max_i = i;
        }
    }
    for(int j = n; j < extended.extent(1); j++) {
        temp = extended(max_i, j);
        extended(max_i, j) = extended(n, j);
        extended(n, j) = temp;
    }
}
*/

//
//不带列主元的高斯消元法，主要用于解以阻抗矩阵为主体的矩量方程ZI=V。
//参数：I 表面电流因子。
//参数：Z 阻抗矩阵。
//参数：V 边元元馈电电压。
//
void gauss_elimination(Array<Complex, 1> & I, Array<Complex, 2> & Z, Array<Complex, 1> & V) {
	int m = Z.extent(0);
	int n = Z.extent(1) + 1;

	if(Z.extent(0) != Z.extent(1))
		cout << "Not a square matrix, linear system cannot be solved..." << endl;

	//Array<Complex, 2> extend(m, n);

	//extend(Range::all(), Range(0, n - 2)) = Z;
	//extend(Range::all(), Range(n - 1, n - 1)) = transpose(V);

	for(int i = 0; i < m; i++) {
        //gauss_elimination_pivot(extended, i);
		for(int j = i + 1; j < m; j++) {
			Complex  scala = Z(j, i) / Z(i, i);
			Z(Range(j, j), Range(i, m - 1))
				= Z(Range(j, j), Range(i, m - 1)) - scala * Z(Range(i, i), Range(i, m - 1));
            V(j) = V(j) - scala * V(i);
		}
	}

	// cout << extend << endl;

	for(int i = m - 1; i >= 0; i--) {
		Complex value = Complex(0);
		for(int j = i + 1; j < m; j++)
			value += I(j) * Z(i, j);
		I(i) = (V(i) - value) / Z(i, i);
	}
}

//求两点连线的中点。
Point edge_center(const Point & point1, const Point & point2) {
	return Point((point1 + point2) / 2);
}


//针对复数的自然指数运算。
Complex exp(Complex comp) {
	double real = comp.real();
	double image = comp.imag();
	return Complex(std::exp(real) * cos(image), std::exp(real) * sin(image));
}

//针对复向量的自然指数运算，返回为复向量。
Vector_c exp(Vector_c comp) {
	Vector_c result;
	for(int i = 0; i < 3; i++)
		result(i) = exp(comp(i));
	return result;
}

Complex operator/(const Complex & c1, const Complex & c2) {
	Complex divid = c1 * Complex(c2.real(), - c2.imag());
	return divid / (square(c2.real()) + square(c2.imag()));
}

Complex transv(const Vector_c & v1, const Vector & v2) {
	return (v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2));
}

Complex transv(const Vector_c & v1, const Vector_c & v2) {
	return (v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2));
}

//将实数向量拓展为复向量，虚步初始化为0。
Vector_c expand(const Vector &v) {
	Vector_c result;
	for(int i = 0; i < 3; i++) {
		result(i) = Complex(v(i), 0);
	}
	return result;
}

Vector_c conj(const Vector_c &vect_c) {
	Vector_c cj = vect_c;
	cj(0) = std::conj(cj(0));
	cj(1) = std::conj(cj(1));
	cj(2) = std::conj(cj(2));
	return cj;
}

Vector real(const Vector_c &vect_c) {
    Vector rl;
    rl(0) = std::real(vect_c(0));
    rl(1) = std::real(vect_c(1));
    rl(2) = std::real(vect_c(2));
    return rl;
}

