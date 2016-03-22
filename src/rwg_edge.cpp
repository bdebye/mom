
#include "rwg_edge.h"
#include "mathlib.h"
#include "gmsh.h"
//
//构造函数，初始化RWG边元数据。
//
//参数： index RWG边元的索引值。
//参数： edge 边元的公共边。
//参数： vertex1 边元中正三角形的自由顶点。
//参数： vertex2 边元中负三角形的自由顶点。
//参数： trian_index1 边元中正三角形的索引值。
//参数： trian_index2 边元中负三角形的索引值。
//参数： data 边元集合所在天线对象的索引，使单个边元也享有对所有数据的访问权。
//
RWGedge::RWGedge(int index,
			const EdgeIndex & edge,
			int vertex1,
			int vertex2,
			int trian_index1,
			int trian_index2
		):
	index(index),
	eg_share(edge),
	ve_plus(vertex1),
	ve_minus(vertex2),
	t_plus_index(trian_index1),
	t_minus_index(trian_index2) {

}

//得到边元中正三角形的自由顶点坐标。
Point RWGedge::vertex_plus() const {
	return mesh.node[ve_plus];
}

//得到边元中正三角形的自由顶点索引值。
int RWGedge::vertex_plus_index() const {
	return ve_plus;
}

//得到边元中负三角形的自由顶点坐标。
Point RWGedge::vertex_minus() const {
	return mesh.node[ve_minus];
}

//得到边元中负三角形的自由顶点索引值。
int RWGedge::vertex_minus_index() const {
	return ve_minus;
}

//得到边元中正三角形的面积。
double RWGedge::plus_area() const {
	return triangle_area(t_plus_index);
}

//得到边元中负三角形的面积。
double RWGedge::minus_area() const {
	return triangle_area(t_minus_index);
}

//得到边元中正三角形中心点的坐标。
Point RWGedge::plus_center() const {
	return trian_center(get_triangle(t_plus_index));
}

//得到边元中负三角形中心点的坐标。
Point RWGedge::minus_center() const {
	return trian_center(get_triangle(t_minus_index));
}

//得到边元公共边的索引，实际上是两个顶点的索引。
EdgeIndex RWGedge::edge_share() const {
	return eg_share;
}

//得到边元正三角形的索引值。
int RWGedge::trian_plus_index() const {
	return t_plus_index;
}

//得到边元负三角形的索引值。
int RWGedge::trian_minus_index() const {
	return t_minus_index;
}

//得到边元公共边的长度。
double RWGedge::share_edge_length() const {
	return distance(get_point(eg_share(0)), get_point(eg_share(1)));
}

const Vector RWGedge::rho_c_plus() const {
	return Vector(plus_center() - vertex_plus());
}

const Vector RWGedge::rho_c_minus() const {
	return Vector(vertex_minus() - minus_center());
}

const Vector RWGedge::rho_plus(const Point & point) const {
	return Vector(point - vertex_plus());
}

const Vector RWGedge::rho_minus(const Point & point) const {
	return Vector(vertex_minus() - point);
}


//边元正三角形的坐标数据
Triangle RWGedge::trian_plus() const {
	return get_triangle(this->trian_plus_index());
}

//边元负三角形的坐标数据
Triangle RWGedge::trian_minus() const {
	return get_triangle(this->trian_minus_index());
}

