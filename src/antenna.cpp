
#include "antenna.h"
#include "type.h"
#include "mathlib.h"
#include "gmsh.h"
#include "rwg_edge.h"
#include "clock.h"
#include "graphic.h"

#include <map>
#include <vector>
#include <fstream>
#include <complex>
#include <algorithm>

//Array<int, 2> trian_data;

//空间点坐标数据，为3乘n的矩阵，每列数据为空间一个点的坐标
//sArray<double, 2> point_data;

//电流扩展系数矢量，由矩量方程ZI=V解出
Array<Complex, 1> Im;

//用于在三角形上进行面积分的小三角形中心
Array<Point, 2> integral_center;

//阻抗矩阵
Array<Complex, 2> impedance;

//电压激励向量
Array<Complex, 1> Vm;

//天线工作在散射模式时的电压激励向量
Vector_c E_inc;

//天线工作在辐射模式时的边元馈电策略
vector<Feed> feed;

//天线结构内部所有的RWG边元
vector<RWGedge> edge_element;

//三角形在哪些RWG边元中是正三角形
vector<vector<int>> trianP_toRWG;

//三角形在哪些RWG边元中是负三角形
vector<vector<int>> trianM_toRWG;

//天线网格结构三角形总数
int trian_num;

//天线网格结构三角形顶点总数
int point_num;

//RWG边元个数
int rwg_num;

double omega_;
//天线工作馈电的相速

const RWGedge& get_RWGedge(int n) {
	return edge_element.at(n);
}

int get_RWGedge_num() {
	return edge_element.size();
}

int find_feed_edge(int trian_index_A, int trian_index_B) {
	for(int i = 0; i < feed.size(); i++) {
		if(trian_index_A == feed[i].trian_index_A && trian_index_B == feed[i].trian_index_B)
			return i;
		if(trian_index_A == feed[i].trian_index_B && trian_index_B == feed[i].trian_index_A)
			return i;
	}
	return -1;
}

bool share_edge_impl(int m, int n, Index2 &edge, int &vertex1, int &vertex2) {

	Index3 trian1 = get_triangle_index(m);
	Index3 trian2 = get_triangle_index(n);

	int count = 0;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++) {
			if(trian1(i) == trian2(j)) {
				edge(count) = trian1(i);
				count++;
			}
		}

	if(count < 2)
		return false;

	for(int i = 0; i < 3; i++) {
		if(trian1(i) != edge(0) && trian1(i) != edge(1))
			vertex1 = trian1(i);
		if(trian2(i) != edge(0) && trian2(i) != edge(1))
			vertex2 = trian2(i);
	}
	//cout << "--------------------------------" << endl;
	//cout << m << " " << n << endl;
	//cout << trian1 << " " << trian2 << endl;
	//cout << edge << endl;
	return true;
}


bool edge_feed_considered(int m, int n, Index2 edge_share) {
	if(find_feed_edge(m, n) >= 0)
		return false;
	for(auto iter = feed.begin(); iter != feed.end(); ++iter) {
		Index2 share1 = edge_share, share2;
		int v1, v2;
		share_edge_impl(iter->trian_index_A, iter->trian_index_B, share2, v1, v2);
		if((share1(0) == share2(0) && share1(1) == share2(1)) || (share1(0) = share2(1) && share1(1) == share2(0)))
			return true;
	}
	return false;
}

bool share_edge(int m, int n, Index2 &edge, int &vertex1, int &vertex2) {

	Index3 trian1 = get_triangle_index(m);
	Index3 trian2 = get_triangle_index(n);

	int count = 0;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++) {
			if(trian1(i) == trian2(j)) {
				edge(count) = trian1(i);
				count++;
			}
		}

	if(count < 2)
		return false;
	if(edge_feed_considered(m, n, edge))
		return false;

	for(int i = 0; i < 3; i++) {
		if(trian1(i) != edge(0) && trian1(i) != edge(1))
			vertex1 = trian1(i);
		if(trian2(i) != edge(0) && trian2(i) != edge(1))
			vertex2 = trian2(i);
	}
	//cout << "--------------------------------" << endl;
	//cout << m << " " << n << endl;
	//cout << trian1 << " " << trian2 << endl;
	//cout << edge << endl;
	return true;
}



void search_RWGedge() {
	trian_num = mesh.triang_num;
	trianP_toRWG.resize(trian_num + 1);
	trianM_toRWG.resize(trian_num + 1);
	cout << "generating the RWG edge elements..." << endl;
	for(int i = 1; i <= trian_num; i++)
		for(int j = i + 1; j <= trian_num; j++) {
			EdgeIndex edge;
			int vertex1, vertex2;
			if(share_edge(i, j, edge, vertex1, vertex2)) {
				edge_element.push_back(RWGedge(edge_element.size(), edge, vertex1, vertex2, i, j));
				trianP_toRWG.at(i).push_back(edge_element.size() - 1);
				trianM_toRWG.at(j).push_back(edge_element.size() - 1);
				if(find_feed_edge(i, j) >= 0) {
					feed[find_feed_edge(i, j)].rwg_index = edge_element.size() - 1;
				}
			}
		}

	rwg_num = edge_element.size();
	cout << "finished gerenation, " << get_RWGedge_num()
		<< " edge elements found..." << endl
		<< endl;
}


void generate_impedance_matrix(double omega) {
	//预先计算一些常量
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double k = omega / c_;
	double Constant1 = mu_ / (4 * pi_);
	double Factor = 1.0 / 9;

	Clock local_clock; //启动计时器
	Complex K = k * j;
	Complex Constant2 = 1 / (j * 4.0 * pi_ * omega * epsilon_);
	Complex FactorA = Factor * (j * omega / 4) * Constant1;
	Complex FactorFi = Factor * Constant2;

	omega_ = omega;
	int dim = get_RWGedge_num();
	impedance.resize(dim);

	//cout << Factor << endl;
	//cout << (j * omega / 4) << endl;
	//cout << Constant1 << endl;
	//cout << Constant2 << endl;
	//cout << k << endl;
	//cout << FactorA << endl;
	//cout << "-------------------------------------------------------------" << endl;

	//初始化阻抗矩阵
	for(int i = 0; i < rwg_num; i++)
		for(int j = 0; j < rwg_num; j++) {
			impedance(i, j) = Complex(0, 0);
		}

	for(int p = 1; p <= trian_num; p++) {
		vector<int> Plus = trianP_toRWG.at(p);
		vector<int> Minus = trianM_toRWG.at(p);
		Point t_center = get_triangle_center(p);
		if(p) {
			printf("generating impedance matrix, finishing %d%%\n",
				(int)((double)(p) / trian_num * 100));
		}

		for(int n = 0; n < rwg_num; n++) {
			const RWGedge &edge_n = get_RWGedge(n);
			TinyVector<Complex, 9> gP;
			TinyVector<Complex, 9> gM;
			Complex ZF;
			Complex Fi;
			int trian_p = edge_n.trian_plus_index();
			int trian_m = edge_n.trian_minus_index();

			for(int h = 0; h < 9; h++) {
				double R = norm(Vector(integral_center(trian_p, h) - t_center));
				gP(h) = exp(-K * R) / R;
				R = norm(Vector(integral_center(trian_m, h) - t_center));
				gM(h) = exp(-K * R) / R;
			}
			Fi = sum(gP) - sum(gM);
			ZF = FactorFi * edge_n.share_edge_length() * Fi;

			for(auto iter = Plus.begin(); iter != Plus.end(); ++iter) {
				int m = *iter;
				const RWGedge &edge_m = get_RWGedge(m);
				int trian_index = edge_m.trian_plus_index();
				Complex A(0, 0);
				for(int h = 0; h < 9; h++) {
					A += gP(h) * transv(edge_n.rho_c_plus(), edge_m.rho_plus(integral_center(trian_index, h))) +
						gM(h) * transv(edge_n.rho_c_minus(), edge_m.rho_plus(integral_center(trian_index, h)));
				}
				Complex Zl = FactorA * A * edge_n.share_edge_length();

				//cout << Zl << endl;
				impedance(m, n) += edge_m.share_edge_length() * (Zl + ZF);
			}

			for(auto iter = Minus.begin(); iter != Minus.end(); ++iter) {
				int m = *iter;
				const RWGedge &edge_m = get_RWGedge(m);
				int trian_index = edge_m.trian_minus_index();
				Complex A(0, 0);
				for(int h = 0; h < 9; h++) {
					A += gP(h) * transv(edge_n.rho_c_plus(), edge_m.rho_minus(integral_center(trian_index, h))) +
						gM(h) * transv(edge_n.rho_c_minus(), edge_m.rho_minus(integral_center(trian_index, h)));
				}
				Complex Zl = FactorA * A * edge_n.share_edge_length();
				impedance(m, n) += edge_m.share_edge_length() * (Zl - ZF);
			}
		}
	}
	cout << endl;
	cout << "finished generation, time elapsed: " << local_clock.time_elapse_format() << endl;
	cout << endl;
}

const Array<complex<double>, 2> & get_impedance_matrix() {
	return impedance;
}

Complex scatter_Vm_value(int m) {
	RWGedge edge_m = get_RWGedge(m);
	Vector rmc_p = edge_m.rho_c_plus();
	Vector rmc_m = edge_m.rho_c_minus();
	double lm    = edge_m.share_edge_length();
	return lm * (transv(E_inc, rmc_p) / 2 + transv(E_inc, rmc_m) / 2);
}

void generate_integral_pieces() {
	Triangle trian;
	TinyVector<Point, 10> pseq;
	Array<Point, 2> &cseq = integral_center;
	cseq.resize(trian_num + 1, 9);
	for(int i = 1; i <= trian_num; i++) {
		trian = get_triangle(i);
		pseq(0) = trian(0);
		pseq(1) = trian(0) * 2 / 3 + trian(1) / 3;
		pseq(2) = trian(0) * 2 / 3 + trian(2) / 3;
		pseq(3) = trian(0) / 3 + trian(1) * 2 / 3;
		pseq(5) = trian(0) / 3 + trian(2) * 2 / 3;
		pseq(4) = pseq(3) / 2 + pseq(5) / 2;
		pseq(6) = trian(1);
		pseq(7) = trian(1) * 2 / 3 + trian(2) / 3;
		pseq(8) = trian(1) / 3 + trian(2) * 2 / 3;
		pseq(9) = trian(2);

		cseq(i, 0) = trian_center(pseq(0), pseq(1), pseq(2));
		cseq(i, 1) = trian_center(pseq(1), pseq(3), pseq(4));
		cseq(i, 2) = trian_center(pseq(1), pseq(2), pseq(4));
		cseq(i, 3) = trian_center(pseq(2), pseq(4), pseq(5));
		cseq(i, 4) = trian_center(pseq(3), pseq(6), pseq(7));
		cseq(i, 5) = trian_center(pseq(3), pseq(4), pseq(7));
		cseq(i, 6) = trian_center(pseq(4), pseq(7), pseq(8));
		cseq(i, 7) = trian_center(pseq(4), pseq(5), pseq(8));
		cseq(i, 8) = trian_center(pseq(5), pseq(8), pseq(9));
	}
}


void generate_edge_elements() {
	search_RWGedge();
	generate_integral_pieces();
}

void set_scatter_model(Vector_c Ei, double omega) {
	feed.clear();
	omega_ = omega;
	cout << "using omega: " << omega_ << endl;
	generate_edge_elements();
	generate_impedance_matrix(omega);
	E_inc = Ei;
	Vm.resize(rwg_num);
	for(int i = 0; i < rwg_num; i++) {
		Vm(i) = scatter_Vm_value(i);
	}
	gauss_elimination(Im, impedance, Vm);
}

void solve_moment_equation() {
	generate_impedance_matrix(omega_);
	Im.resize(rwg_num);
	cout << "Solving the linear system..." << endl;
	gauss_elimination(Im, impedance, Vm);
	cout << "Finished solving the linear system..." << endl;
}

void set_radiate_model(const std::vector<Feed> &fd, double omega) {
	E_inc = Complex(0, 0);
	omega_ = omega;
	feed = fd;
	cout << "using omega: " << omega_ << endl;
	generate_edge_elements();
	Vm.resize(rwg_num);
	for(int i = 0; i < rwg_num; i++) {
		Vm(i) = 0;
	}
	for(vector<Feed>::iterator iter = feed.begin(); iter != feed.end(); ++iter)
		Vm(iter->rwg_index) = iter->V * get_RWGedge(iter->rwg_index).share_edge_length();

}

void print_shared_edge() {
	for(unsigned i = 0; i < edge_element.size(); i++) {
		Index2 node_index = edge_element.at(i).edge_share();
		if(get_point(node_index(0))(1) == 0.01 / 100.0) {
			cout << i << ": " << get_point(node_index(0)) << "----" << get_point(node_index(1)) << endl;
		}
	}
}

Complex get_input_impedance(int feed_edge) {
	bool is_feed = false;
	if(feed.empty()) {
		cout << "Sorry, not working in radiate mode... no feed edge specified..." << endl;
		return Complex(0, 0);
	}
	for(auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if(iter->rwg_index == feed_edge)
			is_feed = true;
	}
	if(!is_feed) {
		cout << "Soory, the RWG edge specified is not an feeding port..." << endl;
	}
	double ln = get_RWGedge(feed_edge).share_edge_length();
	return Vm(feed_edge) / (pow(ln, 2) * Im(feed_edge));
}

Complex get_input_power(int feed_edge) {
	bool is_feed = false;
	if(feed.empty()) {
		cout << "Sorry, not working in radiate mode... no feed edge specified..." << endl;
		return Complex(0, 0);
	}
	for(auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if(iter->rwg_index == feed_edge)
			is_feed = true;
	}
	if(!is_feed) {
		cout << "Soory, the RWG edge specified is not an feeding port..." << endl;
	}
	Complex V = 0.0;
	for(auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if(iter->rwg_index == feed_edge)
			V = iter->V;
	}
	return pow(V, 2) / get_input_impedance(feed_edge) / 2;
}

double near_field_radius() {
	double lambda = c_ / omega_ * 2 * pi_;
	return 2 * pow(geometry_size(), 2) /  lambda;
}

Vector_c RWG_H_strength(const RWGedge &edge, const Point &r) {
	Vector_c result;
	Point r_ = r - (edge.minus_center() + edge.plus_center()) / 2;
	Vector_c m = edge.share_edge_length() * Im(edge.index) * (edge.minus_center() - edge.plus_center());
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double rd = norm(r_);
	double k = omega_ / c_;
	Complex C = 1/ pow(rd, 2) * (Complex(1, 0) + 1 / (j * k * rd));
	result = (j * k) / (4 * pi_) * cross(m, expand(r_)) * C * exp(-j * k * rd);

	return result;
}

Vector_c RWG_E_strength(const RWGedge &edge, const Point &r) {
	Vector_c result;
	Point r_ = r - (edge.minus_center() + edge.plus_center()) / 2;
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double rd = norm(r_);
	double k = omega_ / c_;
	Complex C = 1/ pow(rd, 2) * (Complex(1, 0) + 1 / (j * k * rd));
	Vector_c m = edge.share_edge_length() * Im(edge.index) * (edge.minus_center() - edge.plus_center());
	Vector_c M = transv(r_, m) * r_ / (rd * rd);
	result = eta_ / (4.0 * pi_) * ((M - m) * (j * k / rd + C) + 2.0 * C * M) * exp(-j * k * rd);

	return result;
}

void save_impedance_matrix() {
	const Array<Complex, 2> &data = get_impedance_matrix();
	int m = data.extent(0);
	int n = data.extent(1);
	ofstream fout("Z.txt");
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			fout << std::abs(data(i, j));
			if(j != n - 1)
				fout << ",\t";
		}
		fout << endl;
	}
	fout.close();
}

Vector_c E_strength(const Point &r) {
	Vector_c E(Complex(0, 0), Complex(0, 0), Complex(0, 0));
	for(int i = 0; i < rwg_num; i++) {
		const RWGedge &edge = get_RWGedge(i);
		E = E + RWG_E_strength(edge, r);
	}
	return E;
}

Vector_c H_strength(const Point &r) {
	Vector_c H(Complex(0, 0), Complex(0, 0), Complex(0, 0));
	for(int i = 0; i < rwg_num; i++) {
		const RWGedge &edge = get_RWGedge(i);
		H = H + RWG_H_strength(edge, r);
	}
	return H;
}

Vector W_strength(const Point &r) {
	Vector result;
	Vector_c temp = cross(E_strength(r), conj(H_strength(r))) / 2.0;
	for(int i = 0; i < 3; i++) {
		result(i) = temp(i).real();
	}
	return result;
}

double U_density(const Point &r) {
	return norm(W_strength(r)) * pow(norm(r), 2);
}

void directional_radiate_pattern_XZ() {
    int feed_edge = feed[0].rwg_index;
	double R =  near_field_radius() * 100;
	Array<double, 1> X(1001), Y(1001);
	X = linspace(-0.5 * pi_, 1.5 * pi_, 1000);
	transform(X.begin(), X.end(), Y.begin(), [R] (double rad) -> double {
		return U_density(Point(R * cos(rad), 0, R * sin(rad)));
	});

	double scalar = get_input_power(feed_edge).real() / (4 * pi_);
	for(int i = 0; i < 1001; i++) {
		Y(i) = 10 * log10(Y(i) / scalar) + 40;
	}
	polar(X, Y);
}

void directional_radiate_pattern_XY(int feed_edge) {
	double R =  near_field_radius() * 10;
	Array<double, 1> X(1001), Y(1001);
	X = linspace(0, 2 * pi_, 1000);
	transform(X.begin(), X.end(), Y.begin(), [R] (double rad) -> double {
		return U_density(Point(R * cos(rad), R * sin(rad), 0));
	});

	double scalar = get_input_power(feed_edge).real() / (4 * pi_);
	for(int i = 0; i < 1001; i++) {
		Y(i) = 10 * log10(Y(i) / scalar);
	}
	polar(X, Y);
}

void print_feed_information() {
	for(auto iter = feed.begin(); iter != feed.end(); ++iter) {
		cout << iter->rwg_index << endl;
	}
}

/*
void special_output() {
    int nx = 100, ny = 100;
    double dx = 0.05 / nx, dy = 0.05 / ny;
    double wave_len = c_ / omega_ * 2 * pi_;
    double z_coord = 8.3 / 100 + wave_len / 2;
    Array<double, 2> field(nx, ny);
    cout << "Save field distribution..." << endl;
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++) {
            double x = -2.5 / 100 + i * dx;
            double y = -2.5 / 100 + j * dy;
            double z = z_coord;
            double Ex = E_strength(Point(x, y, z))(0).real();
            double Ey = E_strength(Point(x, y, z))(1).real();
            field(i, j) =  Ex * Ex + Ey * Ey;// norm(real(E_strength(Point(x, y, z))));
        }
    save_csv(field, "field.txt");
}
 */

void special_output() {
    int nx = 100, nz = 240;
    double dx = 0.05 / nx, dz = 0.12 / nz;
    double wave_len = c_ / omega_ * 2 * pi_;
    double y_coord = 0;
    Array<double, 2> field(nx, nz);
    cout << "Save field distribution..." << endl;
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < nz; j++) {
            double x = -2.5 / 100 + i * dx;
            double y = 0;
            double z = -3.2 / 100 + j * dz;
            double Ex = E_strength(Point(x, y, z))(0).real();
            double Ey = E_strength(Point(x, y, z))(1).real();
            field(i, j) =  Ex * Ex + Ey * Ey;// norm(real(E_strength(Point(x, y, z))));
        }
    save_csv(field, "field.txt");
}
