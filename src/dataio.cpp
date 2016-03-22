
#include "dataio.h"

#include <fstream>

using namespace std;
//
//判断文本文件所记录的矩阵的行数和列数。
//参数： filename 文本文件名。
//参数： m 引用参数，用于返回矩阵行数。
//参数： n 引用参数，用于返回矩阵列数。
//调用： load， Antenna::load_point, Antenna::load_trian
//
bool detect_file(const string filename, int& m, int& n) {
	ifstream fin(filename);
    string line;
	int count = 1, wcount = 0;
	double value;

	if(fin.bad()) {
		cout << "cannot open file "
			 << filename << " ..." <<  endl;
		return false;
	}

	getline(fin, line);
	wcount = 0;
	stringstream lins(line);
	while(lins) {
		lins >> value;
		wcount++;
	}
	n = wcount - 1;
    while (getline(fin, line)) {
        if (line.length() > 0) {
            count++;
        }
    }
	m = count;
	fin.close();
	return true;
}

//
//载入文本文件的矩阵数据，矩阵元素在文本中必须规范排列，且不能为复数，这是其中的一个限制。
//
template <typename T>
Array<T, 2> load(const string filename) {
	int m, n;
	detect_file(filename, m, n);
	ifstream fin(filename);
	Array<T, 2> data(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++) {
			fin >> data(i, j);
		}
	fin.close();
	return data;
}

//template <typename T>
//void load(const string filename, Array<T, 2> & matrix) {
//	int m, n;
//	detect_file(filename, m, n);
//	ifstream fin(filename);
//	for(int i = 0; i < m; i++)
//		for(int j = 0; j < n; j++) {
//			fin >> matrix(i, j);
//		}
//	fin.close();
//	return data;
//}

//template <typename T>
//void pe_cemlab::save(const Array<T, 2> & data, string filename = "temp") {
//	int m = data.extent(0);
//	int n = data.extent(1);
//	ofstream fout(filename);
//	for(int i = 0; i < m; i++) {
//		for(int j = 0; j < n; j++)
//			fout << data(i, j) << "\t";
//		fout << endl;
//	}
//	fout.close();
//}

//template <>
//void pe_cemlab::save(const Array<Complex, 2> & data, string filename) {
//	int m = data.extent(0);
//	int n = data.extent(1);
//	ofstream fout(filename);
//	for(int i = 0; i < m; i++) {
//		for(int j = 0; j < n; j++)
//			fout << abs(data(i, j)) << "\t";
//		fout << endl;
//	}
//	fout.close();
//}

//
//将矩阵数据导出到csv格式文本文件中，每行以逗号分割存储矩阵的一行元素。
//参数： data 二维矩阵数据，元素不能是复数。
//参数： filename 要保存的文本文件名称。
//
void save_csv(const Array<double, 2> & data, string filename) {
	int m = data.extent(0);
	int n = data.extent(1);
	ofstream fout(filename);
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			fout << data(i, j);
			if(j != n - 1)
				fout << ",";
		}
		fout << endl;
	}
	fout.close();
}

//
//保存两列数据到csv文件中，每行有两个元素，以逗号分割。
//参数：vect_base 第一列数据。
//参数：vect_value 第二列数据。
//参数：filename 要保存的文件名。
//
void save_csv(const Array<double, 1> vect_base,
			  const Array<double, 1> vect_value,
			  const string filename = "temp") {

	int leng = vect_base.extent(0);
	if(leng != vect_value.extent(0)) {
		cout << "Error, lengths of two vectors are not same..." << endl;
		return;
	}
	ofstream fout(filename);
	for(int i = 0; i < leng; i++) {
		fout << vect_base(i) << ",\t"
			 << vect_value(i)
			 << endl;
	}
	fout.close();
}

//
//保存两列数据到csv文件中，每行有两个元素，以TAB分割。
//参数：vect_base 第一列数据。
//参数：vect_value 第二列数据。
//参数：filename 要保存的文件名。
//
void save(const Array<double, 1> vect_base,
			  const Array<double, 1> vect_value,
			  const string filename = "temp") {

	int leng = vect_base.extent(0);
	if(leng != vect_value.extent(0)) {
		cout << "Error, lengths of two vectors are not same..." << endl;
		return;
	}
	ofstream fout(filename);
	for(int i = 0; i < leng; i++) {
		fout << vect_base(i) << "\t"
			 << vect_value(i)
			 << endl;
	}
	fout.close();
}

//
//将矩阵数据导出到文本文件中，每行以TAB分割存储矩阵的一行元素。
//参数： data 二维矩阵数据，元素不能是复数。
//参数： filename 要保存的文本文件名称。
//
void save(const Array<double, 1> & data, const string filename) {
	int n = data.extent(0);
	ofstream fout(filename);
	for(int i = 0; i < n; i++) {
		fout << data(i) << endl;
	}
	fout.close();
}
