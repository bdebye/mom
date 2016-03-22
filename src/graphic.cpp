
#include "graphic.h"
#include "dataio.h"
//
//用于画天线的方向性图，将数据保存在文本文件中后使用python脚本读取。
//参数： X 极坐标作图中的角参量。
//参数： Y 极坐标作图中的径参量。
//
void polar(const Array<double, 1> &X, const Array<double, 1> &Y) {
	save(X, Y, "graph");
	system("python pplot.py");
}

//
//用于一般性二维作图，将数据保存在文本文件中后使用python脚本读取。
//参数： X 自变量数据。
//参数： Y 因变量数据。
//
void plot(const Array<double, 1> &X, const Array<double, 1> &Y) {
	save(X, Y, "graph");
	system("python plot.py");
}
