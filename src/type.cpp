#include "type.h"

const double epsilon_ = 8.854e-012; //真空介电常数
const double mu_ = 1.257e-006; // 真空磁导率
const double c_ = 2.998e008;	// 光速
const double eta_ = 376.7887; // 自由空间阻抗
//const double epsilon_R = 1.0;
const double pi_ = 3.1415926535897931159979634685441851615906; //圆周率

const Complex j(0, 1); // 纯虚数单位

const std::string GRAGH_DATA_FILE = "graph"; //作图数据暂存位置

const int PLUS = 1;
const int MINUS = 0;
const int OTHER = -1;
