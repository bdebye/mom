
#include "clock.h"

//构造函数，对象产生即开始计时。
Clock::Clock() {
	this->start = clock();
}

//重置计时器，并马上开始新的计时。
void Clock::reset() {
	this->start = clock();
}

//计时器从开始计时到调用此方法所走过的秒数。
int Clock::second() {
	return (clock() - start) / CLOCKS_PER_SEC;
}

//计时器从开始计时到调用此方法所走过的分钟数。
int Clock::minute() {
	return second() / 60;
}

//计时器从开始计时到调用此方法所走过的小时数。
int Clock::hour() {
	return minute() / 60;
}

//计时器从开始计时到调用此方法所走过的毫秒数。
int Clock::millisecond() {
	return (clock() - start);
}

//计时器从开始计时到调用此方法所走过的时间，按通俗表示，分为时分秒各段。
time_ Clock::time_elapse() {
	time_ elp;
	elp.second = this->second() % 60;
	elp.minute = this->minute() % 60;
	elp.hour = this->hour() % 60;
	elp.millisecond = this->millisecond() % 1000;
	return elp;
}

//字符串格式的计时器从开始计时到调用此方法所走过的时间，按通俗表示，分为时分秒各段。
std::string Clock::time_elapse_format() {
	time_ elp = this->time_elapse();
	char format[100];
	sprintf(format, "%d hour %d min %d second %d", elp.hour, elp.minute, elp.second, elp.millisecond);
	return std::string(format);
}
