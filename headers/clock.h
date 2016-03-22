
#ifndef CLOCK_H_
#define CLOCK_H_

#include <string>

typedef struct {
	int hour;
	int minute;
	int second;
	int millisecond;
} time_;

class Clock {
public:
	Clock();

	void reset();

	int hour();
	int minute();
	int second();
	int millisecond();

	time_ time_elapse();
	std::string time_elapse_format();

private:
	int start;
};


#endif /* CLOCK_H_ */
