
#ifndef FEED_H_
#define FEED_H_

#include "type.h"

struct Feed {
	Feed(int trian_index_A, int trian_index_B, Complex V);

	int trian_index_A;
	int trian_index_B;
	Complex V;

	int rwg_index;
};

#endif /* FEED_H_ */
