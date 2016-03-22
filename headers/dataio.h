/*
 * dataio.h
 *
 *  Created on: Sep 28, 2015
 *      Author: lambda
 */

#ifndef DATAIO_H_
#define DATAIO_H_


#include "type.h"


bool detect_file(const string filename, int& m, int& n);

template <typename T>
Array<T, 2> load(const string filename);

//template <typename T>
//void load(const string filename, Array<T, 2> & matrix);

//template <typename T>
//void save(const Array<T, 2> & data, string filename);


void save_csv(const Array<double, 2> & data, string filename);

void save_csv(const Array<double, 1> vect_base,
			  const Array<double, 1> vect_value,
			  const string filename);

void save(const Array<double, 1> vect_base,
			  const Array<double, 1> vect_value,
			  const string filename);

void save(const Array<double, 1> & data, const string filename);


#endif /* DATAIO_H_ */
