/*
 * Mountain.h
 *
 *  Created on: May 11, 2012
 *      Author: Robert
 */

#include <ctime>

#ifndef MOUNTAIN_H_
#define MOUNTAIN_H_

class Mountain
{
private:
	double RandomNumber();

public:
	static const int height = 129;
	static const int width = 513;
	Mountain();
	double heights[width][height];
};

#endif /* MOUNTAIN_H_ */
