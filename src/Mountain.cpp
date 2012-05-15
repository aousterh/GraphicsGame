/*
 * Mountain.cpp
 *
 *  Created on: May 11, 2012
 *      Author: Robert
 */

#include <cstdlib>
#include <cmath>
#include <iostream>

#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#   define TRUE true
#   define FALSE false
#endif

#include "Mountain.h"


////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double Mountain::RandomNumber(void)
{
#ifdef _WIN32
  // Seed random number generator
  static int first = 1;
  if (first) {
    srand(GetTickCount());
    first = 0;
  }

  // Return random number
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) RAND_MAX);
  return (r1 + r2) / ((double) RAND_MAX);
#else
  // Seed random number generator
  static int first = 1;
  if (first) {
    struct timeval timevalue;
    gettimeofday(&timevalue, 0);
    srand48(timevalue.tv_usec);
    first = 0;
  }

  // Return random number
  return drand48();
#endif
}

static const double h = .95;

Mountain::Mountain()
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			heights[i][j] = 0;
		}
	}

	heights[0][0] = 1;
	heights[0][height - 1] = 1;
	heights[width - 1][0] = 1;
	heights[width - 1][height - 1] = 1;

	for (int k = 0; k < 6; k++)
	{
		double range = 50.0;
		int radius = width > height ? height : width;
		radius--;
		while (radius > 1)
		{
			//diamond step
			for (int i = 0; i < width - 1; i += radius)
			{
				for (int j = 0; j < height - 1; j += radius)
				{
					double tl = heights[i][j];
					double tr = heights[i + radius][j];
					double bl = heights[i][j + radius];
					double br = heights[i + radius][j + radius];

					double cent = tl + tr + bl + br;
					cent /= 4.0;
					double r = 2.0 * RandomNumber() - 1.0;
					r *= range;
					cent += r;
					heights[i + radius / 2][j + radius / 2] = cent;
				}
			}

			//square step
			int colCount = 0;
			for (int i = 0; i < width; i += radius / 2)
			{
				for (int j = colCount % 2 == 0 ? radius / 2 : 0; j < height; j += radius)
				{
					int count = 0;
					double avg = 0;
					if (i >= radius / 2)
					{
						avg += heights[i - radius / 2][j];
						count++;
					}
					if (j >= radius / 2)
					{
						avg += heights[i][j - radius / 2];
						count++;
					}
					if (i + radius / 2 < width)
					{
						avg += heights[i + radius / 2][j];
						count++;
					}
					if (j + radius / 2 < height)
					{
						avg += heights[i][j + radius / 2];
						count++;
					}
					avg /= (double) count;
					double cent = avg;
					double r = 2.0 * RandomNumber() - 1.0;
					r *= range;
					cent += r;
					heights[i][j] = cent;
				}
				colCount++;
			}

			radius /= 2;
			range *= pow(2, -h);
		}
	}
    
    double min = 10000;
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (heights[i][j] < min)
            {
                min = heights[i][j];
            }
        }
    }
    
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            heights[i][j] += abs(min);
        }
    }
	printf("done with mountain\n");
}

