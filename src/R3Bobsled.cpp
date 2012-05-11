// Source file for the R3 bobsled class 



// Include files 

#include "R3/R3.h"
#include "R3Bobsled.h"



R3Bobsled::
R3Bobsled(void)
  : position(R3null_point),
    velocity(R3null_vector),
    mass(0)
{
  // Initialize material
	material = NULL;
  
  // Initialize mesh
  mesh = NULL;
}

int R3Bobsled::
Read(const char *filename)
{
  // do stuff?
  // read in mesh?
}

