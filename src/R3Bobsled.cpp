// Source file for the R3 bobsled class 



// Include files 

#include "R3/R3.h"
#include "R3Scene.h"
#include "R3Bobsled.h"



R3Bobsled::
R3Bobsled(void)
  : position(R3null_point),
    velocity(R3null_vector),
    mass(0)
{
  node = NULL;
}

int R3Bobsled::
Read(const char *filename)
{
  // do stuff?
  // read in mesh?
}

////////////////////////////////////////////////////////////
// find force acting on bobsled
////////////////////////////////////////////////////////////

/*
R3Vector Force(R3Bobsled bobsled) {
    R3Vector force(R3null_vector);
    // force of gravity
    R3Vector fg(R3null_vector);
    R3Vector gravity(0, 0, -9.8);
    fg = bobsled.mass * gravity;
    // normal force
    R3Vector fn(R3null_vector);
    fn = fg.Dot(bobsled.track->normal);
    if (bobsled.track->isCurved) {
        fn += bobsled.mass * bobsled.velocity * bobsled.velocity / bobsled.track->R();
    }
    //force of friction
    R3Vector fk(R3null_vector);
    fk = bobsled.track->Cof() * fg.Length() * -1 * bobsled.velocity;
    // total force
    force = fg + fn + fk;
    return force;
}*/

void R3Bobsled::
UpdateBobsled(R3Node *node, double current_time, double delta_time,
              bool force_left, bool force_right)
{
 // fprintf(stderr, "UpdateBobsled unimplemented!\n");
}