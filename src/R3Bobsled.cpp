// Source file for the R3 bobsled class 



// Include files 

// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "R3Bobsled.h"
#include <cmath>
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#   define TRUE true
#   define FALSE false
#endif


////////////////////////////////////////////////////////////
// find force acting on bobsled
////////////////////////////////////////////////////////////

R3Vector Force(R3Bobsled *bobsled) {
    R3Vector force(R3null_vector);

    // force of gravity
    R3Vector fg(R3null_vector);
    R3Vector gravity(0, 0, -9.8);
    fg = bobsled->mass * gravity;
    
	// normal force
    R3Vector fn(fg);
	fn.Dot(bobsled->track->startNormal);
    if (bobsled->track->type != TRACK_STRAIGHT) {
        fn += bobsled->mass * bobsled->velocity * bobsled->velocity / bobsled->track->radius;
    }
    
	//force of friction
    R3Vector fk(R3null_vector);
    fk = bobsled->track->cof * fg.Length() * -1 * bobsled->velocity;
    
	// total force
    force = fg + fn + fk;
    return force;
}

void
UpdateBobsled(R3Scene *scene, double current_time, double delta_time,
			  bool force_left, bool force_right)
{
	return;
}