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

<<<<<<< HEAD
enum angle_shift 0.016

R3Bobsled::
R3Bobsled(void)
  : position(R3null_point),
    velocity(R3null_vector),
    mass(0),
    main_theta(0),
    inside_theta(0)
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
=======
>>>>>>> 1326cd87873ec13084d934af14a7f58bf700b7f2

////////////////////////////////////////////////////////////
// Updating Bobsled
////////////////////////////////////////////////////////////

void UpdateBobsled(R3Scene scene, R3Node *node, double current_time, double delta_time, bool force_left, bool force_right) {
  /*  R3Bobsled bobsled(this);
    double half_height = .5 * node->bbox.YLength();
    //R3Vector along(this.track->rotate_vector);
    R3Vector force(R3null_vector);
    force = Force(bobsled, half_height);
    R3Vector velocity(R3null_vector);
    velocity = this.velocity + delta_time * force/this.mass;
    // do forward translation
    R3Vector v_along(R3null_vector);
    if(this.track->isCurved) {
        v_along = this.velocity.Dot(this.track->along) * this.track->along;
        double delta_dist = delta_time * v_along;
        double delta_theta = delta_dist / (this.track->R() + half_height);
        this.position.Rotate(this.track->rotate_line, delta_theta);
    }
    else {
        v_along = this.velocity.Dot(this.track->along) * this.track->along;
        this.position.Tranlate(v_along);
    }
                                         
    // do side to side translation
    R3Vector v_side(R3null_vector);
    R3Vector rotate_vector;
    if (this.track->isCurved) {
    }
    else {
        v_side = this.velocity.Dot(this.track->side) * this.track->side;
        rotate_vector = this.track->along;
    }
    double delta_dist = delta_time * v_side;
    double delta_theta = delta_dist / (this.track->R() + half_height);
    if (force_left)
        delta_theta -= angle_shift;
    if (force_right)
        delta_theta += angle_shift;
    this.position.Rotate(rotate_vector, delta_theta);

    
    // check if still on track
    R3Vector to_plane(this.track->end_plane.Point() - this.position);
    double dist_plane = to_plane.Dot(this.track->end_plane.Vector());
    if (dist_plane <= 0) {
        this.track = this.track->next;
    }*/
}

////////////////////////////////////////////////////////////
// find force acting on bobsled
////////////////////////////////////////////////////////////
/*
R3Vector Force(R3Bobsled bobsled, double half_height) {
    R3Vector force(R3null_vector);
    // force of gravity
    R3Vector fg(R3null_vector);
    R3Vector gravity(0, 0, -9.8);
    fg = bobsled.mass * gravity;
    // normal force
    R3Vector fn(R3null_vector);
    fn = fg.Dot(bobsled.track->normal);
    if (bobsled.track->isCurved) {
        fn += bobsled.mass * bobsled.velocity * bobsled.velocity / (bobsled.track->R() + half_height);
>>>>>>> 2bc1fea96058d704e9d707c1cc18bf5fb81d8604
    }
    
	//force of friction
    R3Vector fk(R3null_vector);
<<<<<<< HEAD
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
=======
    fk = bobsled.track->cof * fg.Length() * -1 * bobsled.velocity;
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