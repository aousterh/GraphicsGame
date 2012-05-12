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


double ANGLE_SHIFT = 0.016;


////////////////////////////////////////////////////////////
// Updating Bobsled
////////////////////////////////////////////////////////////

void UpdateBobsled(R3Scene *scene, double current_time, double delta_time, 
				   bool force_left, bool force_right)
{
    printf("delta_time = %f\n", delta_time);
	// update each sled in the scene
	for (int i = 0; i < scene->NBobsleds(); i++) {
		
		// get the current sled and its track segment
		R3Bobsled *bobsled = scene->Bobsled(i);
		R3Track *track = bobsled->track;

		// find the closest point on the along vector of the track
        R3Point position(bobsled->position);
		R3Vector temp(bobsled->position - track->start);
        temp.Project(track->along);
        R3Point center_point(track->start);
        center_point += temp;
		double r = R3Distance(center_point, bobsled->position);
        printf("r = %f", r);
		R3Vector force(R3null_vector);
		force = Force(bobsled, r);
		R3Vector velocity(R3null_vector);
		velocity = bobsled->velocity + delta_time * force/bobsled->mass;
		//velocity.Print();
        printf("\n");
		// Forward translation on a straight track
		R3Vector v_along(R3null_vector);
		if (track->type == TRACK_STRAIGHT) {
			v_along = bobsled->velocity.Dot(track->along) * track->along;
            //printf("\n");
            //v_along.Print();
            //printf("track along = ");
            //track->along.Print();
			bobsled->position.Translate(v_along);
            bobsled->sled->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->skates->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->helmets->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->masks->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
			bobsled->camera->eye += v_along;
		}
		
		// Side rotation on a straight track
		R3Vector v_side(R3null_vector);
        R3Vector v_down(R3null_vector);
        R3Vector down(0, -1, 0);
		R3Vector rotate_vector(R3null_vector);
		R3Line rotate_line(track->start, track->along);
        double sign;
        //printf("\n rotate_line = ");
        //rotate_line.Print();
		if (track->type == TRACK_STRAIGHT) {
            R3Vector temp(position - track->start);
            temp.Project(track->along);
            R3Point center_point(track->start);
            center_point += temp;
            center_point.Print();
            printf("\n");
            R3Vector normal(center_point - position);
            normal.Normalize();
            printf("normal = ");
            normal.Print();
            temp = track->along;
            temp.Cross(normal);
            temp.Normalize();
            sign = bobsled->velocity.Dot(temp);
			v_side = sign * temp;
			rotate_vector = track->along;
            /*v_side = bobsled->velocity.Dot(track->side) * track-side;
			rotate_vector = track->along;
            v_down = bobsled->velocity.Dot(down) * down;*/
			
		}
        v_side.Print();
        R3Point new_point(position + delta_time * v_side);
        R3Vector dist_vect(position - new_point);
        double delta_dist = dist_vect.Length();
        if (sign > 0)
            delta_dist *= -1;
        double delta_theta = delta_dist / r;
        //double delta_theta = delta_dist / track->radius;
        if (force_right)
            delta_theta += ANGLE_SHIFT;
        if (force_left)
            delta_theta -= ANGLE_SHIFT;
        bobsled->little_theta += delta_theta;
        printf("\n delta theta = %f", delta_theta);

        bobsled->position.Rotate(rotate_line, delta_theta);
		bobsled->sled->mesh->Rotate(delta_theta, rotate_line);
		bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
		bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
		bobsled->masks->mesh->Rotate(delta_theta, rotate_line);
        
        // check if over the edge
        if (bobsled->little_theta > M_PI/4) {
            // fall of edge or something
        }
        
        // check if bobsled is on new track
        R3Vector to_plane(track->endPlane.Point() - bobsled->position);
        double dist_plane = to_plane.Dot(track->endPlane.Normal());
        if (dist_plane <= 0) {
            bobsled->track = track->next;
            bobsled->big_theta = 0;
        }

		bobsled->velocity = velocity;
	}
}

////////////////////////////////////////////////////////////
// find force acting on bobsled
////////////////////////////////////////////////////////////

R3Vector Force(R3Bobsled *bobsled, double r) {
    R3Vector force(R3null_vector);
    R3Track *track = bobsled->track;
    // force of gravity
    R3Vector fg(R3null_vector);
    R3Vector gravity(0, -9.8, 0);
    fg = bobsled->mass * gravity;
    // normal force
    R3Vector fn(R3null_vector);
    if (track->type == TRACK_STRAIGHT) {
        R3Vector temp(bobsled->position - track->start);
        temp.Project(track->along);
        R3Point center_point(track->start);
        center_point += temp;
        R3Vector normal(center_point - bobsled->position);
        normal.Normalize();
        double dot = fg.Dot(normal);
        printf("\ndot = %f\n", dot);
        R3Vector centripetal(R3null_vector);
        centripetal = bobsled->mass * bobsled->velocity * bobsled->velocity / r;
        centripetal = centripetal.Length() * normal;
        fn = fg.Dot(normal) * normal + centripetal;
    }
    // force of friction
    R3Vector fk(R3null_vector);
    fk = track->cof * fn.Length() * -1 * bobsled->velocity;
    
    // total force
    force = fg + fn + fk;
    //force = fk + fn;
    return force;
}
    
