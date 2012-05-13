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


double ANGLE_SHIFT = 5;



////////////////////////////////////////////////////////////
// Updating Bobsled
////////////////////////////////////////////////////////////

void UpdateBobsled(R3Scene *scene, double current_time, double delta_time, 
				   bool force_left, bool force_right)
{
    //printf("delta_time = %f\n", delta_time);
	// update each sled in the scene
	for (int i = 0; i < 1/*scene->NBobsleds()*/; i++) {
		// get the current sled and its track segment
		R3Bobsled *bobsled = scene->Bobsled(i);
		R3Track *track = bobsled->track;
        double big_r;
        R3Vector new_along(track->along);
        
        //printf("z0: %f\n", bobsled->velocity.Z());

		// find the closest point on the along vector of the track
        R3Point position(bobsled->position);
		R3Vector ve_along(bobsled->position - track->start);
        ve_along.Project(track->along);
        R3Point center_point(track->start);
        center_point += ve_along;
		double r = R3Distance(center_point, bobsled->position);
		R3Vector force(R3null_vector);
		force = Force(bobsled, r);
		R3Vector velocity(R3null_vector);
		velocity = bobsled->velocity + delta_time * force/bobsled->mass;
        
		// Forward translation on a straight track
		R3Vector v_along(R3null_vector);
		if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT) {
			v_along = bobsled->velocity.Dot(track->along) * track->along;
			bobsled->position.Translate(v_along);
            bobsled->sled->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->skates->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->helmets->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
            bobsled->masks->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
			bobsled->camera->eye += v_along;
		}
        
        // Forward translation on a curved track
        if (track->type == TRACK_TURN_LEFT || track->type == TRACK_TURN_RIGHT) {
            R3Vector vel_along = (bobsled->velocity.Dot(track->along)) * track->along;
            R3Point new_point(position + delta_time * vel_along);
            R3Vector dist_vect(new_point - position);
            double dist =  dist_vect.Length();
            R3Vector vect_radius(track->center_point - position);
            double delta_theta = 100 * dist/vect_radius.Length();
            printf("delta_theta = %f\n", delta_theta);
            double percent = delta_theta/(M_PI/2);
            bobsled->big_percent += percent;
            //dist_to_start.Project(bobsled->track->startNormal);
            //double percent = dist/(2*M_PI*bobsled->big_r);
            //new_along = bobsled->big_percent * track->endPlane->Normal() + (1 - bobsled->big_percent) * track->along;
            
            //new_along.Normalize();
            R3Line rotate_line(track->center_pivot);
            new_along.Rotate(rotate_line.Vector(), delta_theta);
            bobsled->position.Rotate(rotate_line, delta_theta);
            bobsled->sled->mesh->Rotate(delta_theta, rotate_line);
            bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
            bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
            bobsled->masks->mesh->Rotate(delta_theta, rotate_line);
			bobsled->camera->eye.Rotate(rotate_line, delta_theta);
			bobsled->camera->right.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->up.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->towards.Rotate(rotate_line.Vector(), delta_theta);
        }
		
		// Side rotation on a straight track
		R3Vector v_side(R3null_vector);
        R3Vector v_down(R3null_vector);
        R3Vector down(0, -1, 0);
		R3Vector rotate_vector(R3null_vector);
		R3Line rotate_line(track->start, track->along);
        double sign;
        R3Vector temp(R3null_vector);
		if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT) {
            temp = position - track->start;
            temp.Project(track->along);
            R3Point center_point(track->start);
            center_point += temp;
            R3Vector normal(center_point - position);
            normal.Normalize();
            temp = track->along;
            temp.Cross(normal);
            temp.Normalize();
            sign = bobsled->velocity.Dot(temp);
			v_side = sign * temp;
			rotate_vector = track->along;
			
		

            R3Point new_point(position + delta_time * v_side);
            R3Vector dist_vect(position - new_point);
            //R3Vector normal(center_point - position);
            //normal.Normalize();
            temp = R3null_vector;
            temp = track->along;
            temp.Cross(normal);
            temp.Normalize();
            temp.SetZ(0);
                double delta_dist = dist_vect.Length();
            if (sign > 0)
                delta_dist *= -1;
            double delta_theta = delta_dist / r;
            if (force_right) {
                double v_change = (ANGLE_SHIFT * r) * delta_time;
                //printf("z change: %f\n", temp.Z());
                velocity += v_change * temp;
                //delta_theta -= ANGLE_SHIFT;
            }
            if (force_left) {
                double v_change = (ANGLE_SHIFT * r) * delta_time;
                velocity -= v_change * temp;
                //delta_theta += ANGLE_SHIFT;
            }
            bobsled->little_theta += delta_theta;

            bobsled->position.Rotate(rotate_line, delta_theta);
            bobsled->sled->mesh->Rotate(delta_theta, rotate_line);
            bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
            bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
            bobsled->masks->mesh->Rotate(delta_theta, rotate_line);
			bobsled->camera->eye.Rotate(rotate_line, delta_theta);
			bobsled->camera->right.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->up.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->towards.Rotate(rotate_line.Vector(), delta_theta);
        }
        
        // check if over the edge
        if (bobsled->little_theta > M_PI/4) {
            // fall of edge or something
        }
        
        track->along = new_along;
        // check if bobsled is on new track
        R3Vector to_plane(track->end - bobsled->position);
        double dist_plane = to_plane.Dot(track->endPlane.Normal());
        if (dist_plane <= 0) {
            bobsled->track = track->next;
            printf("went to new track %d\n", track->type);
            
            if (bobsled->track->type == TRACK_TURN_RIGHT || bobsled->track->type == TRACK_TURN_LEFT) {
                printf("in curved track\n");
                R3Vector dist_to_start(bobsled->position - bobsled->track->start);
                dist_to_start.Project(bobsled->track->along);
                R3Vector vect_radius(track->center_point - position);
                double delta_theta = dist_to_start.Length()/vect_radius.Length();
                //double percent = dist_to_start.Length()/(2*M_PI*big_r);
                //bobsled->big_percent = percent;
                R3Line rotate_line(bobsled->track->center_pivot);
                bobsled->track->along.Rotate(bobsled->track->center_pivot.Vector(), delta_theta);
                bobsled->little_theta -= M_PI/2;
            }
            else
                bobsled->big_percent = 0;
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
    if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT) {
        R3Vector temp(bobsled->position - track->start);
        temp.Project(track->along);
        R3Point center_point(track->start);
        center_point += temp;
        R3Vector normal(center_point - bobsled->position);
        normal.Normalize();
        double dot = fg.Dot(normal);
        //printf("\ndot = %f\n", dot);
        R3Vector centripetal(R3null_vector);
        centripetal = bobsled->mass * bobsled->velocity * bobsled->velocity / r;
        centripetal = centripetal.Length() * normal;
        fn = fg.Dot(normal) * normal + centripetal;
    }
    
    else {
        R3Vector normal(track->center_point - bobsled->position);
        //R3Vector normal(0, 1, 0);
        /*R3Vector init_normal(R3null_vector);
        R3Track *track(bobsled->track);
        init_normal = bobsled->big_percent * track->endNormal + (1 - bobsled->big_percent) * track->startNormal;
        R3Point little_center(track->center_point);
        little_center += -1 * init_normal * (track->big_radius - track->radius);
        R3Vector normal(little_center - bobsled->position); */
        normal.Normalize();
        double dot = fg.Dot(normal);
        R3Vector centripetal(R3null_vector);
        centripetal = bobsled->mass * bobsled->velocity * bobsled->velocity / r;
        centripetal = centripetal.Length() * normal;
        //R3Vector big_normal(track->center_point - bobsled->position);
        //R3Vector big_centripetal(R3null_vector);
        //double R = R3Distance(bobsled->position, bobsled->track->center_point);
        //big_centripetal = bobsled->mass * bobsled->velocity * bobsled->velocity / R;
        //big_centripetal = big_centripetal.Length() * big_normal;
        fn = fg.Dot(normal) * normal;// + centripetal + big_centripetal;
    }
    // force of friction
    R3Vector fk(R3null_vector);
    fk = track->cof * fn.Length() * -1 * bobsled->velocity;
    
    // total force
    force = fg + fn + fk;
    //force = fk + fn;
    return force;
}

    
