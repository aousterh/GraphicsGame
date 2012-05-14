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



double ANGLE_SHIFT = 1;




double Rand(void)
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




////////////////////////////////////////////////////////////
// Updating Bobsled
////////////////////////////////////////////////////////////
void UpdateBobsled(R3Scene *scene, double current_time, double delta_time, 
				   bool force_left, bool force_right)
{
	// update each sled in the scene
	for (int i = 0; i < 1/*scene->NBobsleds()*/; i++) {
		// get the current sled and its track segment
		R3Bobsled *bobsled = scene->Bobsled(i);
		R3Track *track = bobsled->track;

        //if (!bobsled->isFalling) {
            double big_r;
            R3Vector new_along(track->along);
            R3Vector new_normal(track->startNormal);
        
            // find the closest point on the along vector of the track
            R3Point position(bobsled->position);
            R3Vector ve_along(bobsled->position - track->start);
            ve_along.Project(track->along);
            R3Point little_center(track->start);
            little_center += ve_along;
            double r = R3Distance(little_center, bobsled->position);
            R3Vector force(R3null_vector);
            force = Force(bobsled, r, delta_time);
            R3Vector velocity(R3null_vector);
            velocity = bobsled->velocity + delta_time * force/bobsled->mass;
            //velocity.Print();
            //printf("\n");
            // Forward translation on a straight track
            R3Vector v_along(R3null_vector);
            if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT || track->type == TRACK_EXIT_RIGHT || track->type == TRACK_EXIT_LEFT) {
                v_along = bobsled->velocity.Dot(track->along) * track->along * delta_time;
                bobsled->position.Translate(v_along);
                for (int j = 0; j < NUM_SLEDS; j++)
                {
                    bobsled->sleds[j]->mesh->Translate(v_along.X(), v_along.Y(), v_along.Z());
                }
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
                double delta_theta =  dist/vect_radius.Length();
                //printf("delta_theta = %f\n", delta_theta);
                double percent = delta_theta/(M_PI/2);
                //bobsled->big_percent += percent;
                //dist_to_start.Project(bobsled->track->startNormal);
                //double percent = dist/(2*M_PI*bobsled->big_r);
                //new_along = bobsled->big_percent * track->endPlane->Normal() + (1 - bobsled->big_percent) * track->along;
            
                //new_along.Normalize();
                R3Line rotate_line(track->center_pivot);
                new_along.Rotate(rotate_line.Vector(), delta_theta);
                new_normal.Rotate(rotate_line.Vector(), delta_theta);
                R3Point old_position(bobsled->position);
                bobsled->position.Rotate(rotate_line, delta_theta);
                for (int j = 0; j < NUM_SLEDS; j++)
                {
                    bobsled->sleds[j]->mesh->Rotate(delta_theta, rotate_line);
                }
                bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
                bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
                bobsled->masks->mesh->Rotate(delta_theta, rotate_line);


                bobsled->camera->eye.Rotate(rotate_line, delta_theta);
                bobsled->camera->right.Rotate(rotate_line.Vector(), delta_theta);
                bobsled->camera->up.Rotate(rotate_line.Vector(), delta_theta);
                bobsled->camera->towards.Rotate(rotate_line.Vector(), delta_theta);
            }
		
		
            if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT || track->type == TRACK_EXIT_RIGHT || track->type == TRACK_EXIT_LEFT) {
                // Side rotation on a straight track
                R3Vector v_side(R3null_vector);
                R3Vector v_down(R3null_vector);
                R3Vector down(0, -1, 0);
                R3Vector rotate_vector(R3null_vector);
                R3Line rotate_line(track->start, track->along);
                double sign;
                R3Vector temp(R3null_vector);
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
                //temp.SetZ(0);
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
            
                // add vibration
                if (bobsled->x_vibration != 0)
                    delta_theta += bobsled->x_vibration;
      
                bobsled->little_theta += delta_theta;
            
                bobsled->position.Rotate(rotate_line, delta_theta);
                for (int j = 0; j < NUM_SLEDS; j++)
                {
                    bobsled->sleds[j]->mesh->Rotate(delta_theta, rotate_line);
                }
                bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
                bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
                bobsled->masks->mesh->Rotate(delta_theta, rotate_line);
                bobsled->camera->eye.Rotate(rotate_line, delta_theta);
                bobsled->camera->right.Rotate(rotate_line.Vector(), delta_theta);
                bobsled->camera->up.Rotate(rotate_line.Vector(), delta_theta);
                bobsled->camera->towards.Rotate(rotate_line.Vector(), delta_theta);
            }
        
        // Side rotation on a curved track
        
        else {
            R3Vector init_normal(track->startNormal);
            init_normal.Normalize();
            R3Point little_center(track->center_point);
            little_center +=  -1 * init_normal * (track->big_radius);
            R3Line rotate_line(little_center, track->along); 
            R3Vector normal(little_center - bobsled->position); 
            normal.Normalize();
            R3Vector side(track->along);
            side.Cross(normal);
            side.Normalize();
            double R = R3Distance(bobsled->position, little_center);
            double side_dist = -1 * bobsled->velocity.Dot(side) * delta_time;
            double delta_theta = side_dist/R;
            if (force_right) {
                double v_change = (ANGLE_SHIFT * R) * delta_time;
                velocity += v_change * side;
            }
            if (force_left) {
                double v_change = (ANGLE_SHIFT * R) * delta_time;
                velocity -= v_change * side;
            }
            bobsled->little_theta += delta_theta;
            
            bobsled->position.Rotate(rotate_line, delta_theta);
			for (int j = 0; j < NUM_SLEDS; j++)
			{
				bobsled->sleds[j]->mesh->Rotate(delta_theta, rotate_line);
			}
            bobsled->skates->mesh->Rotate(delta_theta, rotate_line);
            bobsled->helmets->mesh->Rotate(delta_theta, rotate_line);
            bobsled->masks->mesh->Rotate(delta_theta, rotate_line);
			bobsled->camera->eye.Rotate(rotate_line, delta_theta);
			bobsled->camera->right.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->up.Rotate(rotate_line.Vector(), delta_theta);
			bobsled->camera->towards.Rotate(rotate_line.Vector(), delta_theta);
            
            
        }
        
        // check if over the edge
        if (!track->isCovered) {
            R3Vector startNormal(track->startNormal);
            if (track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT || track->type == TRACK_EXIT_RIGHT || track->type == TRACK_EXIT_LEFT) {
                double percent = R3Distance(little_center, track->start);
                percent /= R3Distance(track->start, track->end);
                startNormal = percent * track->endNormal + (1 - percent) * track->startNormal;
                startNormal.Normalize();
            }
            R3Vector vec1 = bobsled->position - little_center;
            R3Vector vec2 = startNormal;
            double over_edge = vec1.Dot(vec2);
            if (over_edge >= 0)
            {
                fprintf(stderr, "over edge\n");
                exit(-1);
            }
        }
        
        track->along = new_along;
        track->startNormal = new_normal;
        // check if bobsled is on new track
        R3Vector to_plane(track->end - bobsled->position);
        double dist_plane = to_plane.Dot(track->endPlane.Normal());
        if (dist_plane <= 0) {
            bobsled->track = track->next;
            printf("went to new track %d\n", bobsled->track->type);
            printf("track along = ");
            bobsled->track->along.Print();
            printf("\n");
            printf("track normal = ");
            bobsled->track->startNormal.Print();
            printf("\n");
            printf("velocity = ");
            bobsled->velocity.Print();
            printf("\n");
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
                //bobsled->little_theta -= M_PI/2;
            }
            
        }
		bobsled->velocity = velocity;
	}
}


////////////////////////////////////////////////////////////
// find force acting on bobsled
////////////////////////////////////////////////////////////

R3Vector Force(R3Bobsled *bobsled, double r, double delta_time) {
    R3Vector force(R3null_vector);
    R3Track *track = bobsled->track;
    // force of gravity
    R3Vector fg(R3null_vector);
    R3Vector gravity(0, -9.8, 0);
    fg = bobsled->mass * gravity;
    // normal force
    R3Vector fn(R3null_vector);
    if (track->type == TRACK_STRAIGHT || track->type == TRACK_APPROACH_LEFT || track->type == TRACK_APPROACH_RIGHT || track->type == TRACK_EXIT_RIGHT || track->type == TRACK_EXIT_LEFT) {
        R3Vector temp(bobsled->position - track->start);
        temp.Project(track->along);
        R3Point center_point(track->start);
        center_point += temp;
        R3Vector normal(center_point - bobsled->position);
        normal.Normalize();
        double dot = fg.Dot(normal);
        //printf("\ndot = %f\n", dot);
        R3Vector centripetal_vect(R3null_vector);
        R3Vector side(track->along);
        side.Cross(normal);
        side.Normalize();
         
        double centripetal = bobsled->mass * bobsled->velocity.Dot(side) * bobsled->velocity.Dot(side)/ r;
        centripetal_vect = centripetal * normal;
        R3Vector fg_other(fg);
        fg_other.Flip();
        double normal_fg = fg_other.Dot(normal);
        printf("normal_fg = %f\n", normal_fg);
        fn = centripetal_vect;
        if (normal_fg > 0)
            fn += fg_other.Dot(normal) * normal;
        //fn = fg_other.Dot(normal) * normal + centripetal_vect;
        //printf("fn = ");
        //fn.Print();
        //printf("\n");
    }
    
    else {
        //R3Vector normal(track->center_point - bobsled->position);
        //R3Vector normal(0, 1, 0);
        R3Vector init_normal(track->startNormal);
        init_normal.Normalize();
        R3Point little_center(track->center_point);
        little_center +=  -1 * init_normal * (track->big_radius);
        R3Vector normal(little_center - bobsled->position); 
        normal.Normalize();
        double dot = fg.Dot(normal);
        R3Vector centripetal_vect(R3null_vector);
        double centripetal;
        R3Vector side(track->along);
        side.Cross(normal);
        side.Normalize();
        centripetal = bobsled->mass * bobsled->velocity.Dot(side) * bobsled->velocity.Dot(side)/ R3Distance(little_center, bobsled->position);
        centripetal_vect = centripetal * normal;
        R3Vector big_normal(track->center_point - /*little_center);*/bobsled->position);
        big_normal.Normalize();
        R3Vector big_centripetal_vect(R3null_vector);
        double big_centripetal;
        double R = R3Distance(bobsled->position, bobsled->track->center_point);
        big_centripetal = bobsled->mass * bobsled->velocity.Dot(track->along) * bobsled->velocity.Dot(track->along)  / R;
        big_centripetal_vect = big_centripetal * big_normal;
        R3Vector fg_other(fg);
        fg_other.Flip();
        double normal_fg = fg_other.Dot(normal);
        fn = centripetal_vect + big_centripetal_vect;
        if (normal_fg > 0)
            fn += fg_other.Dot(normal) * normal;
        //fn = fg_other.Dot(normal) * normal + centripetal_vect + big_centripetal_vect;
    }
    // force of friction
    R3Vector fk(R3null_vector);
    fk = track->cof * fn.Length() * -1 * bobsled->velocity;
    
    // total force
    force = (fg + fn + fk);
    //force = fk + fn;
    return force;
}



////////////////////////////////////////////////////////////
// Check for collisions with rocks
////////////////////////////////////////////////////////////
void CheckCollisions(R3Scene *scene)
{
  const double MOVEMENT_WEIGHT = 0.08;

  // check each bobsled - //TODO CHANGE THIS WHEN WE HAVE MULTIPLE BOBSLEDS
  for (unsigned int i = 0; i < 1; i++)
  {
    R3Bobsled *bobsled = scene->bobsleds[i];
    bobsled->x_vibration = 0.0;
    R3Box &bbox = bobsled->sleds[0]->mesh->bbox;
 //   printf("bobsled: %f %f %f %f %f %f\n", bbox.XMin(), bbox.XMax(), bbox.YMin(), bbox.YMax(), bbox.ZMin(), bbox.ZMax());
    
    // check each rock for a collision
    for (unsigned int j = 0; j < scene->obstacles.size(); j++)
    {
      R3Obstacle *obstacle = scene->obstacles[j];
      R3Box intersection = bbox;
      intersection.Intersect(obstacle->bbox);
      bbox = obstacle->bbox;
    //  printf("rock box: %f %f %f %f %f %f\n", bbox.XMin(), bbox.XMax(), bbox.YMin(), bbox.YMax(), bbox.ZMin(), bbox.ZMax());
      if (intersection.XMin() < intersection.XMax() &&
          intersection.YMin() < intersection.YMax() &&
          intersection.ZMin() < intersection.ZMax())
      {
        double current_z = bobsled->velocity.Z();
        printf("intersection!\n");
        // slow down
        if (obstacle->hit_count == 0)
          bobsled->velocity.SetZ(current_z * 0.5);
      
        // add left-to-right vibration
        if (obstacle->hit_count - 2 * ((int) obstacle->hit_count / 2) == 0)
          bobsled->x_vibration = MOVEMENT_WEIGHT * (1 + Rand());
        else
          bobsled->x_vibration = - MOVEMENT_WEIGHT * (1 + Rand());
        
        obstacle->hit_count++;
      }
    }
  }
  
}

