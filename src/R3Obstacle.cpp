// Source file for the R3 bobsled class 



// Include files 

// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "R3Bobsled.h"
#include "R3Obstacle.h"
#include <cmath>
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#   define TRUE true
#   define FALSE false
#endif



double RandNum(void)
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



// Create a new snowball that is a copy of a previous one
R3Obstacle *CopySnowball(R3Obstacle *obstacle)
{
  R3Obstacle *new_obstacle = new R3Obstacle();
  *new_obstacle = *obstacle;
  new_obstacle->obstacle_shape = new R3Shape(*obstacle->obstacle_shape);
  new_obstacle->obstacle_shape->sphere = new R3Sphere(*obstacle->obstacle_shape->sphere);
  return new_obstacle;
}

// create snowball
void CreateSnowballs(R3Scene *scene)
{
  // TODO CHANGE WHEN MORE BOBSLEDS
  for (unsigned int i = 0; i < 1; i++)
  {
    R3Bobsled *bobsled = scene->bobsleds[i];
    R3Track *track = NULL;
    R3Obstacle *track_obstacle = NULL;
    
    // generate a snow ball with some probability
    if (RandNum() > 0.8)
    {
      if (bobsled->track->obstacle != NULL)
      {
        track = bobsled->track;
        track_obstacle = track->obstacle;
      }
      else if (bobsled->track->next != NULL && bobsled->track->next->obstacle != NULL)
      {
        track = bobsled->track->next;
        track_obstacle = track->obstacle;
      }
      if (track_obstacle != NULL)
      {
        // copy this track's obstacle
        R3Obstacle *obstacle = CopySnowball(track_obstacle);
        scene->obstacles.push_back(obstacle);
        printf("x: %f\n", bobsled->position.X());
        double rand_x = 10 - 20 * RandNum();
        R3Point new_center = track->end + 5 * track->endNormal + 
          0.5 * R3Vector(rand_x + bobsled->position.X(), 0, 0);
        obstacle->obstacle_shape->sphere->Reposition(new_center);
        obstacle->velocity = R3Vector(0, -20, 0);
      }
    }
  }
}

// Update Obstacle Positions
void UpdateObstacles(R3Scene *scene, double delta_time)
{
  int i = 0;
  while (i < scene->obstacles.size())
  {
    R3Obstacle *obstacle = scene->obstacles[i];
    // only moves spheres right now!!!
    if (obstacle->obstacle_shape->type == R3_SPHERE_SHAPE)
    {
      R3Sphere *sphere = obstacle->obstacle_shape->sphere;
      sphere->Reposition(sphere->Center() + obstacle->velocity * delta_time);
    }
    
    // delete obstacle if it is out of visible range of all bobsleds
    bool remove = true;
    for (unsigned int j = 0; j < scene->bobsleds.size(); j++)
    {
      R3Track *track = scene->bobsleds[j]->track;
      R3Box *bbox = ObstacleBBox(obstacle);
      if (bbox->ZMin() < track->bbox.ZMax() ||
          bbox->YMax() > track->bbox.YMin())
        remove = false;
    }
    if (remove)
      scene->obstacles.erase(scene->obstacles.begin() + i);
    else
      i++;
  }
}

R3Box *ObstacleBBox(R3Obstacle *obstacle)
{
  R3Box *obstacle_box = new R3Box();
  
  switch (obstacle->obstacle_shape->type)
  {
    case R3_MESH_SHAPE:
    {
      *obstacle_box = obstacle->obstacle_shape->mesh->bbox;
      obstacle_box->Transform(obstacle->transformation);
      break;
    }
    case R3_SPHERE_SHAPE:
      *obstacle_box = obstacle->obstacle_shape->sphere->BBox();
      break;
  }
  return obstacle_box;
}








  
  
  
  
  
  

