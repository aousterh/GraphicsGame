// Source file for the scene file viewer



////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

#include "R3/R3.h"
#include "R3Scene.h"
#include "R3Bobsled.h"
#include "particle.h"
#include "R3Bobsled.h"
#include "cos426_opengl.h"
#include <cmath>
#include "Mountain.h"
#include "R3Obstacle.h"
//#include "glfont.h"

#include <OpenAL/al.h>
#include <OpenAL/alc.h>
#include "../AL/alut.h"

//#include <al.h>
//#include <alc.h>
//#include <alut.h>


////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS
////////////////////////////////////////////////////////////

static const double VIDEO_FRAME_DELAY = 1./25.; // 25 FPS 


////////////////////////////////////////////////////////////
// OPEN AL STUFF
////////////////////////////////////////////////////////////
/*
#define NUM_BUFFERS 1
#define NUM_SOURCES 1
#define NUM_ENVIRONMENTS 1

ALfloat listenerPos[]={0.0,0.0,4.0};
ALfloat listenerVel[]={0.0,0.0,0.0};
ALfloat listenerOri[]={0.0,0.0,1.0, 0.0,1.0,0.0};

ALfloat source0Pos[]={ -2.0, 0.0, 0.0};
ALfloat source0Vel[]={ 0.0, 0.0, 0.0};

ALuint  buffer[NUM_BUFFERS];
ALuint  source[NUM_SOURCES];
ALuint  environment[NUM_ENVIRONMENTS];

ALsizei size,freq;
ALenum  format;
ALvoid  *data;
ALboolean al_bool;
*/


////////////////////////////////////////////////////////////
// GLOBAL VARIABLES
////////////////////////////////////////////////////////////

// Program arguments

static char *input_scene_name = NULL;
static int integration_type = EULER_INTEGRATION;


//Mountain
Mountain * m = new Mountain();


// Display variables

static R3Scene *scene = NULL;
static R3Camera camera;
static R3Camera map_camera;
static R3Box *map_bbox;
static int show_faces = 1;
static int show_edges = 0;
static int show_bboxes = 0;
static int show_lights = 0;
static int show_camera = 0;
static int show_particles = 1;
static int show_particle_springs = 1;
static int show_particle_sources_and_sinks = 0;
static int save_image = 0;
static int save_video = 0;
static int quit = 0;
static int THIRD_PERSON = 1;
static int FIRST_PERSON = 0;
// forces left and right for players
// p1 is A-S-D-W, p2 is left-down-right-up
static bool force_left;
static bool force_right;

static int levDetail = 0;
static bool deadSound = false;

void playDeadSound();




// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 512;
static int GLUTwindow_width = 1024;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// GLUT command list

enum {
  DISPLAY_FACE_TOGGLE_COMMAND,
  DISPLAY_EDGE_TOGGLE_COMMAND,
  DISPLAY_BBOXES_TOGGLE_COMMAND,
  DISPLAY_LIGHTS_TOGGLE_COMMAND,
  DISPLAY_CAMERA_TOGGLE_COMMAND,
  DISPLAY_PARTICLES_TOGGLE_COMMAND,
  DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND,
  DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND,
  SAVE_IMAGE_COMMAND,
  SAVE_VIDEO_COMMAND,
  QUIT_COMMAND,
};



////////////////////////////////////////////////////////////
// TIMER CODE
////////////////////////////////////////////////////////////

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>
#endif

static double GetTime(void)
{
#ifdef _WIN32
  // Return number of seconds since start of execution
  static int first = 1;
  static LARGE_INTEGER timefreq;
  static LARGE_INTEGER start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    QueryPerformanceFrequency(&timefreq);
    QueryPerformanceCounter(&start_timevalue);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    LARGE_INTEGER current_timevalue;
    QueryPerformanceCounter(&current_timevalue);
    return ((double) current_timevalue.QuadPart - 
            (double) start_timevalue.QuadPart) / 
            (double) timefreq.QuadPart;
  }
#else
  // Return number of seconds since start of execution
  static int first = 1;
  static struct timeval start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    gettimeofday(&start_timevalue, NULL);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    struct timeval current_timevalue;
    gettimeofday(&current_timevalue, NULL);
    int secs = current_timevalue.tv_sec - start_timevalue.tv_sec;
    int usecs = current_timevalue.tv_usec - start_timevalue.tv_usec;
    return (double) (secs + 1.0E-6F * usecs);
  }
#endif
}



////////////////////////////////////////////////////////////
// SCENE DRAWING CODE
////////////////////////////////////////////////////////////

void DrawShape(R3Shape *shape)
{
  // Check shape type
  if (shape->type == R3_BOX_SHAPE) shape->box->Draw();
  else if (shape->type == R3_SPHERE_SHAPE) shape->sphere->Draw();
  else if (shape->type == R3_CYLINDER_SHAPE) shape->cylinder->Draw();
  else if (shape->type == R3_CONE_SHAPE) shape->cone->Draw();
  else if (shape->type == R3_MESH_SHAPE) shape->mesh->Draw();
  else if (shape->type == R3_SEGMENT_SHAPE) shape->segment->Draw();
  else if (shape->type == R3_CIRCLE_SHAPE) shape->circle->Draw();
  else fprintf(stderr, "Unrecognized shape type: %d\n", shape->type);
}



void LoadMatrix(R3Matrix *matrix)
{
  // Multiply matrix by top of stack
  // Take transpose of matrix because OpenGL represents vectors with 
  // column-vectors and R3 represents them with row-vectors
  R3Matrix m = matrix->Transpose();
  glMultMatrixd((double *) &m);
}


void LoadMaterial(R3Material *material, bool transparent) 
{
  GLfloat c[4];
  
  // Check if same as current
  static R3Material *current_material = NULL;
  //if (material == current_material) return;
  current_material = material;

  // Compute "opacity"
  double opacity = 1 - material->kt.Luminance();
  if (transparent)
    opacity *= 0.5;

  // Load ambient
  c[0] = material->ka[0];
  c[1] = material->ka[1];
  c[2] = material->ka[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);

  // Load diffuse
  c[0] = material->kd[0];
  c[1] = material->kd[1];
  c[2] = material->kd[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);

  // Load specular
  c[0] = material->ks[0];
  c[1] = material->ks[1];
  c[2] = material->ks[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c);

  // Load emission
  c[0] = material->emission.Red();
  c[1] = material->emission.Green();
  c[2] = material->emission.Blue();
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, c);

  // Load shininess
  c[0] = material->shininess;
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, c[0]);

  
  // Load texture
  if (material->texture && !transparent) {
    if (material->texture_index <= 0) {
      // Create texture in OpenGL
      GLuint texture_index;
      glGenTextures(1, &texture_index);
      material->texture_index = (int) texture_index;
      glBindTexture(GL_TEXTURE_2D, material->texture_index); 
      R2Image *image = material->texture;
      int npixels = image->NPixels();
      R2Pixel *pixels = image->Pixels();
      GLfloat *buffer = new GLfloat [ 4 * npixels ];
      R2Pixel *pixelsp = pixels;
      GLfloat *bufferp = buffer;
      for (int j = 0; j < npixels; j++) { 
        *(bufferp++) = pixelsp->Red();
        *(bufferp++) = pixelsp->Green();
        *(bufferp++) = pixelsp->Blue();
        *(bufferp++) = pixelsp->Alpha();
        pixelsp++;
      }
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
      //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
      glTexImage2D(GL_TEXTURE_2D, 0, 4, image->Width(), image->Height(), 0, GL_RGBA, GL_FLOAT, buffer);
      delete [] buffer;
    }

    // Select texture
    glBindTexture(GL_TEXTURE_2D, material->texture_index); 
    glEnable(GL_TEXTURE_2D);
  }
  else {
    glDisable(GL_TEXTURE_2D);
  }

  // Enable blending for transparent surfaces
  if (opacity < 1) {
    glDepthMask(false);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
  }
  else {
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDepthMask(true);
  }
}


void LoadCamera(R3Camera *camera)
{
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
 
  
  gluPerspective(2*180.0*camera->yfov/M_PI, (GLdouble) GLUTwindow_width /(GLdouble) GLUTwindow_height, camera->neardist, camera->fardist);

  // Set camera transformation
  R3Vector t = -(camera->towards);
  R3Vector& u = camera->up;
  R3Vector& r = camera->right;
  GLdouble camera_matrix[16] = { r[0], u[0], t[0], 0, r[1], u[1], t[1], 0, r[2], u[2], t[2], 0, 0, 0, 0, 1 };
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixd(camera_matrix);
  glTranslated(-(camera->eye[0]), -(camera->eye[1]), -(camera->eye[2]));
}

/* http://www.cs.sonoma.edu/other/CS375_OpenGL_Slides_pdf/Perspective_gluFrustum.pdf
 http://www.opengl.org/archives/resources/faq/technical/viewing.htm */


void LoadMapCamera(R3Camera *map_camera, R3Box *bbox, double ratio)
{
  const int A = -300;
  const int B = 1350;
  const int C = 800;
  
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double x_min = min(bbox->XMin(), bbox->XMax());
  double x_max = max(bbox->XMin(), bbox->XMax());
  double z_min = min(bbox->ZMin(), bbox->ZMax());
  double z_max = max(bbox->ZMin(), bbox->ZMax());
  double x_avg = (bbox->XMin() + bbox->XMax()) / 2;
  double z_avg = (bbox->ZMin() + bbox->ZMax()) / 2;
  
  double dnear = min(map_camera->eye.Y() - bbox->YMin(), map_camera->eye.Y() - bbox->YMax());
  double dfar = max(map_camera->eye.Y() - bbox->YMin(), map_camera->eye.Y() - bbox->YMax());
  dfar += 1000;
  
  double x_dim = x_max - x_avg;
  double z_dim = z_max - z_avg;
  
  if (x_dim >= z_dim * ratio)
  {
    // limited by width
    glOrtho(x_min, x_max + C, z_avg - x_dim / ratio - A, z_avg + x_dim / ratio + B, dnear, dfar); 
  }
  else
  {
    // limited by height
    glOrtho(x_avg - z_dim * ratio, x_avg + z_dim * ratio, z_min - A, z_max + B, dnear, dfar);
  }

  
  // Set camera transformation
  R3Vector t = -(map_camera->towards);
  R3Vector& u = map_camera->up;
  R3Vector& r = map_camera->right;
  GLdouble camera_matrix[16] = { r[0], u[0], t[0], 0, r[1], u[1], t[1], 0, r[2], u[2], t[2], 0, 0, 0, 0, 1 };
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixd(camera_matrix);
  glTranslated(-(map_camera->eye[0]), -(map_camera->eye[1]), -(map_camera->eye[2]));
}


void LoadLights(R3Scene *scene)
{
  GLfloat buffer[4];

  // Load ambient light
  static GLfloat ambient[4];
  ambient[0] = scene->ambient[0];
  ambient[1] = scene->ambient[1];
  ambient[2] = scene->ambient[2];
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Load scene lights
  for (int i = 0; i < (int) scene->lights.size(); i++) {
    R3Light *light = scene->lights[i];
    int index = GL_LIGHT0 + i;

    // Temporarily disable light
    glDisable(index);

    // Load color
    buffer[0] = light->color[0];
    buffer[1] = light->color[1];
    buffer[2] = light->color[2];
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);

    // Load attenuation with distance
    buffer[0] = light->constant_attenuation;
    buffer[1] = light->linear_attenuation;
    buffer[2] = light->quadratic_attenuation;
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);

    // Load spot light behavior
    buffer[0] = 180.0 * light->angle_cutoff / M_PI;
    buffer[1] = light->angle_attenuation;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glLightf(index, GL_SPOT_EXPONENT, buffer[1]);

    // Load positions/directions
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Load direction
      buffer[0] = -(light->direction.X());
      buffer[1] = -(light->direction.Y());
      buffer[2] = -(light->direction.Z());
      buffer[3] = 0.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_POINT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_SPOT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
    else if (light->type == R3_AREA_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
     else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->type);
      return;
    }

    // Enable light
    glEnable(index);
  }
}



void DrawNode(R3Scene *scene, R3Node *node)
{
  // Push transformation onto stack
  glPushMatrix();
  LoadMatrix(&node->transformation);

  // Load material
  if (node->material) LoadMaterial(node->material, false);

  // Draw shape
  if (node->shape) DrawShape(node->shape);

  // Draw children nodes
  for (int i = 0; i < (int) node->children.size(); i++) 
    DrawNode(scene, node->children[i]);

  // Restore previous transformation
  glPopMatrix();

  // Show bounding box
  if (show_bboxes) {
    GLboolean lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    node->bbox.Outline();
    if (lighting) glEnable(GL_LIGHTING);
  }
}


double IntersectGroundPlane(R3Plane p, R3Ray r)
{
	double vDotN = r.Vector().Dot(p.Normal());
	if (vDotN == 0) return -1.0;

	double t = -(p.Normal().Dot(r.Start().Vector()) + p.D()) / vDotN;
	return t;
}

void DrawMountain(R3Scene * scene, R3Camera * cam)
{
	R3Material * mat = new R3Material();
	mat->emission = R3Rgb(0, 0, 0, 0);
	mat->ka = R3Rgb(1, 1, 1, 1);
	mat->kd = R3Rgb(1, 1, 1, 1);
	mat->ks = R3Rgb(1, 1, 1, 1);
	mat->kt = R3Rgb(0, 0, 0, 0);
	mat->shininess = 10;
	mat->texture = NULL;
	LoadMaterial(mat, false);
	delete mat;

	//ground plane
	R3Vector ground_normal(0, 1, 0);
	R3Point ground_pt(0, -60, 0);
	R3Plane ground(ground_pt, ground_normal);


    /*
	double d = cam->neardist;
	double tanTheta = tan(cam->xfov);

    R3Vector side = ground.Normal();
    side = cam->right - cam->right.Dot(side) * ground.Normal();
    side.Normalize();
    
	R3Point bl = cam->eye + d * cam->towards
			- d * tanTheta * cam->up
			//- d * tanTheta * cam->right;
    - d * tanTheta * side;
    
	R3Point br = cam->eye + d * cam->towards
			- d * tanTheta * cam->up
			//+ d * tanTheta * cam->right;
    + d * tanTheta * side;

	R3Ray rray(br, br - cam->eye);
	R3Ray lray(bl, bl - cam->eye);

	//intersect with ground plane
	double t1 = IntersectGroundPlane(ground, rray);
	double t2 = IntersectGroundPlane(ground, lray);
    printf("t1 %f t2 %f\n", t1, t2);
	if (t1 < 0 || t2 < 0) return;

	R3Point r = rray.Point(t1);
	R3Point l = lray.Point(t2);
	R3Vector axis(1, 0, 0);
	R3Vector v1(r[0], 0, r[2]);
	R3Vector v2(l[0], 0, l[2]);
	v1.Normalize();
	v2.Normalize();

	double theta1 = acos(v1.Dot(axis));
	double theta2 = acos(v2.Dot(axis));
	if (v1[2] < 0) theta1 += 3.14159;
	if (v2[2] < 0) theta2 += 3.14159;

	theta1 /= (2.0 * 3.14159);
	theta2 /= (2.0 * 3.14159);
     
     */
    R3Ray r(cam->eye, cam->towards);
    double t = IntersectGroundPlane(ground, r);
    R3Point pt = r.Point(t);
    //pt += cam->eye;
    
    R3Vector towardsInGround = ground.Normal();
    towardsInGround = cam->towards - ground.Normal() * ground.Normal().Dot(cam->towards);
    towardsInGround.Normalize();
    towardsInGround = cam->towards;
    towardsInGround[1] = 0;
    R3Vector v1 = towardsInGround;
    v1.Rotate(ground.Normal(), 2 * cam->xfov);
    R3Vector v2 = towardsInGround;
    v2.Rotate(ground.Normal(), 2 * -cam->xfov);
    v1.Normalize();
    v2.Normalize();
    
	R3Vector axis(1, 0, 0);
	double theta1 = acos(v1.Dot(axis));
	if (v1[2] < 0) theta1 += 3.14159;
	theta1 /= (2.0 * 3.14159);

	const double mountainDist = 2000.0;

	R3Point startPt = (v1 * mountainDist).Point();
	R3Point endPt = (v2 * mountainDist).Point();

	///////TEMP//////////////////////////////////////
	startPt += cam->eye;
	endPt += cam->eye;
	startPt[1] = -ground.D();
	endPt[1] = -ground.D();
	///////TEMP//////////////////////////////////////


	int dist = R3Distance(startPt, endPt);

	int index = (double) m->width * theta1;
	//	printf("start index %f %d\n", theta1, dist);
	//	printf("towards ");
	//	cam->towards.Print();
	//	printf("\n");

	R3Point cur = startPt;
	R3Vector next = endPt - startPt;
	R3Vector back = ground.Normal();
	back.Cross(next);
	back.Flip();

	back.Normalize();
	next.Normalize();

	///////TEMP//////////////////////////////////////
	dist /= 3;
	next *= 3;
	/////END TEMP//////////////////////////////////


	glDisable(GL_LIGHTING);
	//front polygon
	for (int i = 0; i < dist; i++)
	{
		glBegin(GL_POLYGON);
		R3Point nextPt = cur + next;

		glNormal3d(-cam->towards[0], -cam->towards[1], -cam->towards[2]);
		glColor3d(1, 1, 1);
		glVertex3d(cur[0], -ground.D(), cur[2]);
		glVertex3d(cur[0], m->heights[(i + index) % m->width][0] - ground.D(), cur[2]);
		glVertex3d(nextPt[0], m->heights[(i + 1 + index) % m->width][0] - ground.D(), nextPt[2]);
		glVertex3d(nextPt[0], -ground.D(), nextPt[2]);
		glEnd();

		cur = nextPt;
	}
	glEnable(GL_LIGHTING);

	cur = startPt;
	for (int i = 0; i < dist; i++)
	{
		R3Point nextPt = cur + next;

		R3Point anchor = cur;
		R3Point over = nextPt;
		R3Point overUp = over + back;
		for (int j = 0; j < m->height; j++)
		{
			R3Point nextDepth =	anchor + back;

			R3Point p1(anchor[0], m->heights[(i + index) % m->width][j] - ground.D(), anchor[2]);
			R3Point p2(nextDepth[0], m->heights[(i + index) % m->width][j+1] - ground.D(), nextDepth[2]);
			R3Point p3(overUp[0], m->heights[(i + 1 + index) % m->width][j+1] - ground.D(), overUp[2]);
			R3Point p4(over[0], m->heights[(i + 1 + index) % m->width][j] - ground.D(), over[2]);

			glBegin(GL_POLYGON);

			R3Vector v = p2 - p1;
			R3Vector u = p3 - p1;
			R3Vector norm = v;
			norm.Cross(u);
			norm.Normalize();

			glNormal3d(norm[0], norm[1], norm[2]);
			glVertex3d(p1[0], p1[1], p1[2]);
			glVertex3d(p2[0], p2[1], p2[2]);
			glVertex3d(p3[0], p3[1], p3[2]);
			glEnd();

			glBegin(GL_POLYGON);

			v = p4 - p1;
			u = p3 - p1;
			norm = u;
			norm.Cross(v);
			norm.Normalize();

			glNormal3d(norm[0], norm[1], norm[2]);
			glVertex3d(p1[0], p1[1], p1[2]);
			glVertex3d(p3[0], p3[1], p3[2]);
			glVertex3d(p4[0], p4[1], p4[2]);
			glEnd();

			anchor = nextDepth;
			over += back;
			overUp += back;
		}
		cur = nextPt;
	}
}

void DrawBobsleds(R3Scene *scene, bool transparent)
{
    glEnable(GL_LIGHTING);
    // Draw all bobsleds
    for (int i = 0; i < scene->NBobsleds(); i++) {
        R3Bobsled *bobsled = scene->Bobsled(i);
        
        // Push transformation onto stack
        glPushMatrix();
        LoadMatrix(&bobsled->transformation);
        
        // Load sled material
        LoadMaterial(bobsled->sled_material, transparent);
        DrawShape(bobsled->sleds[levDetail]);
        
        // Load sled material
        LoadMaterial(bobsled->skates_material, transparent);
        DrawShape(bobsled->skates);
        
        // Load sled material
        LoadMaterial(bobsled->helmets_material, transparent);
        DrawShape(bobsled->helmets);
        
        // Load sled material
        LoadMaterial(bobsled->masks_material, transparent);
        DrawShape(bobsled->masks);
        
        // Restore previous transformation
        glPopMatrix();
    }
}

void DrawTracks(R3Scene *scene, bool transparent)
{
  glEnable(GL_LIGHTING);
  
  // Draw all tracks
  for (int i = 0; i < scene->NTracks(); i++) {
    R3Track *track = scene->Track(i);
    // Push transformation onto stack
    glPushMatrix();
    LoadMatrix(&track->transformation);

    // Load track material
    LoadMaterial(track->material, transparent);
    DrawShape(track->track_shape);

    // Restore previous transformation
    glPopMatrix();
  }
}

void DrawObstacles(R3Scene *scene, bool transparent)
{
  glEnable(GL_LIGHTING);
  
  // Draw all obstacles
  for (int i = 0; i < scene->NObstacles(); i++) {
    R3Obstacle *obstacle = scene->Obstacle(i);
    
    // Push transformation onto stack
    glPushMatrix();
    LoadMatrix(&obstacle->transformation);
    
    // Load track material
    LoadMaterial(obstacle->material, transparent);
    DrawShape(obstacle->obstacle_shape);
    
    // Restore previous transformation
    glPopMatrix();
  }
}

void DrawGround(R3Scene * scene)
{
    R3Material * mat = new R3Material();
	mat->emission = R3Rgb(0, 0, 0, 0);
	mat->ka = R3Rgb(1, 1, 1, 1);
	mat->kd = R3Rgb(1, 1, 1, 1);
	mat->ks = R3Rgb(1, 1, 1, 1);
	mat->kt = R3Rgb(0, 0, 0, 0);
	mat->shininess = 10;
	mat->texture = NULL;
	LoadMaterial(mat, false);
	delete mat;
    scene->ground->mesh->Draw();
}



void EnableFog(int fog_level)
{
  // Turn on fog
  if (fog_level > 0)
  {
    float fog_color[3] = {0.9f, 0.9f, 0.9f};
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogfv(GL_FOG_COLOR, fog_color);
    glFogi(GL_FOG_START, 750 - 150 * fog_level);
    glFogi(GL_FOG_END, 750);
  }
}

void DisableFog()
{
  glDisable(GL_FOG);
}



void DrawScene(R3Scene *scene, R3Camera * cam)
{
  EnableFog(scene->Bobsled(0)->track->fog);
  DrawGround(scene);
  DrawMountain(scene, cam);
  DrawNode(scene, scene->root);
  DrawObstacles(scene, false);
  DrawBobsleds(scene, false);
  DrawTracks(scene, false);
  DisableFog();
}



void DrawParticles(R3Scene *scene)
{
  // Get current time (in seconds) since start of execution
  double current_time = GetTime();
  static double previous_time = 0;


  static double time_lost_taking_videos = 0; // for switching back and forth
					     // between recording and not
					     // recording smoothly

  // program just started up?
  if (previous_time == 0) previous_time = current_time;

  // time passed since starting
  double delta_time = current_time - previous_time;

// Update particles
  UpdateParticles(scene, current_time - time_lost_taking_videos, delta_time, integration_type);

  // Generate new particles
  GenerateParticles(scene, current_time - time_lost_taking_videos, delta_time);

  // Render particles
  if (show_particles) RenderParticles(scene, current_time - time_lost_taking_videos, delta_time);

  // Remember previous time
  previous_time = current_time;
}


// Overlays the Map on top of existing content
void DrawMap(double width, double height)
{
  static double time;
  // draw another transparent image in center
  // twice as tall as wide for now
  double viewport_width = width * 0.2;
  double viewport_height = width * 0.4;
  //glViewport((width - viewport_width) * 0.5, (height - viewport_height) * 0.5,
  //           viewport_width, viewport_height);
    glViewport(width - viewport_width, 0, viewport_width, viewport_height);
  double ratio = viewport_width / viewport_height;  // ratio of width to height
  
  // Load map camera
  LoadMapCamera(&map_camera, map_bbox, ratio);
  
  // Load scene lights
  LoadLights(scene);
  
  // Do not draw scene camera
  
  // Do not draw scene lights
  
  glDepthMask(false);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  
  DrawTracks(scene, true);
  DrawBobsleds(scene, true);
  
  glDisable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDepthMask(true);
  
  
  // Draw Time
  char time_string [6] = {'0', '0', ':', '0', '0', '\n'};
  char characters[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
  
  // draw the current time
  if (!scene->Bobsled(0)->hasWon && !scene->Bobsled(0)->isFalling)
    time = GetTime();
  
  
  int tens = time / 10;
  int zero = time - tens * 10;
  int tenth = (time - tens * 10 - zero) / 0.1;
  int hundredth = (time - tens * 10 - zero - tenth * 0.1) / 0.01;
  
  
  time_string[0] = characters[tens];
  time_string[1] = characters[zero];
  time_string[3] = characters[tenth];
  time_string[4] = characters[hundredth];
  
  
  char *s = time_string;
  glRasterPos3d(-250, 0, -3800);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *(s++));
}



////////////////////////////////////////////////////////////
// GLUT USER INTERFACE CODE
////////////////////////////////////////////////////////////

void GLUTMainLoop(void)
{
  // Run main loop -- never returns 
  glutMainLoop();
}



void GLUTDrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}


void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Delete scene
  delete scene;

  // Exit
  exit(0);
}



void GLUTIdle(void)
{
  // Set current window
  if ( glutGetWindow() != GLUTwindow ) 
    glutSetWindow(GLUTwindow);  

  // Redraw
  glutPostRedisplay();
}



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize camera vertical field of view to match aspect ratio of viewport
  camera.yfov = atan(tan(camera.xfov) * (double) h/ (double) w); 

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}


void GLUTRedraw(void)
{
  static double old_time;
  
  // Initialize OpenGL drawing modes
  glEnable(GL_LIGHTING);
  glDisable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDepthMask(true);

  // Clear window 
  R3Rgb background = scene->background;
  glClearColor(background[0], background[1], background[2], background[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  R3Bobsled *bobsled = scene->Bobsled(0);
  if (bobsled->hasWon) {
    printf("in hasWon\n");
    R3Point point = R3Point(0,0,0);
    char buffer[40];
    sprintf(buffer, "PLAYER HAS WON!!");
    GLUTDrawText(point, buffer);
  }

  glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);

  // Load camera
  if (THIRD_PERSON)
	  LoadCamera(bobsled->camera3);
  else
	  LoadCamera(bobsled->camera1);

  // Load scene lights 
  LoadLights(scene);

  // Draw particles
  DrawParticles(scene);

  // Get current time (in seconds) since start of execution
  double current_time = GetTime();
  static double previous_time = 0;
  static double time_lost_taking_videos = 0; // for switching back and forth
      // between recording and not
      // recording smoothly
      
  // program just started up?
  if (previous_time == 0) previous_time = current_time;
      
  // time passed since starting
  double delta_time = current_time - previous_time;

  // Create snowballs
  CreateSnowballs(scene);
    
  // Update obstacle
  UpdateObstacles(scene, delta_time);
    
  // Check for collisions
  CheckCollisions(scene);

      
  // Update bobsleds
  UpdateBobsled(scene, current_time - time_lost_taking_videos, delta_time, force_left, force_right);
    if ((bobsled->isFalling == true && bobsled->hasWon == false) && deadSound == false)
    {
        deadSound = true;
        playDeadSound();
    }
  force_left = false;
  force_right = false;
      
  // Remember previous time
  previous_time = current_time;
      

  // Draw scene surfaces
  /*if (show_faces) {
    glEnable(GL_LIGHTING);
    DrawScene(scene, bobsled->camera);
  }*/

    // Draw scene surfaces
    if (show_faces) {
      glEnable(GL_LIGHTING);
      double current_time = GetTime();
    //  printf("delta not drawing: %f\n", current_time - old_time);
      current_time = old_time;
        if (THIRD_PERSON)
        {
            DrawScene(scene, bobsled->camera3);
        }
        else
        {
            DrawScene(scene, bobsled->camera1);
        }
      current_time = GetTime();
   //   printf("delta drawing: %f\n", current_time - old_time);
      current_time = old_time;
    }

    
  // Draw scene edges
  if (show_edges) {
    glDisable(GL_LIGHTING);
    glColor3d(1 - background[0], 1 - background[1], 1 - background[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      if (THIRD_PERSON)
      {
          DrawScene(scene, bobsled->camera3);
      }
      else
      {
          DrawScene(scene, bobsled->camera1);
      }    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // draw another transparent image in center, on top
  glDisable(GL_LIGHTING);
  DrawMap(GLUTwindow_width, GLUTwindow_height);

  // Quit here so that can save image before exit
  if (quit) {
    GLUTStop();
  }
  // Swap buffers 
  glutSwapBuffers();
}    

void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse button event
  if (state == GLUT_DOWN) {
    if (button == GLUT_LEFT_BUTTON) {
    }
    else if (button == GLUT_MIDDLE_BUTTON) {
    }
    else if (button == GLUT_RIGHT_BUTTON) {
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
    case GLUT_KEY_F1:
      save_image = 1;
      break;
    case GLUT_KEY_F2:
      save_video = save_video ^ 1;
      break;
    case GLUT_KEY_LEFT:
      force_left = true;
      break;
    case GLUT_KEY_RIGHT:
      force_right = true;
      break;
  }
  

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
    // Invert y coordinate
    y = GLUTwindow_height - y;
    
    // Process keyboard button event 
    switch (key) {
        case 'B':
        case 'b':
            show_bboxes = !show_bboxes;
            break;
            
        case 'E':
        case 'e':
            show_edges = !show_edges;
            break;
            
        case 'F':
        case 'f':
            show_faces = !show_faces;
            break;
         
        case 'P':
        case 'p':
            show_particles = !show_particles;
            break;
         
        case 'Q':
        case 'q':
        case 27: // ESCAPE
            quit = 1;
            break;
            
        case '=':
        	if (levDetail < NUM_SLEDS - 1)
        		levDetail++;
        	break;

        case '-':
        	if (levDetail > 0)
        		levDetail--;
        	break;

        case '3':
			THIRD_PERSON = true;
			FIRST_PERSON = false;
        	break;

        case '1':
			THIRD_PERSON = false;
			FIRST_PERSON = true;
        	break;

        case ' ': {
            R3Camera cam = *scene->bobsleds[0]->camera3;
            printf("camera %g %g %g  %g %g %g  %g %g %g  %g  %g %g \n",
                   camera.eye[0], camera.eye[1], camera.eye[2], 
                   camera.towards[0], camera.towards[1], camera.towards[2], 
                   camera.up[0], camera.up[1], camera.up[2], 
                   camera.xfov, camera.neardist, camera.fardist);
            /*printf("camera %g %g %g  %g %g %g  %g %g %g  %g  %g %g \n",
                   cam.eye[0], cam.eye[1], cam.eye[2], 
                   cam.towards[0], cam.towards[1], cam.towards[2], 
                   cam.up[0], cam.up[1], cam.up[2], 
                   cam.xfov, cam.neardist, cam.fardist);*/ 
            break; }
    }
    
    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
    
    // Remember modifiers 
    GLUTmodifiers = glutGetModifiers();
    
    // Redraw
    glutPostRedisplay();
}



void GLUTCommand(int cmd)
{
  // Execute command
  switch (cmd) {
  case DISPLAY_PARTICLES_TOGGLE_COMMAND: show_particles = !show_particles; break;
  case DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND: show_particle_springs = !show_particle_springs; break;
  case DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND: show_particle_sources_and_sinks = !show_particle_sources_and_sinks; break;
  case DISPLAY_FACE_TOGGLE_COMMAND: show_faces = !show_faces; break;
  case DISPLAY_EDGE_TOGGLE_COMMAND: show_edges = !show_edges; break;
  case DISPLAY_BBOXES_TOGGLE_COMMAND: show_bboxes = !show_bboxes; break;
  case DISPLAY_LIGHTS_TOGGLE_COMMAND: show_lights = !show_lights; break;
  case DISPLAY_CAMERA_TOGGLE_COMMAND: show_camera = !show_camera; break;
  case SAVE_IMAGE_COMMAND: save_image = 1; break;
  case SAVE_VIDEO_COMMAND: save_video = save_video ^ 1; break;
  case QUIT_COMMAND: quit = 1; break;
  }

  // Mark window for redraw
  glutPostRedisplay();
}



void GLUTCreateMenu(void)
{
  // Display sub-menu
  int display_menu = glutCreateMenu(GLUTCommand);
  glutAddMenuEntry("Particles (P)", DISPLAY_PARTICLES_TOGGLE_COMMAND);
  glutAddMenuEntry("Particle springs (R)", DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND);
  glutAddMenuEntry("Particle sources and sinks (S)", DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND);
  glutAddMenuEntry("Faces (F)", DISPLAY_FACE_TOGGLE_COMMAND);
  glutAddMenuEntry("Edges (E)", DISPLAY_EDGE_TOGGLE_COMMAND);
  glutAddMenuEntry("Bounding boxes (B)", DISPLAY_BBOXES_TOGGLE_COMMAND);
  glutAddMenuEntry("Lights (L)", DISPLAY_LIGHTS_TOGGLE_COMMAND);
  glutAddMenuEntry("Camera (C)", DISPLAY_CAMERA_TOGGLE_COMMAND);

  // Main menu
  glutCreateMenu(GLUTCommand);
  glutAddSubMenu("Display", display_menu);
  glutAddMenuEntry("Save Image (F1)", SAVE_IMAGE_COMMAND); 
  glutAddMenuEntry("Capture Video (F2)", SAVE_VIDEO_COMMAND);
 glutAddMenuEntry("Quit", QUIT_COMMAND);

  // Attach main menu to right mouse button
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}


ALCcontext *context;
ALCdevice *device;


void playDeadSound()
{
    ALuint alSource;

    alSourcef(alSource, AL_PITCH, 1);
	alSourcef(alSource, AL_GAIN, 1);
	alSource3f(alSource, AL_POSITION, 0, 0, 0);
	alSource3f(alSource, AL_VELOCITY, 0, 0, 0);
	alSourcei(alSource, AL_LOOPING, AL_FALSE);
    
	alGenSources(1, &alSource);

	alGetError();
	char*     alBuffer;         //data for the buffer
	ALenum alFormatBuffer;		//buffer format
	ALsizei   alFreqBuffer;     //frequency
	long       alBufferLen;     //bit depth
	unsigned int alSampleSet;
    
	alutLoadWAVFile("../sanka.wav",&alFormatBuffer, (void **) &alBuffer,(ALsizei *)&alBufferLen, &alFreqBuffer);
    
	//create  buffer
	alGenBuffers(1, &alSampleSet);
    
	//put the data into our sampleset buffer
	alBufferData(alSampleSet, alFormatBuffer, alBuffer, alBufferLen, alFreqBuffer);
    
	//assign the buffer to this source
	alSourcei(alSource, AL_BUFFER, alSampleSet);
    
	//release the data
	alutUnloadWAV(alFormatBuffer, alBuffer, alBufferLen, alFreqBuffer);
    
	alSourcei(alSource,AL_LOOPING,AL_FALSE);
    
	//play
	alSourcePlay(alSource);
}

//
void ALinit(int *argc, char **argv)
{
	//alutInit(argc, argv);

	ALCcontext *context;
	ALCdevice *device;

	device = alcOpenDevice(NULL);
	if (device == NULL)
	{
		printf("shit");
	}

	//Create a context
	context=alcCreateContext(device,NULL);

	//Set active context
	alcMakeContextCurrent(context);

	alListener3f(AL_POSITION, 0, 0, 0);
	alListener3f(AL_VELOCITY, 0, 0, 0);
	alListener3f(AL_ORIENTATION, 0, 0, -1);

	ALuint alSource;

	alSourcef(alSource, AL_PITCH, 1);
	alSourcef(alSource, AL_GAIN, 1);
	alSource3f(alSource, AL_POSITION, 0, 0, 0);
	alSource3f(alSource, AL_VELOCITY, 0, 0, 0);
	alSourcei(alSource, AL_LOOPING, AL_FALSE);

	alGenSources(1, &alSource);


	alGetError();
	char*     alBuffer;         //data for the buffer
	ALenum alFormatBuffer;		//buffer format
	ALsizei   alFreqBuffer;     //frequency
	long       alBufferLen;     //bit depth
	ALboolean    alLoop;        //loop
	unsigned int alSampleSet;

	alutLoadWAVFile("../music.wav",&alFormatBuffer, (void **) &alBuffer,(ALsizei *)&alBufferLen, &alFreqBuffer);

	//create  buffer
	alGenBuffers(1, &alSampleSet);

	//put the data into our sampleset buffer
	alBufferData(alSampleSet, alFormatBuffer, alBuffer, alBufferLen, alFreqBuffer);

	//assign the buffer to this source
	alSourcei(alSource, AL_BUFFER, alSampleSet);

	//release the data
	alutUnloadWAV(alFormatBuffer, alBuffer, alBufferLen, alFreqBuffer);

	alSourcei(alSource,AL_LOOPING,AL_TRUE);

	//play
	alSourcePlay(alSource);

	/*
	//to stop
	alSourceStop(alSource);
	alDeleteSources(1,&alSource);

	//delete our buffer
	alDeleteBuffers(1,&alSampleSet);

	context=alcGetCurrentContext();

	//Get device for active context
	device=alcGetContextsDevice(context);

	//Disable context
	alcMakeContextCurrent(NULL);

	//Release context(s)
	alcDestroyContext(context);

	//Close device

	alcCloseDevice(device);*/

	
}


void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("Video Game");

  // Initialize GLUT callback functions 
  glutIdleFunc(GLUTIdle);
  glutReshapeFunc(GLUTResize);
  glutDisplayFunc(GLUTRedraw);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  //glutMotionFunc(GLUTMotion);

  // Initialize graphics modes 
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
 
  // Create menus
  GLUTCreateMenu();

  // make full screen
  glutFullScreen();
}


////////////////////////////////////////////////////////////
// SCENE READING
////////////////////////////////////////////////////////////


R3Scene *
ReadScene(const char *filename)
{
  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Read file
  if (!scene->Read(filename)) {
    fprintf(stderr, "Unable to read scene from %s\n", filename);
    return NULL;
  }

  // Remember initial camera
  camera = scene->camera;
  
  // set track pointers to obstacles and remove track-associated obstacles
  // from the list
  int i = 0;
  while (i < scene->obstacles.size())
  {
    int track_num = scene->obstacles[i]->track_num;
    if (track_num >= 0 && track_num < scene->track_segments.size())
    {
      scene->track_segments[track_num]->obstacle = scene->obstacles[i];
      scene->obstacles.erase(scene->obstacles.begin() + i);
    }
    else
      i++;
  }
  
  // Return scene
  return scene;
}

void SetMapCamera(R3Scene *scene)
{
  // determine bounding box of tracks
  R3Box *bbox = new R3Box(R3null_box);
  for (unsigned int i = 0; i < scene->track_segments.size(); i++)
  {
    R3Track *track = scene->track_segments[i];
    bbox->Union(track->bbox);
  }
  map_bbox = bbox;
  
  
  // determine camera looking down from above (-y direction)
  double x_avg = (bbox->XMax() + bbox->XMin()) / 2;
  double z_avg = (bbox->ZMax() + bbox->ZMin()) / 2;
  double x_width = abs(bbox->XMax() - bbox->XMin());
  double z_width = abs(bbox->ZMax() - bbox->ZMin());
  double y_eye = max(bbox->YMax(), bbox->YMin()) + max(x_width, z_width);
  map_camera.eye = R3Point(x_avg, y_eye, z_avg);
//  printf("eye: %f %f %f\n", map_camera.eye.X(), map_camera.eye.Y(), map_camera.eye.Z());
  map_camera.towards = R3Vector(0, -1, 0);  // looking down in -Y
  map_camera.up = R3Vector(0, 0, -1);       // bobsleds move in -Z direction
  map_camera.right = map_camera.towards;
  map_camera.right.Cross(map_camera.up);
  map_camera.towards.Normalize();
  map_camera.up.Normalize();
  map_camera.right.Normalize();
}



////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////

int 
ParseArgs(int argc, char **argv)
{
  // Innocent until proven guilty
  int print_usage = 0;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if (!input_scene_name) {
      input_scene_name = *argv;
      argc--; argv++;
    }
    else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
  }

  // Check input_scene_name
  if (!input_scene_name || print_usage) {
    printf("Usage: game <input.scn> [-map <map.scn>]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

int 
main(int argc, char **argv)
{
    printf("entering main\n");
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(1);

    printf("about to glut init\n");
  // Initialize GLUT
  GLUTInit(&argc, argv);

    printf("about to al init\n");
  // Initialize AL

    ALinit(&argc, argv);

    printf("about to read\n");
  // Read scene
  scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);
  
    printf("camera\n");
  // Make map camera
  SetMapCamera(scene);
  
  // Set fog
  scene->track_segments[5]->fog = 1;
  scene->track_segments[6]->fog = 2;
  scene->track_segments[7]->fog = 3;
  scene->track_segments[8]->fog = 4;
  scene->track_segments[9]->fog = 5;
  scene->track_segments[10]->fog = 4;
  scene->track_segments[11]->fog = 3;
  scene->track_segments[12]->fog = 2;
  scene->track_segments[13]->fog = 1;
  
    printf("main loop\n");
  // Run GLUT interface
  GLUTMainLoop();

  // Return success 
  return 0;
}





// TODO: add animation to beginning: http://www.swiftless.com/tutorials/opengl/texture_animation.html



