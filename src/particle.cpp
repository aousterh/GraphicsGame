// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include <cmath>
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#   define TRUE true
#   define FALSE false
#endif



////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////

// image properties
R3Matrix current_transformation = R3identity_matrix;


////////////////////////////////////////////////////////////////////////
// FUNCTION DECLARATIONS
////////////////////////////////////////////////////////////////////////

// primitive intersections
R3Intersection *IntersectSphere(R3Ray *ray, R3Sphere *sphere);
R3Intersection *IntersectBox(R3Ray *ray, R3Box *box);
R3Intersection *IntersectMesh(R3Ray *ray, R3Mesh *mesh);
R3Intersection *IntersectCylinder(R3Ray *ray, R3Cylinder *cylinder);
R3Intersection *IntersectCone(R3Ray *ray, R3Cone *cone);

// scene intersections
R3Intersection *IntersectScene(R3Ray *ray, R3Scene *scene);
R3Intersection *IntersectNode(R3Ray *ray, R3Node *node); 

////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double RandomNumber(void)
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


////////////////////////////////////////////////////////////////////////
// RAY-PRIMITIVE INTERSECTION
////////////////////////////////////////////////////////////////////////

R3Intersection
*IntersectSphere(R3Ray *ray, R3Sphere *sphere)
{
	R3Intersection *intersection = new R3Intersection();
	intersection->intersect = FALSE;

	// get required scene info
	R3Point O = sphere->Center();
	R3Point P_zero = ray->Start();
	R3Vector V = ray->Vector();
	double r = sphere->Radius();
	double r2 = pow(r,2);

	// transform info according to current transformation matrix
	O.Transform(current_transformation);
	R3Vector L = O - P_zero;
	double error = pow(10.0, -7.0);

	// assure intersection is not behind camera
	double Tca = L.Dot(V);
	if (Tca <= error)
		return intersection;

	// assure ray does not miss
	double d2 = L.Dot(L) - pow(Tca,2);
	if (d2 > r2)
		return intersection;

	// calculate points of intersection
	double Thc = sqrt(r2 - d2);
	double t1 = Tca - Thc;
	double t2 = Tca + Thc;

	// return nearest intersection
	R3Point P;
	double t;
	
	if (t1 <= t2 && t1 > error)
		t = t1;
	else
		t = t2;
	P = P_zero + t*ray->Vector();
		
	// find normal
	R3Vector N = (P - O) / (P - O).Length();

	// set return values
	intersection->intersection_point = P;
	intersection->intersection_normal = N;
	intersection->parametric_rayValue = t;
	intersection->intersect = TRUE;
	return intersection;
}

R3Intersection
*IntersectBox(R3Ray *ray, R3Box *box)
{
	R3Intersection *intersection = new R3Intersection;
	intersection->intersect = FALSE;

	// get min/max values
	R3Point min = box->Min();
	R3Point max = box->Max();

	// get ray information
	R3Point P_zero = ray->Start();
	R3Vector V = ray->Vector();

	// define 6 faces planes
	R3Plane Face1 = R3Plane(max, R3Vector( 0, 0, 1));
	R3Plane Face2 = R3Plane(max, R3Vector( 0, 1, 0));
	R3Plane Face3 = R3Plane(max, R3Vector( 1, 0, 0));
	R3Plane Face4 = R3Plane(min, R3Vector( 0, 0,-1));
	R3Plane Face5 = R3Plane(min, R3Vector( 0,-1, 0));
	R3Plane Face6 = R3Plane(min, R3Vector(-1, 0, 0));

	// transform points and planes
	min.Transform(current_transformation);
	max.Transform(current_transformation);
	Face1.Transform(current_transformation);
	Face2.Transform(current_transformation);
	Face3.Transform(current_transformation);
	Face4.Transform(current_transformation);
	Face5.Transform(current_transformation);
	Face6.Transform(current_transformation);

	// store face planes in array for easy traversal
	vector<R3Plane> faces;
    if (Face1.Normal().Dot(V) <= 0)
    	faces.push_back(Face1);
	if (Face2.Normal().Dot(V) <= 0)
		faces.push_back(Face2);
	if (Face3.Normal().Dot(V) <= 0)
		faces.push_back(Face3);
	if (Face4.Normal().Dot(V) <= 0)
		faces.push_back(Face4);
	if (Face5.Normal().Dot(V) <= 0)
		faces.push_back(Face5);
	if (Face6.Normal().Dot(V) <= 0)
		faces.push_back(Face6);


	// find intersection with each face plane
	// if one exists and store it as the box intersection
	// if it is the new closest point
	for (unsigned int i = 0; i < faces.size(); i++) {

		// find intersection of ray and plane of box
		R3Plane plane = faces[i];
		R3Vector N = plane.Normal();
		R3Point D = plane.Point();
		double t = (D - P_zero).Dot(N) / (V.Dot(N));

		// continue if intersection is behind camera
		if (t < 0)
			continue;
		
		// if t is greater than the closest found intersection
		// skip it
		if ((intersection->intersect) && (t >= intersection->parametric_rayValue))
			continue;
		
		// check to see if the point found is within
		// the edges defined by the face
        double error = pow(10.0, -7.0);
		R3Point P = P_zero + t*V;
		if ((P.X() >= (min.X() - error)) && (P.X() <= (max.X() + error)) && (P.Y() >= (min.Y() - error)) &&
			(P.Y() <= (max.Y() + error)) && (P.Z() >= (min.Z() - error)) && (P.Z() <= (max.Z() + error)))
		{
			intersection->intersection_normal = N;
			intersection->intersection_point = P;
			intersection->parametric_rayValue = t;
			intersection->intersect = TRUE;
		}
	}
	return intersection;
}

R3Intersection
*IntersectMesh(R3Ray *ray, R3Mesh *mesh)
{
	R3Intersection *intersection = new R3Intersection();
	intersection->intersect = FALSE;

	// get ray information
	R3Point P_zero = ray->Start();
	R3Vector V = ray->Vector();

	for (int i = 0; i < mesh->NFaces(); i++) {
		
		R3MeshFace *face = mesh->Face(i);
		R3Plane plane = face->plane;
		plane.Transform(current_transformation);

		// find intersection of ray and plane of triangle
		R3Vector N = plane.Normal();
		R3Point D = plane.Point();
		double t = (D - P_zero).Dot(N) / (V.Dot(N));

		// continue if intersection is behind camera
		if (t <= 0)
			continue;

		// if t is greater than the closest found intersection
		// skip it
		if ((intersection->intersect) && (t > intersection->parametric_rayValue))
			continue;
		
		// check to see if the point found is within
		// the edges defined by the face
		R3Point P = P_zero + t*V;
		R3Point T1 = face->vertices[0]->position;
		R3Point T2 = face->vertices[1]->position;
		R3Point T3 = face->vertices[2]->position;
		T1.Transform(current_transformation);
		T2.Transform(current_transformation);
		T3.Transform(current_transformation);
		R3Vector V1, V2;

		// check first edge
		V1 = T1 - P;
		V2 = T2 - P;
		V2.Cross(V1);
		V2.Normalize();
		if (V.Dot(V2) < 0)
			continue;

		// check second edge
		V1 = T2 - P;
		V2 = T3 - P;
		V2.Cross(V1);
		V2.Normalize();
		if (V.Dot(V2) < 0)
			continue;

		// check third edge
		V1 = T3 - P;
		V2 = T1 - P;
		V2.Cross(V1);
		V2.Normalize();
		if (V.Dot(V2) < 0)
			continue;

		intersection->intersection_point = P;
		intersection->intersection_normal = N;
		intersection->parametric_rayValue = t;
		intersection->intersect = TRUE;
	}
	return intersection;
}

R3Intersection
*IntersectCylinder(R3Ray *ray, R3Cylinder *cylinder)
{
	R3Intersection *intersection = new R3Intersection();
	intersection->intersect = FALSE;
	double error = pow(10.0, -7.0);

	// get ray information
    ray->Transform(current_transformation.Inverse());
	R3Point P_zero = ray->Start();
	R3Vector V = ray->Vector();

	// check top and bottom planes
	double r = cylinder->Radius();
    double r2 = pow(r, 2.0);
	vector<R3Plane> end_planes;
	R3Point end1 = cylinder->Axis().Start();
	R3Point end2 = cylinder->Axis().End();
	R3Point center = cylinder->Center();
	R3Vector AxisVector = end1 - end2;
	AxisVector.Normalize();
	R3Ray *AxisRay = new R3Ray(center, AxisVector);
	R3Plane top(end1, AxisVector);
	R3Plane bottom(end2, -AxisVector);
	end_planes.push_back(top);
	end_planes.push_back(bottom);

	// if it is the new closest point
	for (unsigned int i = 0; i < end_planes.size(); i++) {

		// find intersection of ray and plane of box
		R3Plane plane = end_planes[i];
		R3Vector N = plane.Normal();
		R3Point D = plane.Point();
		double t = (D - P_zero).Dot(N) / (V.Dot(N));

		// continue if intersection is behind camera
		if (t <= 0)
			continue;
		
		// if t is greater than the closest found intersection
		// skip it
		if ((intersection->intersect) && (t >= intersection->parametric_rayValue))
			continue;
		
		// check to see if the point found is within
		// the edges defined by the face
		R3Point P = P_zero + t*V;

		R3Point O;
		if (i == 0) O = end1;
		else		O = end2;
		double dist = (P - O).Length();

		if (dist <= r + error)
		{
			intersection->intersection_normal = N;
			intersection->intersection_point = P;
			intersection->parametric_rayValue = t;
			intersection->intersect = TRUE;
		}
	}

	// check rounded face
    double a = pow(V.X(), 2.0) + pow(V.Z(), 2.0);
    double b = 2 * P_zero.X() * V.X() + 2 * P_zero.Z() * V.Z();
    double c = pow(P_zero.X(), 2.0) + pow(P_zero.Z(), 2.0) - r2;
	
    double t1 = (-b + sqrt(pow(b,2.0) - 4*a*c))/(2 * a);
    double t2 = (-b - sqrt(pow(b,2.0) - 4*a*c))/(2 * a);
    
    // return nearest intersection
	R3Point P;
	double t;
	if (t1 < error) {
		t = t2;
		P = P_zero + t2*ray->Vector();
	}
	else if (t2 < error) {
		t = t1;
		P = P_zero + t1*ray->Vector();
	}
	else if (t1 <= t2) {
		t = t1;
		P = P_zero + t1*ray->Vector();
	}
	else {
		t = t2;
		P = P_zero + t2*ray->Vector();
	}

    double Axis_T = AxisRay->T(P);
    if (abs(Axis_T) >= cylinder->Height()/2)
        return intersection;

    double radial_dist = (AxisRay->Point(Axis_T) - P).Length();
    if (radial_dist > r+error)
        return intersection;

    else if ((intersection->intersect) && (t >= intersection->parametric_rayValue))
		return intersection;

    else
	{
        R3Vector N = P - AxisRay->Point(AxisRay->T(P));
		N.Normalize(); N.Transform(current_transformation);
		intersection->intersection_normal = N;
		P.Transform(current_transformation);
		intersection->intersection_point = P;
		intersection->parametric_rayValue = t;
		intersection->intersect = TRUE;
	}
    ray->Transform(current_transformation);
    return intersection;
}

R3Intersection
*IntersectCone(R3Ray *ray, R3Cone *cone)
{
	R3Intersection *intersection = new R3Intersection();
	intersection->intersect = FALSE;
	double error = pow(10.0, -7.0);

	// get ray information
    ray->Transform(current_transformation.Inverse());
	R3Point P_zero = ray->Start();
	R3Vector V = ray->Vector();
	V.Normalize();

	// check top and bottom planes
	vector<R3Plane> end_planes;
	R3Point top = cone->Axis().End();
	R3Point base = cone->Axis().Start();
	double h = (top - R3zero_point).Length();
	double h2 = pow(h, 2.0);
	double r = cone->Radius() * (1 - h/cone->Height());
    double r2 = pow(r, 2.0);
	R3Vector AxisVector = top - base;
	AxisVector.Normalize();
	R3Ray *AxisRay = new R3Ray(base, AxisVector);
	R3Plane bottom(base, -AxisVector);
	end_planes.push_back(bottom);

	// if it is the new closest point
	for (unsigned int i = 0; i < end_planes.size(); i++) {

		// find intersection of ray and plane of box
		R3Plane plane = end_planes[i];
		R3Vector N = plane.Normal();
		R3Point D = plane.Point();
		double t = (D - P_zero).Dot(N) / (V.Dot(N));

		// continue if intersection is behind camera
		if (t <= 0)
			continue;
		
		// if t is greater than the closest found intersection
		// skip it
		if ((intersection->intersect) && (t >= intersection->parametric_rayValue))
			continue;
		
		// check to see if the point found is within
		// the edges defined by the face
		R3Point P = P_zero + t*V;

		R3Point O = base;
		double radial_dist = (P - O).Length();

		if (radial_dist <= cone->Radius() + error)
		{
			N.Transform(current_transformation);
			intersection->intersection_normal = N;
			P.Transform(current_transformation);
			intersection->intersection_point = P;
			intersection->parametric_rayValue = t;
			intersection->intersect = TRUE;
		}
	}

	// check rounded face
	double a = pow(V.X(), 2.0) + pow(V.Z(), 2.0) - ((r2/h2) * pow(V.Y(), 2.0));
    double b = (2.0*P_zero.X()*V.X()) + (2.0*P_zero.Z()*V.Z()) - (2.0*(r2/h2)*(P_zero.Y()*V.Y() - V.Y()*h));
    double c = pow(P_zero.X(), 2.0) + pow(P_zero.Z(), 2.0) - (r2/h2)*(pow(P_zero.Y(), 2.0) - 2.0*P_zero.Y()*h + h2);
	
    double t1 = (-b + sqrt(pow(b,2.0) - 4.0*a*c))/(2.0 * a);
    double t2 = (-b - sqrt(pow(b,2.0) - 4.0*a*c))/(2.0 * a);
    
    // return nearest intersection
	R3Point P;
	double t;
	if (t1 < error) {
		if (t2 < error)
			return intersection;
		t = t2;
		P = P_zero + t2*ray->Vector();
	}
	else if (t2 < error) {
		t = t1;
		P = P_zero + t1*ray->Vector();
	}
	else if (t1 <= t2) {
		t = t1;
		P = P_zero + t1*ray->Vector();
	}
	else {
		t = t2;
		P = P_zero + t2*ray->Vector();
	}

    double Axis_T = AxisRay->T(P);
    R3Point AxisPoint = AxisRay->Point(Axis_T);

    if ((AxisPoint.Y() >= top.Y()) || (AxisPoint.Y() <= base.Y()))
        return intersection;

    else if ((intersection->intersect) && (t >= intersection->parametric_rayValue))
		return intersection;

    else
	{
        R3Vector N = (P - AxisPoint) + (AxisVector * r/h);
        N.Normalize();
		N.Transform(current_transformation);
		intersection->intersection_normal = N;
		P.Transform(current_transformation);
		intersection->intersection_point = P;
		intersection->parametric_rayValue = t;
		intersection->intersect = TRUE;
	}
    ray->Transform(current_transformation);
    return intersection;
}

////////////////////////////////////////////////////////////////////////
// RAY-SCENE INTERSECTION
////////////////////////////////////////////////////////////////////////

R3Intersection
*IntersectScene(R3Ray *ray, R3Scene *scene)
{
	return IntersectNode(ray, scene->Root());
}

R3Intersection
*IntersectNode(R3Ray *ray, R3Node *node) 
{
	// update nested transformation
	current_transformation *= node->transformation;

	// closest intersection
	R3Intersection *intersection = new R3Intersection;
	intersection->intersect = FALSE;

	// if the node has children, find the
	// intersection with the child node
	if (node->children.size() != 0) {
		for (unsigned int i = 0; i < node->children.size(); i++) {
            R3Node *node_i = node->children[i];
			R3Intersection *current = IntersectNode(ray, node_i);
			if (current->intersect) {
				if ((! intersection->intersect) || (current->parametric_rayValue < intersection->parametric_rayValue)) {
					R3Intersection *temp = intersection;
					intersection = current;
					delete temp;
				}
			}
			else
				delete current;
        }   
	}

	if (node->parent && node->shape)// && IntersectBBox(ray, node->bbox))
	{
		// return the intersection of this node
		R3ShapeType type = node->shape->type;
		R3Intersection *current;

		if (type == R3_BOX_SHAPE)
			current = IntersectBox(ray, node->shape->box);
		else if (type == R3_SPHERE_SHAPE)
			current = IntersectSphere(ray, node->shape->sphere);
		else if (type == R3_CYLINDER_SHAPE) 
			current = IntersectCylinder(ray, node->shape->cylinder);
		else if (type == R3_CONE_SHAPE)
			current = IntersectCone(ray, node->shape->cone);
		else if (type == R3_MESH_SHAPE)
			current = IntersectMesh(ray, node->shape->mesh);
		else {
			fprintf(stderr, "Error: Undefined shape type\n");
			return NULL;
		}

		if (current->intersect) {
			if ((! intersection->intersect) || (current->parametric_rayValue < intersection->parametric_rayValue)) {
				current->node = node;
				R3Intersection *temp = intersection;
				intersection = current;
				delete temp;
			}
		}
		else
			delete current;
	}
	current_transformation *= node->transformation.Inverse();
	return intersection;
}


////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
	// Generate new particles for every source
	int NSources = scene->NParticleSources();
	for (int i = 0; i < NSources; i++)
	{
		R3ParticleSource *source = scene->ParticleSource(i);
		bool fixed = source->fixed;
		double rate = source->rate;
		double mass = source->mass;
		double drag = source->drag;
		double velocity = source->velocity;
		double lifetime = source->lifetime;
		double elasticity = source->elasticity;
		double angle_cutoff = source->angle_cutoff;
		R3Material *material = source->material;

		double NParticles = ceil(delta_time * rate);
		for (int i = 0; i < NParticles; i++)
		{
			// generate new pixel
			R3Particle *particle = new R3Particle();
			particle->drag = drag;
			particle->mass = mass;
			particle->fixed = fixed;
			particle->lifetime = lifetime;
			particle->material = material;
			particle->elasticity = elasticity;

			// generate random point on source surface
			// for initial particle position
			R3Shape *shape = source->shape;

			// spherical source
			if (shape->type == R3_BOX_SHAPE) {
			
				// get the box
				R3Box *box = shape->box;
				R3Point center = box->Centroid();
				R3Point max = box->Max();
				R3Point min = box->Min();
				double xmax = max.X();
				double ymax = max.Y();
				double zmax = max.Z();
				double xmin = min.X();
				double ymin = min.Y();
				double zmin = min.Z();

				//choose face randomly
				double face = RandomNumber();
				while (face == 1.0)
					face = RandomNumber();
				face *= 6;
				face = floor(face);

				double x, y, z;
				R3Vector N;
				if (face == 0) {
					x = xmin;
					y = (RandomNumber() - 0.5) * (ymax - ymin) + center.Y();
					z = (RandomNumber() - 0.5) * (zmax - zmin) + center.Z();
					N = R3Vector(-1,0,0);
				}
				if (face == 1) {
					x = xmax;
					y = (RandomNumber() - 0.5) * (ymax - ymin) + center.Y();
					z = (RandomNumber() - 0.5) * (zmax - zmin) + center.Z();
					N = R3Vector(-1,0,0);
				}
				if (face == 2) {
					x = (RandomNumber() - 0.5) * (xmax - xmin) + center.X();
					y = ymin;
					z = (RandomNumber() - 0.5) * (zmax - zmin) + center.Z();
					N = R3Vector(0,-1,0);
				}
				if (face == 3) {
					x = (RandomNumber() - 0.5) * (xmax - xmin) + center.X();
					y = ymax;
					z = (RandomNumber() - 0.5) * (zmax - zmin) + center.Z();
					N = R3Vector(0,1,0);
				}
				if (face == 4) {
					x = (RandomNumber() - 0.5) * (xmax - xmin) + center.X();
					y = (RandomNumber() - 0.5) * (ymax - ymin) + center.Y();
					z = zmin;
					N = R3Vector(0,0,-1);
				}
				if (face == 5) {
					x = (RandomNumber() - 0.5) * (xmax - xmin) + center.X();
					y = (RandomNumber() - 0.5) * (ymax - ymin) + center.Y();
					z = zmax;
					N = R3Vector(0,0,1);
				}
				
				R3Point position(x,y,z);
				N.Normalize();
				R3Vector A(0,0,1);
				A.Cross(N);
				double t1 = RandomNumber() * 2 * 3.14159;
				double t2 = RandomNumber() * sin(angle_cutoff);
				A.Rotate(N, t1);
				R3Vector V(A);
				A.Cross(N);
				V.Rotate(A, acos(t2));
                V.Normalize();
				V *= velocity;

				// set fields
				particle->position = position;
				particle->velocity = V;
				scene->particles.push_back(particle);
			}

			// spherical source
			else if (shape->type == R3_SPHERE_SHAPE) {

				// get the sphere
				R3Sphere *sphere = shape->sphere;

				// choose random z in [-1,1]
				double z = RandomNumber();
				double p = RandomNumber();
				if (p > 0.5) {
					z = -z;
				}

				// choose random 1 in [0,2PI]
				double t = RandomNumber();
				t *= 2*3.14159;

				// get x and y values
				double r = sqrt(1 - pow(z, 2.0));
				double x = r * cos(t);
				double y = r * sin(t);

				// translate relative to sphere center
				x += sphere->Center().X();
				y += sphere->Center().Y();
				z += sphere->Center().Z();
				
				// get values
				R3Point position(x,y,z);
				R3Vector N = (position - sphere->Center());
				N.Normalize();
				R3Vector A(0,0,1);
				A.Cross(N);
				double t1 = RandomNumber() * 2 * 3.14159;
				double t2 = RandomNumber() * sin(angle_cutoff);
				A.Rotate(N, t1);
				R3Vector V(A);
				A.Cross(N);
				V.Rotate(A, acos(t2));
				V *= velocity;

				// set fields
				particle->position = position;
				particle->velocity = V;
				scene->particles.push_back(particle);
			}

			// circular source
			else if (shape->type == R3_CIRCLE_SHAPE) {
				
				// get circle
				R3Circle *circle = shape->circle;
				R3Vector N = circle->Normal();
				R3Point  C = circle->Center();
				double   R = circle->Radius();

				// get uv coordinate system for circle
				double u1 = RandomNumber();
				double u2 = RandomNumber();
				double u3 = RandomNumber();
				R3Vector u(u1,u2,u3);
				u.Cross(N);
				R3Vector v(u);
				v.Cross(N);
				u.Normalize(); v.Normalize();

				// get random point on circle
				R3Point position;
				bool found = FALSE;
				while (!found) {
					double Urand = RandomNumber() * R;
					double Vrand = RandomNumber() * R;
					if (pow(Urand, 2.0) + pow(Vrand, 2.0) > pow(R, 2.0))
						continue;
					double p1 = RandomNumber();
					double p2 = RandomNumber();
					if (p1 > 0.5)
						Urand = -Urand;
					if (p2 > 0.5)
						Vrand = -Vrand;
					position = C + Urand * u + Vrand * v;
					found = TRUE;
				}

				// set fields
				R3Vector A(0,0,1);
				A.Cross(N);
				double t1 = RandomNumber() * 2 * 3.14159;
				double t2 = RandomNumber() * sin(angle_cutoff);
				A.Rotate(N, t1);
				R3Vector V(A);
				A.Cross(N);
				V.Rotate(A, acos(t2));
				V *= velocity;

				// set fields
				particle->position = position;
				particle->velocity = V;
				scene->particles.push_back(particle);
			}

			// mesh source
			else if (shape->type == R3_MESH_SHAPE) {
				
				// get mesh
				R3Mesh *mesh = shape->mesh;
				int NFaces = mesh->NFaces();

				// pick a random face
				int f = (int)(RandomNumber() * NFaces);
				R3MeshFace *face = mesh->Face(f);
				R3Vector N = face->plane.Normal();
				R3Point V1 = face->vertices[0]->position;
				R3Point V2 = face->vertices[1]->position;
				R3Point V3 = face->vertices[2]->position;

				R3Vector first = (V3 - V2);
				R3Ray first_ray(V2, first);
				R3Point P1 = first_ray.Point(RandomNumber());

				R3Vector second = (P1 - V1);
				R3Ray second_ray(V1, second);
				R3Point position = second_ray.Point(RandomNumber());

				// set fields
				R3Vector A(0,0,1);
				A.Cross(N);
				double t1 = RandomNumber() * 2 * 3.14159;
				double t2 = RandomNumber() * sin(angle_cutoff);
				A.Rotate(N, t1);
				R3Vector V(A);
				A.Cross(N);
				V.Rotate(A, acos(t2));
				V *= velocity;

				// set fields
				particle->position = position;
				particle->velocity = V;
				scene->particles.push_back(particle);
			}

			// other sources
			else {
				delete particle;
				//fprintf(stderr, "Unsupported particle source shape\n");
			}
		}
	}
}


////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
	// check for expired particles
	int NParticles = scene->NParticles();
	for (int i = NParticles-1; i >= 0; i--)
	{
		R3Particle *particle = scene->Particle(i);
		if (particle->lifetime != 0) {
			particle->lifetime -= delta_time;
			if (particle->lifetime <= 0) 
			{
				NParticles--;
				scene->particles[i] = scene->particles[NParticles];
				scene->particles.pop_back();
				delete particle;
			}
		}
	}
	
	// Update position for every particle
	vector<R3Vector> forces;
	for (int i = 0; i < scene->NParticles(); i++)
	{
		R3Particle *particle = scene->Particle(i);
		double m = particle->mass;


		R3Vector f(0,0,0);
        if (particle->fixed == 1) {
            forces.push_back(f);
            continue;
        }
		f += -particle->drag * particle->velocity;

		// account for spring forces
		int NSprings = particle->springs.size();
		for (int j = 0; j < NSprings; j++)
		{
			// get spring info
			R3ParticleSpring *spring = particle->springs[j];
			double ks = spring->ks;
			double kd = spring->kd;
			double s  = spring->rest_length;

			// get opposite particle
            R3Particle *other = spring->particles[0];
			if (other == particle)
				other = spring->particles[1];

			// compute distance D between particles
			R3Vector D = other->position - particle->position;
			double D_length = D.Length();
			D.Normalize();

			f += (((ks * (D_length - s)) + (kd * (other->velocity - particle->velocity).Dot(D))) * D);
		}

		// take sink forces into account
		bool exists = true;
		int NSinks = scene->NParticleSinks();
		for (int j = 0; j < NSinks; j++)
		{
			R3ParticleSink *sink = scene->ParticleSink(j);
			double intensity = sink->intensity;
			double ca = sink->constant_attenuation;
			double la = sink->linear_attenuation;
			double qa = sink->quadratic_attenuation;
			R3Shape *shape = sink->shape;
			R3Point center;

			// box sink
			if (shape->type == R3_BOX_SHAPE) {
				R3Box *box = shape->box;
				R3Point closest = box->ClosestPoint(particle->position);
				R3Vector d = (closest - particle->position);
				double d_length = d.Length();
				
				if ((particle->position - box->Centroid()).Length() <= (closest - box->Centroid()).Length()) {
					scene->particles[i] = scene->particles[scene->NParticles()-1];
					scene->particles.pop_back();
					exists = false;
					i--;
				}
				else {
					d.Normalize();
					f += (intensity / (ca + la*d_length + qa*d_length*d_length)) * d;
				}
			}

			// spherical sink
			else if (shape->type == R3_SPHERE_SHAPE) {
				R3Sphere *sphere = shape->sphere;
				center = sphere->Center();
				R3Vector d = (center - particle->position);
				double d_length = d.Length() - sphere->Radius();

				if (d_length <= 0) {
					scene->particles[i] = scene->particles[scene->NParticles()-1];
					scene->particles.pop_back();
					exists = false;
					i--;
				}
				else {
					d.Normalize();
					f += (intensity / (ca + la*d_length + qa*d_length*d_length)) * d;
				}
			}

			// cylindrical sink
			else if (shape->type == R3_CYLINDER_SHAPE)
				center = shape->cylinder->Center();
			
			// conical sink
			else if (shape->type == R3_CONE_SHAPE)
				center = shape->cone->Center();
			
			// mesh sink
			else if (shape->type == R3_MESH_SHAPE)
				center = shape->mesh->Center();
			
			// linear segment sink
			else if (shape->type == R3_SEGMENT_SHAPE)
				center = shape->segment->Centroid();
			
			// circular sink
			else {
				R3Circle *circle = shape->circle;
				R3Plane plane = circle->Plane();
				R3Vector N = plane.Normal();
				R3Point D = plane.Point();
				R3Ray perpendicular_ray(particle->position, -N);
				double t = (D - particle->position).Dot(N) / ((-N).Dot(N));
				R3Point closest = perpendicular_ray.Point(t);
				
				R3Vector radial_vector = closest - circle->Center();
				if (radial_vector.Length() > circle->Radius()) {
					radial_vector.Normalize();
					closest = circle->Center() + circle->Radius() * radial_vector;
				}

				R3Vector d = (closest - particle->position);
				double d_length = d.Length();

				if (d_length <= 0.1) { // ARBITRARY
					scene->particles[i] = scene->particles[scene->NParticles()-1];
					scene->particles.pop_back();
					exists = false;
					i--;
				}
				else {
					d.Normalize();
					f += (intensity / (ca + la*d_length + qa*d_length*d_length)) * d;
				}
			}
		}
		if (!exists) {
			delete particle;
			continue;
		}

		//// inter-particle forces
		//for (int j = 0; j < scene->NParticles(); j++)
		//{
		//	if (j == i) continue;
		//	R3Particle *other = scene->Particle(j);
		//	double other_mass = other->mass;
		//	R3Point other_pos = other->position;

		//	R3Vector r = (other_pos - particle->position);
		//	double r_length = r.Length();
		//	r.Normalize();
		//	f += (6.67428 * pow(10.0, -11.0) * m * other_mass) / (r_length * r_length) * r;
		//}

		forces.push_back(f / m);
	}

    // update position while checking for collisions
	for (int i = 0; i < scene->NParticles(); i++)
	{
		R3Particle *particle = scene->Particle(i);

        //double step = 0;
        //while (step < delta_time)
        //{
        //    R3Ray *velocity_ray = new R3Ray(particle->position, particle->velocity);
        //    R3Intersection *intersection = IntersectScene(velocity_ray, scene);
        //    double step_t = intersection->parametric_rayValue;
        //    if ((step_t >= 0.0) && ((step + step_t) <= delta_time))
        //    {
        //        //fprintf(stderr, "here\n");

        //        R3Vector N = intersection->intersection_normal;
        //        //fprintf(stderr, "particle velocity   = (%f,%f,%f)\n", particle->velocity.X(), particle->velocity.Y(), particle->velocity.Z());
        //        //fprintf(stderr, "intersection normal = (%f,%f,%f)\n", N.X(), N.Y(), N.Z());

        //        // reverse normal component of velocity
        //        R3Vector normal_component  = particle->velocity.Dot(N) * N;
        //        R3Vector tangent_component = particle->velocity - normal_component;
        //        //fprintf(stderr, "normal component  = (%f,%f,%f)\n", normal_component.X(), normal_component.Y(), normal_component.Z());
        //        //fprintf(stderr, "tangent component = (%f,%f,%f)\n", tangent_component.X(), tangent_component.Y(), tangent_component.Z());
        //        
		      //  particle->position += ((0.7 * step_t) * particle->velocity);
        //        particle->velocity = tangent_component - particle->elasticity * normal_component;
        //        //fprintf(stderr, "new particle velocity = (%f,%f,%f)\n", particle->velocity.X(), particle->velocity.Y(), particle->velocity.Z());
        //        step += (0.95 * step_t);
        //    }
        //    else {
    		    particle->position += delta_time * particle->velocity;
                //step += delta_time;
           /* }
            delete velocity_ray;
            delete intersection;*/
        //}

		// update velocity
        if (!particle->fixed) {
    		particle->velocity += (delta_time * scene->gravity);
		    particle->velocity += (delta_time * forces[i]);
        }
	}
}



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Draw every particle

  // REPLACE CODE HERE
  glDisable(GL_LIGHTING);
  glPointSize(2);
  glBegin(GL_POINTS);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
    const R3Point& position = particle->position;
    glVertex3d(position[0], position[1], position[2]);
  }   
  glEnd();
}



