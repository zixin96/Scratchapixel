//[header]
// A simple program that uses ray-tracing to render a scene made out of spheres
//[/header]
//[compile]
// Download the simpleshapes.cpp and geometry.h files to a folder.
// Open a shell/terminal, and run the following command where the files is saved:
//
// c++ -o simpleshapes simpleshapes.cpp -O3 -std=c++11 -DMAYA_STYLE
//
// Run with: ./simpleshapes. Open the file ./out.png in Photoshop or any program
// reading PPM files.
//[/compile]
//[ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//[/ignore]

#define _USE_MATH_DEFINES // for M_PI 
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <memory>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>

#include "geometry.h"

const float kInfinity = std::numeric_limits<float>::max();

// create a random number engine (for randomized seeds)
std::random_device rd;

// create a random number engine
std::mt19937 gen(rd());

// create a random number distribution
std::uniform_real_distribution<float> floatBetween0And1(0, 1);

inline float clamp(const float& lo, const float& hi, const float& v)
{
	return std::max(lo, std::min(hi, v));
}

inline float deg2rad(const float& deg)
{
	return deg * (float)M_PI / 180.f;
}

inline Vec3f mix(const Vec3f& a, const Vec3f& b, const float& mixValue)
{
	return a * (1 - mixValue) + b * mixValue;
}

struct Options
{
	uint32_t  width;
	uint32_t  height;
	float     fov;
	Matrix44f cameraToWorld;
};

// [comment]
// Object base class that provides a generic definition of the concept of geometry in the program
// [/comment]
class Object
{
public:
	Object() : color(floatBetween0And1(gen), floatBetween0And1(gen), floatBetween0And1(gen))
	{
	}

	virtual ~Object()
	{
	}

	// Method to compute the intersection of the object with a ray
	// Returns true if an intersection was found, false otherwise
	// See method implementation in children class for details
	virtual bool intersect(const Vec3f&, const Vec3f&, float&) const = 0;

	// Method to compute the surface data such as normal and texture coordinates at the intersection point.
	// See method implementation in children class for details
	virtual void getSurfaceData(const Vec3f&, Vec3f&, Vec2f&) const = 0;

	Vec3f color;
};

// [comment]
// Compute the roots of a quadratic equation
// [/comment]
bool solveQuadratic(const float& a,
                    const float& b,
                    const float& c,
                    float&       x0,
                    float&       x1)
{
	float discr = b * b - 4.f * a * c;

	if (discr < 0.f)
	{
		return false;
	}
	else if (discr == 0.f)
	{
		x0 = -0.5f * b / a;
		x1 = x0;
	}
	else
	{
		//!? use a different formula for solving quadratic equations
		//!? this method avoids "catastrophic cancellation" and is more stable

		float q = (b > 0.f) ? -0.5f * (b + sqrt(discr)) : -0.5f * (b - sqrt(discr));

		x0 = q / a;
		x1 = c / q;
	}

	return true;
}

// [comment]
// Sphere class. A sphere type object
// [/comment]
class Sphere : public Object
{
public:
	Sphere(const Vec3f& c, const float& r)
		: radius(r),
		  radiusSquare(r * r),
		  center(c)
	{
	}

	// [comment]
	// Ray-sphere intersection test
	//
	// \param orig is the ray origin
	//
	// \param dir is the ray (normalized)) direction
	//
	// \param[out] is the distance from the ray origin to the intersection point
	//
	// [/comment]
	bool intersect(const Vec3f& orig, const Vec3f& dir, float& t) const
	{
		float t0, t1; // solutions for t if the ray intersects

		#if 0

		// geometric solution
		Vec3f L = center - orig;

		float t_ca = L.dotProduct(dir);

		//!? if tca is negative, the intersection could potentially be behind the ray's origin
		//!? , which means there are no valid intersections
		if (t_ca < 0) { return false; }

		float dSquare = L.dotProduct(L) - t_ca * t_ca;

		//!? if d is larger than the radius (we compare their squared values, but the underlying resullt is the same)
		//!? , the ray doesn't intersect with the sphere
		if (dSquare > radiusSquare) { return false; }

		float t_hc = sqrt(radiusSquare - dSquare);

		t0 = t_ca - t_hc;
		t1 = t_ca + t_hc;

		#else

		// analytic solution
		Vec3f L = orig - center;

		float a = dir.dotProduct(dir);

		float b = 2 * dir.dotProduct(L);

		float c = L.dotProduct(L) - radiusSquare;

		if (!solveQuadratic(a, b, c, t0, t1)) { return false; }

		#endif

		//!? always make sure t0 < t1
		if (t0 > t1) { std::swap(t0, t1); }

		//!? if t0 < 0, the first intersection point is behind the ray
		if (t0 < 0)
		{
			//!? try the second intersection point: 
			t0 = t1;
			//!? if the section intersection point is also behind the ray, then there are no intersections
			if (t0 < 0) { return false; } // both t0 and t1 are negative
		}

		//!? here, we've found a valid intersection point
		t = t0;

		return true;
	}

	// [comment]
	// Set surface data such as normal and texture coordinates at a given point on the surface
	//
	// \param Phit is the point ont the surface we want to get data on
	//
	// \param[out] Nhit is the normal at Phit
	//
	// \param[out] tex are the texture coordinates at Phit
	//
	// [/comment]
	void getSurfaceData(const Vec3f& Phit, Vec3f& Nhit, Vec2f& tex) const
	{
		// the normal of a point on a sphere can simply be computed as the point position minus the sphere center
		Nhit = Phit - center;
		// don't forget to normalize it! 
		Nhit.normalize();

		//!? the texture coordinates at the intersection point are just the spherical coordinates of the point on
		//!? the sphere remapped to the range [0, 1]

		// In this particular case, the normal can be considered a point on a unit sphere
		// centered around the origin.

		// We can thus use the normal coordinates to compute the spherical coordinates of Phit, and then remap them to [0, 1]

		tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5; // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
		tex.y = acosf(Nhit.y) / M_PI;                     // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
	}

	float radius, radiusSquare;
	Vec3f center;
};

// [comment]
// Ray-plane intersection test
//
// \param n is the plane's normal
//
// \param p0 is a point representing how far the plane is from the world origin
//
// \param l0 is the ray's origin
//
// \param l is the ray's (normalized) direction
//
// \param t [out] is the parametric distance from the ray origin to the intersection point
//
// [/comment]
bool intersectPlane(const Vec3f& n,
                    const Vec3f& p0,
                    const Vec3f& l0,
                    const Vec3f& l,
                    float&       t)
{
	// assuming vectors are all normalized
	float denom = n.dotProduct(l);

	//!? if n and l gets close to 0, then either
	//!? the plane and the ray perfectly coincide
	//!? OR the ray is away from the plane

	if (denom > 1e-6)
	{
		//!? use equation from scratchapixel: 
		Vec3f p0l0 = p0 - l0;

		t = p0l0.dotProduct(n) / denom;

		return (t >= 0);
	}

	//!? if the denominator is lower than a very small value, there are no valid intersections
	return false;
}

// [comment]
// Ray-disk intersection test
//
// \param n is the disk's normal
//
// \param p0 is the disk's center position
//
// \param radius is the disk's radius
//
// \param l0 is the ray's origin
//
// \param l is the ray's (normalized) direction
//
// \param t [out] is the parametric distance from the ray origin to the intersection point
//
// [/comment]
bool intersectDisk(const Vec3f& n,
                   const Vec3f& p0,
                   const float& radius,
                   const Vec3f& l0,
                   const Vec3f& l)
{
	float t = 0;

	//!? does the ray intersect the plane in which the disk lies? 
	if (intersectPlane(n, p0, l0, l, t))
	{
		// compute the intersection point
		Vec3f p = l0 + l * t;

		// compute the distance from the point to the disk's center
		Vec3f v              = p - p0;
		float distanceSquare = v.dotProduct(v);

		// If this distance is lower or equal to the disk radius, then the ray intersects the disk
		float radiusSquare = radius * radius;
		return radiusSquare <= radiusSquare;
	}

	return false;
}


// [comment]
// Returns true if the ray intersects an object. The variable tNear is set to the closest intersection distance and hitObject
// is a pointer to the intersected object. The variable tNear is set to infinity and hitObject is set null if no intersection
// was found.
// [/comment]
bool trace(const Vec3f& orig, const Vec3f& dir, const std::vector<std::unique_ptr<Object>>& objects, float& tNear, const Object*& hitObject)
{
	tNear = kInfinity;

	std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();

	for (; iter != objects.end(); ++iter)
	{
		float t = kInfinity;
		if ((*iter)->intersect(orig, dir, t) && t < tNear)
		{
			hitObject = iter->get();
			tNear     = t;
		}
	}

	return (hitObject != nullptr);
}

// [comment]
// Compute the color at the intersection point if any (returns background color otherwise)
// [/comment]
Vec3f castRay(const Vec3f&                                orig,
              const Vec3f&                                dir,
              const std::vector<std::unique_ptr<Object>>& objects)
{
	Vec3f         hitColor  = 0;
	const Object* hitObject = nullptr; // this is a pointer to the hit object
	float         t;                   // this is the intersection distance from the ray origin to the hit point

	if (trace(orig, dir, objects, t, hitObject))
	{
		Vec3f Phit = orig + dir * t;
		Vec3f Nhit;
		Vec2f tex;
		hitObject->getSurfaceData(Phit, Nhit, tex);
		// Use the normal and texture coordinates to shade the hit point.
		// The normal is used to compute a simple facing ratio and the texture coordinate
		// to compute a basic checker board pattern
		float scale   = 4;
		float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
		hitColor      = std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
	}

	return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void render(const Options&                              options,
            const std::vector<std::unique_ptr<Object>>& objects)
{
	Vec3f* framebuffer = new Vec3f[options.width * options.height];
	Vec3f* pix         = framebuffer;

	float scale            = tan(deg2rad(options.fov * 0.5f));
	float imageAspectRatio = options.width / (float)options.height;

	// [comment]
	// Don't forget to transform the ray origin (which is also the camera origin
	// by transforming the point with coordinates (0,0,0) to world-space using the
	// camera-to-world matrix.
	// [/comment]
	Vec3f orig;
	options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
	for (uint32_t j = 0; j < options.height; ++j)
	{
		for (uint32_t i = 0; i < options.width; ++i)
		{
			// [comment]
			// Generate primary ray direction. Compute the x and y position
			// of the ray in screen space. This gives a point on the image plane
			// at z=1. From there, we simply compute the direction by normalized
			// the resulting vec3f variable. This is similar to taking the vector
			// between the point on the image plane and the camera origin, which
			// in camera space is (0,0,0):
			//
			// ray.dir = normalize(Vec3f(x,y,-1) - Vec3f(0));
			// [/comment]
			#ifdef MAYA_STYLE
            float x = (2 * (i + 0.5) / (float)options.width - 1) * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale * 1 / imageAspectRatio;
			#else

			float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
			#endif
			// [comment]
			// Don't forget to transform the ray direction using the camera-to-world matrix.
			// [/comment]
			Vec3f dir;
			options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
			dir.normalize();
			*(pix++) = castRay(orig, dir, objects);
		}
	}

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i)
	{
		char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}

	ofs.close();

	delete [] framebuffer;
}

// [comment]
// In the main function of the program, we create the scene (create objects)
// as well as set the options for the render (image width and height etc.).
// We then call the render() function.
// [/comment]
int main(int argc, char** argv)
{
	// creating the scene (adding objects and lights)
	std::vector<std::unique_ptr<Object>> objects;

	// generate a scene made of random spheres
	uint32_t numSpheres = 32;

	// fix a seed for debugging purposes
	gen.seed(0);

	for (uint32_t i = 0; i < numSpheres; ++i)
	{
		Vec3f randPos((0.5f - floatBetween0And1(gen)) * 10.f,
		              (0.5f - floatBetween0And1(gen)) * 10.f,
		              (0.5f + floatBetween0And1(gen) * 10.f));
		float randRadius = (0.5f + floatBetween0And1(gen) * 0.5f);
		// objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius))); // before C++ 14
		objects.push_back(std::make_unique<Sphere>(randPos, randRadius)); // make_unique since C++ 14
	}

	// setting up options
	Options options;
	options.width         = 640;
	options.height        = 480;
	options.fov           = 51.52f;
	options.cameraToWorld = Matrix44f(0.945519f, 0.f, -0.325569f, 0.f,
	                                  -0.179534f, 0.834209f, -0.521403f, 0.f,
	                                  0.271593f, 0.551447f, 0.78876f, 0.f,
	                                  4.208271f, 8.374532f, 17.932925f, 1.f);

	// finally, render
	render(options, objects);

	return 0;
}
