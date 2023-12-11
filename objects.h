#pragma once
#include <cmath>
#include "geometry.h"
#include "ray.h"
#include "material.h"
#include <iostream>
#include <limits>

using namespace std;

class Object 
{
  public:
    Material material;
    virtual bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const = 0;  
};

class Sphere : public Object
{
    float r;
    Vec3f c;
      
  public:
    Sphere(const Vec3f &c, const float &r) : c(c), r(r) { }
    Sphere(const Vec3f &c, const float &r, const Material &mat) : c(c), r(r)
    { 
        Object::material = mat;
    }

    bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const
    {
        Vec3f e_minus_c = ray.origin - c;
        float d2 = ray.direction * ray.direction;
        
        float disc = pow(ray.direction * e_minus_c, 2) - d2 * (e_minus_c * e_minus_c - r * r);
        
        if (disc < 0)
        {
            return false;
        }
        else
        {
            bool ray_inside_sphere = e_minus_c * e_minus_c <= r * r;
            if (ray_inside_sphere)
            {
                float t1 = (-ray.direction * e_minus_c + sqrt(disc)) / d2;
                t = t1;
            }
            else
            {
                float t2 = (-ray.direction * e_minus_c - sqrt(disc)) / d2;
                t = t2;
            }
            
            Vec3f hit_point = ray.origin + ray.direction * t;
            normal = (hit_point - c).normalize();
            return true;
        }
    }
};


struct Cuboid : Object
{
	Vec3f p_min;
    Vec3f p_max;
	
    Cuboid(const Vec3f& min, const Vec3f& max, const Material& mat) 
        : p_min{min}, p_max{max} {
            Object::material = mat;
        }
	
bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const
{
    // Initialize minimum and maximum intersection distances along the ray
    float tMin = numeric_limits<float>::min();
    float tMax = numeric_limits<float>::max();
    
    // Define boundaries of the cuboid along each axis
    float minX = min(p_min.x, p_max.x);
    float minY = min(p_min.y, p_max.y);
    float minZ = min(p_min.z, p_max.z);
    float maxX = max(p_min.x, p_max.x);
    float maxY = max(p_min.y, p_max.y);
    float maxZ = max(p_min.z, p_max.z);
    
    // Extract ray origin and direction
    Vec3f e = ray.origin;
    Vec3f d = ray.direction;
    
    // Intersection with X-axis planes
    if (d.x == 0) {
        if (e.x < minX || e.x > maxX)
            return false;
    } else {
        // Calculate intersection distances along X-axis
        float txmin = (minX - e.x) / d.x;
        float txmax = (maxX - e.x) / d.x;

        if (txmin > txmax)
            swap(txmin, txmax);

        tMin = max(tMin, txmin);
        tMax = min(tMax, txmax);

        if (tMin > tMax || tMax < 0)
            return false;
    }

    t = tMin;

    // Intersection with Y-axis planes (similar logic as X-axis)
    if (d.y == 0) {
        if (e.y < minY || e.y > maxY)
            return false;
    } else {
        float tymin = (minY - e.y) / d.y;
        float tymax = (maxY - e.y) / d.y;

        if (tymin > tymax)
            swap(tymin, tymax);

        tMin = max(tMin, tymin);
        tMax = min(tMax, tymax);

        if (tMin > tMax || tMax < 0)
            return false;
    }

    t = tMin;

    // Intersection with Z-axis planes (similar logic as X and Y axes)
    if (d.z == 0) {
        if (e.z < minZ || e.z > maxZ)
            return false;
    } else {
        float tzmin = (minZ - e.z) / d.z;
        float tzmax = (maxZ - e.z) / d.z;

        if (tzmin > tzmax)
            swap(tzmin, tzmax);

        tMin = max(tMin, tzmin);
        tMax = min(tMax, tzmax);

        if (tMin > tMax || tMax < 0)
            return false;
    }

    t = tMin;

    // Calculate the intersection point
    Vec3f hit_point = e + d * t;

    // Define a small epsilon value for accuracy
    const float eps = 0.01;

    // Determine the surface normal at the intersection point based on cuboid faces
    if (abs(hit_point.x - minX) < eps)
        normal = Vec3f(-1, 0, 0);
    else if (abs(hit_point.x - maxX) < eps)
        normal = Vec3f(1, 0, 0);
    else if (abs(hit_point.y - minY) < eps)
        normal = Vec3f(0, -1, 0);
    else if (abs(hit_point.y - maxY) < eps)
        normal = Vec3f(0, 1, 0);
    else if (abs(hit_point.z - minZ) < eps)
        normal = Vec3f(0, 0, -1);
    else if (abs(hit_point.z - maxZ) < eps)
        normal = Vec3f(0, 0, 1);

    
    return true;
}
	
	
};