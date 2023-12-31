#include "image.h"
#include "geometry.h"
#include "ray.h"
#include "objects.h"
#include "scene.h"
#include "light.h"
#include <cmath>
using namespace std;

// compile with: g++ -std=c++17 raytrace.cpp -o out.exe

Ray ray_to_pixel(Vec3f origin, const Image &img, int i, int j)
{
    Ray ray = Ray();

    float aspect_ratio = (float)img.width / img.height;
 
    float x = (-1 + 2 * (i + 0.5) / (float)img.width) * aspect_ratio;
    float y = (1 - 2 * (j + 0.5) / (float)img.height);
    float z = -1;

    ray.origin = origin;
    ray.direction = Vec3f(x, y, z).normalize();
    return ray;
}

void render(const Scene &scene)
{
    // dimenzije slike
    const int width = 1024;
    const int height = 768;
    
    Vec3f origin = Vec3f(0, 0, 0);
    
    Image img = Image(width, height);
    
    // crtanje slike
    for (size_t j = 0; j < height; ++j)
    {
        for (size_t i = 0; i < width; ++i)
        {
            Ray ray = ray_to_pixel(origin, img, i, j);
            Vec3f colour = scene.cast_ray(ray);
            img.set_pixel(i, j, colour);
        }
    }
    
    img.save("./img.ppm");
}

int main()
{
    Scene scene;
    
    Material red(Vec3f(1, 0, 0));
    Material green(Vec3f(0, 1, 0));
    Material blue(Vec3f(0, 0, 1));
    Material grey(Vec3f(0.5, 0.5, 0.5));
	
    Cuboid base(Vec3f(-30, -5, -30), Vec3f(30, -4.5, 9), grey);
    Sphere s1(Vec3f(-3, 0, -16), 2, red);
    Sphere s2(Vec3f(-1, -1.5, -12), 2.1, green);
	Cuboid c1(Vec3f(7, 0,  -15), Vec3f(10, -7,  -10), green);
    Cuboid c2(Vec3f(-7, 0,  -15), Vec3f(-10, -7,  -10), blue);
  
    scene.add_object(&base);
    scene.add_object(&s1);
    scene.add_object(&s2);
    scene.add_object(&c1);
    scene.add_object(&c2);
	
    
    Light l1(Vec3f(-20, 20, 20), 1500);
    Light l2(Vec3f(20, 30, 20), 1800);
    
    scene.add_light(&l1);
    scene.add_light(&l2);
    
    render(scene);
    return 0;
}