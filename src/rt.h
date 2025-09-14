#ifndef RT_H
#define RT_H

#include "util.h"

#include <cglm/cglm.h>
#include <stdatomic.h>
#include <stdint.h>

typedef struct ray {
    vec3 orig;
    vec3 dir;
} ray;

void ray_at(ray r, double t, vec3 out_pos);

typedef struct hit_record hit_record;

typedef struct material {
    bool (*scatter)(ray r, hit_record h, color out_attenuation,
                    ray *out_scattered, void *self);
    void *self;
} material;

typedef struct hit_record {
    vec3 p;
    vec3 norm;
    material mat;
    double t;
    bool front_face;
} hit_record;

typedef struct sphere {
    vec3 center;
    double radius;
    material mat;
} sphere;

typedef struct world {
    sphere *spheres;
    size_t num_spheres;
} world;

typedef struct lambertian {
    color albedo;
} lambertian;

material lambertian_material(lambertian *l);

typedef struct metal {
    color albedo;
    double fuzziness;
} metal;

material metal_material(metal *l);

typedef struct dielectric {
    double refraction_idx;
} dielectric;

material dielectric_material(dielectric *l);

typedef struct pixel {
    uint8_t r, g, b;
} pixel;

typedef struct camera {
    int image_width;
    int image_height;

    vec3 lookfrom;
    vec3 lookat;
    vec3 center;
    vec3 pixel00;
    vec3 pixel_du;
    vec3 pixel_dv;

    vec3 defocus_disk_u;
    vec3 defocus_disk_v;
    double defocus_angle;

    int samples_per_px;
    int max_depth;
} camera;

pixel render_pixel(camera cam, world w, int x, int y);

typedef struct camera_desc {
    int image_width;
    double aspect_ratio;

    vec3 lookfrom;
    vec3 lookat;

    double defocus_angle;
    double focus_dist;

    double fovy;

    int samples_per_px;
    int max_depth;
} camera_desc;

camera camera_init(camera_desc d);

#endif
