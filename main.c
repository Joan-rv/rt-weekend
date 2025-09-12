#include <cglm/cglm.h>

#include <float.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846

typedef vec3 color;

const double infinity = DBL_MAX;

static int maxi(int a, int b) { return a < b ? b : a; }

double random_double(void) { return rand() / (RAND_MAX + 1.0); }

static void random_vec3(vec3 out) {
    out[0] = 2.0 * random_double() - 1.0;
    out[1] = 2.0 * random_double() - 1.0;
    out[2] = 2.0 * random_double() - 1.0;
}

static void random_unit_vec3(vec3 out) {
    while (true) {
        random_vec3(out);
        double lensq = glm_vec3_norm2(out);
        // TODO: better value for rejecting small lensq
        if (1e-15 < lensq && lensq <= 1) {
            glm_vec3_scale(out, 1.0 / sqrt(lensq), out);
            return;
        }
    }
}

static void random_in_unit_disc_vec3(vec3 out) {
    while (true) {
        random_vec3(out);
        if (glm_vec3_norm2(out) < 1)
            return;
    }
}

static bool vec3_near_zero(vec3 v) {
    double s = 1e-8;
    return fabs(v[0]) < s && fabs(v[1]) < s && fabs(v[2]) < s;
}

typedef struct interval {
    double min, max;
} interval;

interval empty = {infinity, -infinity};
interval universe = {-infinity, infinity};

double interval_size(interval i) { return i.max - i.min; }
bool interval_contains(interval i, double x) {
    return i.min <= x && x <= i.max;
}
bool interval_sorrounds(interval i, double x) { return i.min < x && x < i.max; }
double interval_clamp(interval i, double x) {
    if (x < i.min)
        return i.min;
    if (x > i.max)
        return i.max;
    return x;
}

typedef struct ray {
    vec3 orig;
    vec3 dir;
} ray;

void ray_at(double t, ray r, vec3 out_pos) {
    glm_vec3_scale(r.dir, t, out_pos);
    glm_vec3_add(r.orig, out_pos, out_pos);
}

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

void hit_set_front_face(ray r, hit_record *rec) {
    rec->front_face = glm_vec3_dot(r.dir, rec->norm) < 0.0;
    if (!rec->front_face)
        glm_vec3_scale(rec->norm, -1.0, rec->norm);
}

typedef struct hittable {
    bool (*hit)(ray r, interval valid_t, hit_record *out, void *self);
    void *self;
} hittable;

typedef struct sphere {
    vec3 center;
    double radius;
    material mat;
} sphere;

bool sphere_hit(ray r, interval t, hit_record *out, void *self) {
    sphere *s = self;
    vec3 oc;
    glm_vec3_sub(s->center, r.orig, oc);
    double a = glm_vec3_norm2(r.dir);
    double h = glm_vec3_dot(r.dir, oc);
    double c = glm_vec3_norm2(oc) - s->radius * s->radius;
    double discriminant = h * h - a * c;
    if (discriminant < 0.0)
        return false;
    double d_sqrt = sqrt(discriminant);
    double root = (h - d_sqrt) / a;
    if (!interval_sorrounds(t, root)) {
        root = (h + d_sqrt) / a;
        if (!interval_sorrounds(t, root))
            return false;
    }

    out->mat = s->mat;
    ray_at(root, r, out->p);
    glm_vec3_sub(out->p, s->center, out->norm);
    glm_vec3_scale(out->norm, 1.0 / s->radius, out->norm);
    out->t = root;
    hit_set_front_face(r, out);

    return true;
}

hittable sphere_hittable(sphere *s) {
    return (hittable){
        .hit = sphere_hit,
        .self = s,
    };
}

typedef struct hittable_list {
    hittable *objects;
    size_t len;
} hittable_list;

bool hittable_list_hit(ray r, interval t, hit_record *out, void *self) {
    hittable_list *s = self;
    hit_record tmp_rec;
    bool hit_anything = false;
    double closest = t.max;

    for (size_t i = 0; i < s->len; ++i) {
        if (s->objects[i].hit(r, (interval){t.min, closest}, &tmp_rec,
                              s->objects[i].self)) {
            hit_anything = true;
            closest = tmp_rec.t;
            *out = tmp_rec;
        }
    }

    return hit_anything;
}

hittable hittable_list_hittable(hittable_list *h) {
    return (hittable){
        .hit = hittable_list_hit,
        .self = h,
    };
}

typedef struct lambertian {
    color albedo;
} lambertian;

bool lambertian_scatter(ray r, hit_record h, color out_attenuation,
                        ray *out_scattered, void *self) {
    lambertian *s = self;
    random_unit_vec3(out_scattered->dir);
    glm_vec3_add(r.dir, h.norm, r.dir);
    if (vec3_near_zero(r.dir))
        glm_vec3_copy(h.norm, r.dir);
    glm_vec3_copy(h.p, out_scattered->orig);
    glm_vec3_copy(s->albedo, out_attenuation);
    return true;
}

material lambertian_material(lambertian *l) {
    return (material){
        .scatter = lambertian_scatter,
        .self = l,
    };
}

typedef struct metal {
    color albedo;
    double fuzziness;
} metal;

bool metal_scatter(ray r, hit_record h, color out_attenuation,
                   ray *out_scattered, void *self) {
    metal *s = self;
    glm_vec3_reflect(r.dir, h.norm, out_scattered->dir);
    vec3 fuzz;
    random_unit_vec3(fuzz);
    glm_vec3_scale(fuzz, s->fuzziness, fuzz);
    glm_vec3_add(out_scattered->dir, fuzz, out_scattered->dir);
    glm_vec3_copy(h.p, out_scattered->orig);
    glm_vec3_copy(s->albedo, out_attenuation);
    return true;
}

material metal_material(metal *m) {
    return (material){
        .scatter = metal_scatter,
        .self = m,
    };
}

typedef struct dielectric {
    double refraction_idx;
} dielectric;

static double reflectance(double cosine, double refraction_idx) {
    // Schlick's approximation for reflectance
    double r0 = (1 - refraction_idx) / (1 + refraction_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow(1 - cosine, 5);
}

bool dielectric_scatter(ray r, hit_record h, color out_attenuation,
                        ray *out_scattered, void *self) {
    dielectric *s = self;
    glm_vec3_copy((color){1.0, 1.0, 1.0}, out_attenuation);
    double refraction_idx =
        h.front_face ? 1.0 / s->refraction_idx : s->refraction_idx;
    glm_vec3_normalize(r.dir);

    vec3 tmp;
    glm_vec3_negate_to(r.dir, tmp);
    double cos_theta = fmin(glm_vec3_dot(tmp, h.norm), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    bool cannot_refract = refraction_idx * sin_theta > 1.0;
    if (cannot_refract ||
        reflectance(cos_theta, refraction_idx) > random_double())
        glm_vec3_reflect(r.dir, h.norm, out_scattered->dir);
    else
        glm_vec3_refract(r.dir, h.norm, refraction_idx, out_scattered->dir);

    glm_vec3_copy(h.p, out_scattered->orig);
    return true;
}

material dielectric_material(dielectric *d) {
    return (material){
        .scatter = dielectric_scatter,
        .self = d,
    };
}

void ray_color(ray r, int depth, hittable world, color out) {
    if (depth <= 0) {
        glm_vec3_copy((color){0.0, 0.0, 0.0}, out);
        return;
    }
    hit_record h;
    if (world.hit(r, (interval){0.001, infinity}, &h, world.self)) {
        ray scattered;
        color attenuation;
        if (h.mat.scatter(r, h, attenuation, &scattered, h.mat.self)) {
            ray_color(scattered, depth - 1, world, out);
            glm_vec3_mul(out, attenuation, out);
            return;
        }
        glm_vec3_copy((color){0.0, 0.0, 0.0}, out);
        return;
    }

    glm_vec3_normalize(r.dir);
    double a = 0.5 * (r.dir[1] + 1.0);
    color start = {1.0, 1.0, 1.0};
    color end = {0.5, 0.7, 1.0};
    glm_vec3_scale(start, 1.0 - a, start);
    glm_vec3_scale(end, a, end);
    glm_vec3_add(start, end, out);
}

static double linear_to_gamma(double linear) {
    if (linear > 0)
        return sqrt(linear);
    return 0.0;
}

void write_color(FILE *out, const color color) {
    interval i = {0.0, 0.999};
    int r = 256.0 * interval_clamp(i, linear_to_gamma(color[0]));
    int g = 256.0 * interval_clamp(i, linear_to_gamma(color[1]));
    int b = 256.0 * interval_clamp(i, linear_to_gamma(color[2]));
    fprintf(out, "%d %d %d\n", r, g, b);
}

typedef struct camera {
    double aspect_ratio;
    double fovy;
    double defocus_angle;
    double focus_dist;
    int image_width;
    int image_height;
    int samples_per_px;
    int max_depth;
    vec3 lookfrom;
    vec3 lookat;
    vec3 vup;
    vec3 center;
    vec3 pixel00;
    vec3 pixel_du;
    vec3 pixel_dv;
    vec3 u, v, w;
    vec3 defocus_disk_u;
    vec3 defocus_disk_v;
} camera;

camera camera_init(vec3 lookfrom, vec3 lookat, double defocus_angle,
                   double focus_dist, int image_width, double aspect_ratio,
                   double fovy, int samples_per_px, int max_depth) {
    camera ret;
    vec3 tmp;
    ret.aspect_ratio = aspect_ratio;
    ret.fovy = fovy;
    ret.image_width = image_width;
    ret.samples_per_px = samples_per_px;
    ret.max_depth = max_depth;
    ret.defocus_angle = defocus_angle;
    ret.focus_dist = focus_dist;
    ret.image_height = maxi(1, ret.image_width / ret.aspect_ratio);
    glm_vec3_copy(lookfrom, ret.center);

    double h = tan(fovy / 2);
    double viewport_height = 2.0 * h * focus_dist;
    double viewport_width =
        viewport_height * (double)ret.image_width / ret.image_height;

    glm_vec3_copy((vec3){0.0, 1.0, 0.0}, ret.vup);
    glm_vec3_sub(lookfrom, lookat, ret.w);
    glm_normalize(ret.w);
    glm_vec3_cross(ret.vup, ret.w, ret.u);
    glm_vec3_normalize(ret.u);
    glm_vec3_cross(ret.w, ret.u, ret.v);

    vec3 viewport_u, viewport_v;
    glm_vec3_scale(ret.u, viewport_width, viewport_u);
    glm_vec3_scale(ret.v, viewport_height, viewport_v);
    glm_vec3_negate(viewport_v);
    glm_vec3_scale(viewport_u, 1.0 / ret.image_width, ret.pixel_du);
    glm_vec3_scale(viewport_v, 1.0 / ret.image_height, ret.pixel_dv);

    vec3 viewport_upper_left;
    glm_vec3_copy(ret.center, viewport_upper_left);
    glm_vec3_scale(ret.w, focus_dist, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);
    glm_vec3_scale(viewport_u, 0.5, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);
    glm_vec3_scale(viewport_v, 0.5, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);

    glm_vec3_add(ret.pixel_du, ret.pixel_dv, tmp);
    glm_vec3_scale(tmp, 0.5, tmp);
    glm_vec3_add(viewport_upper_left, tmp, ret.pixel00);

    double defocus_radius = focus_dist * tan(defocus_angle / 2.0);
    glm_vec3_scale(ret.u, defocus_radius, ret.defocus_disk_u);
    glm_vec3_scale(ret.v, defocus_radius, ret.defocus_disk_v);
    return ret;
}

static void random_square_offset(vec3 out) {
    out[0] = random_double() - 0.5;
    out[1] = random_double() - 0.5;
    out[2] = 0.0;
}

static ray camera_get_ray(camera cam, int i, int j) {
    ray ret;
    vec3 offset;
    random_square_offset(offset);
    glm_vec3_copy(cam.center, ret.orig);
    if (cam.defocus_angle > 0.0) {
        vec3 tmp, tmp2;
        random_in_unit_disc_vec3(tmp);
        glm_vec3_scale(cam.defocus_disk_u, tmp[0], tmp2);
        glm_vec3_add(ret.orig, tmp2, ret.orig);
        glm_vec3_scale(cam.defocus_disk_v, tmp[1], tmp2);
        glm_vec3_add(ret.orig, tmp2, ret.orig);
    }

    vec3 du;
    glm_vec3_scale(cam.pixel_du, j + offset[0], du);
    vec3 dv;
    glm_vec3_scale(cam.pixel_dv, i + offset[1], dv);
    vec3 px_center;
    glm_vec3_add(cam.pixel00, du, px_center);
    glm_vec3_add(px_center, dv, px_center);
    glm_vec3_sub(px_center, ret.orig, ret.dir);
    return ret;
}

void camera_render(camera cam, hittable world) {
    printf("P3\n%d %d\n255\n", cam.image_width, cam.image_height);
    for (int i = 0; i < cam.image_height; ++i) {
        fprintf(stderr, "\rScanlines remaining: %5d", cam.image_height - i);
        fflush(stderr);
        for (int j = 0; j < cam.image_width; ++j) {
            color color = {0.0, 0.0, 0.0};
            for (int sample = 0; sample < cam.samples_per_px; ++sample) {
                vec3 sample_color;
                ray ray = camera_get_ray(cam, i, j);
                ray_color(ray, cam.max_depth, world, sample_color);
                glm_vec3_add(color, sample_color, color);
            }
            glm_vec3_scale(color, 1.0 / cam.samples_per_px, color);
            write_color(stdout, color);
        }
    }
    fprintf(stderr, "\rDone                              \n");
}

int sample_scene(void) {
    // NOTE: we don't seed `rand` for more predictable results.

    camera cam = camera_init((vec3){-2.0, 2.0, 1.0}, (vec3){0.0, 0.0, -1.0},
                             PI * 10.0 / 180.0, 3.4, 400, 16.0 / 9.0,
                             PI * 20.0 / 180.0, 100, 50);

    lambertian mat_ground = {.albedo = {0.8, 0.8, 0.0}};
    lambertian mat_center = {.albedo = {0.1, 0.2, 0.5}};
    dielectric mat_left = {1.5};
    dielectric mat_bubble = {1.0 / 1.5};
    metal mat_right = {.albedo = {0.8, 0.6, 0.2}, .fuzziness = 1.0};

    sphere sphere1 = {
        .center = {0.0, 0.0, -1.0},
        .radius = 0.5,
        .mat = lambertian_material(&mat_center),
    };
    sphere sphere2 = {
        .center = {0.0, -100.5, -1.0},
        .radius = 100.0,
        .mat = lambertian_material(&mat_ground),
    };
    sphere sphere3 = {
        .center = {-1.0, 0.0, -1.0},
        .radius = 0.5,
        .mat = dielectric_material(&mat_left),
    };
    sphere sphere4 = {
        .center = {-1.0, 0.0, -1.0},
        .radius = 0.4,
        .mat = dielectric_material(&mat_bubble),
    };
    sphere sphere5 = {
        .center = {1.0, 0.0, -1.0},
        .radius = 0.5,
        .mat = metal_material(&mat_right),
    };
    hittable objects[] = {sphere_hittable(&sphere1), sphere_hittable(&sphere2),
                          sphere_hittable(&sphere3), sphere_hittable(&sphere4),
                          sphere_hittable(&sphere5)};
    hittable_list world = {
        .objects = objects,
        .len = sizeof(objects) / sizeof(*objects),
    };

    camera_render(cam, hittable_list_hittable(&world));
    return 0;
}

int cover_scene(void) {
    sphere spheres[22 * 22 + 4];
    int sphere_idx = 0;
    hittable objects[22 * 22 + 4];
    int object_idx = 0;
    lambertian lambertians[22 * 22 + 4];
    int lambertian_idx = 0;
    metal metals[22 * 22 + 4];
    int metal_idx = 0;
    dielectric dielectrics[22 * 22 + 4];
    int dielectric_idx = 0;

    lambertians[lambertian_idx++] = (lambertian){.albedo = {0.5, 0.5, 0.5}};
    spheres[sphere_idx++] = (sphere){
        .center = {0.0, -1000.0, 0.0},
        1000.0,
        lambertian_material(&lambertians[lambertian_idx - 1]),
    };

    for (int a = -11; a < 11; ++a) {
        for (int b = -11; b < 11; ++b) {
            double choose_mat = random_double();
            vec3 center = {a + 0.9 * random_double(), 0.2,
                           b + 0.9 * random_double()};
            if (glm_vec3_distance(center, (vec3){4.0, 0.2, 0.0}) <= 0.9) {
                continue;
            }

            if (choose_mat < 0.8) {
                color albedo;
                random_vec3(albedo);
                glm_vec3_mul(albedo, albedo, albedo);
                glm_vec3_copy(albedo, lambertians[lambertian_idx++].albedo);
                glm_vec3_copy(center, spheres[sphere_idx].center);
                spheres[sphere_idx].radius = 0.2;
                spheres[sphere_idx].mat =
                    lambertian_material(&lambertians[lambertian_idx - 1]);
                sphere_idx++;
            } else if (choose_mat < 0.95) {
                color albedo;
                random_vec3(albedo);
                glm_vec3_scale(albedo, 0.5, albedo);
                glm_vec3_add(albedo, (vec3){0.5, 0.5, 0.5}, albedo);
                glm_vec3_copy(albedo, metals[metal_idx].albedo);
                metals[metal_idx++].fuzziness = random_double() / 2.0;
                glm_vec3_copy(center, spheres[sphere_idx].center);
                spheres[sphere_idx].radius = 0.2;
                spheres[sphere_idx].mat =
                    metal_material(&metals[metal_idx - 1]);
                sphere_idx++;
            } else {
                dielectrics[dielectric_idx++].refraction_idx = 1.5;
                glm_vec3_copy(center, spheres[sphere_idx].center);
                spheres[sphere_idx].radius = 0.2;
                spheres[sphere_idx].mat =
                    dielectric_material(&dielectrics[dielectric_idx - 1]);
                sphere_idx++;
            }
        }
    }

    dielectrics[dielectric_idx++].refraction_idx = 1.5;
    spheres[sphere_idx++] = (sphere){
        .center = {0.0, 1.0, 0.0},
        .radius = 1.0,
        .mat = dielectric_material(&dielectrics[dielectric_idx - 1]),
    };

    glm_vec3_copy((vec3){0.4, 0.2, 0.1}, lambertians[lambertian_idx++].albedo);
    spheres[sphere_idx++] = (sphere){
        .center = {-4.0, 1.0, 0.0},
        .radius = 1.0,
        .mat = lambertian_material(&lambertians[lambertian_idx - 1]),
    };

    glm_vec3_copy((vec3){0.7, 0.6, 0.5}, metals[metal_idx].albedo);
    metals[metal_idx++].fuzziness = 0.0;
    spheres[sphere_idx++] = (sphere){
        .center = {4.0, 1.0, 0.0},
        .radius = 1.0,
        .mat = metal_material(&metals[metal_idx - 1]),
    };

    for (int i = 0; i < sphere_idx; ++i) {
        objects[object_idx++] = sphere_hittable(&spheres[i]);
    }
    hittable_list world = {.objects = objects, .len = object_idx};

    camera cam = camera_init((vec3){13.0, 2.0, 3.0}, (vec3){0.0, 0.0, 0.0},
                             PI * 0.6 / 180.0, 10.0, 1200, 16.0 / 9.0,
                             PI * 20.0 / 180.0, 500, 50);
    camera_render(cam, hittable_list_hittable(&world));

    return 0;
}

int main(void) { return cover_scene(); }
