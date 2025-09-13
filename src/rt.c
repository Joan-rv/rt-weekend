#include "rt.h"

void ray_at(ray r, double t, vec3 out_pos) {
    glm_vec3_scale(r.dir, t, out_pos);
    glm_vec3_add(r.orig, out_pos, out_pos);
}

static void hit_set_front_face(ray r, hit_record *rec) {
    rec->front_face = glm_vec3_dot(r.dir, rec->norm) < 0.0;
    if (!rec->front_face)
        glm_vec3_negate(rec->norm);
}

static bool sphere_hit(ray r, interval t, hit_record *out, void *self) {
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
    ray_at(r, root, out->p);
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

static bool hittable_list_hit(ray r, interval t, hit_record *out, void *self) {
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

static bool lambertian_scatter(ray r, hit_record h, color out_attenuation,
                               ray *out_scattered, void *self) {
    lambertian *s = self;
    random_normalized_vec3(out_scattered->dir);
    glm_vec3_add(r.dir, h.norm, r.dir);
    if (vec3_is_near_zero(r.dir))
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

static bool metal_scatter(ray r, hit_record h, color out_attenuation,
                          ray *out_scattered, void *self) {
    metal *s = self;
    glm_vec3_reflect(r.dir, h.norm, out_scattered->dir);
    vec3 fuzz;
    random_normalized_vec3(fuzz);
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

static double reflectance(double cosine, double refraction_idx) {
    // Schlick's approximation for reflectance
    double r0 = (1 - refraction_idx) / (1 + refraction_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow(1 - cosine, 5);
}

static bool dielectric_scatter(ray r, hit_record h, color out_attenuation,
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

static void ray_color(ray r, int depth, hittable world, color out) {
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

static void write_color(int width, int x, int y, const color color,
                        pixel out_buf[]) {
    interval i = {0.0, 0.999};
    pixel p;
    p.r = 256.0 * interval_clamp(i, linear_to_gamma(color[0]));
    p.g = 256.0 * interval_clamp(i, linear_to_gamma(color[1]));
    p.b = 256.0 * interval_clamp(i, linear_to_gamma(color[2]));
    out_buf[width * y + x] = p;
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
        vec3 tmp;
        random_in_unit_disc_vec3(tmp);
        glm_vec3_muladds(cam.defocus_disk_u, tmp[0], ret.orig);
        glm_vec3_muladds(cam.defocus_disk_v, tmp[1], ret.orig);
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

void camera_render(camera cam, hittable world, pixel out_buf[]) {
    // printf("P3\n%d %d\n255\n", cam.image_width, cam.image_height);
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
            write_color(cam.image_width, j, i, color, out_buf);
        }
    }
    fprintf(stderr, "\rDone                              \n");
}
camera camera_init(camera_desc d) {
    camera ret;
    vec3 tmp;
    ret.image_width = d.image_width;
    ret.samples_per_px = d.samples_per_px;
    ret.max_depth = d.max_depth;
    ret.defocus_angle = d.defocus_angle;
    ret.image_height = maxi(1, ret.image_width / d.aspect_ratio);
    glm_vec3_copy(d.lookfrom, ret.center);

    double h = tan(d.fovy / 2);
    double viewport_height = 2.0 * h * d.focus_dist;
    double viewport_width =
        viewport_height * (double)ret.image_width / ret.image_height;

    vec3 u, v, w, vup;
    glm_vec3_copy((vec3){0.0, 1.0, 0.0}, vup);
    glm_vec3_sub(d.lookfrom, d.lookat, w);
    glm_normalize(w);
    glm_vec3_cross(vup, w, u);
    glm_vec3_normalize(u);
    glm_vec3_cross(w, u, v);

    vec3 viewport_u, viewport_v;
    glm_vec3_scale(u, viewport_width, viewport_u);
    glm_vec3_scale(v, viewport_height, viewport_v);
    glm_vec3_negate(viewport_v);
    glm_vec3_scale(viewport_u, 1.0 / ret.image_width, ret.pixel_du);
    glm_vec3_scale(viewport_v, 1.0 / ret.image_height, ret.pixel_dv);

    vec3 viewport_upper_left;
    glm_vec3_copy(ret.center, viewport_upper_left);
    glm_vec3_scale(w, d.focus_dist, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);
    glm_vec3_scale(viewport_u, 0.5, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);
    glm_vec3_scale(viewport_v, 0.5, tmp);
    glm_vec3_sub(viewport_upper_left, tmp, viewport_upper_left);

    glm_vec3_add(ret.pixel_du, ret.pixel_dv, tmp);
    glm_vec3_scale(tmp, 0.5, tmp);
    glm_vec3_add(viewport_upper_left, tmp, ret.pixel00);

    double defocus_radius = d.focus_dist * tan(d.defocus_angle / 2.0);
    glm_vec3_scale(u, defocus_radius, ret.defocus_disk_u);
    glm_vec3_scale(v, defocus_radius, ret.defocus_disk_v);
    return ret;
}
