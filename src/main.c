#include "rt.h"
#include "util.h"

#include <SDL3/SDL.h>
#include <stb_image_write.h>

#include <pthread.h>
#include <stdatomic.h>
#include <sys/sysinfo.h>

const int width = 2000;
const double aspect_ratio = 16.0 / 9.0;

atomic_bool rt_running = false;

typedef struct render_range {
    int px_start, px_end;
} render_range;

void render_scene_range(camera cam, hittable world, pixel out_buf[],
                        render_range r) {
    for (int px = r.px_start; px < r.px_end && rt_running; ++px) {
        int x = px % cam.image_width;
        int y = px / cam.image_width;
        out_buf[cam.image_width * y + x] = render_pixel(cam, world, x, y);
    }
}

typedef struct render_thread_args {
    camera cam;
    hittable world;
    pixel *out_buf;
    render_range range;
} render_thread_args;

void *render_thread(void *data) {
    render_thread_args *args = data;
    seed_rng(1);
    render_scene_range(args->cam, args->world, args->out_buf, args->range);
    return NULL;
}

void render_scene(camera cam, hittable world, pixel out_buf[]) {
    int nproc = get_nprocs();
    int num_pixels = cam.image_width * cam.image_height;
    int px_per_thread = num_pixels / nproc;
    render_thread_args thread_args[nproc];
    pthread_t threads[nproc];
    for (int i = 0; i < nproc; ++i) {
        int px_start = i * px_per_thread;
        int px_end = (i + 1) * px_per_thread;
        thread_args[i] = (render_thread_args){
            .cam = cam,
            .world = world,
            .out_buf = out_buf,
            .range =
                (render_range){
                    .px_start = px_start,
                    .px_end = (i + 1 == nproc ? num_pixels : px_end),
                },
        };
        assert(pthread_create(&threads[i], NULL, render_thread,
                              &thread_args[i]) == 0);
    }

    for (int i = 0; i < nproc; ++i) {
        assert(pthread_join(threads[i], NULL) == 0);
    }
}

void render_scene2(camera cam, hittable world, pixel out_buf[]) {
    for (int y = 0; y < cam.image_height && rt_running; ++y) {
        for (int x = 0; x < cam.image_width && rt_running; ++x) {
            out_buf[cam.image_width * y + x] = render_pixel(cam, world, x, y);
        }
    }
}

int sample_scene(pixel out_buf[]) {
    srand(1);

    camera cam = camera_init((camera_desc){
        .lookfrom = {-2.0, 2.0, 1.0},
        .lookat = {0.0, 0.0, -1.0},
        .defocus_angle = PI * 10.0 / 180.0,
        .focus_dist = 3.4,
        .image_width = width,
        .aspect_ratio = aspect_ratio,
        .fovy = PI * 20.0 / 180.0,
        .samples_per_px = 100,
        .max_depth = 50,
    });

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

    render_scene(cam, hittable_list_hittable(&world), out_buf);
    return 0;
}

int cover_scene(pixel out_buf[]) {
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

    camera cam = camera_init((camera_desc){
        .lookfrom = {13.0, 2.0, 3.0},
        .lookat = {0.0, 0.0, 0.0},
        .defocus_angle = PI * 0.6 / 180.0,
        .focus_dist = 10.0,
        .image_width = width,
        .aspect_ratio = aspect_ratio,
        .fovy = PI * 20.0 / 180.0,
        .samples_per_px = 500,
        .max_depth = 50,
    });
    render_scene(cam, hittable_list_hittable(&world), out_buf);

    return 0;
}

void *rt_run(void *data) {
    sample_scene(data);
    rt_running = false;
    return NULL;
}

int main(int argc, char **argv) {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init: %s", SDL_GetError());
        return 1;
    }
    SDL_Window *window;
    if (!(window = SDL_CreateWindow("rt", width, width / aspect_ratio, 0))) {
        SDL_Log("SDL_CreateWindow: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_Surface *surface = SDL_GetWindowSurface(window);
    assert(surface);
    assert(!SDL_MUSTLOCK(surface));
    assert(SDL_FillSurfaceRect(surface, NULL,
                               SDL_MapSurfaceRGB(surface, 0x00, 0x00, 0x00)));

    pixel *out_buf = calloc(width * width / aspect_ratio, sizeof(pixel));

    rt_running = true;
    pthread_t rt_thrd;
    assert(pthread_create(&rt_thrd, NULL, rt_run, out_buf) == 0);

    while (rt_running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT)
                rt_running = false;
        }

        Uint32 *pixels = surface->pixels;
        for (int i = 0; i < surface->h; ++i) {
            for (int j = 0; j < surface->w; ++j) {
                pixel p = out_buf[i * surface->w + j];
                pixels[i * surface->w + j] =
                    SDL_MapSurfaceRGB(surface, p.r, p.g, p.b);
            }
        }

        assert(SDL_UpdateWindowSurface(window));
    }

    if (argc >= 2) {
        if (!stbi_write_png(argv[1], width, (int)(width / aspect_ratio),
                            sizeof(pixel), out_buf, width * sizeof(pixel)))
            fprintf(stderr, "failed to write image to %s\n", argv[1]);
    }

    assert(pthread_join(rt_thrd, NULL) == 0);

    SDL_DestroyWindow(window);
    SDL_Quit();

    free(out_buf);

    return 0;
}
