#ifndef UTIL_H
#define UTIL_H

#include <SDL3/SDL_surface.h>
#include <cglm/cglm.h>

#define PI 3.14159265358979323846

typedef vec3 color;

extern const double infinity;

int maxi(int a, int b);

void seed_rng(uint32_t seed);
uint32_t random_u32(void);
double random_double(void);
void random_vec3(vec3 out);
void random_normalized_vec3(vec3 out);
void random_in_unit_disc_vec3(vec3 out);
bool vec3_is_near_zero(vec3 v);

typedef struct interval {
    double min, max;
} interval;

extern interval empty;
extern interval universe;

double interval_size(interval i);
bool interval_contains(interval i, double x);
bool interval_sorrounds(interval i, double x);
double interval_clamp(interval i, double x);

#endif
