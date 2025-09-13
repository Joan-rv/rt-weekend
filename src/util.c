#include "util.h"

#include <float.h>

const double infinity = DBL_MAX;

int maxi(int a, int b) { return a < b ? b : a; }

double random_double(void) { return rand() / (RAND_MAX + 1.0); }

void random_vec3(vec3 out) {
    out[0] = 2.0 * random_double() - 1.0;
    out[1] = 2.0 * random_double() - 1.0;
    out[2] = 2.0 * random_double() - 1.0;
}

void random_normalized_vec3(vec3 out) {
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

void random_in_unit_disc_vec3(vec3 out) {
    while (true) {
        random_vec3(out);
        if (glm_vec3_norm2(out) < 1)
            return;
    }
}

bool vec3_is_near_zero(vec3 v) {
    double s = 1e-8;
    return fabs(v[0]) < s && fabs(v[1]) < s && fabs(v[2]) < s;
}

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
