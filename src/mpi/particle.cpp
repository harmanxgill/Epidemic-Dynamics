#include "particle.h"
#include <stdlib.h>

void init_particle(Particle *p, float x, float y) {
    p->x = x;
    p->y = y;
    p->vx = p->vy = 0.0;
    p->fx = p->fy = 0.0;
    p->rho = 0.0;
    p->p = 0.0;
}

float randab(float a, float b) {
    return a + (b - a) * rand() / (float)(RAND_MAX);
}
