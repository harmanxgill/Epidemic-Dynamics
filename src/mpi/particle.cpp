#include "particle.h"
#include <stdlib.h>

/**
 * Initialize a particle with given position and type.
 * - `p`: Pointer to the particle to initialize.
 * - `x`, `y`: Position of the particle.
 * - `type`: The type of the particle (HUMAN, ZOMBIE, IMMUNE).
 */
void init_particle(Particle *p, float x, float y, ParticleType type) {
    p->x = x;
    p->y = y;
    p->vx = 0.0f;
    p->vy = 0.0f;
    p->fx = 0.0f;
    p->fy = 0.0f;
    p->rho = 0.0f;
    p->p = 0.0f;
    p->type = type;

    if (type == ZOMBIE) {
        p->health = 100; // Initial health for zombies
    } else {
        p->health = 0;   // No health for other particles
    }
}

/**
 * Generate a random float in the range [a, b].
 * - `a`: Lower bound.
 * - `b`: Upper bound.
 * - Returns: A random float between `a` and `b`.
 */
float randab(float a, float b) {
    return a + (b - a) * rand() / (float)RAND_MAX;
}
