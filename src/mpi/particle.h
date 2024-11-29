#ifndef PARTICLE_H
#define PARTICLE_H

typedef struct {
    float x, y;         // position
    float vx, vy;       // velocity
    float fx, fy;       // force
    float rho, p;       // density, pressure
} Particle;

/* Initialize a particle at position (x, y) with default attributes */
void init_particle(Particle *p, float x, float y);

/* Utility function to generate a random float in range [a, b] */
float randab(float a, float b);

#endif
