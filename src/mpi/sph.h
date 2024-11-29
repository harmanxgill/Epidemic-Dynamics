#ifndef SPH_H
#define SPH_H

#ifdef GUI
#define MAX_PARTICLES 5000
#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 768
#else
#define MAX_PARTICLES 20000
#define WINDOW_WIDTH 3000
#define WINDOW_HEIGHT 2000
#endif

#define DAM_PARTICLES 500

#define VIEW_WIDTH (1.5 * WINDOW_WIDTH)
#define VIEW_HEIGHT (1.5 * WINDOW_HEIGHT)
#define SPH_KERNEL_RADIUS 16.0f

#include "particle.h"

// Declare SPH-related functions and variables
extern Particle *sph_particles;
extern int sph_num_particles;

void init_sph(int num_particles);
void compute_density_and_pressure(void);
void compute_sph_forces(void);
void integrate_sph(void);
float calculate_avg_velocity(void);
void update_sph(void);
int is_in_domain(float x, float y);

#endif