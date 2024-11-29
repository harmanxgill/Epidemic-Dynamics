#include "sph.h"
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <stddef.h>

/* Global Variables */
Particle *sph_particles = NULL;
int sph_num_particles = 0;

/**
 * Initialize the SPH system with random distribution of humans, zombies, and immune particles.
 */
/**
 * Initialize the SPH system with random distribution of humans, zombies, and immune particles.
 */
void init_zombie_simulation(int num_humans, int num_zombies, int num_immune, int rank, int size) {
    sph_num_particles = 0;

    // Total particles per process
    int total_particles = num_humans + num_zombies + num_immune;
    int particles_per_rank = total_particles / size;
    int remainder = total_particles % size;

    int start_index = rank * particles_per_rank + (rank < remainder ? rank : remainder);
    int end_index = start_index + particles_per_rank + (rank < remainder ? 1 : 0);

    for (int i = start_index; i < end_index; i++) {
        float x = (rand() / (float)RAND_MAX) * (VIEW_WIDTH - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
        float y = (rand() / (float)RAND_MAX) * (VIEW_HEIGHT - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;

        if (i < num_humans) {
            init_particle(&sph_particles[sph_num_particles++], x, y, HUMAN);
        } else if (i < num_humans + num_zombies) {
            init_particle(&sph_particles[sph_num_particles++], x, y, ZOMBIE);
        } else {
            init_particle(&sph_particles[sph_num_particles++], x, y, IMMUNE);
        }
    }

    printf("Rank %d initialized %d particles (start: %d, end: %d).\n", rank, sph_num_particles, start_index, end_index);
}

/**
 * Synchronize boundary particles between ranks.
 */
void synchronize_particles(int rank, int size) {
    int counts[size], displs[size];

    MPI_Allgatherv(
        &sph_particles[0], sph_num_particles, MPI_PARTICLES,
        sph_particles, counts, displs, MPI_PARTICLES,
        MPI_COMM_WORLD
    );
}


void move_zombies() {
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == ZOMBIE) {
            Particle *zombie = &sph_particles[i];
            float closest_dist = VIEW_WIDTH * VIEW_HEIGHT;
            Particle *closest_human = NULL;

            // Find the nearest human
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - zombie->x;
                    float dy = sph_particles[j].y - zombie->y;
                    float dist = dx * dx + dy * dy;
                    if (dist < closest_dist) {
                        closest_dist = dist;
                        closest_human = &sph_particles[j];
                    }
                }
            }

            // Move toward the nearest human
            if (closest_human) {
                float dx = closest_human->x - zombie->x;
                float dy = closest_human->y - zombie->y;
                float length = hypot(dx, dy);
                zombie->vx = (dx / length) * 230.0f; // Adjust speed
                zombie->vy = (dy / length) * 230.0f;
            }
        }
    }
}

void move_immune() {
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == IMMUNE) {
            Particle *immune = &sph_particles[i];
            float closest_dist = VIEW_WIDTH * VIEW_HEIGHT;
            Particle *closest_human = NULL;

            // Find the nearest human
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - immune->x;
                    float dy = sph_particles[j].y - immune->y;
                    float dist = dx * dx + dy * dy;
                    if (dist < closest_dist) {
                        closest_dist = dist;
                        closest_human = &sph_particles[j];
                    }
                }
            }

            // Move toward the nearest human
            if (closest_human) {
                float dx = closest_human->x - immune->x;
                float dy = closest_human->y - immune->y;
                float length = hypot(dx, dy);
                immune->vx = (dx / length) * 230.0f; // Adjust speed (slower than zombies)
                immune->vy = (dy / length) * 230.0f;
            }
        }
    }
}

void move_humans() {
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == HUMAN) {
            Particle *human = &sph_particles[i];
            float dx = 0, dy = 0;

            // Move away from nearby zombies
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == ZOMBIE) {
                    float dist_x = human->x - sph_particles[j].x;
                    float dist_y = human->y - sph_particles[j].y;
                    float dist = hypot(dist_x, dist_y);
                    if (dist < 100.0f) { // Threshold distance
                        dx += dist_x / dist;
                        dy += dist_y / dist;
                    }
                }
            }

            // Normalize and apply velocity
            float length = hypot(dx, dy);
            if (length > 0) {
                human->vx = (dx / length) * 200.0f; // Adjust speed (faster than zombies)
                human->vy = (dy / length) * 200.0f;
            }
        }
    }
}


/**
 * Compute density and pressure for each particle.
 */
/**
 * Compute density and pressure for each particle.
 */
void compute_density_and_pressure(int rank, int size) {
    const float h_squared = SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS;
    const float poly6_constant = 4.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 8.0f));

    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p_i = &sph_particles[i];
        p_i->rho = 0.0f;

        for (int j = 0; j < sph_num_particles; j++) {
            Particle *p_j = &sph_particles[j];

            float dx = p_j->x - p_i->x;
            float dy = p_j->y - p_i->y;
            float distance_squared = dx * dx + dy * dy;

            if (distance_squared < h_squared) {
                p_i->rho += SPH_MASS * poly6_constant * pow(h_squared - distance_squared, 3.0f);
            }
        }

        p_i->p = SPH_GAS_CONSTANT * (p_i->rho - SPH_REST_DENSITY);
    }

    // Synchronize particles after density update
    synchronize_particles(rank, size);
}



int is_in_domain(float x, float y) {
    return ((x < VIEW_WIDTH - SPH_KERNEL_RADIUS) &&
            (x > SPH_KERNEL_RADIUS) &&
            (y < VIEW_HEIGHT - SPH_KERNEL_RADIUS) &&
            (y > SPH_KERNEL_RADIUS));
}

/**
 * Compute forces for each particle.
 */
vvoid compute_sph_forces(int rank, int size) {
    const float spiky_gradient = -10.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float visc_laplacian = 40.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float small_eps = 1e-6f;

    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p_i = &sph_particles[i];
        float pressure_force_x = 0.0f, pressure_force_y = 0.0f;
        float viscosity_force_x = 0.0f, viscosity_force_y = 0.0f;

        for (int j = 0; j < sph_num_particles; j++) {
            Particle *p_j = &sph_particles[j];
            if (p_i == p_j) continue;

            float dx = p_j->x - p_i->x;
            float dy = p_j->y - p_i->y;
            float dist = hypotf(dx, dy) + small_eps;

            if (dist < SPH_KERNEL_RADIUS) {
                float weight = SPH_MASS / p_j->rho;
                float spiky_term = pow(SPH_KERNEL_RADIUS - dist, 3.0f);

                // Pressure force
                float pressure_term = (p_i->p + p_j->p) * spiky_gradient * spiky_term;
                pressure_force_x += -weight * pressure_term * (dx / dist);
                pressure_force_y += -weight * pressure_term * (dy / dist);

                // Viscosity force
                float visc_term = (SPH_KERNEL_RADIUS - dist) * visc_laplacian;
                viscosity_force_x += visc_term * (p_j->vx - p_i->vx) * weight;
                viscosity_force_y += visc_term * (p_j->vy - p_i->vy) * weight;
            }
        }

        // Set total forces
        p_i->fx = pressure_force_x + viscosity_force_x;
        p_i->fy = pressure_force_y + viscosity_force_y;
    }

    // Synchronize particles after forces update
    synchronize_particles(rank, size);
}


/**
 * Integrate particle positions and velocities over time.
 */
void integrate_sph(int rank, int size) {
    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p = &sph_particles[i];

        // Update velocity
        p->vx += SPH_TIME_STEP * p->fx / SPH_MASS;
        p->vy += SPH_TIME_STEP * p->fy / SPH_MASS;

        // Update position
        p->x += SPH_TIME_STEP * p->vx;
        p->y += SPH_TIME_STEP * p->vy;

        // Enforce boundary conditions
        if (p->x < SPH_KERNEL_RADIUS) p->x = SPH_KERNEL_RADIUS;
        if (p->x > VIEW_WIDTH - SPH_KERNEL_RADIUS) p->x = VIEW_WIDTH - SPH_KERNEL_RADIUS;
        if (p->y < SPH_KERNEL_RADIUS) p->y = SPH_KERNEL_RADIUS;
        if (p->y > VIEW_HEIGHT - SPH_KERNEL_RADIUS) p->y = VIEW_HEIGHT - SPH_KERNEL_RADIUS;
    }

    // Synchronize particles after integration
    synchronize_particles(rank, size);
}



void handle_collisions(void) {
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == ZOMBIE) {
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - sph_particles[i].x;
                    float dy = sph_particles[j].y - sph_particles[i].y;
                    float dist_squared = dx * dx + dy * dy;

                    if (dist_squared < SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS) {
                        // Turn human into zombie
                        sph_particles[j].type = ZOMBIE;
                    }
                }
            }
        } else if (sph_particles[i].type == IMMUNE) {
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - sph_particles[i].x;
                    float dy = sph_particles[j].y - sph_particles[i].y;
                    float dist_squared = dx * dx + dy * dy;

                    if (dist_squared < SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS) {
                        // Turn human into immune
                        sph_particles[j].type = IMMUNE;
                    }
                }
            }
        }
    }
}



/**
 * Update the SPH system by computing densities, forces, and integrating.
 */
void update_sph(int rank, int size) {
    compute_density_and_pressure(rank, size);
    compute_sph_forces(rank, size);
    move_humans();
    move_zombies();
    move_immune();
    handle_collisions();
    integrate_sph(rank, size);
}
