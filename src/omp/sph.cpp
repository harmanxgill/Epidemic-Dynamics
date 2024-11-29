#include "sph.h"
#include <math.h>
#include <stdio.h>
#include <omp.h>

/* Global Variables */
Particle *sph_particles = NULL;
int sph_num_particles = 0;

/**
 * Initialize the SPH system with random distribution of humans, zombies, and immune particles.
 */
void init_zombie_simulation(int num_humans, int num_zombies, int num_immune) {
    sph_num_particles = 0;
    int total_particles = num_humans + num_zombies + num_immune;

    // Allocate memory for particles (ensure enough space)
    sph_particles = (Particle *)malloc(total_particles * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        int humans_start = (num_humans * thread_num) / num_threads;
        int humans_end = (num_humans * (thread_num + 1)) / num_threads;
        for (int i = humans_start; i < humans_end; i++) {
            float x = (rand() / (float)RAND_MAX) * (VIEW_WIDTH - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            float y = (rand() / (float)RAND_MAX) * (VIEW_HEIGHT - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            init_particle(&sph_particles[i], x, y, HUMAN);
            sph_particles[i].type = HUMAN;
        }

        int zombies_start = num_humans + (num_zombies * thread_num) / num_threads;
        int zombies_end = num_humans + (num_zombies * (thread_num + 1)) / num_threads;
        for (int i = zombies_start; i < zombies_end; i++) {
            float x = (rand() / (float)RAND_MAX) * (VIEW_WIDTH - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            float y = (rand() / (float)RAND_MAX) * (VIEW_HEIGHT - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            init_particle(&sph_particles[i], x, y, ZOMBIE);
            sph_particles[i].type = ZOMBIE;
        }

        int immune_start = num_humans + num_zombies + (num_immune * thread_num) / num_threads;
        int immune_end = num_humans + num_zombies + (num_immune * (thread_num + 1)) / num_threads;
        for (int i = immune_start; i < immune_end; i++) {
            float x = (rand() / (float)RAND_MAX) * (VIEW_WIDTH - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            float y = (rand() / (float)RAND_MAX) * (VIEW_HEIGHT - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
            init_particle(&sph_particles[i], x, y, IMMUNE);
            sph_particles[i].type = IMMUNE;
        }
    }

    sph_num_particles = total_particles;
    printf("Initialized with %d humans, %d zombies, %d immune particles.\n", num_humans, num_zombies, num_immune);
}



void move_zombies() {
#pragma omp parallel for
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
#pragma omp parallel for
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
                immune->vx = (dx / length) * 230.0f; // Adjust speed
                immune->vy = (dy / length) * 230.0f;
            }
        }
    }
}


void move_humans() {
#pragma omp parallel for
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
void compute_density_and_pressure(void) {
    const float h_squared = SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS; // Kernel radius squared
    const float poly6_constant = 4.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 8.0f));

    // Allocate an array to store densities
    float *rho = (float *)calloc(sph_num_particles, sizeof(*rho));
    if (rho == NULL) {
        fprintf(stderr, "Memory allocation for density array failed.\n");
        exit(EXIT_FAILURE);
    }

    // Parallelize density computation using reduction
#pragma omp parallel default(none) shared(sph_num_particles, sph_particles, rho, h_squared, poly6_constant)
    {
        #pragma omp for reduction(+:rho[:sph_num_particles]) collapse(2)
        for (int i = 0; i < sph_num_particles; i++) {
            for (int j = 0; j < sph_num_particles; j++) {
                float dx = sph_particles[j].x - sph_particles[i].x;
                float dy = sph_particles[j].y - sph_particles[i].y;
                float distance_squared = dx * dx + dy * dy;

                if (distance_squared < h_squared) {
                    rho[i] += SPH_MASS * poly6_constant * pow(h_squared - distance_squared, 3.0f);
                }
            }
        }

        // Parallelize pressure calculation and particle updates
        #pragma omp for
        for (int i = 0; i < sph_num_particles; i++) {
            Particle *p_i = &sph_particles[i];
            p_i->rho = rho[i];
            p_i->p = SPH_GAS_CONSTANT * (p_i->rho - SPH_REST_DENSITY);

            // Special rule: Humans flee high-density areas
            if (p_i->type == HUMAN && p_i->rho > 2 * SPH_REST_DENSITY) {
                p_i->vx -= 0.1f * p_i->p; // Adjust velocity based on pressure
                p_i->vy -= 0.1f * p_i->p;
            }
        }
    }

    free(rho);
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
void compute_sph_forces(void) {
    const float spiky_gradient = -10.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float visc_laplacian = 40.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float small_eps = 1e-6f;

    // Arrays for reduction pattern
    float *pressure_force_x = (float *)calloc(sph_num_particles, sizeof(*pressure_force_x));
    float *pressure_force_y = (float *)calloc(sph_num_particles, sizeof(*pressure_force_y));
    float *viscosity_force_x = (float *)calloc(sph_num_particles, sizeof(*viscosity_force_x));
    float *viscosity_force_y = (float *)calloc(sph_num_particles, sizeof(*viscosity_force_y));

    if (!pressure_force_x || !pressure_force_y || !viscosity_force_x || !viscosity_force_y) {
        fprintf(stderr, "Memory allocation failed for force arrays.\n");
        exit(EXIT_FAILURE);
    }

#pragma omp parallel default(none) shared(sph_particles, sph_num_particles, pressure_force_x, pressure_force_y, viscosity_force_x, viscosity_force_y, spiky_gradient, visc_laplacian, small_eps)
    {
#pragma omp for reduction(+:pressure_force_x[:sph_num_particles], pressure_force_y[:sph_num_particles], viscosity_force_x[:sph_num_particles], viscosity_force_y[:sph_num_particles]) collapse(2)
        for (int i = 0; i < sph_num_particles; i++) {
            for (int j = 0; j < sph_num_particles; j++) {
                if (i == j) continue;

                Particle *p_i = &sph_particles[i];
                Particle *p_j = &sph_particles[j];

                float dx = p_j->x - p_i->x;
                float dy = p_j->y - p_i->y;
                float dist = hypotf(dx, dy) + small_eps;

                if (dist < SPH_KERNEL_RADIUS) {
                    float weight = SPH_MASS / p_j->rho;
                    float spiky_term = pow(SPH_KERNEL_RADIUS - dist, 3.0f);

                    // Pressure force contribution
                    float pressure_term = (p_i->p + p_j->p) * spiky_gradient * spiky_term;
                    pressure_force_x[i] += -weight * pressure_term * (dx / dist);
                    pressure_force_y[i] += -weight * pressure_term * (dy / dist);

                    // Viscosity force contribution
                    float velocity_difference_x = p_j->vx - p_i->vx;
                    float velocity_difference_y = p_j->vy - p_i->vy;
                    float visc_term = (SPH_KERNEL_RADIUS - dist) * visc_laplacian;
                    viscosity_force_x[i] += visc_term * velocity_difference_x * weight;
                    viscosity_force_y[i] += visc_term * velocity_difference_y * weight;

                    // Special rule: Zombies repel each other based on pressure
                    if (p_i->type == ZOMBIE && p_j->type == ZOMBIE) {
                        pressure_force_x[i] += 0.5f * dx / dist;
                        pressure_force_y[i] += 0.5f * dy / dist;
                    }
                }
            }
        }

#pragma omp for
        for (int i = 0; i < sph_num_particles; i++) {
            Particle *p_i = &sph_particles[i];

            // Add gravitational forces
            float gravity_x = SPH_GRAVITY_X * SPH_MASS / p_i->rho;
            float gravity_y = SPH_GRAVITY_Y * SPH_MASS / p_i->rho;

            // Combine forces
            p_i->fx = pressure_force_x[i] + viscosity_force_x[i] + gravity_x;
            p_i->fy = pressure_force_y[i] + viscosity_force_y[i] + gravity_y;
        }
    }

    // Free allocated memory
    free(pressure_force_x);
    free(pressure_force_y);
    free(viscosity_force_x);
    free(viscosity_force_y);
}



/**
 * Integrate particle positions and velocities over time.
 */
void integrate_sph(void) {

    #pragma omp parallel for default(none) shared(sph_num_particles, sph_particles, DT, EPS, SPH_BOUNDARY_DAMPING, SPH_KERNEL_RADIUS, VIEW_WIDTH, VIEW_HEIGHT)
    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p = &sph_particles[i];

        // Update velocity
        p->vx += SPH_TIME_STEP * p->fx / SPH_MASS;
        p->vy += SPH_TIME_STEP * p->fy / SPH_MASS;

        // Update position
        p->x += SPH_TIME_STEP * p->vx;
        p->y += SPH_TIME_STEP * p->vy;

        // Enforce boundary conditions
        if (p->x < SPH_KERNEL_RADIUS) {
            p->vx *= SPH_BOUNDARY_DAMPING;
            p->x = SPH_KERNEL_RADIUS;
        }
        if (p->x > VIEW_WIDTH - SPH_KERNEL_RADIUS) {
            p->vx *= SPH_BOUNDARY_DAMPING;
            p->x = VIEW_WIDTH - SPH_KERNEL_RADIUS;
        }
        if (p->y < SPH_KERNEL_RADIUS) {
            p->vy *= SPH_BOUNDARY_DAMPING;
            p->y = SPH_KERNEL_RADIUS;
        }
        if (p->y > VIEW_HEIGHT - SPH_KERNEL_RADIUS) {
            p->vy *= SPH_BOUNDARY_DAMPING;
            p->y = VIEW_HEIGHT - SPH_KERNEL_RADIUS;
        }
    }
}



void handle_collisions(void) {
    const float kernel_radius_squared = SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS;

#pragma omp parallel for default(none) shared(sph_num_particles, sph_particles, kernel_radius_squared)
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == ZOMBIE) {
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - sph_particles[i].x;
                    float dy = sph_particles[j].y - sph_particles[i].y;
                    float dist_squared = dx * dx + dy * dy;

                    if (dist_squared < kernel_radius_squared) {
                        // Safely turn human into zombie
#pragma omp critical
                        {
                            sph_particles[j].type = ZOMBIE;
                        }
                    }
                }
            }
        } else if (sph_particles[i].type == IMMUNE) {
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - sph_particles[i].x;
                    float dy = sph_particles[j].y - sph_particles[i].y;
                    float dist_squared = dx * dx + dy * dy;

                    if (dist_squared < kernel_radius_squared) {
                        // Safely turn human into immune
#pragma omp critical
                        {
                            sph_particles[j].type = IMMUNE;
                        }
                    }
                }
            }
        }
    }
}




/**
 * Update the SPH system by computing densities, forces, and integrating.
 */
void update_sph(void) {
#pragma omp parallel sections
    {
        #pragma omp section
        compute_density_and_pressure();

        #pragma omp section
        compute_sph_forces();

        #pragma omp section
        move_humans();

        #pragma omp section
        move_zombies();

        #pragma omp section
        move_immune();

        #pragma omp section
        handle_collisions();

        #pragma omp section
        integrate_sph();
    }
}

