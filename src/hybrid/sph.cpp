#include "sph.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>

Particle* sph_particles = NULL;
int sph_num_particles = 0;

void init_zombie_simulation(int num_humans, int num_zombies, int num_immune, int rank, int size) {
    int total_particles = num_humans + num_zombies + num_immune;
    int particles_per_rank = total_particles / size;
    int start_index = rank * particles_per_rank;
    int end_index = (rank == size - 1) ? total_particles : start_index + particles_per_rank;

    sph_num_particles = end_index - start_index;

    #pragma omp parallel for
    for (int i = start_index; i < end_index; i++) {
        float x = (rand() / (float)RAND_MAX) * (VIEW_WIDTH - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;
        float y = (rand() / (float)RAND_MAX) * (VIEW_HEIGHT - 2 * SPH_KERNEL_RADIUS) + SPH_KERNEL_RADIUS;

        if (i < num_humans) {
            sph_particles[i - start_index] = (Particle){x, y, 0, 0, 0, 0, 0, 0, HUMAN, 0};
        } else if (i < num_humans + num_zombies) {
            sph_particles[i - start_index] = (Particle){x, y, 0, 0, 0, 0, 0, 0, ZOMBIE, 100};
        } else {
            sph_particles[i - start_index] = (Particle){x, y, 0, 0, 0, 0, 0, 0, IMMUNE, 0};
        }
    }

    printf("Rank %d initialized %d particles.\n", rank, sph_num_particles);
}

void compute_density_and_pressure(int rank, int size) {
    const float h_squared = SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS;
    const float poly6_constant = 4.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 8.0f));

    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        Particle* p_i = &sph_particles[i];
        p_i->rho = 0.0f;

        for (int j = 0; j < sph_num_particles; j++) {
            Particle* p_j = &sph_particles[j];

            float dx = p_j->x - p_i->x;
            float dy = p_j->y - p_i->y;
            float distance_squared = dx * dx + dy * dy;

            if (distance_squared < h_squared) {
                p_i->rho += SPH_MASS * poly6_constant * pow(h_squared - distance_squared, 3.0f);
            }
        }

        p_i->p = SPH_GAS_CONSTANT * (p_i->rho - SPH_REST_DENSITY);
    }
}

void compute_sph_forces(int rank, int size) {
    const float spiky_gradient = -10.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float visc_laplacian = 40.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float small_eps = 1e-6f;

    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        Particle* p_i = &sph_particles[i];
        float pressure_force_x = 0.0f, pressure_force_y = 0.0f;
        float viscosity_force_x = 0.0f, viscosity_force_y = 0.0f;

        for (int j = 0; j < sph_num_particles; j++) {
            Particle* p_j = &sph_particles[j];
            if (p_i == p_j) continue;

            float dx = p_j->x - p_i->x;
            float dy = p_j->y - p_i->y;
            float dist = hypotf(dx, dy) + small_eps;

            if (dist < SPH_KERNEL_RADIUS) {
                float weight = SPH_MASS / p_j->rho;
                float spiky_term = pow(SPH_KERNEL_RADIUS - dist, 3.0f);

                pressure_force_x += -weight * (p_i->p + p_j->p) * spiky_gradient * spiky_term * (dx / dist);
                pressure_force_y += -weight * (p_i->p + p_j->p) * spiky_gradient * spiky_term * (dy / dist);
            }
        }

        p_i->fx = pressure_force_x + viscosity_force_x;
        p_i->fy = pressure_force_y + viscosity_force_y;
    }
}

void synchronize_particles(int rank, int size) {
    MPI_Allgatherv(
        &sph_particles[0], sph_num_particles, MPI_PARTICLES,
        sph_particles, NULL, NULL, MPI_PARTICLES,
        MPI_COMM_WORLD
    );
}

void move_humans() {
    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == HUMAN) {
            Particle* human = &sph_particles[i];
            float dx = 0.0f, dy = 0.0f;

            // Compute avoidance vector
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == ZOMBIE) {
                    float dist_x = human->x - sph_particles[j].x;
                    float dist_y = human->y - sph_particles[j].y;
                    float dist_squared = dist_x * dist_x + dist_y * dist_y;

                    if (dist_squared < SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS) {
                        dx += dist_x / sqrtf(dist_squared);
                        dy += dist_y / sqrtf(dist_squared);
                    }
                }
            }

            // Normalize and scale
            float length = sqrtf(dx * dx + dy * dy);
            if (length > 0.0f) {
                human->vx = (dx / length) * 200.0f; // Adjust speed
                human->vy = (dy / length) * 200.0f;
            }
        }
    }
}

void move_zombies() {
    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == ZOMBIE) {
            Particle* zombie = &sph_particles[i];
            float closest_dist = VIEW_WIDTH * VIEW_HEIGHT;
            Particle* closest_human = NULL;

            // Find the nearest human
            for (int j = 0; j < sph_num_particles; j++) {
                if (sph_particles[j].type == HUMAN) {
                    float dx = sph_particles[j].x - zombie->x;
                    float dy = sph_particles[j].y - zombie->y;
                    float dist_squared = dx * dx + dy * dy;

                    if (dist_squared < closest_dist) {
                        closest_dist = dist_squared;
                        closest_human = &sph_particles[j];
                    }
                }
            }

            // Move toward the nearest human
            if (closest_human) {
                float dx = closest_human->x - zombie->x;
                float dy = closest_human->y - zombie->y;
                float length = sqrtf(dx * dx + dy * dy);
                zombie->vx = (dx / length) * 150.0f; // Zombies are slower
                zombie->vy = (dy / length) * 150.0f;
            }
        }
    }
}

void move_immune() {
    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        if (sph_particles[i].type == IMMUNE) {
            Particle* immune = &sph_particles[i];

            // Apply random movement
            immune->vx = ((rand() / (float)RAND_MAX) - 0.5f) * 100.0f; // Random velocity in range [-50, 50]
            immune->vy = ((rand() / (float)RAND_MAX) - 0.5f) * 100.0f;
        }
    }
}

void handle_collisions() {
    #pragma omp parallel for
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
                        sph_particles[j].health = 100;
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
                        // Protect human (e.g., by preventing infection)
                        // Placeholder: currently does nothing specific
                    }
                }
            }
        }
    }
}

void integrate_sph(int rank, int size) {
    #pragma omp parallel for
    for (int i = 0; i < sph_num_particles; i++) {
        Particle* p = &sph_particles[i];

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
}



// Implement other functions: move_humans, move_zombies, move_immune, handle_collisions, integrate_sph
