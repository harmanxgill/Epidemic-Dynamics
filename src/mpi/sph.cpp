#include "sph.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

/* Constants for SPH */
const float SPH_GRAVITY_X = 0.0f, SPH_GRAVITY_Y = -10.0f; 
const float SPH_REST_DENSITY = 300.0f;                   
const float SPH_GAS_CONSTANT = 2000.0f;                  
const float SPH_VISCOSITY = 200.0f;                      
const float SPH_TIME_STEP = 0.0007f;                     
const float SPH_BOUNDARY_DAMPING = -0.5f;                
const float SPH_MASS = 2.5;   

#ifdef GUI
const int SPH_MAX_PARTICLES = 5000;                     
#else
const int SPH_MAX_PARTICLES = 20000;                    
#endif

const int SPH_INITIAL_PARTICLES = 500;                  

/* Global variables */
Particle *sph_particles;
int sph_num_particles = 0;
int *counts, *displs;
MPI_Datatype mpi_particle_type;
int comm_rank, comm_size;

/**
 * Initialize the SPH system with `n` particles.
 */
void init_sph(int num_particles) {
    if (comm_rank == 0) {
        sph_num_particles = 0;
        printf("Initializing SPH system with %d particles.\n", num_particles);

        for (float y = SPH_KERNEL_RADIUS; y < VIEW_HEIGHT - SPH_KERNEL_RADIUS; y += SPH_KERNEL_RADIUS) {
            for (float x = SPH_KERNEL_RADIUS; x < VIEW_WIDTH * 0.8f; x += SPH_KERNEL_RADIUS) {
                if (sph_num_particles < SPH_MAX_PARTICLES) {
                    float jitter = rand() / (float)RAND_MAX;
                    init_particle(&sph_particles[sph_num_particles], x + jitter, y);
                    sph_num_particles++;
                } else {
                    return;
                }
            }
        }
        assert(sph_num_particles == num_particles);
    }

    // Broadcast particle data to all processes
    MPI_Bcast(sph_particles, SPH_MAX_PARTICLES, mpi_particle_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sph_num_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * Compute density and pressure for each particle.
 */
void compute_density_and_pressure(void) {
    const float h_squared = SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS; 
    const float poly6_constant = 4.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 8.0f));

    int start = displs[comm_rank];
    int local_n = counts[comm_rank];

    for (int i = start; i < start + local_n; i++) {
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

    // Synchronize particles across processes
    MPI_Allgatherv(&sph_particles[start], local_n, mpi_particle_type, sph_particles, counts, displs, mpi_particle_type, MPI_COMM_WORLD);
}

/**
 * Compute forces for each particle.
 */
void compute_sph_forces(void) {
    const float spiky_gradient = -10.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float visc_laplacian = 40.0f / (M_PI * pow(SPH_KERNEL_RADIUS, 5.0f));
    const float small_eps = 1e-6f;

    int start = displs[comm_rank];
    int local_n = counts[comm_rank];

    for (int i = start; i < start + local_n; i++) {
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
                float velocity_difference_x = p_j->vx - p_i->vx;
                float velocity_difference_y = p_j->vy - p_i->vy;
                float visc_term = (SPH_KERNEL_RADIUS - dist) * visc_laplacian;
                viscosity_force_x += visc_term * velocity_difference_x * weight;
                viscosity_force_y += visc_term * velocity_difference_y * weight;
            }
        }

        // Add gravitational forces
        float gravity_x = SPH_GRAVITY_X * SPH_MASS / p_i->rho;
        float gravity_y = SPH_GRAVITY_Y * SPH_MASS / p_i->rho;

        // Set total forces
        p_i->fx = pressure_force_x + viscosity_force_x + gravity_x;
        p_i->fy = pressure_force_y + viscosity_force_y + gravity_y;
    }

    // Synchronize particles across processes
    MPI_Allgatherv(&sph_particles[start], local_n, mpi_particle_type, sph_particles, counts, displs, mpi_particle_type, MPI_COMM_WORLD);
}

/**
 * Integrate particle positions and velocities over time.
 */
void integrate_sph(void) {
    int start = displs[comm_rank];
    int local_n = counts[comm_rank];

    for (int i = start; i < start + local_n; i++) {
        Particle *p = &sph_particles[i];

        // Update velocity
        p->vx += SPH_TIME_STEP * p->fx / p->rho;
        p->vy += SPH_TIME_STEP * p->fy / p->rho;

        // Update position
        p->x += SPH_TIME_STEP * p->vx;
        p->y += SPH_TIME_STEP * p->vy;

        // Boundary conditions
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

/**
 * Compute the average velocity of particles using MPI.
 */
float calculate_avg_velocity(void) {
    double local_velocity_sum = 0.0;
    int local_particle_count = counts[comm_rank]; // Number of particles handled by this process
    int start = displs[comm_rank];               // Start index for this process

    // Compute local velocity sum
    for (int i = start; i < start + local_particle_count; i++) {
        local_velocity_sum += hypot(sph_particles[i].vx, sph_particles[i].vy);
    }

    // Aggregate velocity sums and particle counts across all processes
    double global_velocity_sum = 0.0;
    MPI_Allreduce(&local_velocity_sum, &global_velocity_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Global average velocity
    return global_velocity_sum / sph_num_particles;
}

/**
 * Update the SPH system by computing densities, forces, and integrating using MPI.
 */
void update_sph(void) {
    // Compute density and pressure for local particles
    compute_density_and_pressure();

    // Synchronize particle data across all processes
    MPI_Allgatherv(
        &sph_particles[displs[comm_rank]],
        counts[comm_rank],
        mpi_particle_type,
        sph_particles,
        counts,
        displs,
        mpi_particle_type,
        MPI_COMM_WORLD
    );

    // Compute forces for local particles
    compute_sph_forces();

    // Synchronize particle data across all processes
    MPI_Allgatherv(
        &sph_particles[displs[comm_rank]],
        counts[comm_rank],
        mpi_particle_type,
        sph_particles,
        counts,
        displs,
        mpi_particle_type,
        MPI_COMM_WORLD
    );

    // Integrate positions and velocities for local particles
    integrate_sph();

    // Synchronize particle data across all processes
    MPI_Allgatherv(
        &sph_particles[displs[comm_rank]],
        counts[comm_rank],
        mpi_particle_type,
        sph_particles,
        counts,
        displs,
        mpi_particle_type,
        MPI_COMM_WORLD
    );
}
