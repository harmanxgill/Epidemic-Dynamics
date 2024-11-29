#ifndef SPH_H
#define SPH_H

#include <mpi.h>
#include <stddef.h>

#define SPH_KERNEL_RADIUS 16.0f
#define SPH_MASS 2.5f
#define SPH_REST_DENSITY 300.0f
#define SPH_GAS_CONSTANT 2000.0f
#define SPH_TIME_STEP 0.0007f

#define VIEW_WIDTH 1000.0f
#define VIEW_HEIGHT 1000.0f

typedef enum { HUMAN, ZOMBIE, IMMUNE } ParticleType;

typedef struct {
    float x, y;         // Position
    float vx, vy;       // Velocity
    float fx, fy;       // Force
    float rho, p;       // Density, pressure
    ParticleType type;  // Particle type
    int health;         // Health (used for zombies)
} Particle;

extern MPI_Datatype MPI_PARTICLES; // MPI data type for Particle

// Global variables
extern Particle* sph_particles;
extern int sph_num_particles;

// Function prototypes
void init_zombie_simulation(int num_humans, int num_zombies, int num_immune, int rank, int size);
void compute_density_and_pressure(int rank, int size);
void compute_sph_forces(int rank, int size);
void move_humans();
void move_zombies();
void move_immune();
void handle_collisions();
void integrate_sph(int rank, int size);
void synchronize_particles(int rank, int size);

#endif
