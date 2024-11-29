#ifndef SPH_H
#define SPH_H

#include "particle.h"

// Constants
const float SPH_KERNEL_RADIUS = 10.0f;                  // Kernel radius
const float SPH_GRAVITY_X = 0.0f, SPH_GRAVITY_Y = 0.0f; // Gravity components
const float SPH_GAS_CONSTANT = 2000.0f;                  // Gas constant
const float SPH_VISCOSITY = 200.0f;                      // Viscosity constant
const float SPH_TIME_STEP = 0.0007f;                     // Time step
const float SPH_BOUNDARY_DAMPING = -0.5f;                // Boundary damping
const float SPH_MASS = 2.5f;                             // Particle mass
const float SPH_REST_DENSITY = 300.0f;                   // Rest density

#ifdef GUI
const int MAX_PARTICLES = 5000;                          // Max particles for GUI
#else
const int MAX_PARTICLES = 1000000;                         // Max particles for non-GUI
#endif

#ifndef GUI
#define VIEW_WIDTH 3000  // Larger width for non-GUI
#define VIEW_HEIGHT 2000 // Larger height for non-GUI
#else
#define VIEW_WIDTH 1024  // Standard width for GUI
#define VIEW_HEIGHT 768  // Standard height for GUI
#endif

// Global Variables
extern Particle *sph_particles;                          // Array of SPH particles
extern int sph_num_particles;                            // Current number of particles

// Function Declarations
void init_zombie_simulation(int humans, int zombies, int immune);
void compute_density_and_pressure(void);                 // SPH density and pressure computation
void compute_sph_forces(void);                           // Compute SPH forces
void integrate_sph(void);                                // Integrate positions and velocities
void update_sph(void);                                   // Update the SPH system
void move_humans(void);                                  // Move humans based on nearby zombies
void move_zombies(void);                                 // Move zombies toward humans
void handle_collisions(void);                            // Handle collisions between particles
int is_in_domain(float x, float y);                      // Check if a particle is in the domain
void immune_attack(void);

#endif // SPH_H
