#ifndef PARTICLE_H
#define PARTICLE_H

#include "sph.h"
// Enum to define different types of particles in the simulation
enum ParticleType {
    HUMAN,       // Human particles
    ZOMBIE,      // Zombie particles
    IMMUNE,       // Immune particles
    INERT
};

// Particle structure
typedef struct {
    float x, y;           // Position
    float vx, vy;         // Velocity
    float fx, fy;         // Force
    float rho, p;         // Density and pressure
    ParticleType type;    // Type of particle (HUMAN, ZOMBIE, IMMUNE)
    int health;         // Health of the particle (for zombies only)
} Particle;

// Function declarations
void init_particle(Particle *p, float x, float y, ParticleType type);
float randab(float a, float b); // Generate a random float in [a, b]

#endif // PARTICLE_H
