#include "sph.h"
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

MPI_Datatype MPI_PARTICLES; // Define MPI datatype for Particle

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create custom MPI datatype for Particle
    MPI_Type_contiguous(sizeof(Particle), MPI_BYTE, &MPI_PARTICLES);
    MPI_Type_commit(&MPI_PARTICLES);

    // Particle setup
    int num_humans = 50, num_zombies = 50, num_immune = 50;
    int total_particles = num_humans + num_zombies + num_immune;

    sph_particles = (Particle*)malloc(total_particles * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Rank %d: Memory allocation failed.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Initialize particles
    init_zombie_simulation(num_humans, num_zombies, num_immune, rank, size);

    // Main simulation loop
    int steps = 100;
    for (int step = 0; step < steps; step++) {
        compute_density_and_pressure(rank, size);
        compute_sph_forces(rank, size);
        move_humans();
        move_zombies();
        move_immune();
        handle_collisions();
        integrate_sph(rank, size);

        synchronize_particles(rank, size);

        if (rank == 0 && step % 10 == 0) {
            printf("Step %d completed.\n", step);
        }
    }

    // Free memory and MPI resources
    free(sph_particles);
    MPI_Type_free(&MPI_PARTICLES);
    MPI_Finalize();
    return 0;
}
