#include <stdio.h>
#include <stdlib.h>
#include "sph.h"
#include <time.h>
#include <chrono>
#include <mpi.h>

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype mpi_particles;
    create_mpi_particle_type(&mpi_particles);

    // Allocate memory for particles
    sph_particles = (Particle*)malloc(MAX_PARTICLES * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Memory allocation failed.\n");
        return EXIT_FAILURE;
    }


        int num_humans = 50;  // Adjust particle counts for scaling analysis
        int num_zombies = 50;
        int num_immune = 50;
        int total_particles = num_humans + num_zombies + num_immune;

        // Broadcast total number of particles
        MPI_Bcast(&total_particles, 1, mpi_particles, 0, MPI_COMM_WORLD);

        // Initialize particles on rank 0
        if (rank == 0) {
            init_zombie_simulation(num_humans, num_zombies, num_immune, rank, size);
        }

        // Broadcast particles to all processes
        MPI_Bcast(sph_particles, total_particles, MPI_PARTICLE, 0, MPI_COMM_WORLD);

        // Start timing
        clock_t start_time = clock();
        int step = 0; // Step counter

    while (1) {
        update_sph(rank, size);

        // Count particle types
        int local_humans = 0, local_zombies = 0, local_immune = 0;
        for (int i = 0; i < total_particles; i++) {
            Particle* p = &sph_particles[i];
            if (p->type == HUMAN) local_humans++;
            else if (p->type == ZOMBIE) local_zombies++;
            else if (p->type == IMMUNE) local_immune++;
        }

        // Reduce counts across all ranks
        int humans = 0, zombies = 0, immune = 0;
        MPI_Reduce(&local_humans, &humans, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_zombies, &zombies, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_immune, &immune, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // Only rank 0 outputs results
        if (rank == 0 && step % 5 == 0) {
            printf("Step %d: Humans: %d, Zombies: %d, Immune: %d\n", step, humans, zombies, immune);
        }

        // End simulation if no humans remain
        if (rank == 0 && humans == 0) {
            printf("No humans remain at step %d.\n", step);
            if (zombies > immune) {
                printf("Zombies win with %d remaining!\n", zombies);
            } else if (immune > zombies) {
                printf("Immune win with %d remaining!\n", immune);
            } else {
                printf("It's a draw between Zombies and Immune!\n");
            }
            break;
        }

        step++;
    }

        // End timing
        clock_t end_time = clock();

        double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        printf("Simulation completed in %.2f seconds for %d particles and %d steps.\n",
            elapsed_time, total_particles, step);

    free(sph_particles);
    MPI_Finalize();
    return 0;
}
