#include <stdio.h>
#include <stdlib.h>
#include "sph.h"
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define custom MPI datatype for Particle
    MPI_Type_contiguous(sizeof(Particle) / sizeof(float), MPI_FLOAT, &MPI_PARTICLES);
    MPI_Type_commit(&MPI_PARTICLES);

    int num_humans = 50;
    int num_zombies = 50;
    int num_immune = 50;
    int total_particles = num_humans + num_zombies + num_immune;

    // Arrays for counts and displacements
    int *counts = (int*)malloc(size * sizeof(int));
    int *displs = (int*)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        int start = total_particles * i / size;
        int end = total_particles * (i + 1) / size;
        counts[i] = end - start;
        displs[i] = start;
    }

    // Allocate memory for particles
    sph_particles = (Particle*)malloc(total_particles * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    // Initialize particles (only rank 0 initializes)
    if (rank == 0) {
        init_zombie_simulation(num_humans, num_zombies, num_immune, rank, size);
    }

    // Broadcast total number of particles and particles array
    MPI_Bcast(&total_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(sph_particles, total_particles, MPI_PARTICLES, 0, MPI_COMM_WORLD);

    // Start timing
    clock_t start_time = clock();
    int step = 0;

    while (1) {
        // Update SPH system
        update_sph(rank, size);

        // Count particle types locally
        int local_humans = 0, local_zombies = 0, local_immune = 0;
        for (int i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
            Particle *p = &sph_particles[i];
            if (p->type == HUMAN) local_humans++;
            else if (p->type == ZOMBIE) local_zombies++;
            else if (p->type == IMMUNE) local_immune++;
        }

        // Reduce counts across ranks
        int humans = 0, zombies = 0, immune = 0;
        MPI_Reduce(&local_humans, &humans, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_zombies, &zombies, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_immune, &immune, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // Gather updated particles
        MPI_Allgatherv(
            &sph_particles[displs[rank]],
            counts[rank],
            MPI_PARTICLES,
            sph_particles,
            counts,
            displs,
            MPI_PARTICLES,
            MPI_COMM_WORLD
        );

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
    if (rank == 0) {
        double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        printf("Simulation completed in %.2f seconds for %d particles and %d steps.\n",
               elapsed_time, total_particles, step);
    }

    // Free memory
    free(sph_particles);
    free(counts);
    free(displs);
    MPI_Type_free(&MPI_PARTICLES);
    MPI_Finalize();
    return 0;
}
