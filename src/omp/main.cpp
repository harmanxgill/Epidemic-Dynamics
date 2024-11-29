#include <stdio.h>
#include <stdlib.h>
#include "sph.h"
#include <time.h>
#include <omp.h>

int main(int argc, char **argv) {
    srand(1234);
    sph_particles = (Particle*)malloc(MAX_PARTICLES * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Memory allocation failed.\n");
        return EXIT_FAILURE;
    }


        int num_humans = 500000;  // Adjust particle counts for scaling analysis
        int num_zombies = 500;
        int num_immune = 500;

        // Initialize simulation
        init_zombie_simulation(num_humans, num_zombies, num_immune);

        // Start timing
        clock_t start_time = clock();

        int step = 0; // Step counter
    while (1) {
        update_sph();

        // Count particle types
        int humans = 0, zombies = 0, immune = 0;
        #pragma omp parallel for reduction(+:humans, zombies, immune)
        for (int i = 0; i < sph_num_particles; i++) {
            if (sph_particles[i].type == HUMAN) humans++;
            else if (sph_particles[i].type == ZOMBIE) zombies++;
            else if (sph_particles[i].type == IMMUNE) immune++;
        }

        // Print particle counts every 50 steps
        if (step % 5 == 0) {
            printf("Step %d: Humans: %d, Zombies: %d, Immune: %d\n", step, humans, zombies, immune);
        }

        // End simulation if no humans remain
        if (humans == 0) {
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
            elapsed_time, sph_num_particles, step);

    free(sph_particles);
    return 0;
}
