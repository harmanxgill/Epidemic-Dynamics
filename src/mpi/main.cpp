#include <stdio.h>
#include <stdlib.h>
#include "sph.h"
#ifdef GUI
#include "visualization.h"
#include <GL/glut.h>
#endif
#include <mpi.h>

int main(int argc, char **argv) {
    int comm_size, rank;
    int n_particles = DAM_PARTICLES;
    int n_steps = 50;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Allocate memory for particles
    sph_particles = (Particle *)malloc(MAX_PARTICLES * sizeof(Particle));
    if (!sph_particles) {
        fprintf(stderr, "Memory allocation failed.\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Initialize particles on the root process
    if (rank == 0) {
        init_sph(n_particles);
    }

    // Broadcast particles to all processes
    MPI_Bcast(sph_particles, MAX_PARTICLES, mpi_particle_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef GUI
    if (rank == 0) { // Only the root process handles GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
        glutInitWindowSize(1024, 768);
        glutCreateWindow("SPH Simulation");
        glutDisplayFunc(render);
        glutIdleFunc(update_sph);
        glutKeyboardFunc(keyboard_handler);
        glutMouseFunc(mouse_handler);
        init_gl();
        glutMainLoop();
    }
#else
    // Non-GUI distributed computation
    for (int step = 0; step < n_steps; step++) {
        update_sph();

        // Gather particles' updates across all processes
        MPI_Allgatherv(
            &sph_particles[displs[rank]],
            counts[rank],
            mpi_particle_type,
            sph_particles,
            counts,
            displs,
            mpi_particle_type,
            MPI_COMM_WORLD
        );

        // Optionally, calculate and print the average velocity
        float local_avg = avg_velocities();
        float global_avg = 0.0;
        MPI_Allreduce(&local_avg, &global_avg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0 && step % 10 == 0) {
            printf("Step %d: Avg Velocity = %f\n", step, global_avg);
        }
    }
#endif

    // Free resources
    free(sph_particles);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
