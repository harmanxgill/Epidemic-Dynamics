#include "visualization.h"
#include "sph.h"
#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h> // For sleep()

// Frame limit for logging
static const int FRAME_LIMIT = 100;
static int current_frame = 0;

/**
 * Initialize OpenGL settings for rendering the simulation.
 */
void init_gl(void) {
    glClearColor(0.1, 0.1, 0.1, 1.0);  // Dark background
    glEnable(GL_POINT_SMOOTH);         // Smooth particle rendering
    glMatrixMode(GL_PROJECTION);       // Set projection matrix
    glPointSize(5.0f);                 // Default particle size
}

void render_text(float x, float y, const char* text) {
    glRasterPos2f(x, y); // Set the position for the text
    while (*text) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *text); // Draw each character
        text++;
    }
}

/**
 * Render the particles on the screen with different colors based on their type.
 */
void render(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

    int humans = 0, zombies = 0, immune = 0;

    // Count and render particles
    glBegin(GL_POINTS);
    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p = &sph_particles[i];
        if (p->type == HUMAN) {
            glColor3f(0.2f, 0.8f, 0.2f); // Green
            humans++;
        } else if (p->type == ZOMBIE) {
            glColor3f(0.8f, 0.2f, 0.2f); // Red
            zombies++;
        } else if (p->type == IMMUNE) {
            glColor3f(0.2f, 0.2f, 0.8f); // Blue
            immune++;
        }
        glVertex2f(p->x, p->y);
    }
    glEnd();

    // Log rendering
    render_text(10, VIEW_HEIGHT - 30, "Zombie Epidemic");
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "Humans: %d", humans);
    render_text(10, VIEW_HEIGHT - 60, buffer);
    snprintf(buffer, sizeof(buffer), "Zombies: %d", zombies);
    render_text(10, VIEW_HEIGHT - 90, buffer);
    snprintf(buffer, sizeof(buffer), "Immune: %d", immune);
    render_text(10, VIEW_HEIGHT - 120, buffer);

    // Check for game end condition
    if (humans == 0) {
        if (zombies > immune) {
            snprintf(buffer, sizeof(buffer), "Zombies win!");
        } else if (immune > zombies) {
            snprintf(buffer, sizeof(buffer), "Immune wins!");
        } else {
            snprintf(buffer, sizeof(buffer), "It's a draw!");
        }
        render_text(VIEW_WIDTH / 2 - 50, VIEW_HEIGHT / 2, buffer);
        glutSwapBuffers();
        sleep(5); // Pause for 5 seconds
        exit(0); // End the simulation
    }

    glutSwapBuffers();
    glutPostRedisplay();
}


/**
 * Handle keyboard input events.
 * - 'r' or 'R': Resets the particle system.
 * - 'q' or 'Q': Exits the simulation.
 */
void keyboard_handler(unsigned char key, int x, int y) {
    if (key == 'r' || key == 'R') {
        init_zombie_simulation(400, 50, 50); // Reset simulation with 400 humans, 50 zombies, 50 immune
        printf("Simulation reset.\n");
    } else if (key == 'q' || key == 'Q') {
        printf("Exiting simulation.\n");
        exit(0);
    }
}

/**
 * Handle mouse input events.
 * Left-click adds an explosion of particles at the mouse location.
 */
void mouse_handler(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        float center_x = (float)x / glutGet(GLUT_WINDOW_WIDTH) * VIEW_WIDTH;
        float center_y = (float)(glutGet(GLUT_WINDOW_HEIGHT) - y) / glutGet(GLUT_WINDOW_HEIGHT) * VIEW_HEIGHT;

        // Create an explosion of humans at the mouse location
        for (int i = 0; i < 10 && sph_num_particles < MAX_PARTICLES; i++) {
            float angle = (rand() / (float)RAND_MAX) * 2.0f * M_PI;
            float speed = 50.0f + (rand() / (float)RAND_MAX) * 50.0f;
            float vx = speed * cos(angle);
            float vy = speed * sin(angle);

            init_particle(&sph_particles[sph_num_particles], center_x, center_y, HUMAN);
            sph_particles[sph_num_particles].vx = vx;
            sph_particles[sph_num_particles].vy = vy;
            sph_num_particles++;
        }

        printf("Explosion of particles at (%.2f, %.2f).\n", center_x, center_y);
    }
}
