#include "visualization.h"
#include "sph.h"
#include <GL/glut.h>
#include <stdio.h>

static const int FRAME_LIMIT = 100;
static int current_frame = 0;

/**
 * Initialize OpenGL settings for rendering the simulation.
 */
void init_gl(void) {
    glClearColor(0.0, 0.0, 0.1, 1.0); // Dark space-like background
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Add transparency
    glEnable(GL_POINT_SMOOTH);        // Smooth particle rendering
    glPointSize(SPH_KERNEL_RADIUS / 2.5); // Point size
    glMatrixMode(GL_PROJECTION);      // Set projection matrix
}

/**
 * Render the particles on the screen.
 */
void render(void) {
    glClear(GL_COLOR_BUFFER_BIT);     // Clear the buffer
    glLoadIdentity();
    glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1); // Set orthogonal projection

    glBegin(GL_POINTS);

    // Draw all particles
    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p = &sph_particles[i];

        // Set color based on density (cosmic glow effect)
        float intensity = p->rho / SPH_REST_DENSITY;
        glColor4f(0.2f + intensity, 0.5f * intensity, 1.0f, 1.0f); // Blue-ish color gradient

        // Set size based on velocity magnitude
        glPointSize(2.0 + hypot(p->vx, p->vy) * 10.0);

        glVertex2f(p->x, p->y);
    }

    glEnd();

    glutSwapBuffers();               // Swap buffers to display the frame
    glutPostRedisplay();             // Trigger next frame rendering

    current_frame++;
    if (current_frame >= FRAME_LIMIT) {
        float avg_velocity = calculate_avg_velocity();
        printf("Average Velocity: %f\n", avg_velocity);
        current_frame = 0;
    }
}

void place_explosion(float cx, float cy, float energy) {
    for (int i = 0; i < 50; i++) { // Create 50 new particles
        if (sph_num_particles < MAX_PARTICLES) {
            float angle = (rand() / (float)RAND_MAX) * 2.0f * M_PI; // Random direction
            float speed = (rand() / (float)RAND_MAX) * energy; // Random speed based on energy
            float vx = speed * cos(angle);
            float vy = speed * sin(angle);

            init_particle(&sph_particles[sph_num_particles], cx, cy);
            sph_particles[sph_num_particles].vx = vx;
            sph_particles[sph_num_particles].vy = vy;
            sph_num_particles++;
        }
    }
    printf("Supernova explosion at (%.2f, %.2f) with energy %.2f. Total particles: %d.\n",
           cx, cy, energy, sph_num_particles);
}


/**
 * Handle keyboard input events.
 * - 'r' or 'R' resets the particle system.
 */
void keyboard_handler(unsigned char key, int x_pos, int y_pos) {
    if (key == 'r' || key == 'R') {
        init_sph(DAM_PARTICLES);
        printf("Simulation reset with %d particles.\n", DAM_PARTICLES);
    }
}

void mouse_handler(int button, int state, int mouse_x, int mouse_y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        float center_x = 1.5 * mouse_x;                   // Scale mouse X to simulation space
        float center_y = VIEW_HEIGHT - 1.5 * mouse_y;     // Scale and flip mouse Y
        place_explosion(center_x, center_y, 50.0f);       // Energy of 50.0
    }
}

void apply_gravity_well(float cx, float cy, float strength) {
    for (int i = 0; i < sph_num_particles; i++) {
        Particle *p = &sph_particles[i];
        float dx = cx - p->x;
        float dy = cy - p->y;
        float distance = hypot(dx, dy) + 1e-5f; // Avoid division by zero

        // Apply gravitational pull
        float force = strength / (distance * distance);
        p->vx += force * dx / distance;
        p->vy += force * dy / distance;
    }
}

void keyboard_handler(unsigned char key, int x_pos, int y_pos) {
    if (key == 'r' || key == 'R') {
        init_sph(DAM_PARTICLES);
        printf("Simulation reset with %d particles.\n", DAM_PARTICLES);
    } else if (key == 'g' || key == 'G') {
        apply_gravity_well(VIEW_WIDTH / 2.0, VIEW_HEIGHT / 2.0, 1000.0f); // Gravity at center
        printf("Gravity well applied.\n");
    }
}


/**
 * Add a ball (a group of particles in a circular shape) into the simulation.
 * - `cx`, `cy`: Center coordinates of the ball.
 * - `r`: Radius of the ball.
 */
void place_ball(float cx, float cy, float r) {
    for (float y = cy - r; y <= cy + r; y += SPH_KERNEL_RADIUS) {
        for (float x = cx - r; x <= cx + r; x += SPH_KERNEL_RADIUS) {
            if (sph_num_particles < MAX_PARTICLES &&
                is_in_domain(x, y) &&
                ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r)) {
                float jitter_x = rand() / (float)RAND_MAX; // Add random jitter
                float jitter_y = rand() / (float)RAND_MAX;

                init_particle(&sph_particles[sph_num_particles], x + jitter_x, y + jitter_y);
                sph_num_particles++;
            }
        }
    }
    printf("Placed a ball with radius %.2f at (%.2f, %.2f). Total particles: %d.\n", 
           r, cx, cy, sph_num_particles);
}
