#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "sph.h"

// Function declarations
void init_gl(void);                       // Initialize OpenGL settings
void render(void);                        // Render particles on the screen
void keyboard_handler(unsigned char key, int x, int y); // Handle keyboard input
void mouse_handler(int button, int state, int x, int y); // Handle mouse input

#endif // VISUALIZATION_H
