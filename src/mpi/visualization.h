#ifndef VISUALIZATION_H
#define VISUALIZATION_H

/* Initialize OpenGL settings */
void init_gl(void);

/* Render the particles */
void render(void);

/* Handle keyboard input */
void keyboard_handler(unsigned char c, int x, int y);

/* Handle mouse input */
void mouse_handler(int button, int state, int x, int y);

/* Add a ball into the simulation */
void place_ball(float cx, float cy, float r);

#endif
