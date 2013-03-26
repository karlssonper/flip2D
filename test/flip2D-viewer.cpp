#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "../src/flip2D.h"

#include <iostream>

Settings::Ptr settings;
Particles::Ptr particles;

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < particles->numParticles(); ++i) {
        glVertex2f(particles->pos(i).x, particles->pos(i).y);
    }
    glEnd();
    
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key) {

    }
}

void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode( GL_PROJECTION);
    glLoadIdentity();

    const float sx = (settings->nx + 1) * settings->dx;
    const float sy = (settings->ny + 1) * settings->dx;
    glOrtho(0.0, sx, 0.0, sy, 0.0, 1.0);
    glMatrixMode( GL_MODELVIEW);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "flip2d-viewer error: wrong number of arguments"
                  << std::endl;
        exit(0);
    } else {
        settings = Settings::create();
        particles = Particles::create();
        FLIP2D::read(argv[1], settings, particles);
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800,800);
    glutCreateWindow("flip2D-viewer");
    glClearColor(0.f,0.f,0.f,1.0f);
        
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutReshapeFunc(reshape);
    glDisable(GL_DEPTH_TEST);
    glutMainLoop();

    return 0; 
}
