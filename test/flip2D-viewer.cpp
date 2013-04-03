#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "../src/flip2D.h"

#include <iostream>
#include <string>
#include <sstream>

Settings::Ptr settings;

int frame = 0;
std::string simfile;

void filename(std::string & input, int frame)
{
    size_t pos = simfile.find("$F");
    if (pos != std::string::npos) {
        std::stringstream ss;
        ss << frame;
        input.replace(pos,2, ss.str());
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    //Read frame
    std::string simfileFrame = simfile;
    filename(simfileFrame,frame);
    Particles::Ptr particles = Particles::create();
    FLIP2D::read(simfileFrame.c_str(), settings, particles);
    
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
        case 's':
            ++frame;
            break;
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
        simfile = std::string(argv[1]);
        std::string simfileFirst = simfile;
        filename(simfileFirst,0);
        settings = Settings::create();
        Particles::Ptr particles = Particles::create();
        FLIP2D::read(simfileFirst.c_str(), settings, particles);
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
