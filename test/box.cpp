#include <iostream>
#include "../src/flip2D.h"

int main(int argc, char *argv[]) {
    std::cout << "<<< Box Test >>>" << std::endl;

    Settings::Ptr s = Settings::create();
    s->nx = 128;
    s->ny = 128;
    s->dx = 1.0 / 129.0;
    s->solidWidth = 3.0f;
    s->initialFluidCenter = Vec2f(0.5,0.25);
    s->initialFluidRadius = 0.33;
    s->particlesPerCell = 4;
    s->R = 1.0 * s->dx;
    s->r = 0.6 * s->dx;
    s->numPhiSweepIterations = 2;
    s->gravity = Vec2f(0.0f, -0.82f);
    s->numVelSweepIterations = 4;
    s->usePCG = true;
    s->tolerance = 1e-5;
    s->maxIterations = 100;

    FLIP2D::Ptr flip = FLIP2D::create(s);
    flip->write("test.flip2D");
    flip->step(5/24.0);
    flip->write("step.flip2D");
}
