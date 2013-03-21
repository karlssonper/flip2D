#include <iostream>

#include "../src/sdf.h"
#include "../src/particles.h"
#include "../src/util.h"

bool printPassed = false;

int test(bool cond, const char * msg)
{
    if (cond) {
        if (printPassed) {
            std::cout << msg << " ... PASSED" << std::endl;
        }
        return 0;
    } else {
        std::cout << msg << " ... FAILED" << std::endl;
        return 1;
    }
}

int main(int argc, char *argv[]) {
    std::cout << "Starting SDF test..." << std::endl;

    Settings::Ptr solidSettings = Settings::create();
    solidSettings->nx = 30;
    solidSettings->ny = 30;
    solidSettings->dx = 1.0;
    SolidSDF::Ptr s = SolidSDF::create(solidSettings);
    s->initBoxBoundary(2);

    int numFailed = 0;

    numFailed += test(s->center(0,0) == (-1.5*3 - 0.5)/4.0, "boundary");
    numFailed += test(s->center(1,1) == (-0.5*3 + 0.5)/4.0, "boundary");
    numFailed += test(s->center(2,2) == (0.5*3 + 1.5)/4.0, "boundary");
    numFailed += test(s->center(3,3) == (1.5*3 + 2.5)/4.0, "boundary");
    numFailed += test(s->center(4,4) == (2.5*3 + 3.5)/4.0, "boundary");
    numFailed += test(s->center(15,15) == 3.5, "boundary");

    numFailed += test(SolidSDF::fractionInside(-1.0,1.0) == 0.5, "fraction");

    FaceArray2Xf u(30,30,1.0);
    FaceArray2Yf v(30,30,1.0);

    s->createWeights(u,v);

    numFailed += test(u.face<RIGHT>(2,0) == 0, "weights");
    numFailed += test(u.face<RIGHT>(2,1) == 0.5, "weights");
    numFailed += test(u.face<RIGHT>(2,2) == 1, "weights");
    numFailed += test(u.face<RIGHT>(0,2) == 0, "weights");
    numFailed += test(u.face<RIGHT>(1,2) == 1, "weights");
    numFailed += test(u.face<RIGHT>(2,2) == 1, "weights");

    numFailed += test(v.face<TOP>(2,0) == 0, "weights");
    numFailed += test(v.face<TOP>(2,1) == 1, "weights");
    numFailed += test(v.face<TOP>(2,2) == 1, "weights");
    numFailed += test(v.face<TOP>(0,2) == 0, "weights");
    numFailed += test(v.face<TOP>(1,2) == 0.5, "weights");
    numFailed += test(v.face<TOP>(2,2) == 1, "weights");

    int res = 128;
    Settings::Ptr fluidSettings = Settings::create();
    fluidSettings->nx = res;
    fluidSettings->ny = res;
    fluidSettings->dx = 1.0;
    FluidSDF::Ptr f = FluidSDF::create(fluidSettings);
    Particles::Ptr p = Particles::create();
    Vec2f mid(res*0.5,res*0.5);
    float radius = 40;
    const int particlesPerCell = 4;
    for (unsigned int y = 0; y < res; ++y) {
        for (unsigned int x = 0; x < res; ++x) {
            for (unsigned int n = 0; n < particlesPerCell; ++n) {
                Vec2f pos(x + 0.5f + random(-0.495, 0.495),
                          y + 0.5f + random(-0.495, 0.495));

                // Inside test
                const float px = (pos.x - mid.x);
                const float py = (pos.y - mid.y);
                if (sqr((pos-mid).length()) - radius * radius < 0.0) {
                    p->addParticle(pos);
                }
            }
        }
    }
    
    f->reconstructSurface(p, 1.0f, 0.6f);
    f->reinitialize(2);
    f->extrapolateIntoSolid(s);
    numFailed += test(f->isFluid(res*0.5,res*0.5), "weights");
    
    std::cout << "Number of failed tests: " << numFailed << std::endl;
}
