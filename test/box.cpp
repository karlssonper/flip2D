#include "../src/flip2D.h"
#include "../src/log.h"

#include <iostream>
#include <sstream>

void filename(std::string & input, int frame)
{
    size_t pos = input.find("$F");
    if (pos != std::string::npos) {
        std::stringstream ss;
        ss << frame;
        input.replace(pos,2, ss.str());
    }
}

std::string frame(int i)
{
    std::stringstream ss;
    ss << std::endl
       << "-------------------------------------------------------------------"
       << std::endl
       << "                              Frame " << i
       << std::endl
       << "-------------------------------------------------------------------";
    return ss.str();
}

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
    std::string simOutput;
    if (argc > 1) {
        simOutput = argv[1];
    } else {
        simOutput = "sim/boxSim.$F.flip2D";
    }
    int nFrames = 24;
    for(int i = 0; i < nFrames; ++i) {
        LOG_OUTPUT_WITHOUT_TIMESTAMPS(frame(i));
        flip->step(1.0/24.0);
        std::string simOutputFrame = simOutput;
        filename(simOutputFrame,i);
        flip->write(simOutputFrame.c_str());
    }
}
