#include "flip2D.h"
#include "util.h"
#include "pcg.h"
#include "multigrid.h"
#include "gaussSeidel.h"
#include "jacobi.h"
#include <fstream>

FLIP2D::FLIP2D(Settings::Ptr s) : _settings(s)
{
    _grid = Grid::create(s);
    _particles = Particles::create();
    _fluid = FluidSDF::create(s);
    _solid = SolidSDF::create(s);

    if (s->usePCG) {
        _pressureSolver = PCG::create(s);
    } else if (s->useMultigrid) {
        _pressureSolver = Multigrid::create(s);
    } else if (s->useGaussSeidel) {
        _pressureSolver = GaussSeidel::create(s);
    } else if (s->useJacobi) {
        _pressureSolver = Jacobi::create(s);
    } else {
        // error
    }

    // Init Solid SDF and spawn fluid particles
    _solid->initBoxBoundary(s->solidWidth);
    _particles->initSphere(_solid->phi(),
                           _settings->initialFluidCenter,
                           _settings->initialFluidRadius,
                           _settings->particlesPerCell,
                           _settings->initialVelocity);
}

void FLIP2D::step(float dt)
{
    float tStep = 0;
    while (tStep < dt) {
        const float t = min(_grid->CFL(), dt - tStep);
        _grid->sampleVelocities(_particles);
        _fluid->reconstructSurface(_particles, _settings->R, _settings->r);
        _fluid->reinitialize(_settings->numPhiSweepIterations);
        _fluid->extrapolateIntoSolid(_solid);
        _grid->applyGravity(_settings->gravity, t);
        _grid->extrapolateVelocities(_fluid, _settings->numVelSweepIterations);
        _pressureSolver->buildLinearSystem(_grid, _solid, _fluid, t);
        _pressureSolver->solveLinearSystem(_fluid, t);
        _grid->pressureProjection(_pressureSolver->pressure(), _fluid, t);
        _grid->extrapolateVelocities(_fluid, _settings->numVelSweepIterations);
        _grid->enforceBoundaryConditions(_solid);
        _particles->advect(_grid->u(), _grid->v(), t);
        _particles->updateVelocities(_grid->u(), _grid->v());
        tStep += t;
    }
}

bool FLIP2D::write(const char * filename) const
{
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if (out.is_open()) {
        _settings->write(out);
        _particles->write(out);    
        out.close();
        return true;
    } else {
        //error
        return false;
    }
}

bool FLIP2D::read(const char * filename, Settings::Ptr s, Particles::Ptr p)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in.is_open()) {
        s->read(in);
        p->read(in);
        in.close();
        return true;
    } else {
        //error
        return false;
    }
}

