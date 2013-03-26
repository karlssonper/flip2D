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
    _initFluid();
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
        _particles->updateVelocities(_grid->u(), _grid->v());
        _particles->advect(t);
        tStep += t;
    }
}

bool FLIP2D::write(const char * filename)
{
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if (out.is_open()) {
        _particles->write(out);    
        out.close();
    } else {
        //error
    }
}

void FLIP2D::_initFluid()
{
    const float dx = _fluid->phi().dx();
    const float r2 = sqr(_settings->initialFluidRadius);
    for (int i = 0; i < _fluid->phi().nx(); ++i) {
        for (int j = 0; j < _fluid->phi().ny(); ++j) {
            for (int n = 0; n < _settings->particlesPerCell; ++n) {
                const Vec2f offset(random(-0.495, 0.495),random(-0.495, 0.495));
                const Vec2f pos = _fluid->phi().pos(i,j) + dx * offset;
                const Vec2f d = pos - _settings->initialFluidCenter;

                // Inside test
                if (sqr(d.x) + sqr(d.y) - r2 < 0 && !_solid->inside(pos)) {
                    _particles->addParticle(pos, _settings->initialVelocity);
                }
            }
        }
    }
}

