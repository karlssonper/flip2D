#include "pcg.h"

PCG::PCG(Settings::Ptr s) :
        PressureSolver(s, PRECONDITIONED_CONJUGATE_GRADIENT),
        _tol(s->tolerance),
        _maxIterations(s->maxIterations)
{
    _z.resize(s->nx,s->ny,s->dx);
    _s.resize(s->nx,s->ny,s->dx);
    _q.resize(s->nx,s->ny,s->dx);
    _precon.resize(s->nx,s->ny,s->dx);
}

void PCG::buildLinearSystem(Grid::Ptr grid,
                            SolidSDF::Ptr solid,
                            FluidSDF::Ptr fluid,
                            float dt)
{
    PressureSolver::buildLinearSystem(grid, solid, fluid,dt);
    _buildIncompleteCholeskyPreconditioner(fluid);
}

void PCG::solveLinearSystem(FluidSDF::Ptr f, float dt)
{
    float tol = _tol * _b.infNorm();
    _pressure.reset();
    if (_b.infNorm() == 0) {
        return;
    }

    _applyPreconditioner(f);
    _s.copy(_z);
    float rho = _z.dot(_b);
    if (rho == 0) {
        return;
    }

    int iter;
    for (iter = 0; iter < _maxIterations; ++iter) {
        _applyLaplace(f, _s, _z);
        float alpha = rho / _s.dot(_z);
        _pressure.add(_s, alpha);
        _b.add(_z, -alpha);
        if (_b.infNorm() <= tol) {
            std::cout << "Pressure converged, |r| =  " << _b.infNorm() 
                      << " in " << iter << " iterations." << std::endl;
            return;
        }
        _applyPreconditioner(f);
        float rhoNew = _z.dot(_b);
        float beta = rhoNew / rho;
        _s.scaleAndAdd(beta, _z);
        rho = rhoNew;
    }
    std::cout << "Didn't converge in PCG solver, tol = " << tol 
              << ", |r|= " << _b.infNorm() << "." << std::endl;
}

void PCG::_applyPreconditioner(FluidSDF::Ptr f)
{
    _q.reset();
    _z.reset();
    float t;
    // Solve Lq = r
    for (int i = 1; i < _q.nx() ; ++i) {
        for (int j = 1; j < _q.ny() ; ++j) {
            if (f->isFluid(i,j)) {
                t = _b(i,j)- _A.value<LEFT>(i,j) * _precon(i-1,j) * _q(i-1,j)
                           - _A.value<BOTTOM>(i,j) * _precon(i,j-1) * _q(i,j-1);
                _q(i,j) = t * _precon(i,j);
            }
        }
    }

    // Solve L^Tz = q
    for (int i = _z.nx() - 2; i != -1; --i) {
        for (int j = _z.ny() - 2; j != -1; --j) {
            if (f->isFluid(i,j)) {
                t = _q(i,j) - _A.value<RIGHT>(i,j) * _precon(i,j) * _z(i+1,j)
                            - _A.value<TOP>(i,j) * _precon(i,j) * _z(i,j+1);
                _z(i,j) = t * _precon(i,j);
            }
        }
    }
}

void PCG::_applyLaplace(FluidSDF::Ptr f, const Array2f x, Array2f & b)
{
    b.reset();
    for (int i = 0; i < b.nx(); ++i) {
        for (int j = 0; j < b.ny(); ++j) {
            if (f->isFluid(i,j)){
                b(i,j) = _A.mult(x,i,j);
            }
        }
    }
}

void PCG::_buildIncompleteCholeskyPreconditioner(FluidSDF::Ptr f)
{
    const float mic = 0.99;
    const float safety = 0.25;
    float e;
    _precon.reset();
    for (int i = 1; i < _precon.nx(); ++i) {
        for (int j = 1; j < _precon.ny(); ++j) {
            if (f->isFluid(i,j)) {
                const float a = _A.value<CENTER>(i,j);
                const float aii = _A.value<LEFT>(i,j);
                const float aij = _A.value<RIGHT>(i,j-1);
                const float aji = _A.value<TOP>(i-1,j);
                const float ajj = _A.value<BOTTOM>(i,j);
                const float pi = _precon(i-1,j);
                const float pj = _precon(i,j-1);       

                e = a - sqr(aii*pi) - sqr(ajj*pj) - 
                        mic*(aii*aji*sqr(pi) + aij*ajj*sqr(pj));

                if (e < safety * _A.value<CENTER>(i,j) ) {
                    e =  _A.value<CENTER>(i,j);
                }
                _precon(i,j) = 1.0 / sqrt(e + 1e-6);
            }
        }
    }
}
