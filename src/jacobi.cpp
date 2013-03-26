#include "jacobi.h"

Jacobi::Jacobi(Settings::Ptr s) : PressureSolver(s, JACOBI)
{
    _pressure.resize(s->nx,s->ny,s->dx);
    _iterations = s->numJacobiIterations;
}

void Jacobi::solveLinearSystem(FluidSDF::Ptr fluid, float dt)
{
    _pressure.reset();
    _pressureFrom.reset();
    for (int i = 0; i < _iterations; ++i) {
        iteration(fluid->phi(), _A, _b, _pressureFrom, _pressure);
        _pressure.swap(_pressureFrom);
    }
}

void Jacobi::iteration(const Array2f & phi,
                       const SparseLaplacianMatrix<float> & A,
                       const Array2f & b,
                       const Array2f & pFrom,
                       Array2f & p)
{    
    for (int i = 0; i < p.nx(); ++i) {
        for (int j = 0; j < p.ny(); ++j) {
            if (phi(i,j) < 0 && A.value<CENTER>(i,j)){
                const float sum = A.multNeighbors(pFrom,i,j);
                p(i,j) = (b(i,j) - sum) / A.value<CENTER>(i,j);
            }
        }
    }
}

