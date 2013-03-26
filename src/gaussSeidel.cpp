#include "gaussSeidel.h"

GaussSeidel::GaussSeidel(Settings::Ptr s) : PressureSolver(s, GAUSS_SEIDEL)
{
    _iterations = s->numGaussSeidelIterations;
}

void GaussSeidel::solveLinearSystem(FluidSDF::Ptr fluid, float dt)
{
    _pressure.reset();
    for (int i = 0; i < _iterations; ++i) {
        redBlackIteration(true, fluid->phi(), _A, _b, _pressure);
        redBlackIteration(false, fluid->phi(),_A,  _b, _pressure);
    }
}

void GaussSeidel::redBlackIteration(bool red,
                                    const Array2f & phi,
                                    const SparseLaplacianMatrix<float> & A,
                                    const Array2f & b,
                                    Array2f & p)
{
    int start = red ? 0 : 1;
    for (int i = 0; i < p.nx(); ++i) {
        for (int j = start + i % 2; j < p.ny(); j+=2) {
            if (phi(i,j) < 0 && A.value<CENTER>(i,j)){
                p(i,j) = (b(i,j)-A.multNeighbors(p,i,j)) / A.value<CENTER>(i,j);
            }
        }
    }
}

