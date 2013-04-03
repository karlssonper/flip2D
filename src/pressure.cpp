#include "pressure.h"
#include "util.h"
#include "log.h"

PressureSolver::PressureSolver(Settings::Ptr s, SolverType type) : _type(type)
{
  _resize(s->nx,s->ny,s->dx);
}

void PressureSolver::buildLinearSystem(Grid::Ptr grid,
                                       SolidSDF::Ptr solid,
                                       FluidSDF::Ptr fluid,
                                       float dt)
{
    LOG_OUTPUT("Building the linear system for the pressure equation.");
    _buildLaplace(grid->uWeights(), grid->vWeights(), fluid->phi(), _A, dt);
    _buildRHS(grid->u(),grid->v(),grid->uWeights(),grid->vWeights(),fluid,_b);
}

void PressureSolver::_buildLaplace(const FaceArray2Xf & uw,
                                   const FaceArray2Yf & vw,
                                   const Array2f & phi,
                                   SparseLaplacianMatrix<float> & A,
                                   float dt)
{
    A.reset();
    for (int i = 0; i < _pressure.nx(); ++i) {
        for (int j = 0; j < _pressure.ny(); ++j) {
            if (phi(i,j) < 0) {
                if (i > 0) {
                    A.value<CENTER>(i,j) += _laplaceCenter(uw.face<LEFT>(i,j),
                                                           phi(i,j),
                                                           phi(i-1,j));
                }
                if (j > 0) {
                    A.value<CENTER>(i,j) += _laplaceCenter(vw.face<BOTTOM>(i,j),
                                                           phi(i,j),
                                                           phi(i,j-1));
                }

                // We only have to do it on two faces because of symmetri.
                if (i < _pressure.nx() - 1) {
                    A.value<RIGHT>(i,j) = -uw.face<RIGHT>(i,j)*(phi(i+1,j) < 0);
                    A.value<CENTER>(i,j) += _laplaceCenter(uw.face<RIGHT>(i,j),
                                                           phi(i,j),
                                                           phi(i+1,j));
                }
                if (j < _pressure.ny() - 1) {
                    A.value<TOP>(i,j) = -vw.face<TOP>(i,j) * (phi(i,j+1) < 0);
                    A.value<CENTER>(i,j) += _laplaceCenter(vw.face<TOP>(i,j),
                                                           phi(i,j),
                                                           phi(i,j+1));
                }
            }
        }
    }
    const float scale = dt / (sqr(_pressure.dx()));
    A.multiply(scale);
}

float PressureSolver::_laplaceCenter(float w, float phiFluid, float phiAir)
{
    if (phiAir >= 0) {
        float theta = SolidSDF::fractionInside(phiFluid, phiAir);
        if (theta < 0.01) {
            theta = 0.01;
        }
        return w / theta;
    } else {
        return w;
    }
}

void PressureSolver::_buildRHS(const FaceArray2Xf & u,
                               const FaceArray2Yf & v,
                               const FaceArray2Xf & uw,
                               const FaceArray2Yf & vw,
                               FluidSDF::Ptr f,
                               Array2f & b)
{
    b.reset();
    const float scale = 1.0 / _pressure.dx();
    for (int i = 0; i < _pressure.nx(); ++i) {
        for (int j = 0; j < _pressure.ny(); ++j) {
            if (f->isFluid(i,j)) {
                b(i,j) = scale * (u.face<LEFT>(i,j) * uw.face<LEFT>(i,j) -
                                  u.face<RIGHT>(i,j) * uw.face<RIGHT>(i,j) +
                                  v.face<BOTTOM>(i,j) * vw.face<BOTTOM>(i,j) -
                                  v.face<TOP>(i,j) * vw.face<TOP>(i,j));
            }
        }
    }
}

void PressureSolver::_computeResidual(const Array2f & phi,
                                      const SparseLaplacianMatrix<float> & A,
                                      const Array2f & pressure,
                                      const Array2f & b,
                                      Array2f & r)
{
    r.reset();
    for (int i = 0; i < pressure.nx(); ++i) {
        for (int j = 0; j < pressure.ny(); ++j) {
            if (phi(i,j) < 0) {
                r(i,j) = b(i,j) - A.mult(pressure,i,j);
            }
        }
    }
}

void PressureSolver::_resize(int nx, int ny, float dx)
{
    _pressure.resize(nx,ny,dx);
    _b.resize(nx,ny,dx);
    _A.resize(nx,ny,dx);
}
