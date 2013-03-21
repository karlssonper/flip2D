#include "pressure.h"
#include "util.h"

PressureSolver::PressureSolver(Settings::Ptr s)
{
    _pressure.resize(s->nx,s->ny,s->dx);
    _b.resize(s->nx,s->ny,s->dx);
    _A.resize(s->nx,s->ny,s->dx);
}

void PressureSolver::buildLinearSystem(const FaceArray2Xf & u,
                                       const FaceArray2Xf & v,
                                       const FaceArray2Xf & uWeights,
                                       const FaceArray2Xf & vWeights,
                                       FluidSDF::Ptr f,
                                       float dt)
{
    _buildLaplace(uWeights, vWeights, f, _A, dt);
    _buildRHS(u, v, uWeights, vWeights, f, _b);
}

void PressureSolver::_buildLaplace(const FaceArray2Xf & uw,
                                   const FaceArray2Xf & vw,
                                   FluidSDF::Ptr f,
                                   SparseLaplacianMatrix<float> & A,
                                   float dt)
{
    A.reset();
    for (int i = 0; i < _pressure.nx(); ++i) {
        for (int j = 0; j < _pressure.ny(); ++j) {
            if (f->isFluid(i,j)) {
                if (i > 0) {
                    A.value<CENTER>(i,j) += _laplaceCenter(uw.face<LEFT>(i,j),
                                                           f->phi(i,j),
                                                           f->phi(i-1,j));
                }
                if (j > 0) {
                    A.value<CENTER>(i,j) += _laplaceCenter(vw.face<BOTTOM>(i,j),
                                                           f->phi(i,j),
                                                           f->phi(i,j-1));
                }

                // We only have to do it on two faces because of symmetri.
                if (i < _pressure.nx() - 1) {
                    A.value<RIGHT>(i,j) =-uw.face<RIGHT>(i,j)*f->isFluid(i+1,j);
                    A.value<CENTER>(i,j) += _laplaceCenter(uw.face<RIGHT>(i,j),
                                                           f->phi(i,j),
                                                           f->phi(i+1,j));
                }
                if (j < _pressure.ny() - 1) {
                    A.value<TOP>(i,j) = -vw.face<TOP>(i,j) * f->isFluid(i,j+1);
                    A.value<CENTER>(i,j) += _laplaceCenter(vw.face<TOP>(i,j),
                                                           f->phi(i,j),
                                                           f->phi(i,j+1));
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
                               const FaceArray2Xf & v,
                               const FaceArray2Xf & uw,
                               const FaceArray2Xf & vw,
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
