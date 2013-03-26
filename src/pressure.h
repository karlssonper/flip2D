#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "ptr.h"
#include "settings.h"
#include "sparse.h"
#include "sdf.h"
#include "grid.h"

class PressureSolver : public SmartPtrInterface<PressureSolver>
{
  public:
    typedef SmartPtr<PressureSolver> Ptr;

    enum SolverType
    {
        JACOBI, GAUSS_SEIDEL, PRECONDITIONED_CONJUGATE_GRADIENT, MULTIGRID
    };

    SolverType type() const { return _type; }
    
    virtual void buildLinearSystem(Grid::Ptr grid,
                                   SolidSDF::Ptr solid,
                                   FluidSDF::Ptr fluid,
                                   float dt);
    
    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt) = 0;
    
    const Array2f & pressure() const { return _pressure; }
    
  protected:
    SolverType _type;
    Array2f _pressure;
    Array2f _b;
    SparseLaplacianMatrix<float> _A;
    
    PressureSolver(Settings::Ptr s, SolverType type);
    PressureSolver();
    PressureSolver(const PressureSolver &);
    void operator=(const PressureSolver&);

    void _buildLaplace(const FaceArray2Xf & uWeights,
                       const FaceArray2Yf & vWeights,
                       const Array2f & fluidPhi,
                       SparseLaplacianMatrix<float> & A,
                       float dt);

    void _buildRHS(const FaceArray2Xf & u,
                   const FaceArray2Yf & v,
                   const FaceArray2Xf & uWeights,
                   const FaceArray2Yf & vWeights,
                   FluidSDF::Ptr f,
                   Array2f & b);

    void _computeResidual(const Array2f & phi,
                          const SparseLaplacianMatrix<float> & A,
                          const Array2f & pressure,
                          const Array2f & b,
                          Array2f & r);

    void _resize(int nx, int ny, float dx = 1.0f);
    
    private:
    float _laplaceCenter(float w, float phiFluid, float phiAir);
};

#endif
